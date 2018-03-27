classdef OptProb
    %OPTPROB This class defines an optimization for a two link arm
    %   Detailed explanation goes here
    
    properties
        arm
        q_init
        q_final
        T           % Number of time steps (0..T)
        num_states
        num_controls
        num_lambdas
        lambda_max
        h_min
        h_max
        g           % Function that is used to compute objective
        g_f       
        k           % Number of timesteps to consider when computing EE dist from final pose
        eps         % Tolerance on EE dist from final pose
    end
    
    methods
        function obj = OptProb(arm, q_init, q_final, T, g, g_f, k, eps)
            %OPTPROB Construct an instance of this class
            %   
            %   g and g_f are used compute objective: 
            %       g_f(q_T, \dot{q}_T) + h \sum^T-1_t=0 g(q_t, u_{t+1})
            obj.arm = arm;
            % TODO right now, q_init, q_final have both joint position and
            % velocities, change name?
            obj.q_init = q_init;
            obj.q_final = q_final;
            obj.T = T;
            
            dof = obj.arm.dof;
            c = obj.arm.c;

            obj.num_states = 2 * dof * (obj.T + 1);
            obj.num_controls = dof * obj.T;
            obj.num_lambdas = c * obj.T;
            
            obj.lambda_max = 1000;
            obj.h_min = 0.01;
            obj.h_max = 0.2;
            
            obj.g = g;
            obj.g_f = g_f;
            
            obj.k = k;
            obj.eps = eps;
        end
        
        function [x, xlow, xupp, F, Flow, Fupp] = generate(obj)
            %GENERATE 
            %   Detailed explanation goes here
            
            % x = [q_0, \dot{q}_0, ..., q_T, \dot{q}_T | u_1, ..., u_T | \lambda_1, ..., \lambda_T | h_1, ..., h_T]
            % where q_i, \dot{q}_i \in \mathbb{R}^{dof}, u_i \in
            % \mathbb{R}^{dof}, \lambda_i \in \mathbb{R}^{c}, h_i \in
            % \mathbb{R}
            % dof = number of degrees of freedom 
            % c = number of contact points
            
            x = zeros(obj.num_states + obj.num_controls + obj.num_lambdas + obj.T, 1);
            xlow = zeros(obj.num_states + obj.num_controls + obj.num_lambdas + obj.T, 1);
            xupp = zeros(obj.num_states + obj.num_controls + obj.num_lambdas + obj.T, 1);
            
            for t = 0:obj.T
                q_interp = obj.lin_interp(obj.q_init(1:obj.arm.dof), obj.q_final(1:obj.arm.dof), t);
                x = obj.set_q(x, t, q_interp);
                
                % TODO interpolate over joint velocities?
                
                if t == 0
                    % Constrain the initial joint position and velocity
                    xlow = obj.set_q(xlow, t, obj.q_init(1:obj.arm.dof));
                    xupp = obj.set_q(xupp, t, obj.q_init(1:obj.arm.dof));
                    xlow = obj.set_dq(xlow, t, obj.q_init(obj.arm.dof + 1:2 * obj.arm.dof));
                    xupp = obj.set_dq(xupp, t, obj.q_init(obj.arm.dof + 1:2 * obj.arm.dof));
                else
                    % Set the joint limits
                	xlow = obj.set_q(xlow, t, obj.arm.q_min);
                    xupp = obj.set_q(xupp, t, obj.arm.q_max);
                    xlow = obj.set_dq(xlow, t, obj.arm.dq_min);
                    xupp = obj.set_dq(xupp, t, obj.arm.dq_max);
                end
                
                if t > 0
                    % Set the torque limits
                    xlow = obj.set_u(xlow, t, obj.arm.u_min);
                    xupp = obj.set_u(xupp, t, obj.arm.u_max);
                    
                    xupp = obj.set_lambda(xupp, t, obj.lambda_max);
                    xlow = obj.set_h(xlow, t, obj.h_min);
                    xupp = obj.set_h(xupp, t, obj.h_max);
                end
            end
            
            function f = func(x)
                % compute final cost
                q_T = obj.get_q(x,obj.T);
                dq_T = obj.get_dq(x,obj.T);
                objective = obj.g_f(q_T, dq_T);
                
                % stores q_{t+1} = q_{t} + h * \dot{q}_{t+1}
                kin_constraints = zeros(obj.arm.dof * obj.T, 1);
                
                % stores equation of motion constraints 
                dyn_constraints = zeros(obj.arm.dof * obj.T, 1);
                
                % signed distance constraints
                phi_constraints = zeros(obj.T * obj.arm.c, 1);
                
                % enforce that phi' * lambda = 0
                phi_lambda_constraints = zeros(obj.T, 1);
                
                % enforce that h_{t-1} + h_{t} = h_{t+1} + h_{t+2}
                num_h_constraints = floor((obj.T-3)/2);
                h_constraints = zeros(num_h_constraints, 1);

                for t=0:obj.T-1
                    
                    h_t1 = obj.get_h(x,t+1);
                    q_t = obj.get_q(x,t);
                    q_t1 = obj.get_q(x,t+1);
                    dq_t = obj.get_dq(x,t);
                    dq_t1 = obj.get_dq(x,t+1);
                    ddq_t1 = (dq_t1 - dq_t)/h_t1;
                    u_t1 = obj.get_u(x,t+1);
                    lambda_t1 = obj.get_lambda(x,t+1);
                    
                    [M,C,N] = obj.arm.dynamics(q_t1,dq_t1);
                    u_contact = obj.arm.u_contact(q_t1,lambda_t1);
                    
                    phi = obj.arm.signed_dist(q_t1);
                    
                    % compute running cost
                    objective = objective + h_t1*obj.g(q_t,dq_t,u_t1);
                    
                    % compute kinematic and dynamics constraints
                    kin_constraints(obj.arm.dof*t+1 : obj.arm.dof*(t+1)) = q_t - q_t1 + h_t1 * dq_t1;
                    dyn_constraints(obj.arm.dof*t+1 : obj.arm.dof*(t+1)) = M*ddq_t1 + C*dq_t1 + N - u_t1 - u_contact;
                    
                    % compute signed distance from contact point to contact
                    % surface
                    % TODO If c > 1 this indexing will not work!
                    phi_constraints(t+1) = phi;
                    phi_lambda_constraints(t+1) = phi' * lambda_t1;
                end
                
                % compute all pairwise h constraints
                for j=1:num_h_constraints
                    h1 = obj.get_h(x, 2*j-1);
                    h2 = obj.get_h(x, 2*j);
                    h3 = obj.get_h(x, 2*j+1);
                    h4 = obj.get_h(x, 2*j+2);

                    h_constraints(j) = h1+h2-h3-h4;
                end
                
                [~,target_ee] = obj.arm.fwd_kinematics(obj.q_final);
                ee_dist = 0;
                % sum up the cost of the EE deviating from final pose
                for i=(obj.T-obj.k):obj.T
                    [~,pos_ee] = obj.arm.fwd_kinematics(obj.get_q(x,i));
                    ee_dist = ee_dist + norm(target_ee - pos_ee);
                end
                
                f = [objective; kin_constraints; dyn_constraints; phi_constraints; phi_lambda_constraints; h_constraints; ee_dist];
            end
            
            F = @func;
            
            % objective                 1
            % kin_constraints           dof * T
            % dyn_constraints           dof * T
            % phi_constraints           T * c
            % phi_lambda_constraints    T
            % h_constraints             floor((T - 3) / 2)
            % ee_dist                   1
            num_constraints = 2 + 2 * obj.arm.dof * obj.T + obj.arm.c * obj.T + obj.T + floor((obj.T-3)/2);
            Flow = zeros(num_constraints, 1);
            Fupp = zeros(num_constraints, 1);
            
            % Bounds on the objective
            Flow(1) = -1000;
            Fupp(1) = 1000;
            
            % Tolerance on EE dist to final pose
            Flow(num_constraints) = 0;
            Fupp(num_constraints) = obj.eps;
            
            % Upper bound on the signed distance 
            Fupp(1 + 2 * obj.arm.dof * obj.T + 1:1 + 2 * obj.arm.dof * obj.T + obj.arm.c * obj.T) = 50;
        end
        
        function p = lin_interp(obj, p_start, p_end, t)
            % Interpolates a single point at time t
            p = p_start + (t / obj.T) * (p_end - p_start);
        end
        
        function q = get_q(obj, x, t)
            q = x(2 * obj.arm.dof * t + 1:2 * obj.arm.dof * t + obj.arm.dof);
        end
        
        function x = set_q(obj, x, t, q)
            x(2 * obj.arm.dof * t + 1:2 * obj.arm.dof * t + obj.arm.dof) = q;
        end
        
        function dq = get_dq(obj, x, t)
            dq = x(2 * obj.arm.dof * t + obj.arm.dof + 1:2 * obj.arm.dof * t + 2 * obj.arm.dof);
        end
        
        function x = set_dq(obj, x, t, dq)
            x(2 * obj.arm.dof * t + obj.arm.dof + 1:2 * obj.arm.dof * t + 2 * obj.arm.dof) = dq;
        end
        
        function u = get_u(obj, x, t)
            u_start = obj.num_states + 1;
            u = x(u_start + obj.arm.dof * (t - 1):u_start + obj.arm.dof * t - 1);
        end
        
        function x = set_u(obj, x, t, u)
            u_start = obj.num_states + 1;
            x(u_start + obj.arm.dof * (t - 1):u_start + obj.arm.dof * t - 1) = u;
        end
        
        function lambda = get_lambda(obj, x, t)
            lambda_start = obj.num_states + obj.num_controls + 1;
            lambda = x(lambda_start + obj.arm.c * (t - 1):lambda_start + obj.arm.c * t - 1);
        end
        
        function x = set_lambda(obj, x, t, lambda)
            lambda_start = obj.num_states + obj.num_controls + 1;
            x(lambda_start + obj.arm.c * (t - 1):lambda_start + obj.arm.c * t - 1) = lambda;
        end        
        
        function h = get_h(obj, x, t)
            h_start = obj.num_states + obj.num_controls + obj.num_lambdas + 1;
            h = x(h_start + t - 1);
        end
        
        function x = set_h(obj, x, t, h)
            h_start = obj.num_states + obj.num_controls + obj.num_lambdas + 1;
            x(h_start + t - 1) = h;
        end
    end
end

