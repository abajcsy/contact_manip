classdef OptProb
    %OPTPROB This class defines an optimization for a two link arm
    %   Detailed explanation goes here
    
    properties
        arm
        q_init
        q_final
        T % Number of time steps (0..T)
        num_states
        num_controls
        num_lambdas
        lambda_max
        h_min
        h_max
    end
    
    methods
        function obj = OptProb(arm, q_init, q_final, T)
            %OPTPROB Construct an instance of this class
            %   Detailed explanation goes here
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
                f = [0];
            end
            
            F = @func;
            
            Flow = [];
            Fupp = [];
        end
        
        function p = lin_interp(obj, p_start, p_end, t)
            % TODO implement linear interpolation
            p = zeros(size(p_start));
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

