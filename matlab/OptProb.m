classdef OptProb
    
    properties
        arm             % Robot arm object
        q_init          % Starting robot configuration (q,dq)
        q_final         % Final robot configuration (q,dq)
        T               % Number of time steps (0..T)
        num_states      % Number of states (q,dq) in the optimization variable
        num_controls    % Number of controls (u) 
        num_lambdas     % Number of force lambda's
        lambda_max      % Upper bound on lambda
        h_min           % Minimum timestep
        h_max           % Maximum timestep
        g           % Function used to compute running cost 
        g_f         % Function used to compute final cost
    end
    
    methods
        function obj = OptProb(arm, q_init, q_final, T, g, g_f)
            %OPTPROB Constructor. 
            % Creates and optimization problem object that generates the
            % initial optimization variable, constraints, and relevant
            % bounds for both. 
            % 
            % Inputs:   arm     - robot arm object
            %           q_init  - initial robot state
            %           q_final - final robot state
            %           T       - number of timesteps (0...T)
            %           g       - used to compute running cost: h*\sum^T-1_t=0 g(q_t, u_{t+1})
            %           g_f     - used to compute final cost: g_f(q_T, \dot{q}_T) 
            % Output:   obj     - optimization problem object
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
        end
        
        function [x, xlow, xupp, F, Flow, Fupp] = generate(obj)
            % Generates the optimization variable x, its bounds, and the
            % constraints F as well as their bounds
            % 
            % Inputs:   n/a
            % Outputs:  x       - optimization variable that consists of:
            %                   [q_0, \dot{q}_0, ..., q_T, \dot{q}_T | u_1, ..., u_T | \lambda_1, ..., \lambda_T | h_1, ..., h_T]
            %           xlow    - lower bound for optimization variable
            %           xupp    - upper bound for optimization variable
            %           F       - function that computes constraint vector
            %           Flow    - lower bounds on constraint vector
            %           Fupp    - upper bounds on constraint vector
            %   where:
            %           q_i \in R^{dof} 
            %           \dot{q}_i \in R^{dof}
            %           u_i \in R^{dof}
            %           \lambda_i \in R^{c}
            %           h_i \in \mathbb{R}
            %   and:
            %           dof     - robot degrees of freedom 
            %           c       - number of contact points with robot
            
            x_length = obj.num_states + obj.num_controls + obj.num_lambdas + obj.T;
            x = zeros(x_length, 1);
            xlow = zeros(x_length, 1);
            xupp = zeros(x_length, 1);
            
            for t = 0:obj.T
                % interpolate robot position of each timestep
                q_interp = obj.lin_interp(obj.q_init(1:obj.arm.dof), obj.q_final(1:obj.arm.dof), t);
                x = obj.set_q(x, t, q_interp);
                
                % TODO interpolate over joint velocities?
                % TODO initialize time steps?
                
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
                    % Set the lambda and timestep limits
                    xupp = obj.set_lambda(xupp, t, obj.lambda_max);
                    xlow = obj.set_h(xlow, t, obj.h_min);
                    xupp = obj.set_h(xupp, t, obj.h_max);
                end
            end
            
            function f = func(x)
                % Generates the constraint vector, F(x).
                %
                % Input:    x   - optimization variable 
                % Output:   f   - constraint vector that contains the
                %                 objective, kinematic and dynamic 
                %                 constraints for all timesteps, the phi 
                %                 (signed-distance to obstacles), the 
                %                 phi-lambda constraints, and the timestep
                %                 constraints
                
                % compute final cost
                q_T = obj.get_q(x,obj.T);
                dq_T = obj.get_dq(x,obj.T);
                objective = obj.g_f(q_T, dq_T, obj.arm);
                
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
                    
%                     h_t1 = obj.get_h(x,t+1);
                    % DEBUGGING (Ellis) -- fix the time step to 0.1
                    h_t1 = 0.1;
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
                    objective = objective + h_t1*obj.g(q_t,dq_t,u_t1,t+1,obj.T);
                    
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
                
                f = [objective; kin_constraints; dyn_constraints; phi_constraints; phi_lambda_constraints; h_constraints];
            end
            
            F = @func;
            
            % Size of the constraint vector: 
            %   objective                 1
            %   kin_constraints           dof * T
            %   dyn_constraints           dof * T
            %   phi_constraints           T * c
            %   phi_lambda_constraints    T
            %   h_constraints             floor((T - 3) / 2)
            num_constraints = 1 + 2 * obj.arm.dof * obj.T + obj.arm.c * obj.T + obj.T + floor((obj.T-3)/2);
            Flow = zeros(num_constraints, 1);
            Fupp = zeros(num_constraints, 1);
            
            % Bounds on the objective
            Flow(1) = -1000;
            Fupp(1) = 1000;
            
            % Upper bound on the signed distance 
            Fupp(1 + 2 * obj.arm.dof * obj.T + 1:1 + 2 * obj.arm.dof * obj.T + obj.arm.c * obj.T) = 50;
        end
        
        function p = lin_interp(obj, p_start, p_end, t)
            % Interpolates a single point at time t
            %
            % Input:    p_start     - start point of line
            %           p_end       - end point of line
            %           t           - desired interpolation time
            % Output:   P           - interpolated point at time t
            p = p_start + (t / obj.T) * (p_end - p_start);
        end
        
        function q = get_q(obj, x, t)
            % Gets the configuration at time t from the optimization variable
            %
            % Input:    x   - optimization vector
            %           t   - time point between [0,...T]
            % Output:   q   - robot configuration at time t
            q = x(2 * obj.arm.dof * t + 1:2 * obj.arm.dof * t + obj.arm.dof);
        end
        
        function x = set_q(obj, x, t, q)            
            % Sets the configuration at time t in vector x
            %
            % Input:    x   - optimization vector
            %           t   - time point between [0,...T]
            %           q   - new configuration 
            % Output:   x   - updated optimization vector with q_t inside
            x(2 * obj.arm.dof * t + 1:2 * obj.arm.dof * t + obj.arm.dof) = q;
        end
        
        function dq = get_dq(obj, x, t)
            % Gets the joint velocity at time t from the optimization variable
            %
            % Input:    x   - optimization vector
            %           t   - time point between [0,...T]
            % Output:   dq  - robot joint velocity at time t
            dq = x(2 * obj.arm.dof * t + obj.arm.dof + 1:2 * obj.arm.dof * t + 2 * obj.arm.dof);
        end
        
        function x = set_dq(obj, x, t, dq)
            % Sets the joint velocity at time t in vector x
            %
            % Input:    x   - optimization vector
            %           t   - time point between [0,...T]
            %           dq  - new joint velocity
            % Output:   x   - updated optimization vector with dq_t inside
            x(2 * obj.arm.dof * t + obj.arm.dof + 1:2 * obj.arm.dof * t + 2 * obj.arm.dof) = dq;
        end
        
        function u = get_u(obj, x, t)
            % Gets the control at time t from the optimization variable
            %
            % Input:    x   - optimization vector
            %           t   - time point between [0,...T]
            % Output:   u   - robot control at time t
            u_start = obj.num_states + 1;
            u = x(u_start + obj.arm.dof * (t - 1):u_start + obj.arm.dof * t - 1);
        end
        
        function x = set_u(obj, x, t, u)
            % Sets the control at time t in vector x
            %
            % Input:    x   - optimization vector
            %           t   - time point between [0,...T]
            %           u   - new control
            % Output:   x   - updated optimization vector with u_t inside
            u_start = obj.num_states + 1;
            x(u_start + obj.arm.dof * (t - 1):u_start + obj.arm.dof * t - 1) = u;
        end
        
        function lambda = get_lambda(obj, x, t)
            % Gets the lambda at time t from the optimization variable
            %
            % Input:    x       - optimization vector
            %           t       - time point between [0,...T]
            % Output:   lambda  - lambda at time t
            lambda_start = obj.num_states + obj.num_controls + 1;
            lambda = x(lambda_start + obj.arm.c * (t - 1):lambda_start + obj.arm.c * t - 1);
        end
        
        function x = set_lambda(obj, x, t, lambda)
            % Sets the lambda at time t in vector x
            %
            % Input:    x       - optimization vector
            %           t       - time point between [0,...T]
            %           lambda  - new lambda
            % Output:   x       - updated optimization vector
            lambda_start = obj.num_states + obj.num_controls + 1;
            x(lambda_start + obj.arm.c * (t - 1):lambda_start + obj.arm.c * t - 1) = lambda;
        end        
        
        function h = get_h(obj, x, t)
            % Gets the timestep at time t from the optimization variable
            %
            % Input:    x   - optimization vector
            %           t   - time point between [0,...T]
            % Output:   h   - current timestep to take
            h_start = obj.num_states + obj.num_controls + obj.num_lambdas + 1;
            h = x(h_start + t - 1);
        end
        
        function x = set_h(obj, x, t, h)
            % Sets the timestep h at time t in vector x
            %
            % Input:    x   - optimization vector
            %           t   - time point between [0,...T]
            %           h   - new timestep
            % Output:   x   - updated optimization vector with h_t inside
            h_start = obj.num_states + obj.num_controls + obj.num_lambdas + 1;
            x(h_start + t - 1) = h;
        end
        
        function traj = get_traj(obj, x)
            % Extract just the joint angles from optimization variable
            % 
            % Input:   x    - optimization vector
            % Output   traj - configs over time [q_0, q_1, ... , q_T]^T
            traj = zeros(obj.T + 1, obj.arm.dof);
            for t = 0:obj.T
               q = obj.get_q(x, t);
               traj(t + 1, :) = q';
            end
        end
        
        function traj_dq = get_traj_dq(obj, x)
            % Extract just the joint velocities from the solution
            traj_dq = zeros(obj.T + 1, obj.arm.dof);
            for t = 0:obj.T
               dq = obj.get_dq(x, t);
               traj_dq(t + 1, :) = dq';
            end
        end
        
        function traj_u = get_traj_u(obj, x)
            % Extract just the joint torques (control inputs) from the solution
            traj_u = zeros(obj.T, obj.arm.dof);
            for t = 1:obj.T
               u = obj.get_u(x, t);
               traj_u(t, :) = u';
            end            
        end
        
        function traj_h = get_traj_h(obj, x)
            % Extract just the timsteps from the solution
            traj_h = zeros(obj.T, 1);
            for t = 1:obj.T
               h = obj.get_h(x, t);
               traj_h(t) = h;
            end          
        end
        

    end
end

