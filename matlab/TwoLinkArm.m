classdef TwoLinkArm
    
    properties
        dof     % Degrees of freedom
        c       % Number of contact points
        m1      % Mass of each link
        m2
        l1      % Length of each link
        l2
        g
        q_min   % Joint limits
        q_max
        dq_min  % Velocity limits
        dq_max
        u_min   % Control torque limits
        u_max
        alpha   % Constants for equations of motion
        beta
        delta
    end
    
    methods
        function obj = TwoLinkArm(dof, c, m1, m2, l1, l2)
            %TWOLINKARM Constructor. 
            % Creates a robot arm object with dof degrees of freedom,
            % c contact points, and desired masses and link lengths
            %
            % Inputs:   dof - degrees of freedom (2 for twolinkarm)
            %           c   - number of contact points
            %           m1  - mass of link 1 (attached to base)
            %           m2  - mass of link 2 (attached to link 1)
            %           l1  - length of link 1
            %           l2  - length of link 2
            % Output:   obj - a robot arm object
            
            obj.dof = dof;
            obj.c = c;
            obj.m1 = m1;
            obj.m2 = m2;
            obj.l1 = l1;
            obj.l2 = l2;
            
            % Standard acceleration due to gravity
            obj.g = 9.81;
            
            obj.q_min = [0; -pi];
            obj.q_max = [pi; pi];
            obj.dq_min = [-20; -20];
            obj.dq_max = [20; 20];
            obj.u_min = [-10; -10];
            obj.u_max = [10; 10];
            
            % moment of inertia for rod of length l and mass m rotating
            % about its center
            Ix1 = (m1*l1^2)/12.0;
            Ix2 = (m2*l2^2)/12.0;

            % compute constants in the matrix
            obj.alpha = Ix1 + Ix2 + m1*(l1/2.)^2 + m2*(l1^2 + (l2/2.)^2);
            obj.beta = m2*l1*(l2/2.);
            obj.delta = Ix2 + m2*(l2/2.)^2;

        end
        
        % -------------------- BEGIN DYNAMICS ------------------- %
        
        function [M,C,N] = dynamics(obj, q_t1, dq_t1)
            % Returns equation of motion matrices satisfying 
            %   M(q_t1)*ddq_t1 + C(q_t1, dq_t1)*dq_t1 + N(q_t1) = u_t1 
            % 
            % Inputs:   q_t1 - current joint angles of robot (\in R^{dof})
            %           dq_t1 - current joint velocity of robot
            % Outputs:  M   - inertial matrix (\in R^{dof, dof})
            %           C   - coriolis matrix (\in R^{dof, dof})
            %           N   - gravity matrix  (\in R^{dof})
            M = obj.get_M(q_t1);        
            C = obj.get_C(q_t1, dq_t1); 
            N = obj.get_N(q_t1);        
        end
        
        function M = get_M(obj, q_t1)
            % Computes the inertial matrix at time t for 2-link arm
            %
            % Inputs:   q_t1    - configuration of robot
            % Outputs:  M       - intertial matrix
            th2 = q_t1(2);
            M = [obj.alpha + 2*obj.beta*cos(th2) obj.delta + obj.beta*cos(th2);
                 obj.delta + obj.beta*cos(th2)   obj.delta];
        end

        function C = get_C(obj, q_t1, dq_t1)
            %  Computes the coriolis and centrifugal forces acting on the joints
            %
            % Inputs:   q_t1    - configuration of robot
            %           dq_t1   - joint velocities of robot
            % Outputs:  C       - coriolis matrix
            dth1 = dq_t1(1);
            th2 = q_t1(2);
            dth2 = dq_t1(2);

            C = [-obj.beta*sin(th2)*dth2 -obj.beta*sin(th2)*(dth1+dth2);
                 obj.beta*sin(th2)*dth1  0];
        end

        function N = get_N(obj, q_t1)
            %  Computes the effect of gravity on the links
            %
            % Inputs:   q_t1    - configuration of robot
            % Outputs:  N       - gravity matrix
            th1 = q_t1(1);
            th2 = q_t1(2);
            N = [(obj.m1+obj.m2)*obj.g*obj.l1*cos(th1)+obj.m2*obj.g*obj.l2*cos(th1 + th2);
                 obj.m2*obj.g*obj.l2*cos(th1 + th2)];
        end
        
        % -------------------- END DYNAMICS ------------------- %

        function [pos_elbow, pos_ee] = fwd_kinematics(obj, q_t1)
            % Returns the (x,y) position of elbow and end-effector
            %
            % Inputs:   q_t1        - current joint angles of robot
            % Outputs:  pos_elbow   - (x,y) position of elbow
            %           pos_ee      - (x,y) position of end-effector
            th1 = q_t1(1);
            th2 = q_t1(2);
            pos_elbow = [obj.l1*cos(th1); obj.l1*sin(th1)];
            pos_ee = [obj.l1*cos(th1) + obj.l2*cos(th1+th2); obj.l1*sin(th1) + obj.l2*sin(th1+th2)];
        end
        
        function u_c = u_contact(obj, q_t1, lambda_t1)
            % Computes joint torques as a result of contact forces by
            %   F_joints = J^T(q) * F_cartesian
            % 
            % Inputs:   q_t1        - current joint angles of robot
            %           lambda_t1   - magnitude of contact force in y-dir 
            % Outputs:  u_c     - joint torque as a result of contact force
            J_c = obj.get_elbow_jacobian(q_t1);
            F_c = [0; lambda_t1];
            u_c = J_c' * F_c;
        end
        
        function phi = signed_dist(obj, q_t1)
            % Computes signed distance from contact point (elbow) to contact
            % surface (ground-plane)
            %
            % Inputs:   q_t1    - current joint angles of robot
            % Outputs:  phi     - y-dir distance of elbow to ground-plane
            [pos_elbow, ~] = obj.fwd_kinematics(q_t1);
            phi = pos_elbow(2);
        end
        
        function J_c = get_elbow_jacobian(obj, q_t1)
            % Jacoboian of the elbow of arm
            %   J(q) = [dx/dq(1) dx/dq(2);
            %           dy/dq(1) dy/dq(2)]
            %   x = l1*cos(q(1)) 
            %   y = l1*sin(q(2))
            % 
            % Inputs:   q_t1    - current joint angles of the robot
            % Outputs:  J_c     - Jacobian of elbow (i.e. contact point) 
            th1 = q_t1(1);
            J_c = [-obj.l1*sin(th1) 0;
                    obj.l1*cos(th1) 0];
        end
        
        function plot(obj, q, alpha)
            % Plots the arm at configuration q
            [pos_elbow, pos_ee] = obj.fwd_kinematics(q);

            p0 = [0 0];

            link1 = plot([p0(1) pos_elbow(1)], [p0(2) pos_elbow(2)], 'k');
            link2 = plot([pos_elbow(1) pos_ee(1)], [pos_elbow(2) pos_ee(2)], 'k');
            set(link1, 'LineWidth', 5);
            set(link2, 'LineWidth', 5);
            link1.Color(4) = alpha;
            link2.Color(4) = alpha;

            scatter([p0(1) pos_elbow(1)], [p0(2) pos_elbow(2)], 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
            scatter(pos_ee(1), pos_ee(2), 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
        end
        
        function plot_traj(obj, traj, filename)
            % Plots the arm following a trajectory
            %
            % Inputs:   traj    - sequence of configurations over time [q_0, q_1, ..., q_T]^T
            % Outputs:  n/a
            clf
            fig = figure(1);
            hold on;

            for t = 0:size(traj, 1) - 1
                                
                floor = rectangle('Position',[-3 -3 6 3]', 'FaceColor',[0.7 0.7 0.7], ... 
                'EdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
                
                % get the current configuration and elbow/ee positions
                q = traj(t + 1, :)';

                a = (t+1)/(size(traj, 1)+1);
    
                obj.plot(q, a);

                axis([-2.5 2.5 -2.5 2.5]);
    
                set(gca,'YTick',[]);
                set(gca,'XTick',[]);
            
                if ~isempty(filename)
                    frame = getframe(fig);
                    im = frame2im(frame);
                    [imind, cm] = rgb2ind(im, 256); 
                    if t == 0
                        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf); 
                    else
                        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
                    end
                end
    
                if t ~= size(traj, 1)
                    pause(0.1);
                end
            end
            hold off;
        end
        
        function plot_dq_u(obj, traj_dq, traj_u, T)
            figure(2);

            subplot(2, 2, 1);
            plot([1:1:T]', traj_u(:, 1));
            axis([1 T min(obj.u_min) - 1 max(obj.u_max) + 1]);
            xlabel('time $t$', 'Interpreter', 'latex');
            ylabel('$u_1(t)$', 'Interpreter', 'latex');

            subplot(2, 2, 3);
            plot([1:1:T]', traj_u(:, 2));
            axis([1 T min(obj.u_min) - 1 max(obj.u_max) + 1]);
            xlabel('time $t$', 'Interpreter', 'latex');
            ylabel('$u_2(t)$', 'Interpreter', 'latex');

            subplot(2, 2, 2);
            plot([1:1:T+1]', traj_dq(:, 1));
            axis([1 T min(obj.u_min) - 1 max(obj.u_max) + 1]);
            xlabel('time $t$', 'Interpreter', 'latex');
            ylabel('$\dot{q}_1(t)$', 'Interpreter', 'latex');

            subplot(2, 2, 4);
            plot([1:1:T+1]', traj_dq(:, 2));
            axis([1 T min(obj.u_min) - 1 max(obj.u_max) + 1]);
            xlabel('time $t$', 'Interpreter', 'latex');
            ylabel('$\dot{q}_2(t)$', 'Interpreter', 'latex');
        end
        
        function simulate_no_contact(obj, q_init, dq_init, dt, T)
            figure(1);
            
            q = q_init;
            dq = dq_init;

            hold on;

            axis([-3 3 -3 3]);

            set(gca,'YTick',[]);
            set(gca,'XTick',[]);

            [pos_elbow, pos_ee] = obj.fwd_kinematics(q);

            p0 = [0 0];
            link1 = plot([p0(1) pos_elbow(1)], [p0(2) pos_elbow(2)], 'k');
            link2 = plot([pos_elbow(1) pos_ee(1)], [pos_elbow(2) pos_ee(2)], 'k');
            joint1 = scatter([p0(1) pos_elbow(1)], [p0(2) pos_elbow(2)], 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
            joint2 = scatter(pos_ee(1), pos_ee(2), 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');

            for t = 2:T + 1
                [M, C, N] = obj.dynamics(q, dq);
                ddq = M \ (-N - C * dq);
                dq = dq + dt * ddq
                q = q + dt * dq

                [pos_elbow, pos_ee] = obj.fwd_kinematics(q)
                set(link1, 'Xdata', [p0(1) pos_elbow(1)]);
                set(link1, 'Ydata', [p0(2) pos_elbow(2)]);
                set(link2, 'Xdata', [pos_elbow(1) pos_ee(1)]);
                set(link2, 'Ydata',  [pos_elbow(2) pos_ee(2)]);
                set(joint1, 'Xdata', [p0(1) pos_elbow(1)]);
                set(joint1, 'Ydata', [p0(2) pos_elbow(2)]);
                set(joint2, 'Xdata', pos_ee(1));
                set(joint2, 'Ydata', pos_ee(2));

                pause(dt);
            end

            hold off;            
        end
        

        %------------ HELPER FUNCTIONS -----------%

        function fwd_simulate(obj, q_traj, dt)
            figure(1);

            q_curr = q_traj(1,:);

            hold on;
            axis([-3 3 -3 3]);

            p0 = [0 0];
            [pos_elbow, pos_ee] = obj.fwd_kinematics(q_curr);

            link1 = plot([p0(1) pos_elbow(1)], [p0(2) pos_elbow(2)], 'k');
            link2 = plot([pos_elbow(1) pos_ee(1)], [pos_elbow(2) pos_ee(2)], 'k');
            joint1 = scatter([p0(1) pos_elbow(1)], [p0(2) pos_elbow(2)], 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
            joint2 = scatter(pos_ee(1), pos_ee(2), 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');

            for t = 1:size(q_traj, 1)
                q_curr = q_traj(t,:);
                [pos_elbow, pos_ee] = obj.fwd_kinematics(q_curr);
                set(link1, 'Xdata', [p0(1) pos_elbow(1)]);
                set(link1, 'Ydata', [p0(2) pos_elbow(2)]);
                set(link2, 'Xdata', [pos_elbow(1) pos_ee(1)]);
                set(link2, 'Ydata',  [pos_elbow(2) pos_ee(2)]);
                set(joint1, 'Xdata', [p0(1) pos_elbow(1)]);
                set(joint1, 'Ydata', [p0(2) pos_elbow(2)]);
                set(joint2, 'Xdata', pos_ee(1));
                set(joint2, 'Ydata', pos_ee(2));

                pause(dt);
            end

            hold off;            
        end
        
    end
end


