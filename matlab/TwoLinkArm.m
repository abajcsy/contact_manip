classdef TwoLinkArm
    %TWOLINKARM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dof % Degrees of freedom
        c   % Number of contact points
        m1
        m2
        l1
        l2
        I1x
        I2x
        g
        q_min
        q_max
        dq_min
        dq_max
        u_min
        u_max
    end
    
    methods
        function obj = TwoLinkArm(dof, c, m1, m2, l1, l2)
            %TWOLINKARM Construct an instance of this class
            %   Detailed explanation goes here
            obj.dof = dof;
            obj.c = c;
            obj.m1 = m1;
            obj.m2 = m2;
            obj.l1 = l1;
            obj.l2 = l2;
            
            % Approximate the links as thin rods
            obj.I1x = (obj.m1 * obj.l1^2) / 12.0;
            obj.I2x = (obj.m2 * obj.l2^2) / 12.0;
            
            % Standard acceleration due to gravity
            obj.g = 9.81;
            
            obj.q_min = [-pi / 2; -pi / 2];
            obj.q_max = [pi / 2; pi / 2];
            obj.dq_min = [-20; -20];
            obj.dq_max = [20; 20];
            obj.u_min = [-5; -5];
            obj.u_max = [5; 5];
        end
        
        function [M, C, N] = dynamics(obj, q, dq)
            % Returns the manipulator equation matrices
            %   M(q)*ddq + C(q, dq)*dq + N(q) = u
            M = two_link_inertia_matrix(obj.I1x, obj.I2x, obj.l1, obj.l2, obj.m1, obj.m2, q(1), q(2));
            C = two_link_coriolis_matrix(dq(1), dq(2), obj.l1, obj.l2, obj.m2, q(1), q(2));
            N = two_link_gravity_matrix(obj.g, obj.l1, obj.l2, obj.m1, obj.m2, q(1), q(2));
        end
        
        function [pos_elbow, pos_ee] = fwd_kinematics(obj, q_t1)
            % Returns the (x,y) position of elbow and end-effector
            th1 = q_t1(1);
            th2 = q_t1(2);
            pos_elbow = [obj.l1*cos(th1); obj.l1*sin(th1)];
            pos_ee = [obj.l1*cos(th1) + obj.l2*cos(th1+th2); obj.l1*sin(th1) + obj.l2*sin(th1+th2)];
        end
        
        function u_c = u_contact(obj, q_t1, lambda_t1)
            % Computes joint torques as a result of contact forces
            J_c = obj.get_elbow_jacobian(q_t1);
            F_c = [0; lambda_t1];
            u_c = J_c' * F_c;
        end
        
        function phi = signed_dist(obj, q_t1)
            % Computes signed distance from contact point (elbow) to contact
            % surface (ground-plane)
            [pos_elbow, ~] = obj.fwd_kinematics(q_t1);
            % get the y-distance
            phi = pos_elbow(2);
        end
        
        function J_c = get_elbow_jacobian(obj, q_t1)
            % Jacoboian of the elbow of arm
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
        
        function plot_traj(obj, traj)
            clf
            fig = figure(1);
            hold on;

            for t = 0:size(traj, 1) - 1
                q = traj(t + 1, :)';

                alpha = (t+1)/(size(traj, 1)+1);
    
                obj.plot(q, alpha);

                axis([-2 2 0 3]);
    
                set(gca,'YTick',[]);
                set(gca,'XTick',[]);
    
                if t ~= size(traj, 1)
                    pause(0.1);
                end
            end
            hold off;
        end
        
        function simulate(obj, q_init, dq_init, dt, T)
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
                dq = dq + dt * ddq;
                q = q + dt * dq;

                [pos_elbow, pos_ee] = obj.fwd_kinematics(q);
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

