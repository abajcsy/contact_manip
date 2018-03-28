classdef TwoLinkArm
    
    properties
        dof     % Degrees of freedom
        c       % Number of contact points
        m1      % Mass of each link
        m2
        l1      % Length of each link
        l2
        q_min   % Joint limits
        q_max
        dq_min  % Velocity limits
        dq_max
        u_min   % Control torque limits
        u_max
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
            
            obj.q_min = [-pi / 2; -pi / 2];
            obj.q_max = [pi / 2; pi / 2];
            obj.dq_min = [-20; -20];
            obj.dq_max = [20; 20];
            obj.u_min = [-5; -5];
            obj.u_max = [5; 5];
        end
        
        function [M,C,N] = dynamics(obj, q_t1, dq_t1)
            % Returns equation of motion matrices satisfying 
            %   M(q_t1)*ddq_t1 + C(q_t1, dq_t1)*dq_t1 + N(q_t1) = u_t1 
            % 
            % Inputs:   q_t1 - current joint angles of robot (\in R^{dof})
            %           dq_t1 - current joint velocity of robot
            % Outputs:  M   - inertial matrix (\in R^{dof, dof})
            %           C   - coriolis matrix (\in R^{dof, dof})
            %           N   - gravity matrix  (\in R^{dof})
            M = zeros(obj.dof);
            C = zeros(obj.dof);
            N = zeros(obj.dof,1);
        end
        
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
        
        function plot_traj(obj, traj)
            % Plots the arm following a trajectory
            %
            % Inputs:   traj    - sequence of configurations over time [q_0, q_1, ..., q_T]^T
            % Outputs:  n/a
            clf
            figure(1);
            hold on;

            for t = 0:size(traj, 1) - 1
                % get the current configuration and elbow/ee positions
                q = traj(t + 1, :)';
                [pos_elbow, pos_ee] = obj.fwd_kinematics(q);

                % robot base position
                p0 = [0 0];
                % fade trajectory from light (t=0) to saturated (t=T) 
                alpha = (t+1)/(size(traj, 1)+1);
    
                % plot the robot links
                link1 = plot([p0(1) pos_elbow(1)], [p0(2) pos_elbow(2)], 'k');
                link2 = plot([pos_elbow(1) pos_ee(1)], [pos_elbow(2) pos_ee(2)], 'k');
                set(link1,'LineWidth',5);
                set(link2,'LineWidth',5);
                link1.Color(4) = alpha;
                link2.Color(4) = alpha;

                % plot the robot joints
                scatter([p0(1) pos_elbow(1)], [p0(2) pos_elbow(2)], 'o', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
                scatter(pos_ee(1), pos_ee(2), 'o', 'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    
                axis([-2 2 0 3]);
    
                set(gca,'YTick',[]);
                set(gca,'XTick',[]);
    
                if t ~= size(traj, 1)
                    pause(0.1);
                end
            end
            hold off;
        end
        
    end
end

