classdef TwoLinkArmWithDoor

    properties
        dof     % Degrees of freedom (2 for arm + 1 for door)
        c       % Number of contact points
        m1      % Mass of each link
        m2
        m3 
        l1      % Length of each link
        l2
        l3 
        g       % Gravity constant
        q_min   % Joint limits
        q_max
        dq_min  % Velocity limits
        dq_max
        u_min   % Control torque limits
        u_max
        alpha   % Constants for equations of motion
        beta
        delta
        I_door  % Moment of inertia of the door
        damping % Damping coefficient of door
        door_hinge 
    end
    
    methods
        function obj = TwoLinkArmWithDoor(dof, c, m1, m2, m3, l1, l2, l3)
            %TWOLINKARMWITHDOOR Constructor. 
            
            obj.dof = dof;
            obj.c = c;
            obj.m1 = m1;
            obj.m2 = m2;
            obj.m3 = m3;
            obj.l1 = l1;
            obj.l2 = l2;
            obj.l3 = l3;
            
            % Standard acceleration due to gravity
            obj.g = 9.81;
            
            obj.q_min = [0; -pi / 2; -pi];
            obj.q_max = [pi; pi / 2; 0];
            obj.dq_min = [-1; -1; -1];
            obj.dq_max = [1; 1; 1];
            obj.u_min = [-5; -5; 0];
            obj.u_max = [5; 5; 0];
            
            % moment of inertia for rod of length l and mass m rotating
            % about its center
            Ix1 = (m1*l1^2)/12.0;
            Ix2 = (m2*l2^2)/12.0;

            % compute constants in the matrix
            obj.alpha = Ix1 + Ix2 + m1*(l1/2.)^2 + m2*(l1^2 + (l2/2.)^2);
            obj.beta = m2*l1*(l2/2.);
            obj.delta = Ix2 + m2*(l2/2.)^2;
            
            obj.I_door = (m3*l3^2)/12.0;
            obj.damping = 2.0;
            
            obj.door_hinge = [-l1/2.0; l1+l2+l3/2.0];
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
            M = [obj.alpha + 2*obj.beta*cos(th2) obj.delta + obj.beta*cos(th2) 0;
                 obj.delta + obj.beta*cos(th2)   obj.delta 0;
                 0 0 obj.I_door];
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

            C = [-obj.beta*sin(th2)*dth2 -obj.beta*sin(th2)*(dth1+dth2) 0;
                 obj.beta*sin(th2)*dth1  0 0;
                 0 0 obj.damping];
        end

        function N = get_N(obj, q_t1)
            %  Computes the effect of gravity on the links
            %
            % Inputs:   q_t1    - configuration of robot
            % Outputs:  N       - gravity matrix
            th1 = q_t1(1);
            th2 = q_t1(2);
            N = [(obj.m1+obj.m2)*obj.g*obj.l1*cos(th1)+obj.m2*obj.g*obj.l2*cos(th1 + th2);
                 obj.m2*obj.g*obj.l2*cos(th1 + th2);
                 0 ];
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
        
        function pos_door = fwd_kin_door(obj, q_t1)
            % Returns (x,y) position of the door.
            
            th3 = q_t1;
            pos_door = obj.door_hinge + [obj.l3*cos(th3); obj.l3*sin(th3)];
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
        
        function door_contact_vec = get_contact_vec(obj, th3)
            door_contact_vec = [(obj.l3/2)*cos(th3); (obj.l3/2)*sin(th3); 0];
        end
        
        function u_c = u_contact(obj, q_t1, lambda_t1)
            % Computes joint torques as a result of contact forces by
            %   F_joints = J^T(q) * lambda_dor
            % 
            % Inputs:   q_t1        - current joint angles of robot
            %           lambda_t1   - magnitude of contact force in y-dir 
            % Outputs:  u_c     - joint torque as a result of contact force
            J_c = obj.get_elbow_jacobian(q_t1);
            th3 = q_t1(3);
            F_c = [-lambda_t1*sin(th3); lambda_t1*cos(th3)];
            u_r = J_c' * F_c; % joint torques experienced by robot
            
            F_door = [lambda_t1*sin(th3); -lambda_t1*cos(th3); 0];
            door_contact_vec = obj.get_contact_vec(th3);
            u_d = norm(cross(F_door, door_contact_vec)); % joint torques experienced by door
            u_c = [u_r; u_d];
        end
        
        function dist_link_door = signed_dist_contact(obj, q)
            % Returns signed distance between the first link (which is
            %   allowed to make contact) and the door
            [pos_elbow, ~] = obj.fwd_kinematics(q);
            pos_base = [0; 0];
            pos_door_tip = obj.fwd_kin_door(q(3));
            % Compute intersection between link 1 and the door
            [s1, s2] = obj.intersection(pos_base, pos_elbow, obj.door_hinge, pos_door_tip);
            % TODO (Ellis) make max signed distance a parameter
            dist_link_door = max([abs(s1 - 0.5) - 0.5, abs(s2 - 0.5) - 0.5, 500]);
        end
        
        function [dist_floor_elbow, dist_ceil_elbow, dist_floor_ee, dist_ceil_ee, dist_link_door] = signed_dist_no_contact_separate(obj, q)
            % Returns signed distances from the floor and ceiling planes to
            %   the elbow and the end effector, and the signed distance
            %   between door and the second link of the arm (which is not
            %   allowed to make contact)
            [pos_elbow, pos_ee] = obj.fwd_kinematics(q);
            
            dist_floor_elbow = pos_elbow(2);
            dist_floor_ee = pos_ee(2);
            
            ceil_y = obj.door_hinge(2);
            dist_ceil_elbow = ceil_y - pos_elbow(2);
            dist_ceil_ee = ceil_y - pos_elbow(2);
            
            pos_door_tip = obj.fwd_kin_door(q(3));
            [link_w_door_s1, link_w_door_s2] = obj.intersection(pos_elbow, pos_ee, obj.door_hinge, pos_door_tip);
            if link_w_door_s1 == Inf || link_w_door_s2 == Inf
                % TODO should set to some large value that is still within
                % the variable bounds for the optimization problem
                dist_link_door = 500;
            else
                dist_link_door = max(abs(link_w_door_s1 - 0.5) - 0.5, abs(link_w_door_s2 - 0.5) - 0.5);
            end
        end
        
        function dist = signed_dist_no_contact(obj, q)
            [dist_floor_elbow, dist_ceil_elbow, dist_floor_ee, dist_ceil_ee, dist_link_door] = obj.signed_dist_no_contact_separate(q);
            % TODO (Ellis) make max signed distance a parameter
            dist = min([dist_floor_elbow, dist_ceil_elbow, dist_floor_ee, dist_ceil_ee, dist_link_door, 500]);
        end
        
        function [s1, s2] = intersection(obj, l1_low, l1_high, l2_low, l2_high)
            % Checks for intersection of line segments (l1_low, l1_high)
            %   and (l2_low, l2_high) (some degenerate cases, such as one
            %   or both of the segments being points, are not handled)
            %
            % Inputs: l1_low  - start point of the first line segment
            %         l1_high - end point of the first line segment
            %         l2_low  - start point of the second line segment
            %         l2_high - end point of the second line segment
            % Outputs: s1 - parameterized point of intersection on the
            %               first line segment, given by (1 - s1) * l1_low
            %               + s1 * l1_high
            %          s2 - parameterized point of intersection on the
            %               second line segment, given by (1 - s2) * l2_low
            %               + s2 * l2_high
            delta_l1 = l1_high - l1_low;
            delta_l2 = l2_high - l2_low;
            diff = l1_low - l2_low;
            
            if abs(delta_l1(1) * delta_l2(2) - delta_l1(2) * delta_l2(1)) == 0
                fprintf("Lines are parallel, no intersection point\n");
                s1 = Inf;
                s2 = Inf;
            else
                s1 = (diff(2) * delta_l2(1) - diff(1) * delta_l2(2)) / (delta_l1(1) * delta_l2(2) - delta_l1(2) * delta_l2(1));
                if abs(delta_l2(1)) > 0
                    s2 = (diff(1) + s1 * delta_l1(1)) / delta_l2(1);
                elseif abs(delta_l2(2)) > 0
                    s2 = (diff(2) + s1 * delta_l1(2)) / delta_l2(2);
                else
                    fprintf("Degenerate case?\n");
                    s2 = Inf;
                end
                
%                 fprintf("s1 = %f, s2 = %f\n", s1, s2);
            end
        end
        
        function simulate(obj, q_init, dq_init, dt, T, filename)
            fig = figure(1);
            
            q = q_init;
            dq = dq_init;

            hold on;

            axis([-3 3 -3 3]);

            set(gca,'YTick',[]);
            set(gca,'XTick',[]);

            [pos_elbow, pos_ee] = obj.fwd_kinematics(q(1:2));
            pos_door = obj.fwd_kin_door(q(3));

            p0 = [0; 0];
            p0_door = obj.door_hinge;
            
            ceiling = rectangle('Position',[-3 obj.door_hinge(2) 6 3]', 'FaceColor',[0.7 0.7 0.7], ... 
                'EdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
            
            floor = rectangle('Position',[-3 -3 6 3]', 'FaceColor',[0.7 0.7 0.7], ... 
                'EdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
            
            link1 = plot([p0(1) pos_elbow(1)], [p0(2) pos_elbow(2)], 'k');
            link2 = plot([pos_elbow(1) pos_ee(1)], [pos_elbow(2) pos_ee(2)], 'k');
            link3 = plot([p0_door(1) pos_door(1)], [p0_door(2) pos_door(2)], 'k');
            joint1 = scatter([p0(1) p0_door(1) pos_elbow(1)], [p0(2) p0_door(2) pos_elbow(2)], 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
            joint2 = scatter([pos_door(1) pos_ee(1)], [pos_door(2) pos_ee(2)], 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
            
%             force_counter = 0;
            for t = 2:T + 1
                [M, C, N] = obj.dynamics(q, dq);
                % simulate constant force on the door
%                 if force_counter < 50
%                     u = [0; 0; -2];
%                 elseif force_counter > 100 && force_counter < 160
%                     u = [0;0;5];
%                 else
%                     u = [0;0;0];
%                 end
%                 force_counter = force_counter + 1;
                u = [0; 0; 0];
                
                fprintf("----- M: \n");
                disp(M);
                fprintf("----- q: \n");
                disp(q);
                fprintf("----- dq: \n");
                disp(dq);
                
                ddq = M \ (u - N - C * dq);
                dq = dq + dt * ddq;
                q = q + dt * dq;

                [pos_elbow, pos_ee] = obj.fwd_kinematics(q(1:2));
                pos_door = obj.fwd_kin_door(q(3));
                
                set(link1, 'Xdata', [p0(1) pos_elbow(1)]);
                set(link1, 'Ydata', [p0(2) pos_elbow(2)]);
                set(link2, 'Xdata', [pos_elbow(1) pos_ee(1)]);
                set(link2, 'Ydata',  [pos_elbow(2) pos_ee(2)]);
                set(link3, 'Xdata', [p0_door(1) pos_door(1)]);
                set(link3, 'Ydata', [p0_door(2) pos_door(2)]);
                set(joint1, 'Xdata', [p0(1) p0_door(1) pos_elbow(1)]);
                set(joint1, 'Ydata', [p0(2) p0_door(2) pos_elbow(2)]);
                set(joint2, 'Xdata', [pos_door(1) pos_ee(1)]);
                set(joint2, 'Ydata', [pos_door(2) pos_ee(2)]);

                % save fig
                if ~isempty(filename)
                    frame = getframe(fig);
                    im = frame2im(frame);
                    [imind, cm] = rgb2ind(im, 256); 
                    if t == 2
                        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf); 
                    else
                        imwrite(imind, cm, filename, 'gif', 'DelayTime',0.01, 'WriteMode', 'append');
                    end
                end
                
                pause(dt);
            end

            hold off;            
        end
        
        function plot_traj(obj, traj, filename)
            % Plots the arm and door following a trajectory
            %
            % Inputs:   traj    - sequence of configurations over time [q_0, q_1, ..., q_T]^T
            % Outputs:  n/a
            fig = figure(5);
            hold on;
            
            axis([-3 3 -3 3]);

            set(gca,'YTick',[]);
            set(gca,'XTick',[]);

            pos_base = [0; 0];
            
            ceiling = rectangle('Position',[-3 obj.door_hinge(2) 6 3]', 'FaceColor',[0.7 0.7 0.7], ... 
                'EdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
            
            floor = rectangle('Position',[-3 -3 6 3]', 'FaceColor',[0.7 0.7 0.7], ... 
                'EdgeColor',[0.5 0.5 0.5], 'LineWidth',1);

            for t = 0:size(traj, 1) - 1
                % get the current configuration and elbow/ee/door positions
                q = traj(t + 1, :)';

                a = (t+1)/(size(traj, 1)+1);
                
                % TODO plot the arm
                [pos_elbow, pos_ee] = obj.fwd_kinematics(q(1:2));
                pos_door = obj.fwd_kin_door(q(3));
                
                link1 = plot([pos_base(1) pos_elbow(1)], [pos_base(2) pos_elbow(2)], 'k');
                link2 = plot([pos_elbow(1) pos_ee(1)], [pos_elbow(2) pos_ee(2)], 'k');
                link3 = plot([obj.door_hinge(1) pos_door(1)], [obj.door_hinge(2) pos_door(2)], 'k');
                joint1 = scatter([pos_base(1) obj.door_hinge(1) pos_elbow(1)], [pos_base(2) obj.door_hinge(2) pos_elbow(2)], 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
                joint2 = scatter([pos_door(1) pos_ee(1)], [pos_door(2) pos_ee(2)], 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
                
                set(link1, 'LineWidth', 5);
                set(link2, 'LineWidth', 5);
                set(link3, 'LineWidth', 5);
                link1.Color(4) = a;
                link2.Color(4) = a;
                link3.Color(4) = a;
                
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
        
        function show_state_debug(obj, q, fig)
            % TODO function should display state of the robot and door,
            % with collision 
            figure(fig);
            
            hold on;

            axis([-3 3 -3 3]);

            set(gca,'YTick',[]);
            set(gca,'XTick',[]);
            
            pos_base = [0; 0];
            [pos_elbow, pos_ee] = obj.fwd_kinematics(q(1:2));
            pos_door_hinge = obj.door_hinge;
            pos_door_tip = obj.fwd_kin_door(q(3));
            
            % Compute intersection between link 1 and the door
            [link1_w_door_s1, link1_w_door_s2] = obj.intersection(pos_base, pos_elbow, pos_door_hinge, pos_door_tip);
            % Compute intersection between link 2 and the door
            [link2_w_door_s1, link2_w_door_s2] = obj.intersection(pos_elbow, pos_ee, pos_door_hinge, pos_door_tip);
            
            % TODO handle parallel case
            
            % Check if link 1 intersects the door
            p1 = (1 - link1_w_door_s1) * pos_base + link1_w_door_s1 * pos_elbow;
            if 0 <= link1_w_door_s1 && link1_w_door_s1 <= 1 && 0 <= link1_w_door_s2 && link1_w_door_s2 <= 1
                fprintf('Link 1 intersects with door!\n');
                scatter(p1(1), p1(2), 'rx');
            end
            
            % Check if link 2 intersects the door
            p2 = (1 - link2_w_door_s1) * pos_elbow + link2_w_door_s1 * pos_ee;
            if 0 <= link2_w_door_s1 && link2_w_door_s1 <= 1 && 0 <= link2_w_door_s2 && link2_w_door_s2 <= 1
                fprintf('Link 2 intersects with door!\n');
                scatter(p2(1), p2(2), 'rx');
            end
            
            % Compute signed distance functions
            fprintf('%%%% signed distances no contact %%%%\n');
            [dist_floor_elbow, dist_ceil_elbow, dist_floor_ee, dist_ceil_ee, dist_link_door] = obj.signed_dist_no_contact_separate(q);
            fprintf('\tdist floor to elbow:\t%f\n', dist_floor_elbow);
            fprintf('\tdist ceil to elbow:\t%f\n', dist_ceil_elbow);
            fprintf('\tdist floor to ee:\t%f\n', dist_floor_ee);
            fprintf('\tdist ceil to ee:\t%f\n', dist_ceil_ee);
            fprintf('\tdist link2 to door:\t%f\n', dist_link_door);
            fprintf('\tmin dist:\t\t%f\n', obj.signed_dist_no_contact(q));
            fprintf('%%%% signed distance contact %%%%\n');
            dist_door = obj.signed_dist_contact(q);
            fprintf('\tdist elbow to door:\t%f\n', dist_door);
            
            ceiling = rectangle('Position',[-3 pos_door_hinge(2) 6 3]', 'FaceColor',[0.7 0.7 0.7], ... 
                'EdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
            
            floor = rectangle('Position',[-3 -3 6 3]', 'FaceColor',[0.7 0.7 0.7], ... 
                'EdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
            
            link1 = plot([pos_base(1) pos_elbow(1)], [pos_base(2) pos_elbow(2)], 'k');
            link2 = plot([pos_elbow(1) pos_ee(1)], [pos_elbow(2) pos_ee(2)], 'k');
            link3 = plot([pos_door_hinge(1) pos_door_tip(1)], [pos_door_hinge(2) pos_door_tip(2)], 'k');
            joint1 = scatter([pos_base(1) pos_door_hinge(1) pos_elbow(1)], [pos_base(2) pos_door_hinge(2) pos_elbow(2)], 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
            joint2 = scatter([pos_door_tip(1) pos_ee(1)], [pos_door_tip(2) pos_ee(2)], 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
        end
        
    end
end

