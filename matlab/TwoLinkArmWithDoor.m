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
            obj.u_min = [-5; -5];
            obj.u_max = [5; 5];
            
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
            door_contact_vec = [(obj.l3/2)*cos(th3); (obj.l3/2)*sin(th3)];
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
            
            F_door = [lambda_t1*sin(th3); -lambda_t1*cos(th3)];
            door_contact_vec = obj.get_contact_vec(th3);
            u_d = norm(cross(F_door, door_contact_vec)); % joint torques experienced by door
            u_c = [u_r; u_d];
        end
        
        function phi = signed_dist(obj, q_t1)
            % TODO DESCRIBE HOW THIS WORKS
            pos_base = [0; 0];
            [pos_elbow, pos_ee] = obj.fwd_kinematics(q_t1(1:2));
            pos_door_hinge = obj.door_hinge;
            pos_door_tip = obj.fwd_kin_door(q_t1);
            
            % intersection between door and link 1
            s1 = obj.intersect(pos_base, pos_elbow, pos_door_hinge, pos_door_tip);
             % intersection between door and link 2
            s2 = obj.intersect(pos_elbow, pos_ee, pos_door_hinge, pos_door_tip);
            
            % link 1 intersects
            if s1 < 1 && s1 > 0
                phi = -1;
            elseif s1 == 1 % at elbow
                phi = 0;
            elseif s1 == 0 % at base
                phi = 1;
            else % not intersecting
                phi = 1;
            end
            
            if phi >= 0
                % link 2 intersects
                if s2 < 1 && s2 > 0
                    phi = -1;
                elseif s2 == 0 % at elbow
                    phi = 0;
                elseif s2 == 1 % at ee
                    phi = 1;
                else % not intersecting
                    phi = 1;
                end
            end
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
        
        
    end
end

