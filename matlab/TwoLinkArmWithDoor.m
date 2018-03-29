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
        door_center
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
            
            obj.q_min = [0; -pi / 2];
            obj.q_max = [pi; pi / 2];
            obj.dq_min = [-20; -20];
            obj.dq_max = [20; 20];
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
            
            obj.door_center = [-l1/2.0; l1+l2+l3/2.0];
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
            pos_door = obj.door_center + [obj.l2*cos(th3); obj.l1*sin(th3)];
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
        
        function u_c = u_contact(obj, q_t1, lambda_t1)
            % TODO IMPLEMENT ME
        end
        
        function phi = signed_dist(obj, q_t1)
            % TODO IMPLEMENT ME
            % TODO This is wrong -- needs to account for door intersection
            %   with the two links in addition to distance from elbow to
            %   contact point
            [pos_elbow, ~] = obj.fwd_kinematics(q_t1(1:2));
            pos_door = obj.fwd_kin_dor(q_t1(3));
            pos_door_center = [pos_door(1)/2. pos_door(2)/2.];
            phi = norm(pos_elbow - pos_door_center);
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
            p0_door = obj.door_center;
            
            ceiling = rectangle('Position',[-3 obj.door_center(2) 6 3]', 'FaceColor',[0.7 0.7 0.7], ... 
                'EdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
            
            link1 = plot([p0(1) pos_elbow(1)], [p0(2) pos_elbow(2)], 'k');
            link2 = plot([pos_elbow(1) pos_ee(1)], [pos_elbow(2) pos_ee(2)], 'k');
            link3 = plot([p0_door(1) pos_door(1)], [p0_door(2) pos_door(2)], 'k');
            joint1 = scatter([p0(1) p0_door(1) pos_elbow(1)], [p0(2) p0_door(2) pos_elbow(2)], 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
            joint2 = scatter([pos_door(1) pos_ee(1)], [pos_door(2) pos_ee(2)], 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
            
            force_counter = 0;
            for t = 2:T + 1
                [M, C, N] = obj.dynamics(q, dq);
                % simulate constant force on the door
                if force_counter < 50
                    u = [0; 0; -2];
                elseif force_counter > 100 && force_counter < 160
                    u = [0;0;5];
                else
                    u = [0;0;0];
                end
                force_counter = force_counter + 1;
                ddq = M \ (u -N - C * dq);
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

