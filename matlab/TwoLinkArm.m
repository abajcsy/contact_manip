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
            
            obj.q_min = [-pi / 2; -pi / 2];
            obj.q_max = [pi / 2; pi / 2];
            obj.dq_min = [-20; -20];
            obj.dq_max = [20; 20];
            obj.u_min = [-5; -5];
            obj.u_max = [5; 5];
        end
        
        function [M,C,N] = dynamics(obj, q_t1,dq_t1)
            % Returns equation of motion matrices for
            %   M(q_t1)*ddq_t1 + C(q_t1, dq_t1)*dq_t1 + N(q_t1) = u_t1
            M = zeros(obj.dof);
            C = zeros(obj.dof);
            N = zeros(obj.dof,1);
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
        
        
    end
end

