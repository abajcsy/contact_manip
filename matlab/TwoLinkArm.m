classdef TwoLinkArm
    %TWOLINKARM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dof % Degrees of freedom
        c % Number of contact points
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
    end
end

