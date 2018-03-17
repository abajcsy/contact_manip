function [g_sl1, g_st] = two_link_kinematics(th1, th2)
%TWO_LINK_KINEMATICS Computes the forward kinematics for a two link arm
%   Detailed explanation goes here
l1 = 1.0; % Length of link 1
l2 = 1.0; % Length of link 2

g_sl1 = [1 0 0 0; 
         0 cos(th1) -sin(th1) -l1 * sin(th1);
         0 sin(th1) cos(th1) l1 * cos(th1);
         0 0 0 1];
    
g_st = [1 0 0 0;
        0 cos(th1 + th2) -sin(th1 + th2) -l1 * sin(th1) - l2 * sin(th1 + th2);
        0 sin(th1 + th2) cos(th1 + th2) l1 * cos(th1) + l2 * cos(th1 + th2);
        0 0 0 1];
end

