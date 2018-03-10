function [dist] = two_link_signed_dist(q)
%TWO_LINK_SIGNED_DIST Summary of this function goes here
%   Detailed explanation goes here
[g_sl1, g_st] = two_link_kinematics(q(1), q(2));
% dist = min(g_sl1(3, 4), g_st(3, 4));
dist = g_sl1(3, 4);
end

