function [g] = twist_angle_to_trans(twist, angle)
%TWIST_ANGLE_TO_TRANS Summary of this function goes here
%   Detailed explanation goes here
R = axis_angle_to_rot(twist(4:6, 1), angle);
g = sym(zeros(4, 4));
g(1:3, 1:3) = R;
g(1:3, 4) = (eye(3) - R) * cross(twist(4:6, 1), twist(1:3, 1)) + twist(4:6, 1) * twist(4:6, 1)' * twist(1:3, 1) * angle;
g(4, 4) = 1;
end

