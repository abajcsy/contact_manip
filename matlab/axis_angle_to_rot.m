function [R] = axis_angle_to_rot(axis, angle)
%AXIS_ANGLE_TO_ROT Summary of this function goes here
%   Detailed explanation goes here
R = eye(3) + sin(angle) * hat(axis) + (1 - cos(angle)) * hat(axis)^2;
end

