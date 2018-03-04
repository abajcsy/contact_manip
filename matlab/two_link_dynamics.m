function [M, C, N] = two_link_dynamics(th1, th2, dth1, dth2)
%TWO_LINK_DYNAMICS Computes the manipulator dynamics matrices of a 2 link arm
%   Detailed explanation goes here
m1 = 1.0; % Mass of link 1
m2 = 1.0; % Mass of link 2
l1 = 1.0; % Length of link 1
l2 = 1.0; % Length of link 2
g = 9.81;

% Approximate the links as thin rods
I1x = (m1 * l1^2) / 12.0;
I2x = (m2 * l2^2) / 12.0;

% Inertia matrix
M = [I1x + I2x + 0.25 * l1^2 * m1 + l1^2 * m2 * sin(th2)^2 + m2 * (0.5 * l1 * cos(th2) + 0.25 * l2)^2 I2x + 0.5 * l2 * m2 * (0.5 * l1 * cos(th2) + 0.25 * l2); 
     I2x + 0.5 * l2 * m2 * (0.5 * l1 * cos(th2) + 0.25 * l2) I2x + 0.25 * l2^2 * m2];

% Coriolis matrix 
C = [dth2 * l1 * m2 * (0.75 * l1 * cos(th2) - 0.125 * l2) * sin(th2) l1 * m2 * (0.75 * dth1 * l1 * cos(th2) - 0.125 * dth1 * l2 - 0.25 * dth2 * l2) * sin(th2);
     l1 * m2 * (-0.75 * dth1 * l1 * cos(th2) + 0.125 * dth1 * l2 - 0.125 * dth2 * l2) * sin(th2) -0.375 * dth1 * l1^2 * m2 * sin(2 * th2)];
 
% Potential forces matrix 
N = [-g * (0.5 * l1 * m1 * sin(th1) + l1 * m2 * sin(th1) + 0.5 * l2 * m2 * sin(th1 + th2));
     -0.5 * g * l2 * m2 * sin(th1 + th2)];
end
