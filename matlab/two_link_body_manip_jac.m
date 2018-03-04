function [J] = two_link_body_manip_jac(q)
%TWO_LINK_BODY_MANIP_JAC Computes the body manipulator Jacobian for a two
%link arm
%   Detailed explanation goes here
l1 = 1.0; % Length of link 1
l2 = 1.0; % Length of link 2

% TODO has the wrong sign (why?)
J = [0 0; 
     (-l1 * sin(q(1)) - l2 * sin(q(1) + q(2))) (-l2 * sin(q(1) + q(2)));
     (l1 * cos(q(1)) + q(2) * cos(q(1) + q(2))) (l1 * cos(q(1) + q(2)));
     1 1;
     0 0;
     0 0];
end

