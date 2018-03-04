function [adj] = adjoint(g)
%ADJOINT Summary of this function goes here
%   Detailed explanation goes here
rot = g(1:3, 1:3);
p = g(1:3, 4);

adj = zeros(6, 6);
adj(1:3, 1:3) = rot;
adj(1:3, 4:6) = hat(p) * rot;
adj(4:6, 4:6) = rot;
end

