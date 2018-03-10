function [adj_inv] = adjoint_inv(g)
%ADJOINT_INV Summary of this function goes here
%   Detailed explanation goes here
rot = g(1:3, 1:3);
p = g(1:3, 4);

adj_inv = sym(zeros(6, 6));
adj_inv(1:3, 1:3) = rot';
adj_inv(1:3, 4:6) = -rot' * hat(p);
adj_inv(4:6, 4:6) = rot';
end

