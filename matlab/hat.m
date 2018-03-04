function [a_hat] = hat(a)
%HAT Summary of this function goes here
%   Detailed explanation goes here
a_hat = [0 -a(3) a(2);
         a(3) 0 -a(1);
         -a(2) a(1) 0];
end

