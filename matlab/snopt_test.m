x = [5; 5]; % Initial value of x
xlow = [0; 0]; % Lower bound on x
xupp = [10; 10]; % Upper bound on x
xmul = [0; 0]; % ?
xstate = [0; 0]; % ?
% xmul = [];
% xstate = [];
Flow = [0; 0; 0]; % Lower bound on F (objective and constraints)
Fupp = [100; 10; 10]; % Upper bound on F (objective and constraints)
Fmul = [0; 0; 0]; % ?
Fstate = [0; 0; 0]; % ?
% Fmul = [];
% Fstate = [];

[x, F, INFO, xmul, Fmul, xstate, Fstate, output] = snopt(x, xlow, xupp, xmul, xstate, ...
                Flow, Fupp, Fmul, Fstate, @userfun, 0, 1);

% userfun(x)
x
F
           
function [F] = userfun(x)
c = [1; 2];
% c = [0; 0];
F = [(x(1) - c(1))^2 + (x(2) + c(2))^2; x];
end
