x = [5; 5];
xlow = [0; 0];
xupp = [10; 10];
xmul = [0; 0];
xstate = [0; 0];
Flow = [0; 0; 0];
Fupp = [100; 10; 10];
Fmul = [0; 0; 0];
Fstate = [0; 0; 0];

[x, F, INFO, xmul, Fmul, xstate, Fstate, output] = snopt(x, xlow, xupp, ...
                    xmul, xstate, Flow, Fupp, Fmul, Fstate, @userfun, 0, 1);
x
F

function [F] = userfun(x)
c = [0; 2];
F = [norm(x+c)^2; x];
end