arm = TwoLinkArm(2, 1, 1, 1, 1, 1);
optProb = OptProb(arm, [0; 0; 0; 0], [0; pi / 2; 0; 0], 50);

x = [1:1:optProb.num_states + optProb.num_controls + optProb.num_lambdas + optProb.T]';

s = [];
for t = 0:optProb.T
    q = optProb.get_q(x, t);
    dq = optProb.get_dq(x, t);
    s = [s; q; dq];
end

us = [];
ls = [];
hs = [];
for t = 1:optProb.T
    u = optProb.get_u(x, t);
    l = optProb.get_lambda(x, t);
    h = optProb.get_h(x, t);
    us = [us; u];
    ls = [ls; l];
    hs = [hs; h];
end

xp = [s; us; ls; hs];
assert(isequal(xp, x))

xp = -1*ones(optProb.num_states + optProb.num_controls + optProb.num_lambdas + optProb.T, 1);

for t = 0:optProb.T
    x = optProb.set_q(x, t, [-1; -1]);
    x = optProb.set_dq(x, t, [-1; -1]);
end

for t = 1:optProb.T
    x = optProb.set_u(x, t, [-1; -1]);
    x = optProb.set_lambda(x, t, -1);
    x = optProb.set_h(x, t, -1);
end

assert(isequal(xp, x));

[x, xlow, xupp, F, ~, ~] = optProb.generate();
F(x)