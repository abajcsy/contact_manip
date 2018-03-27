arm = TwoLinkArm(2, 1, 1, 1, 1, 1);
q_init = [pi / 2; 0; 0; 0];
q_final = [pi / 2; pi / 2; 0; 0];
T = 20;
k = 2;
eps = 0.01;
optProb = OptProb(arm, q_init, q_final, T, @g, @g_f, k, eps);

[x, xlow, xupp, F, Flow, Fupp] = optProb.generate();

xmul = [];
xstate = [];
Fmul = [];
Fstate = [];

snscreen on;
[x_opt, F_opt, INFO, xmul, Fmul, xstate, Fstate, output] = snopt(x, xlow, xupp, xmul, xstate, ...
                Flow, Fupp, Fmul, Fstate, F, 0, 1);

snprint off;
snend;

% x_opt = x;

traj = zeros(T + 1, arm.dof);
for t = 0:T
   q = optProb.get_q(x_opt, t);
   traj(t + 1, :) = q';
end

arm.plot_traj(traj);

function run_cost = g(q_t, dq_t, u_t1)
    run_cost = norm(dq_t);
end

function final_cost = g_f(q_T, dq_T)
    final_cost = 0;
end