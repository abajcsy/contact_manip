
% setup the robot arm parameters
dof = 2;
c = 1;
m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1;

% create a robot arm
arm = TwoLinkArm(dof, c, m1, m2, l1, l2);

% initial and final configurations
q_init = [0; 0; 0; 0];
q_final = [0; pi / 2; 0; 0];

% number of timesteps
T = 5;
% number of timesteps to consider when computing EE dist from final pose
k = 2;
eps = 0.01;

% create the optimization problem
optProb = OptProb(arm, q_init, q_final, T, @g, @g_f, k, eps);

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

[x, xlow, xupp, F, Flow, Fupp] = optProb.generate();
% F(x)
% Flow
% Fupp
% x
traj = zeros(T + 1, arm.dof);
for t = 0:T
   q = optProb.get_q(x, t);
   traj(t + 1, :) = q';
end

arm.plot_traj(traj);

function run_cost = g(q_t, dq_t, u_t1)
    run_cost = 1;
end

function final_cost = g_f(q_T, dq_T)
    final_cost = 0.5;
end