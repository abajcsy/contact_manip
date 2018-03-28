%% Trajectory Optimization with Contact

% setup the robot arm parameters
dof = 2;
c = 1;
m1 = 1.0;
m2 = 1.0;
l1 = 1.0;
l2 = 1.0;

% create a robot arm
arm = TwoLinkArm(dof, c, m1, m2, l1, l2);

% initial and final configurations
q_init = [pi / 2; 0; 0; 0];
q_final = [pi / 2; pi / 2; 0; 0];

% number of timesteps
T = 20;

% create the optimization problem
optProb = OptProb(arm, q_init, q_final, T, @g, @g_f);
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

F_opt(1)

% get just the joint angles from the opt variable
traj = optProb.get_traj(x_opt);

% plot the arm's trajectory
arm.plot_traj(traj);

%% Running cost function g(q,dq,u)
function run_cost = g(q_t, dq_t, u_t1)
    run_cost = 0;
end

%% Final cost function g_f(q,dq)
function final_cost = g_f(q_T, dq_T, arm)
    q_final = [pi / 2; pi / 2; 0; 0];
    % get the desired EE position from the final configuration
    [~,target_ee] = arm.fwd_kinematics(q_final);
    [~, curr_ee] = arm.fwd_kinematics(q_T);
    final_cost = norm(target_ee - curr_ee);
end
