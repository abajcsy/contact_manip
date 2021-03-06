%% Trajectory Optimization with Contact (2 Link Arm)

% setup the robot arm parameters
dof = 2;
c = 1;
m1 = 0.1;
m2 = 0.05;
l1 = 1.0;
l2 = 1.0;

% create a robot arm
arm = TwoLinkArm(dof, c, m1, m2, l1, l2);

global q_final;

% initial and final configurations
q_init = [pi / 2; 0; 0; 0];
% q_final = [pi / 2; pi / 2; 0; 0]; % First way to initialize
q_final = [pi; -pi / 2; 0; 0];  % Second way to initialize
% q_final = [0; 0; 0; 0];         % Should initialize to all zeros 

% number of timesteps
T = 20;

% create the optimization problem
optProb = OptProb(arm, q_init, q_final, T, @g, @g_f);
[x, xlow, xupp, F, Flow, Fupp] = optProb.generate();

            
% DEBUGGING (Ellis) -- constrain the final joint velocities to be zero
xlow = optProb.set_dq(xlow, optProb.T - 1, 0);
xlow = optProb.set_dq(xlow, optProb.T, 0);
xupp = optProb.set_dq(xupp, optProb.T - 1, 0);
xupp = optProb.set_dq(xupp, optProb.T, 0);

xmul = [];
xstate = [];
Fmul = [];
Fstate = [];

snscreen on;
[x_opt, F_opt, INFO, xmul, Fmul, xstate, Fstate, output] = snopt(x, xlow, xupp, xmul, xstate, ...
                Flow, Fupp, Fmul, Fstate, F, 0, 1);

snprint off;
snend;

fprintf('objective value: %d\n', F_opt(1));

% get components from optimization variable
traj = optProb.get_traj(x_opt);
traj_dq = optProb.get_traj_dq(x_opt);
traj_u = optProb.get_traj_u(x_opt);
traj_h = optProb.get_traj_h(x_opt);

% plot the arm's trajectory
arm.plot_traj(traj, 'out.gif');

% plot velocities and controls 
arm.plot_dq_u(traj_dq, traj_u, T);

for t = 1:optProb.T
   h = traj_h(t);
   lambda = optProb.get_lambda(x_opt, t);
   fprintf('h(%d) = %f, lambda(%d) = %f\n', t, h, t, lambda);
end

[~, target_ee] = arm.fwd_kinematics(q_final);
fprintf('ee target position = (%f, %f)\n', target_ee(1), target_ee(2));


%% Running cost function g(q,dq,u)
function run_cost = g(q, dq, u, t, T)
    run_cost = 0;

    if t >= T - 1
        run_cost = 100 * u' * u;
    end
end

%% Final cost function g_f(q,dq)
function final_cost = g_f(q, dq, arm)
    global q_final;
    % get the desired EE position from the final configuration
    [~, target_ee] = arm.fwd_kinematics(q_final);
    [~, curr_ee] = arm.fwd_kinematics(q);
    final_cost = 1000 * norm(target_ee - curr_ee);
%     final_cost = 0;
end
