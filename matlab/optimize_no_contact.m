
snspec('settings.spc');
snscreen on;

% setup the number of timesteps
N = 10;

% setup the start and goal configurations
q_start = [0;0;0;0];
q_goal = [0;-pi/2;0;0];
% get the end-effector (x,y) goal position 
[g_sl1, g_st] = two_link_kinematics(q_goal(1),q_goal(2));
EE_goal = [g_st(2, 4) g_st(3, 4)];

% torque limits
tau_min = [-5,-5];
tau_max = [5,5];

% joint limits
q_min = [-pi/2;-pi/2];
q_max = [pi/2;pi/2];
% velocity limits
dq_min = [-20;-20];
dq_max = [20;20];

% limit on how far the EE can be from the goal EE position
EE_epsilon = 0.1;

% setup the initial guess of the trajectory
x = zeros(4*(N+1)+2*N,1);
% set the start (at t=0) and goal (at t=N)
x(1:4) = q_start;
x(4 * N + 1:4 * N + 4) = q_goal;

% lower bound the arm positions
xlow = zeros(4*(N+1)+2*N,1);
% constrain the start to be at the start config
xlow(1:4) = q_start; 
% constrain the end
% xlow(4 * N + 1:4 * N + 4) = q_goal;
xlow(4 * N + 1:4 * N + 2) = q_min;
xlow(4 * N + 3:4 * N + 4) = dq_min;

% upper bound the arm positions
xupp = ones(4*(N+1)+2*N,1);
% constrain the start 
xupp(1:4) = q_start; 
% constrain the end 
% xupp(4 * N + 1:4 * N + 4) = q_goal; 
xupp(4 * N + 1:4 * N + 2) = q_max; 
xupp(4 * N + 3:4 * N + 4) = dq_max;

% linearly interpolate between the start and goal
for t = 1:N-1
    target = (q_goal - q_start)*(t/N) + q_start;
    x(4 * t + 1:4 * t + 2) = [target(1);target(2)];
    % set the lower and upper bounds on configuration and velocity 
    % during the middle of the traj
    xlow(4 * t + 1:4 * t + 2) = q_min;
    xlow(4 * t + 3:4 * t + 4) = dq_min;
    xupp(4 * t + 1:4 * t + 2) = q_max;
    xupp(4 * t + 3:4 * t + 4) = dq_max;
end

u_start_idx = 4 * (N + 1) + 1;
% set the control bounds
for t=1:N
    xlow(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1) = tau_min;
    xupp(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1) = tau_max;
end
    
xmul = [];
xstate = [];
Flow = zeros(1+4*N + 1,1);
Fupp = zeros(1+4*N + 1,1); % all of the constraints are equality constraints (so keep them at zero)
Fupp(1) = 100; % bound the objective separately 
Fupp(1+4*N + 1) = EE_epsilon % bound the final EE position
Fmul = [];
Fstate = [];
% solve the optimization problem
[x, F, INFO, xmul, Fmul, xstate, Fstate, output] = snopt(x, xlow, xupp, ...
                    xmul, xstate, Flow, Fupp, Fmul, Fstate, @robotFun, 0, 1);
snprint off;
snend;

% extract just the trajectory (no controls) from the solution
traj = x(1:4*N+4);
 
% plot the resulting trajectory
path = '/home/abajcsy/hybrid_ws/src/contact_manip/matlab/gifs/';
filename = strcat(path,'nocontact.gif');
plot_arm(traj,N,filename);

