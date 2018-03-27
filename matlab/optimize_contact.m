
snspec('settings.spc');
snscreen on;

% setup the number of timesteps
N = 20;

% setup the start and goal configurations
q_start = [0;0;0;0];%[-pi/2+pi/16;0;0;0];
q_goal = [0;-pi/2;0;0];%[-pi/2;0;0;0];

% torque limits
tau_min = [-5,-5];
tau_max = [5,5];

% joint limits
q_min = [-pi/2;-pi/2];
q_max = [pi/2;pi/2];
% velocity limits
dq_min = [-20;-20];
dq_max = [20;20];

% limits on what the timesteps could be
h_min = 0.01;
h_max = 0.2;

% limit on how far the EE can be from the goal EE position
EE_epsilon = 1e-2;

% phi limits
phi_min = [0; 0];
phi_max = [1000; 1000];

% lambda limits
lambda_min = [0; 0];
lambda_max = [1000; 1000];

% 4*(N+1)   for configuration and velocity
% 2*N       for control torques
% 2*N       for external lambda forces
% N         for variable timesteps
num_vars = 4*(N+1) + 2*N + 2*N + N;

% setup the initial guess of the trajectory
x = zeros(num_vars,1);
% set the start (at t=0) and goal (at t=N)
x(1:4) = q_start;
x(4 * N + 1:4 * N + 4) = q_goal;
% initialize the timesteps
time_start_idx = 4*(N + 1) + 2*N + 2*N + 1;
for t=1:N
    x(time_start_idx + (t-1)) = 0.1;
end

% lower bound the arm positions, torques, and lambdas
xlow = zeros(num_vars,1);
% constrain the start to be at the start config
xlow(1:4) = q_start; 
% constrain the end
% xlow(4 * N + 1:4 * N + 4) = q_goal;
xlow(4 * N + 1:4 * N + 2) = q_min;
xlow(4 * N + 3:4 * N + 4) = [0;0]; %dq_min;

% upper bound the arm positions
xupp = ones(num_vars,1);
% constrain the start 
xupp(1:4) = q_start; 
% constrain the end 
% xupp(4 * N + 1:4 * N + 4) = q_goal; 
xupp(4 * N + 1:4 * N + 2) = q_max; 
xupp(4 * N + 3:4 * N + 4) = [0;0]; %dq_max;

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
    
    % constrain the N-1th configuration to have zero velocity
    if t == N-1 
        xlow(4*t + 3:4*t + 4) = [0;0];
        xupp(4*t + 3:4*t + 4) = [0;0];
    end
    if t == N-2
        xlow(4*t + 3:4*t + 4) = [0;0];
        xupp(4*t + 3:4*t + 4) = [0;0];
    end
end

u_start_idx = 4 * (N + 1) + 1;
% set the control bounds
for t=1:N
    xlow(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1) = tau_min;
    xupp(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1) = tau_max;
end
    
lambda_start_idx = 4 * (N + 1) + 2 * N + 1;
% set the lambda bounds
for t=1:N
    xlow(lambda_start_idx + 2 * (t - 1):lambda_start_idx + 2 * (t - 1) + 1) = lambda_min;
    xupp(lambda_start_idx + 2 * (t - 1):lambda_start_idx + 2 * (t - 1) + 1) = lambda_max;
end

% set the timestep bounds
time_start_idx = 4*(N + 1) + 2*N + 2*N + 1;
for t=1:N
    xlow(time_start_idx + (t-1)) = h_min;
    xupp(time_start_idx + (t-1)) = h_max;
end

xmul = [];
xstate = [];

% --- denote all the constraints --- %
% 1     for objective
% 2N    for q_t+1 - q_t + h * dq_t = 0
% 2N    for M_t+1* ddq_t+1 _ ... - u_t+1 - J^T(q_t+1)*lambda_t+1
% 1     for ||EE_pos - EE_goal|| <= EE_epsilon
% 2N    for phi signed distance constraint
% N     for phi*lambda equality constraint
% floor((N-3)/2)    for the variable timestep =0 constraint
% N     for EE constraint not going through the ground plane
num_constraints = 1 + 2*N + 2*N + 1 + 2*N + N + floor((N-3)/2) + N;
Flow = zeros(num_constraints,1);
Fupp = zeros(num_constraints,1);        % all of the constraints are equality constraints (so keep them at zero)
Fupp(1) = 100;                          % bound the objective separately 
Fupp(1 + 2*N + 2*N + 1) = EE_epsilon;   % bound the final EE position

ee_start = 1 + 2*N + 2*N + 1 + 2*N + N + floor((N-3)/2) + 1;
% bound the phi collision constraints
phi_start_idx = 1 + 2*N + 2*N + 1 + 1;
phi_lambda_start_idx = 1 + 2*N + 2*N + 1 + 2*N + 1;

for t=1:N
    Fupp(phi_start_idx + 2 * (t - 1):phi_start_idx + 2 * (t - 1) + 1) = phi_max;
    %Fupp(phi_lambda_start_idx + (t-1)) = 1e-2;
    Fupp(ee_start + (t-1)) = 100;
end
Fmul = [];
Fstate = [];

robotFun_contact(x)

% solve the optimization problem
[x, F, INFO, xmul, Fmul, xstate, Fstate, output] = snopt(x, xlow, xupp, ...
                    xmul, xstate, Flow, Fupp, Fmul, Fstate, @robotFun_contact, 0, 1);

snprint off;
snend;



% extract just the trajectory (no controls) from the solution
traj = x(1:4*N+4);
% extract the timesteps
times = x(4*(N+1) + 2*N + 2*N + 1: 4*(N+1) + 2*N + 2*N + N);
% extract the controls
controls = x(4*(N+1) + 1: 4*(N+1) + 2*N);

% print lambdas
lambdas = x(4*(N+1)+2*N+1:4*(N+1)+2*N+2*N)
% print objective value
objective = F(1)
% print state and velocities
state = x(1:4*N+4)

% plot the resulting trajectory
path = '/home/abajcsy/hybrid_ws/src/contact_manip/matlab/gifs/';
arm_filename = strcat(path,'contact.gif');
plot_arm(traj,times,N,arm_filename);
controls_filename = strcat(path,'contact_controls.png');
plot_controls(controls,N,controls_filename);
% plot the final velocity profile
velocity_filename = strcat(path,'contact_velocity.png');
plot_velocities(traj,N,velocity_filename);