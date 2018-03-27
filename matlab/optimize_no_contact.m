
% setup the number of timesteps
N = 30;

% setup the start and goal configurations
q_start = [0;0;0;0];
q_goal = [0;-pi/2;0;0];

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
%h_min = 0.1; %0.05;
%h_max = 0.1; %0.5;

% limit on how far the EE can be from the goal EE position
EE_epsilon = 1e-2;

% 4*(N+1)   for configuration and velocity
% 2*N       for control torques
% N       	for variable timesteps
% num_vars = 4*(N+1) + 2*N + N;
num_vars = 4*(N+1) + 2*N;

% setup the initial guess of the trajectory
x = zeros(num_vars,1);
% set the start (at t=0) and goal (at t=N)
%x(1:4) = q_start;
%x(4 * N + 1:4 * N + 4) = q_goal;
% linearly interpolate between the start and goal
for t = 0:N
    target = (q_goal - q_start)*(t/N) + q_start;
    x(4 * t + 1:4 * t + 2) = [target(1);target(2)];
end

% lower bound the arm positions
%xlow = zeros(num_vars,1);
% constrain the start to be at the start config
%xlow(1:4) = q_start; 
% constrain the end
%xlow(4 * (N-1) + 1:4 * (N-1) + 4) = q_goal; 
%xlow(4 * N + 1:4 * N + 4) = q_goal;
%xlow(4 * N + 1:4 * N + 2) = q_min;
%xlow(4 * N + 3:4 * N + 4) = [0;0]; %dq_min;

% upper bound the arm positions
%xupp = zeros(num_vars,1);
% constrain the start 
%xupp(1:4) = q_start; 
% constrain the end 
%xupp(4 * (N-1) + 1:4 * (N-1) + 4) = q_goal; 
%xupp(4 * N + 1:4 * N + 4) = q_goal; 
%xupp(4 * N + 1:4 * N + 2) = q_max; 
%xupp(4 * N + 3:4 * N + 4) = [0;0]; %dq_max;

% lower and upper bound the arm positions
xlow = zeros(num_vars,1);
xupp = zeros(num_vars,1);
% linearly interpolate between the start and goal
for t = 0:N
    if t == 0
        % initial configuration constraint
        xlow(4 * t + 1:4 * t + 4) = q_start;
        xupp(4 * t + 1:4 * t + 4) = q_start; 
    elseif t == N
        % constraint final config
        %xlow(4 * N + 1:4 * N + 4) = q_goal;
        %xupp(4 * N + 1:4 * N + 4) = q_goal; 
        
        % constrain final position
        xlow(4 * t + 1:4 * t + 4) = [q_min; 0; 0]; 
        xupp(4 * t + 1:4 * t + 4) = [q_max; 0; 0]; 
%     elseif t == N-1
%         % constrain the N-1th configuration to have zero velocity
%         xlow(4 * t + 1:4 * t + 4) = [q_min; 0; 0];
%         xupp(4 * t + 1:4 * t + 4) = [q_max; 0; 0];
    else
        % set the lower and upper bounds on configuration and velocity 
        % during the middle of the traj
        xlow(4 * t + 1:4 * t + 4) = [q_min;dq_min];
        xupp(4 * t + 1:4 * t + 4) = [q_max;dq_max];
    end
end

u_start_idx = 4 * (N + 1) + 1;
% set the control bounds
for t=1:N
    xlow(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1) = tau_min;
    xupp(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1) = tau_max;
end

% set the timestep bounds
% time_start_idx = 4*(N + 1) + 2*N + 1;
% for t=1:N
%     xlow(time_start_idx + (t-1)) = h_min;
%     xupp(time_start_idx + (t-1)) = h_max;
% end

xmul = [];
xstate = [];

% --- denote all the constraints --- %
% 1         for objective
% 2N        for q_t+1 - q_t + h * dq_t = 0
% 2N        for M_t+1* ddq_t+1 _ ... - u_t+1 - J^T(q_t+1)*lambda_t+1
% 1         for ||EE_pos - EE_goal|| <= EE_epsilon
% floor((N-3)/2)    for the timestep =0 constraints
% N         for EE constraint not going through the ground plane
%num_constraints = 1 + 2*N + 2*N + 1 + floor((N-3)/2) + N;
% ----------------------------------- %

% 1         for objective
% 2N        for q_t+1 - q_t + h * dq_t = 0
% 2N        for M_t+1* ddq_t+1 _ ... - u_t+1 - J^T(q_t+1)*lambda_t+1
% 1         for ||EE_pos - EE_goal|| <= EE_epsilon
num_constraints = 1 + 2*N + 2*N + 1;

Flow = zeros(num_constraints,1);
Fupp = zeros(num_constraints,1);    % all of the constraints are equality constraints (so keep them at zero)
Fupp(1) = 100;                      % bound the objective separately 
Fupp(1+2*N+2*N+1) = EE_epsilon;     % bound the final EE position
Fmul = [];
Fstate = [];

% set how high up the EE can be from ground plane
% ee_start = 1 + 2*N + 2*N + 1 + floor((N-3)/2) +1;
% for t=1:N
%     Fupp(ee_start + (t-1)) = 100;
% end

% solve the optimization problem
snprint('no_contact.out');
snspec('settings.spc');
snscreen on;
[x, F, INFO, xmul, Fmul, xstate, Fstate, output] = snopt(x, xlow, xupp, ...
                    xmul, xstate, Flow, Fupp, Fmul, Fstate, @robotFun, 0, 1);
snprint off;
snend;

states = x(1:4*(N+1))     
obj = F(1)

% extract just the trajectory (no controls) from the solution
traj = x(1:4*(N+1));
% extract the timesteps
% times = x(4*(N+1) + 2*N + 1: 4*(N+1) + 2*N + N);
times = ones(N)*0.1;
% extract the controls
controls = x(4*(N+1) + 1: 4*(N+1) + 2*N);

% plot the resulting trajectory
path = '/home/abajcsy/hybrid_ws/src/contact_manip/matlab/gifs/';
arm_filename = strcat(path,'nocontact.gif');
plot_arm(traj,times,N,arm_filename);
% plot the final controls
controls_filename = strcat(path,'nocontact_controls.png');
plot_controls(controls,N,controls_filename);
% plot the final velocity profile
velocity_filename = strcat(path,'nocontact_velocity.png');
plot_velocities(traj,N,velocity_filename);