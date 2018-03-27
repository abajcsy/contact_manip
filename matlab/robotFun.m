%% Setup Robot Optimization Problem (objective and constraints)
% State x is are the joint states (angles and velocities) proceeded by the
% controls (joint torques)
% x = [[th1; th2; dth1; dth2]; ...; [u1; u2]; ...]
function [F] = robotFun(x)

global l1 l2 m1 m2 r1 r2 Ix1 Ix2 alpha beta delta N;

% final time
N = 30;

% length and mass of links
l1 = 1.0;
l2 = 1.0;
m1 = 1.0;
m2 = 1.0;

% length to center of mass
r1 = l1/2.0;
r2 = l2/2.0;

% moment of inertia for rod of length l and mass m rotating about its
% center: https://en.wikipedia.org/wiki/List_of_moments_of_inertia#List_of_3D_inertia_tensors
Ix1 = (m1*l1^2)/12.0;
Ix2 = (m2*l2^2)/12.0;

% compute constants in the matrix
alpha = Ix1 + Ix2 + m1*r1^2 + m2*(l1^2 + r2^2);
beta = m2*l1*r2;
delta = Ix2 + m2*r2^2;

% fixed timestep
h = 0.1;

% get the end-effector (x,y) goal position [THIS IS A HACK] %
q_goal = [0;-pi/2;0;0];
[~, g_st_goal] = two_link_kinematics(q_goal(1),q_goal(2));
EE_goal = [g_st_goal(2, 4); g_st_goal(3, 4)];

% objective to optimize
obj = get_u(x,N)'*get_u(x,N); %norm(get_u(x,N))^2;

% constraints of the form:
%   q_t+1 - q_t - h*dq_t+1 = 0
deriv_constraints = zeros(2*N,1);

% constraints of the form:
%   M_t+1 * (ddq_t+1) + C_t+1(q_t+1) * (dq_t+1) + N_t+1(q_t+1, dq_t+1) - u_t+1 = 0 
% where
%   ddq_t+1 = (dq_t+1-dq_t)/h
dyn_constraints = zeros(2*N,1);

% constraint so that the EE doesn't go through the ground plane
% ee_constraints = zeros(N,1);

% for each timestep, setup the constraints
for t = 0:N-1
    % get the current delta t
    %h  = get_h(x,t+1);
    
    % get current and next state
    q_t = get_q(x,t);
    q_t1 = get_q(x,t+1);
    
    % get the next velocity, acceleration, and control
    dq_t = get_dq(x,t);
    dq_t1 = get_dq(x,t+1);
    ddq_t1 = (dq_t1-dq_t)/h;
    u_t1 = get_u(x,t+1);
    
    % gravity
    g = -9.81;
    
    % compute the inertial, coriolis, and normal force matrices
    %M_t1 = two_link_inertia_matrix(Ix1,Ix2,l1,l2,m1,m2,q_t1(1),q_t1(2)); 
    M_t1 = get_M(q_t1);
    %C_t1 = two_link_coriolis_matrix(dq_t1(1),dq_t1(2),l1,l2,m1,m2,q_t1(1),q_t1(2)); 
    C_t1 = get_C(q_t1, dq_t1);
    %N_t1 = two_link_gravity_matrix(g,l1,l2,m1,m2,q_t1(1),q_t1(2)); 
    N_t1 = get_N(q_t1);
    
    deriv = q_t1 - q_t - h*dq_t1;
    dyn = M_t1*ddq_t1 + C_t1*dq_t1 + N_t1 - u_t1;
    
    % setup the constraints for this timestep
    deriv_constraints(2*t+1:2*t+2) = deriv;
    dyn_constraints(2*t+1:2*t+2) = dyn;
    
    % setup the endeffector ground-plane collision constraint
   % ee_constraints(t+1) = get_ee_dist(q_t1);
end

% ensure that
%   h_2*j-1 + h_2*j = h_2*j+1 + h_2*j+2
% where j = 1, 2, ... floor(N-3/2)
% timestep_constraints = zeros(floor((N-3)/2),1);

% add in the timestep constraints
% for j=1:floor((N-3)/2)
%     h = get_h(x,2*j-1);
%     h1 = get_h(x,2*j);
%     h2 = get_h(x,2*j+1);
%     h3 = get_h(x,2*j+2);
%     timestep_constraints(j) = h+h1-h2-h3;
% end

dist_to_goal_pos = 0;
for t = N-2:N
    % setup the cost of the EE position from goal position
    % Note: we need to constrain the last two positions to make it stabilize
    q = get_q(x,t);
    [~, g_st] = two_link_kinematics(q(1),q(2));
    EE_final = [g_st(2, 4); g_st(3, 4)];
    dist_to_goal_pos = dist_to_goal_pos + norm(EE_final - EE_goal);
end

% setup constraints
F = [obj; deriv_constraints; dyn_constraints; dist_to_goal_pos]; %timestep_constraints]; %ee_constraints];
end

%% Getting the inertial matrix 
%  Computes the inertial matrix at time t for 2-link arm
function [M] = get_M(q_t1)

global alpha beta delta Ix1 Ix2 l1 l2 m1 m2;

% extract the angles from the configuration
theta2 = q_t1(2);
M = [alpha + 2*beta*cos(theta2) delta + beta*cos(theta2);
      delta + beta*cos(theta2)   delta];

% th2 = q_t1(2);
% M = [Ix1 + Ix2 + 0.25 * l1^2 * m1 + l1^2 * m2 * sin(th2)^2 + m2 * (0.5 * l1 * cos(th2) + 0.25 * l2)^2 Ix2 + 0.5 * l2 * m2 * (0.5 * l1 * cos(th2) + 0.25 * l2); 
%     Ix2 + 0.5 * l2 * m2 * (0.5 * l1 * cos(th2) + 0.25 * l2) Ix2 + 0.25 * l2^2 * m2];
end

%% Getting the coriolis matrix
%  Computes the coriolis and centrifugal forces acting on the joints
function [C] = get_C(q_t1, dq_t1)

global beta l1 l2 m1 m2;

% extract the angles and velocities from the configuration
dtheta1 = dq_t1(1);
theta2 = q_t1(2);
dtheta2 = dq_t1(2);

C = [-beta*sin(theta2)*dtheta2 -beta*sin(theta2)*(dtheta1+dtheta2);
     beta*sin(theta2)*dtheta1 0];

% dth1 = dq_t1(1);
% th2 = q_t1(2);
% dth2 = dq_t1(2);
% C = [dth2 * l1 * m2 * (0.75 * l1 * cos(th2) - 0.125 * l2) * sin(th2) l1 * m2 * (0.75 * dth1 * l1 * cos(th2) - 0.125 * dth1 * l2 - 0.25 * dth2 * l2) * sin(th2);
%      l1 * m2 * (-0.75 * dth1 * l1 * cos(th2) + 0.125 * dth1 * l2 - 0.125 * dth2 * l2) * sin(th2) -0.375 * dth1 * l1^2 * m2 * sin(2 * th2)];
end

%% Getting the gravitational matrix
%  Computes the effect of gravity on the links
function [N] = get_N(q_t1)

global m1 m2 l1 l2;

theta1 = q_t1(1);
theta2 = q_t1(2);
g = -9.81;
N = [(m1+m2)*g*l1*cos(theta1)+m2*g*l2*cos(theta1 + theta2);
     m2*g*l2*cos(theta1 + theta2)];

% th1 = q_t1(1);
% th2 = q_t1(2);
% g = 9.81;
% N = [-g * (0.5 * l1 * m1 * sin(th1) + l1 * m2 * sin(th1) + 0.5 * l2 * m2 * sin(th1 + th2));
%      -0.5 * g * l2 * m2 * sin(th1 + th2)];

end

%% End-effector distance constraints
% Computes the distance of the end-effector to the ground.
% function [ee_dist] = get_ee_dist(q_t1)
%     [~, g_st] = two_link_kinematics(q_t1(1), q_t1(2));
%     ee_dist = g_st(3, 4);
% end

%% Getting the timestep duration
% function [h] = get_h(x, t)
%     global N
%     time_start_idx = 4*(N + 1) + 2*N + 1;
%     h = x(time_start_idx + (t-1));
% end

%% Getting control at a specific time
%  Extracts the control vector at time t from the optimization vector x
%  using the equation
function [u] = get_u(x, t)
    global N
    u_start_idx = 4 * (N + 1) + 1;
    u = x(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1);
end

%% Getting a configuration at a specific time
function [q] = get_q(x, t)
    q = x(4 * t + 1:4 * t + 2);
end

%% Getting a velocity at specific time
function [dq] = get_dq(x, t)
    dq = x(4 * t + 3:4 * t + 4);
end