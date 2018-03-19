% State x is are the joint states (angles and velocities) proceeded by the
% controls (joint torques), proceeded by the contact forces lambda
% x = [[th1; th2; dth1; dth2]; ...; [u1; u2]; ...; lmda1; ...; h; ...]

% Constraints are dth(t+1) - th(t) - h*th(t+1) = 0
%                 M(t+1)*(dth(t+1) - dth(t)) + h*(C(t+1)*dth(t+1) + N(t+1)
%                 - u(t+1) - J(t+1)'*lmda(t+1)) = 0
% Objective is u(T)'*u(T)

% T = 10;
% T = 15;
T = 20;
h = 0.1;

% q_init = [0; 0];
q_init = [-pi / 2 + pi / 16; 0];
% q_init = [-pi / 8; 0];
% q_init = [-pi / 2 + pi / 90; pi / 10];
% q_init = [-pi / 2; pi / 10];
% q_init = [-pi / 4; 0];
% q_init = [0; pi / 2]; 
% q_init = [0; pi / 4];
% q_final = [0.5; 0.25];
% q_final = [-0.75; 0.5];
% q_final = [0; 0];
dq_init = [0; 0];
q_final = [-pi / 2; 0];
dq_final = [0; 0];

% Tolerance on the final position of the end effector
% final_pos_thresh = 0.12;
final_pos_thresh = 1e-2;

% Limits on the joint angles
% q_limit_upp = [pi / 2; pi / 2];
% q_limit_low = [-pi / 2; -pi / 2];
q_limit_upp = [pi; pi];
q_limit_low = [-pi; -pi];

% Limits on the joint velocities
dq_limit_upp = [20; 20];
dq_limit_low = [-20; -20];

% Limits on the joint torques
u_limit_upp = [5; 5];
u_limit_low = [-5; -5];
% u_limit_upp = [10; 10];
% u_limit_low = [-10; -10];
% u_limit_upp = [20; 20];
% u_limit_low = [-20; -20];

% Limits on the contact force magnitudes
lmda_limit_upp = 1000;
lmda_limit_low = 0;

% Limits on the signed distances
phi_limit_upp = 1000;
% phi_limit_low = -1000;
phi_limit_low = 0;

% Limits on size of time step 
% h_limit_upp = 5.0;
h_limit_upp = 0.2;
h_limit_low = 0.01;
n_h_constraints = floor((T - 3) / 2);

% Construct an initial value of x
% x = zeros(4 * (T + 1) + 2 * T + T, 1);
x = zeros(4 * (T + 1) + 4 * T, 1);
% qs = zeros(2, T);
% for t = 0:T
%     k = t / T;
%     x(4 * t + 1:4 * t + 2, 1) = (1 - k) * q_init + k * q_final;
%     qs(:, t + 1) = get_q(x, t);
% end

% Initialize the timesteps
for t = 1:T
    h_start_idx = 4 * (T + 1) + 3 * T + 1;
    x(h_start_idx + t - 1, 1) = 0.1;
end

% two_link_draw(qs, h, '');

% Set upper and lower bounds on x
% xupp = zeros(4 * (T + 1) + 2 * T + T, 1);
% xlow = zeros(4 * (T + 1) + 2 * T + T, 1);
xupp = zeros(4 * (T + 1) + 4 * T, 1);
xlow = zeros(4 * (T + 1) + 4 * T, 1);
for t = 0:T
    if t == 0
        % Set initial config constraint
        xupp(4 * t + 1:4 * t + 2, 1) = q_init;
        xlow(4 * t + 1:4 * t + 2, 1) = q_init;
        xupp(4 * t + 3:4 * t + 4, 1) = dq_init;
        xlow(4 * t + 3:4 * t + 4, 1) = dq_init;
    elseif t == T
        % Set final config constraint
        xupp(4 * t + 1:4 * t + 2, 1) = q_final;
        xlow(4 * t + 1:4 * t + 2, 1) = q_final;
%         xupp(4 * t + 1:4 * t + 2, 1) = q_limit_upp;
%         xlow(4 * t + 1:4 * t + 2, 1) = q_limit_low;       
%         xupp(4 * t + 3:4 * t + 4, 1) = dq_limit_upp;
%         xlow(4 * t + 3:4 * t + 4, 1) = dq_limit_low; 
        xupp(4 * t + 3:4 * t + 4, 1) = dq_final;
        xlow(4 * t + 3:4 * t + 4, 1) = dq_final;
    elseif t == T - 1
        xupp(4 * t + 1:4 * t + 2, 1) = q_final;
        xlow(4 * t + 1:4 * t + 2, 1) = q_final;  
%         xupp(4 * t + 1:4 * t + 2, 1) = q_limit_upp;
%         xlow(4 * t + 1:4 * t + 2, 1) = q_limit_low;       
%         xupp(4 * t + 3:4 * t + 4, 1) = dq_limit_upp;
%         xlow(4 * t + 3:4 * t + 4, 1) = dq_limit_low; 
        xupp(4 * t + 3:4 * t + 4, 1) = dq_final;
        xlow(4 * t + 3:4 * t + 4, 1) = dq_final;
    else
        % Set joint limits and joint velocity limits
        xupp(4 * t + 1:4 * t + 2, 1) = q_limit_upp;
        xlow(4 * t + 1:4 * t + 2, 1) = q_limit_low;
        xupp(4 * t + 3:4 * t + 4, 1) = dq_limit_upp;
        xlow(4 * t + 3:4 * t + 4, 1) = dq_limit_low;
    end
    
    if t > 0
        % Set torque limits
        u_start_idx = 4 * (T + 1) + 1;
        xupp(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1, 1) = u_limit_upp;
        xlow(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1, 1) = u_limit_low;
        
        % Set the contact force magnitude limits
        lmda_start_idx = 4 * (T + 1) + 2 * T + 1;
        xupp(lmda_start_idx + t - 1, 1) = lmda_limit_upp;
        xlow(lmda_start_idx + t - 1, 1) = lmda_limit_low;
        
        % Set the time step limits
        h_start_idx = 4 * (T + 1) + 3 * T + 1;
        xupp(h_start_idx + t - 1, 1) = h_limit_upp;
        xlow(h_start_idx + t - 1, 1) = h_limit_low;
    end
end

% TODO what do these do?
% xmul = zeros(5 * (T + 1) + 2 * T + T, 1);
% xstate = zeros(5 * (T + 1) + 2 * T + T, 1);
xmul = [];
xstate = [];

% Set upper and lower bounds on the objective and constraints
Flow = zeros(7 * T + 2 + n_h_constraints, 1);
Fupp = zeros(7 * T + 2 + n_h_constraints, 1);
for k = 1:7 * T + 2 + n_h_constraints
    if k == 1
        % Bounds on the objective function
        Flow(k, 1) = 0;
        Fupp(k, 1) = 100;
%     elseif k >= 5 * T + 1 && k <= 6 * T + 1
    elseif k > 4 * T + 1 && k <= 5 * T + 1
        Flow(k, 1) = 0;
        Fupp(k, 1) = 1e-2;
    elseif k > 5 * T + 1 && k <= 7 * T + 1
        % Signed distance function bounds
        Flow(k, 1) = phi_limit_low;
        Fupp(k, 1) = phi_limit_upp;
    elseif k == 7 * T + 2
        % Bounds on final end-effector position constraint
        Flow(k, 1) = 0;
        Fupp(k, 1) = final_pos_thresh;
    else
        % All equality constraints, so upper and lower bounds both 0
        Flow(k, 1) = 0;
        Fupp(k, 1) = 0;
    end
end

% TODO what do these do?
% Fmul = zeros(6 * T + 2, 1);
% Fstate = zeros(6 * T + 2, 1);
Fmul = [];
Fstate = [];

% userfun(x)

xlow
xupp
Flow
Fupp

snprint('contact.out');
snspec('settings.spc');
snscreen on;
[x, F, INFO, xmul, Fmul, xstate, Fstate, output] = snopt(x, xlow, xupp, xmul, xstate, ...
                Flow, Fupp, Fmul, Fstate, @userfun, 0, 1);

snprint off;
snend;

F
x

qs = zeros(2, T);
for t = 0:T
    qs(:, t + 1) = get_q(x, t);
end

two_link_draw(qs, h, 'test2.gif');

us_1 = zeros(1, T);
us_2 = zeros(1, T);
dqs_1 = zeros(1, T);
dqs_2 = zeros(1, T);
ts = linspace(1, T, T);
u_opp = zeros(3, T);
for t = 1:T
    u = get_u(x, t, T);
    us_1(1, t) = u(1);
    us_2(1, t) = u(2);
    
    dq = get_dq(x, t);
    dqs_1(1, t) = dq(1);
    dqs_2(1, t) = dq(2);
    
    q = get_q(x, t);
    J = two_link_body_manip_jac(q);
    [g_s1, g_st] = two_link_kinematics(q(1), q(2));
%     g_s1(1:3, 4) = zeros(3, 1);
%     adj = adjoint(g_s1);
%     lmda = adj' * [0; 0; 1; 0; 0; 0];
    lmda = [0; 0; 1; 0; 0; 0];

%     disp("---");
%     disp(J');
%     disp(lmda);
    
    u_opp(1:2, t) = J' * lmda;
    u_opp(3, t) = two_link_signed_dist(q);
    
    fprintf('elbow z: %f\n', g_s1(3, 4));
    fprintf('lambda %d: %f\n', t, get_lmda(x, t, T));
end

u_opp

% for t = 1:T
%     q = get_q(x, t);
%     dq = get_dq(x, t);
%     fprintf('q(%d) = (%f, %f)\ndq(%d) = (%f, %f)\n', t, q(1), q(2), t, dq(1), dq(2));
% end
% for t = 1:T
%     u = get_u(x, t, T);
%     fprintf('u(%d) = (%f, %f)\n', t, u(1), u(2));
% end
% for t = 1:T
%     lmda = get_lmda(x, t, T);
%     fprintf('lmda(%d) = %f\n', t, lmda);
% end

% Plot controls (joint torques)
figure(2);

subplot(2, 1, 1);
plot(ts, us_1);
axis([1 T u_limit_low(1) - 1 u_limit_upp(1) + 1]);
xlabel('time $t$', 'Interpreter', 'latex');
ylabel('$u_1(t)$', 'Interpreter', 'latex');

subplot(2, 1, 2);
plot(ts, us_2);
axis([1 T u_limit_low(2) - 1 u_limit_upp(2) + 1]);
xlabel('time $t$', 'Interpreter', 'latex');
ylabel('$u_2(t)$', 'Interpreter', 'latex');

% Plot joint velocities
figure(3);

subplot(2, 1, 1);
plot(ts, dqs_1);
axis([1 T min(dqs_1) - 1 max(dqs_1) + 1]);
xlabel('time $t$', 'Interpreter', 'latex');
ylabel('$\dot{q}_1(t)$', 'Interpreter', 'latex');

subplot(2, 1, 2);
plot(ts, dqs_2);
axis([1 T min(dqs_2) - 1 max(dqs_2) + 1]);
xlabel('time $t$', 'Interpreter', 'latex');
ylabel('$\dot{q}_2(t)$', 'Interpreter', 'latex');

function [F] = userfun(x)
% T = 10;
% T = 15;
T = 20;
% h = 0.1;
pos_final = [0; 1; 1];
% pos_final = [0; 1.75; 0];
% pos_final = [0; 2.0; 0];
% pos_final = [0; 1.55; 0.2];

q_cons = zeros(2 * T, 1);
dynam_cons = zeros(2 * T, 1);
lmda_cons = zeros(T, 1);
signed_dist_cons = zeros(T, 1);
signed_dist_eef_cons = zeros(T, 1);

for t = 0:T - 1
    h = get_h(x, t + 1, T);
    
    % Constraints give the relationship between joint angles and velocities
    q_cons(2 * t + 1:2 * t + 2, 1) = get_q(x, t + 1) - get_q(x, t) - h * get_dq(x, t + 1);    
    % Constraints give the relationship between torques, angles, and
    % velocities
    q = get_q(x, t + 1);
    dq = get_dq(x, t + 1);
    [M, C, N] = two_link_dynamics(q(1), q(2), dq(1), dq(2));
    J = two_link_body_manip_jac(q);
%     J = zeros(6, 2);
    [g_s1, g_st] = two_link_kinematics(q(1), q(2));
    
    % Transform the wrench applied at the origin of a frame translated from 
    % the inertial from the the point of contact to a wrench applied at the
    % origin of the end-effector frame 
%     g_st(1:3, 4) = zeros(3, 1);
%     adj = adjoint(g_st);
%     g_s1(1:3, 4) = zeros(3, 1);
%     adj = adjoint(g_s1);
%     lmda = adj' * [0; 0; get_lmda(x, t + 1, T); 0; 0; 0];
    lmda = [0; 0; get_lmda(x, t + 1, T); 0; 0; 0];

%     lmda_force = g_st * [0; 0; get_lmda(x, t + 1, T); 0];
%     lmda = [lmda_force(1:3); 0; 0; 0];

%     lmda = [0; 0; get_lmda(x, t + 1, T); 0; 0; 0];
    
%     dynam_cons(2 * t + 1:2 * t + 2, 1) = dq - get_dq(x, t) - h * get_u(x, t + 1, T);

    dynam_cons(2 * t + 1:2 * t + 2, 1) = M * (dq - get_dq(x, t)) + h * (C * dq + N - get_u(x, t + 1, T) - J' * lmda); 

    dist = two_link_signed_dist(q);
    signed_dist_cons(t + 1, 1) = dist;
    signed_dist_eef_cons(t + 1, 1) = two_link_signed_dist_eef(q);
    lmda_cons(t + 1, 1) = dist * get_lmda(x, t + 1, T);
end

% Constrain the final position of the end effector
% Note: we actually constrain the last k positions of the end effector
% because we want the manipulator to stabilize at that position
d = 0;
% k = 0;
k = 1;
% k = 2;
% for t = T - k:T
%     q = get_q(x, t);
%     [~, g_st] = two_link_kinematics(q(1), q(2));
%     d = d + norm(g_st(1:3, 4) - pos_final);
% end

% Constrain time steps
% TODO change this so that h(1) + h(2) = c, h(3) + h(4) = c, ..., where c
% is constrained to be within some small interval
n_h_constraints = floor((T - 3) / 2);
h_cons = zeros(n_h_constraints, 1);
for j = 1:n_h_constraints
    h_cons(j) = get_h(x, 2 * j - 1, T) + get_h(x, 2 * j, T) - get_h(x, 2 * j + 1, T) - get_h(x, 2 * j + 2, T);
end

% obj = 0;

% Define objective as squared joint torques at the end
% obj = get_u(x, T, T)' * get_u(x, T, T);
% obj = get_u(x, T - 1, T)' * get_u(x, T - 1, T) + get_u(x, T, T)' * get_u(x, T, T) + get_dq(x, T)' * get_dq(x, T);

% obj = get_u(x, T, T)' * get_u(x, T, T) + get_dq(x, T)' * get_dq(x, T);

% Define objective as the sum of squared joint torques
obj = 0;
% for t = 1:T
%     obj = obj + get_h(x, t, T) * get_u(x, t, T)' * get_u(x, t, T);
%     obj = obj + norm(get_u(x, t, T)) + norm(get_dq(x, t));
% end

obj = obj + get_u(x, t, T)' * get_u(x, t, T);

% obj = obj + get_dq(x, T)' * get_dq(x, T);

F = [obj; q_cons; dynam_cons; lmda_cons; signed_dist_cons; signed_dist_eef_cons; d; h_cons];
end

function [u] = get_u(x, t, T)
u_start_idx = 4 * (T + 1) + 1;
u = x(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1);
end

function [q] = get_q(x, t)
q = x(4 * t + 1:4 * t + 2);
end

function [dq] = get_dq(x, t)
dq = x(4 * t + 3:4 * t + 4);
end

function [lmda] = get_lmda(x, t, T)
lmda_start_idx = 4 * (T + 1) + 2 * T + 1;
lmda = x(lmda_start_idx + t - 1, 1);
end

function [h] = get_h(x, t, T)
h_start_idx = 4 * (T + 1) + 3 * T + 1;
h = x(h_start_idx + t - 1, 1);
end
