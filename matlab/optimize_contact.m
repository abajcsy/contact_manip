% State x is are the joint states (angles and velocities) proceeded by the
% controls (joint torques), proceeded by the contact forces lambda
% x = [[th1; th2; dth1; dth2]; ...; [u1; u2]; ...; lmda1; ...]

% Constraints are dth(t+1) - th(t) - h*th(t+1) = 0
%                 M(t+1)*(dth(t+1) - dth(t)) + h*(C(t+1)*dth(t+1) + N(t+1)
%                 - u(t+1) - J(t+1)'*lmda(t+1)) = 0
% Objective is u(T)'*u(T)

T = 15;
h = 0.1;

q_init = [0; 0];
% q_init = [-pi / 2; 0];
% q_init = [0; pi / 2]; 
% q_init = [0; pi / 4];
% q_final = [0.5; 0.25];
% q_final = [-0.75; 0.5];
% q_final = [0; 0];
q_final = [-pi / 2; 0];

% Tolerance on the final position of the end effector
final_pos_thresh = 0.12;

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
% u_limit_upp = [20; 20];
% u_limit_low = [-20; -20];

% Limit on the contact force magnitudes
lmda_limit_upp = 100;
lmda_limit_low = 0;

% Construct an initial value of x
x = zeros(4 * (T + 1) + 2 * T + T, 1);
% qs = zeros(2, T);
% for t = 0:T
%     k = t / T;
%     x(4 * t + 1:4 * t + 2, 1) = (1 - k) * q_init + k * q_final;
%     qs(:, t + 1) = get_q(x, t);
% end

% two_link_draw(qs, h, '');

% Set upper and lower bounds on x
xupp = zeros(4 * (T + 1) + 2 * T + T, 1);
xlow = zeros(4 * (T + 1) + 2 * T + T, 1);
for t = 0:T
    if t == 0
        % Set initial config constraint
        xupp(4 * t + 1:4 * t + 2, 1) = q_init;
        xlow(4 * t + 1:4 * t + 2, 1) = q_init;
%     elseif t == T
%         % Set final config constraint
%         xupp(4 * t + 1:4 * t + 2, 1) = q_final;
%         xlow(4 * t + 1:4 * t + 2, 1) = q_final;
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
    end
end

% TODO what do these do?
% xmul = zeros(5 * (T + 1) + 2 * T + T, 1);
% xstate = zeros(5 * (T + 1) + 2 * T + T, 1);
xmul = [];
xstate = [];

% Set upper and lower bounds on the objective and constraints
Flow = zeros(6 * T + 2, 1);
Fupp = zeros(6 * T + 2, 1);
for k = 1:6 * T + 2
    if k == 1
        % Bounds on the objective function
        Flow(k, 1) = 0;
        Fupp(k, 1) = 100;
    elseif k >= 5 * T + 1 && k <= 6 * T + 1
        Flow(k, 1) = -100;
        Fupp(k, 1) = 100;
    elseif k == 6 * T + 2
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

% xlow
% xupp
% Flow
% Fupp

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
ts = linspace(1, T, T);
u_opp = zeros(3, T);
for t = 1:T
    u = get_u(x, t, T);
    us_1(1, t) = u(1);
    us_2(1, t) = u(2);
    
    q = get_q(x, t);
    [J] = two_link_body_manip_jac(q);
    
    [~, g_st] = two_link_kinematics(q(1), q(2));
    
%     tmp = g_st * [0; 0; 1; 1];
    
%     lmda = [tmp(1:3, 1); 0; 0; 0];
    g_st(1:3, 4) = zeros(3, 1);
    adj = adjoint(g_st);
    disp(adj);
    lmda = adj' * [0; 0; 1; 0; 0; 0];
    
%     u_opp(1:2, t) = J' * [0; 0; 1; 0; 0; 0];
    u_opp(1:2, t) = J' * lmda;
    u_opp(3, t) = two_link_signed_dist(q);
    
    fprintf('lambda %d: %f\n', t, get_lmda(x, t, T));
end

u_opp

for t = 1:T
    q = get_q(x, t);
    dq = get_dq(x, t);
    fprintf('q(%d) = (%f, %f)\ndq(%d) = (%f, %f)\n', t, q(1), q(2), t, dq(1), dq(2));
end
for t = 1:T
    u = get_u(x, t, T);
    fprintf('u(%d) = (%f, %f)\n', t, u(1), u(2));
end
for t = 1:T
    lmda = get_lmda(x, t, T);
    fprintf('lmda(%d) = %f\n', t, lmda);
end

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

function [F] = userfun(x)
T = 15;
h = 0.1;
% pos_final = [0; 1; 1];
pos_final = [0; 1.75; 0];
% pos_final = [0; 2.0; 0];
% pos_final = [0; 1.55; 0.2];

q_cons = zeros(2 * T, 1);
dynam_cons = zeros(2 * T, 1);
lmda_cons = zeros(T, 1);
signed_dist_cons = zeros(T, 1);

for t = 0:T - 1
    % Constraints give the relationship between joint angles and velocities
    q_cons(2 * t + 1:2 * t + 2, 1) = get_q(x, t + 1) - get_q(x, t) - h * get_dq(x, t + 1);    
    % Constraints give the relationship between torques, angles, and
    % velocities
    q = get_q(x, t + 1);
    dq = get_dq(x, t + 1);
    [M, C, N] = two_link_dynamics(q(1), q(2), dq(1), dq(2));
    % TODO looks like the manipulator Jacobian has the wrong sign (why?)
    J = -two_link_body_manip_jac(q);
%     J = zeros(6, 2);
    [~, g_st] = two_link_kinematics(q(1), q(2));
    
    % Transform the wrench applied at the origin of a frame translated from 
    % the inertial from the the point of contact to a wrench applied at the
    % origin of the end-effector frame 
    g_st(1:3, 4) = zeros(3, 1);
    adj = adjoint(g_st);
    lmda = adj' * [0; 0; get_lmda(x, t + 1, T); 0; 0; 0];
    
%     lmda_force = g_st * [0; 0; get_lmda(x, t + 1, T); 0];
%     lmda = [lmda_force(1:3); 0; 0; 0];

%     lmda = [0; 0; get_lmda(x, t + 1, T); 0; 0; 0];
    
    dynam_cons(2 * t + 1:2 * t + 2, 1) = M * (dq - get_dq(x, t)) + h * (C * dq + N - get_u(x, t + 1, T) - J' * lmda); 
    
    dist = two_link_signed_dist(q);
    signed_dist_cons(t + 1, 1) = dist;
    lmda_cons(t + 1, 1) = dist * get_lmda(x, t + 1, T);
    % TODO add in T constraints that dist(t) >= 0 (i.e. prevent manipulator
    % from going through the contact surface)
end

% Constrain the final position of the end effector
% Note: we actually constrain the last k positions of the end effector
% because we want the manipulator to stabilize at that position
d = 0;
% k = 1;
k = 2;
for t = T - k:T
    q = get_q(x, t);
    [~, g_st] = two_link_kinematics(q(1), q(2));
    d = d + norm(g_st(1:3, 4) - pos_final);
end

obj = get_u(x, T, T)' * get_u(x, T, T);
% obj = 0;
% for t = 1:T
%     obj = obj + get_u(x, t, T)' * get_u(x, t, T);
% end

F = [obj; q_cons; dynam_cons; lmda_cons; signed_dist_cons; d];
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
