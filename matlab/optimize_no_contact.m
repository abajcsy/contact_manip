% State x is are the joint states (angles and velocities) proceeded by the
% controls (joint torques)
% x = [[th1; th2; dth1; dth2]; ...; [u1; u2]; ...]

% Constraints are dth(t+1) - th(t) - h*th(t+1) = 0
%                 M(t+1)*(dth(t+1) - dth(t)) + h*(C(t+1)*dth(t+1) + N(t+1)
%                 - u(t+1)) = 0
% Objective is u(T)'*u(T)

T = 10;
h = 0.1;

% q_init = [-1.5690; 1.5685];
% q_init = [0; pi / 2];
% q_init = [0; pi / 4];
q_init = [0; 0];
dq_init = [0; 0];
% q_final = [0.5; 0.25];
% q_final = [-0.75; 0.5];
q_final = [0; pi / 2];
dq_final = [0; 0];

% Tolerance on the final position of the end effector
% final_pos_thresh = 0.12;
final_pos_thresh = 0.01;

% Limits on the joint angles
q_limit_upp = [pi / 2; pi / 2];
q_limit_low = [-pi / 2; -pi / 2];

% Limits on the joint velocities
dq_limit_upp = [20; 20];
dq_limit_low = [-20; -20];

% Limits on the joint torques
u_limit_upp = [5; 5];
u_limit_low = [-5; -5];

% Construct an initial value of x
x = zeros(4 * (T + 1) + 2 * T, 1);
qs = zeros(2, T);
for t = 0:T
    k = t / T;
    x(4 * t + 1:4 * t + 2, 1) = (1 - k) * q_init + k * q_final;
    qs(:, t + 1) = get_q(x, t);
end

% two_link_draw(qs, h, '');

% Set upper and lower bounds on x
xupp = zeros(4 * (T + 1) + 2 * T, 1);
xlow = zeros(4 * (T + 1) + 2 * T, 1);
for t = 0:T
    if t == 0
        % Set initial config constraint
        xupp(4 * t + 1:4 * t + 2, 1) = q_init;
        xlow(4 * t + 1:4 * t + 2, 1) = q_init;
        xupp(4 * t + 3:4 * t + 4, 1) = dq_init;
        xlow(4 * t + 3:4 * t + 4, 1) = dq_init;
%     elseif t == T
%         % Set final config constraint
%         xupp(4 * t + 1:4 * t + 2, 1) = q_final;
%         xlow(4 * t + 1:4 * t + 2, 1) = q_final;
%         xupp(4 * t + 3:4 * t + 4, 1) = dq_final;
%         xlow(4 * t + 3:4 * t + 4, 1) = dq_final;
%     elseif t == T - 1
%         % Set final config constraint
%         xupp(4 * t + 1:4 * t + 2, 1) = q_final;
%         xlow(4 * t + 1:4 * t + 2, 1) = q_final;
%         xupp(4 * t + 3:4 * t + 4, 1) = dq_final;
%         xlow(4 * t + 3:4 * t + 4, 1) = dq_final;        
    else
        % Set joint limits and joint velocity limits
        xupp(4 * t + 1:4 * t + 2, 1) = q_limit_upp;
        xlow(4 * t + 1:4 * t + 2, 1) = q_limit_low;
        xupp(4 * t + 3:4 * t + 4, 1) = dq_limit_upp;
        xlow(4 * t + 3:4 * t + 4, 1) = dq_limit_low;
    end
    
    % Set torque limits
    if t > 0
        u_start_idx = 4 * (T + 1) + 1;
        xupp(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1, 1) = u_limit_upp;
        xlow(u_start_idx + 2 * (t - 1):u_start_idx + 2 * (t - 1) + 1, 1) = u_limit_low;      
    end
end

% TODO what do these do?
% xmul = zeros(4 * (T + 1) + 2 * T, 1);
% xstate = zeros(4 * (T + 1) + 2 * T, 1);
xmul = [];
xstate = [];

% Set upper and lower bounds on the objective and constraints
Flow = zeros(4 * T + 2, 1);
Fupp = zeros(4 * T + 2, 1);
for k = 1:4 * T + 2
    if k == 1
        % Bounds on the objective function
        Flow(k, 1) = 0;
        Fupp(k, 1) = 100;
    elseif k == 4 * T + 2
        Flow(k, 1) = 0;
        Fupp(k, 1) = final_pos_thresh;
    else
        % All equality constraints, so upper and lower bounds both 0
        Flow(k, 1) = 0;
        Fupp(k, 1) = 0;
    end
end

% TODO what do these do?
% Fmul = zeros(4 * T + 2, 1);
% Fstate = zeros(4 * T + 2, 1);
Fmul = [];
Fstate = [];

% userfun(x)

snprint('no_contact.out');
snspec('settings.spc');
snscreen on;
[x, F, INFO, xmul, Fmul, xstate, Fstate, output] = snopt(x, xlow, xupp, xmul, xstate, ...
                Flow, Fupp, Fmul, Fstate, @userfun, 0, 1);

snprint off;
snend;

F

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
for t = 1:T
    u = get_u(x, t, T);
    us_1(1, t) = u(1);
    us_2(1, t) = u(2);
    dq = get_dq(x, t);
    dqs_1(1, t) = dq(1);
    dqs_2(1, t) = dq(2);
end

% Plot joint torques
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
T = 10;
h = 0.1;
pos_final = [0; 1; 1];

q_cons = zeros(2 * T, 1);
dynam_cons = zeros(2 * T, 1);

for t = 0:T - 1
    % Constraints give the relationship between joint angles and velocities
    q_cons(2 * t + 1:2 * t + 2, 1) = get_q(x, t + 1) - get_q(x, t) - h * get_dq(x, t + 1);    
    % Constraints give the relationship between torques, angles, and
    % velocities
    q = get_q(x, t + 1);
    dq = get_dq(x, t + 1);
    [M, C, N] = two_link_dynamics(q(1), q(2), dq(1), dq(2));
    dynam_cons(2 * t + 1:2 * t + 2, 1) = M * (dq - get_dq(x, t)) + h * (C * dq + N - get_u(x, t + 1, T));    
end

% Constrain the final position of the end effector
% Note: we actually constrain the last two positions of the end effector
% because we want the manipulator to stabilize at that position
d = 0;
for t = T-1:T
    q = get_q(x, t);
    [~, g_st] = two_link_kinematics(q(1), q(2));
    d = d + norm(g_st(1:3, 4) - pos_final);
end

F = [get_u(x, T, T)' * get_u(x, T, T); q_cons; dynam_cons; d];
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
    