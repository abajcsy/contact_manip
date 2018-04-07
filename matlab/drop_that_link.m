function drop_that_link()

l = 1.0;
m = 0.1;
g = 9.81;
I = (m*l^2)/12.0;

M = [m 0 0; 0 m 0; 0 0 I];
F = [0; m*g; 0];
T = 5;
dt = 0.001;
t = 0;
q_init = [0; 0.5; pi/2]; % x, y, theta
dq_init = [0; 0; 0];

q = q_init;
dq = dq_init;

q_traj = zeros(T/dt, 3);

while t < 10

    idx = floor(t/dt)+1;
    q_traj(idx,:) = q';
    
    theta = q(3);
    dtheta = dq(3);
    
    A_nu = (1+3*(cos(theta))^2)/m;
    w_u = (l/2)*(dtheta^2)*sin(theta)-g;
    
    % x = LCP(M,q) solves the LCP
    % 
    % x >= 0 
    % Mx + q >= 0 
    % x'(Mx + q) = 0
    lambda = LCP(A_nu, w_u);
    
    gradh = compute_gradh(theta);
    ddq = M\(gradh*lambda - F);
    dq = dq + dt * ddq;
    q = q + dt * dq;
    
    t = t + dt;
end

simulate(q_traj);

%------------ HELPER FUNCTIONS -----------%

function simulate(q_traj)
    figure(1);

    q_curr = q_traj(1,:);

    hold on;
    axis([-3 3 -1 5]);
    
    [p0, p1, p2] = get_points(q_curr);
    link1 = plot([p1(1) p2(1)], [p1(2) p2(2)], 'k', 'LineWidth', 2);
    joint1 = plot(p0(1), p0(2), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

    for t = 1:size(q_traj, 1)
        q_curr = q_traj(t,:);
        [p0, p1, p2] = get_points(q_curr);
        set(link1, 'Xdata', [p1(1) p2(1)]);
        set(link1, 'Ydata', [p1(2) p2(2)]);
        set(joint1, 'Xdata', p0(1));
        set(joint1, 'Ydata', p0(2));
        
        pause(0.01);
    end

    hold off;            
end

function [p0, p1, p2] = get_points(q)
    % gets endpoints of link
    x = q(1); y = q(2); th = q(3);
    u = [cos(th); sin(th)];
    p0 = [x; y];
    p1 = (l/2)*u + p0;
    p2 = -(l/2)*u + p0;
end

function gradh = compute_gradh(theta)
    gradh = [0; 1; -(l/2)*cos(theta)];
end

function dtgradh = compute_dtgradh(theta, dtheta)
    dtgradh = [0; 0; (l/2)*sin(theta)*dtheta]
end

end