function drop_that_link()
% m*ddx     = 0
% m*ddy     = -m*g + lambda
% I*ddtheta = -(l/2)*cos(theta)*lambda

l = 1.0;
m = 0.1;
g = 9.81;
I = (m*l^2)/12.0;


T = 10;
dt = 0.001;
t = 0;

q_init = [0; 1; pi/6; 0; 0; 0]; % x, y, theta, dx, dy, dtheta
q = q_init;

q_traj = zeros(T/dt, 6);
lambda_traj = zeros(T/dt, 1);

fprintf("Simulating for %f s...\n", T);

while t < T

    idx = floor(t/dt)+1;
    q_traj(idx,:) = q';
    
    % positions
    x = q(1);
    y = q(2);
    theta = q(3);
    
    % velocities
    dx = q(4);
    dy = q(5);
    dtheta = q(6);
    
    A = dt^2/m; %- (l/2)*sin(theta + dt*dtheta - ((dt^2*l)/2*I)*cos(theta)*lambda);
    b = y + dt*dy - dt^2*g;

    % x = LCP(M,q) solves the LCP
    % 
    % x >= 0 
    % Mx + q >= 0 
    % x'(Mx + q) = 0
    [w, lambda, retcode] = LCPSolve(A, b);
    lambda
    lambda_traj(idx,1) = lambda; 
    
    % update x
    x = x + dt*dx;
    
    % update y, dx
    ddy = -g + lambda/m;    
    dy = dy + dt * ddy;
    y = y + dt * dy;
    
    % update theta, dtheta
    ddtheta = -(l/(2*I))*cos(theta)*lambda;
    dtheta = dtheta + dt * ddtheta;
    theta = theta + dt * dtheta;
    
    q = [x; y; theta; dx; dy; dtheta];
    
    t = t + dt;
end

fprintf("Showing trajectory...\n");

simulate(q_traj, lambda_traj, dt);

%------------ HELPER FUNCTIONS -----------%

function simulate(q_traj, lambda_traj, dt)
    figure(1);

    q_curr = q_traj(1,:);
    lamb_curr = lambda_traj(1,1);

    hold on;
    axis([-3 3 -3 3]);
    
    [p0, p1, p2] = get_points(q_curr);
    link1 = plot([p1(1) p2(1)], [p1(2) p2(2)], 'k', 'LineWidth', 2);
    joint1 = plot(p0(1), p0(2), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    dir = 0;
    len = lamb_curr;
    q = quiver(0,0,dir,len);
    set(q,'LineWidth',1);
    
    time = 0;
    h = text(1.1,1.8,strcat('t=',num2str(time)));
    la = text(1.1,1.6,strcat('lambda=',num2str(lamb_curr)));
    
	for t = 1:size(q_traj, 1)
        q_curr = q_traj(t,:);
        lamb_curr = lambda_traj(t,1);
        
        [p0, p1, p2] = get_points(q_curr);
        set(link1, 'Xdata', [p1(1) p2(1)]);
        set(link1, 'Ydata', [p1(2) p2(2)]);
        set(joint1, 'Xdata', p0(1));
        set(joint1, 'Ydata', p0(2));
        
        set(h, 'String', strcat('t=',num2str(time)));
        set(la, 'String', strcat('lambda=',num2str(lamb_curr)));
        
        if lamb_curr ~= 0
            set(q, 'Vdata', lamb_curr);
            set(q, 'Color', 'red');
        end
        
        pause(dt);
        time = time + dt;
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
    dtgradh = [0; 0; (l/2)*sin(theta)*dtheta];
end
end