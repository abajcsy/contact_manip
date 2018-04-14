function drop_that_ball()

m = 0.1;
g = 9.81;

k = 0; % spring constant

T = 5;
dt = 0.001;
t = 0;
q_init = [1; 0]; % y, dy

q = q_init;

q_traj = zeros(T/dt, 2);
lambda_traj = zeros(T/dt, 1);

fprintf("Simulating for %f s...\n", T);

while t < T

    idx = floor(t/dt)+1;
    q_traj(idx,:) = q';
    
    y = q(1);
    dy = q(2);
    
    A = dt^2/m;
    b = y + dt*dy + dt^2*(-g - k*y/m);
    
    % x = LCPSolve(A,b) solves the LCP
    % 
    % x >= 0 
    % Ax + b >= 0 
    % x'(Ax + b) = 0
    [w, lambda, retcode] = LCPSolve(A, b);
    lambda
    lambda_traj(idx,1) = lambda; 
    
    ddy = -g + lambda/m - k*y/m;
    dy = dy + dt * ddy;
    y = y + dt * dy;
    
    q = [y; dy];
    
    t = t + dt;
end

fprintf("Showing trajectory...\n");

simulate(q_traj, lambda_traj, dt);

%------------ HELPER FUNCTIONS -----------%

function simulate(q_traj, lambda_traj, dt)
    fig = figure(1);

    q_curr = q_traj(1,:);
    lamb_curr = lambda_traj(1,1);
    
    hold on;
    axis([-2 2 -1 2]);
    ax = gca;
    ax.Visible = 'off';
    
    floor = rectangle('Position',[-3 -3 6 3]', 'FaceColor',[0.9 0.9 0.9], ... 
                'EdgeColor',[0 0 0], 'LineWidth',1);
    
    p0 = [0; q_curr(1)];
    joint1 = plot(p0(1), p0(2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    %lambvec = plot([0, 0], [0, lamb_curr]);
    dir = 0;
    len = lamb_curr;
    q = quiver(0,0,dir,len);
    set(q,'LineWidth',1);
   
    time = 0;
    h = text(1.1,1.8,strcat('t=',num2str(time)));
    l = text(1.1,1.6,strcat('lambda=',num2str(lamb_curr)));
    
    for t = 1:size(q_traj, 1)
        q_curr = q_traj(t,:);
        lamb_curr = lambda_traj(t,1);
        
        p0 = [0; q_curr(1)];
        set(joint1, 'Xdata', p0(1));
        set(joint1, 'Ydata', p0(2));
        
        set(h, 'String', strcat('t=',num2str(time)));
        set(l, 'String', strcat('lambda=',num2str(lamb_curr)));
        
        if lamb_curr ~= 0
            set(q, 'Vdata', lamb_curr);
            set(q, 'Color', 'red');
        end
        
        %filename = '/home/abajcsy/hybrid_ws/src/contact_manip/matlab/balldrop.gif';
        %make_gif(t, fig, filename);
        
        pause(dt);
        time = time + dt;
    end

    hold off;            
end

function make_gif(t, fig, filename)
    if ~isempty(filename)
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256); 
        if t == 1
            imwrite(imind, cm, filename, 'gif', 'DelayTime',0.001, 'Loopcount', inf); 
        else
            imwrite(imind, cm, filename, 'gif', 'DelayTime',0.001, 'WriteMode', 'append');
        end
    end
end

end