function plot_velocities(traj,N,filename)
clf

fig = figure(3);
subplot(2,1,1);
time = linspace(1,N,N);
dq1 = zeros(N,1);
dq2 = zeros(N,1);
idx = 1;
for t=1:N
    dq1(t) = traj(idx+2);
    dq2(t) = traj(idx+3);
    idx = idx+4;
end
% first joint controls
plot(time, dq1, 'k');
% title('$$ Controls~Over~Time $$', 'interpreter', 'latex')
ylabel('$$ \dot{q}_1 $$', 'interpreter', 'latex'); % x-axis label
% set the viewing axes
axis([1 N -5 5]);
    
subplot(2,1,2);
plot(time, dq2, 'k');
% set the viewing axes
axis([1 N -5 5]);
ylabel('$$ \dot{q}_2 $$', 'interpreter', 'latex'); % x-axis label
xlabel('time (s)', 'interpreter', 'latex'); % y-axis label

if ~isempty(filename)
    saveas(fig,filename);
end
end