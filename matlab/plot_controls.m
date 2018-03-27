%% Plots the controls over time
function plot_controls(u,N,filename)

clf
fig = figure(2);
subplot(2,1,1);
time = linspace(1,N,N);
u1 = zeros(N,1);
u2 = zeros(N,1);
idx = 1;
for t=1:N
    u1(t) = u(idx);
    u2(t) = u(idx+1);
    idx = idx+2;
end
% first joint controls
plot(time, u1, 'k');
% title('$$ Controls~Over~Time $$', 'interpreter', 'latex')
ylabel('$$ u_1 $$', 'interpreter', 'latex'); % x-axis label
% set the viewing axes
axis([1 N -8 8]);
    
subplot(2,1,2);
plot(time, u2, 'k');
% set the viewing axes
axis([1 N -8 8]);
ylabel('$$ u_2 $$', 'interpreter', 'latex'); % x-axis label
xlabel('time (s)', 'interpreter', 'latex'); % y-axis label

if ~isempty(filename)
    saveas(fig,filename);
end