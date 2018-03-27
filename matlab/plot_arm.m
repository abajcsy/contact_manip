%% Plot arm trajectory over time
function plot_arm(traj,times,N,filename)
clf
fig = figure(1);
hold on;

for t = 0:N
    q = get_q(traj, t);
    [g_sl1, g_st] = two_link_kinematics(q(1), q(2));

    p0 = [0 0];
    p1 = [g_sl1(2, 4) g_sl1(3, 4)];
    p2 = [g_st(2, 4) g_st(3, 4)];

    alpha = (t+1)/(N+1);
    
    link1 = plot([p0(1) p1(1)], [p0(2) p1(2)], 'k');
    link2 = plot([p1(1) p2(1)], [p1(2) p2(2)], 'k');
    set(link1,'LineWidth',5);
    set(link2,'LineWidth',5);
    link1.Color(4) = alpha;
    link2.Color(4) = alpha;
    % plot the delta t
    if t ~= 0
        dt_text = text(p2(1)+0.05,p2(2)+0.1,strcat('dt=',num2str(times(t))));
        dt_text.Color = [0.5 0.5 0.5];
    end
    % plot the sequence in time
    %text(p2(1)-0.2,p2(2)+0.2,strcat('t=',num2str(t)));
    % color the two active joints in red
    scatter([p0(1) p1(1)], [p0(2) p1(2)], 'o', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    % color the EE separately
    scatter([p2(1)], [p2(2)], 'o', 'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    
    % set the viewing axes
    axis([-2 2 0 3]);
    
    % remove tick labels
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    
    if ~isempty(filename)
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256); 
        if t == 0 
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf); 
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
        end
    end

    if t ~= N
        pause(times(t+1));
    end
end
hold off;
end

%% Getting a configuration at a specific time
function [q] = get_q(x, t)
    q = x(4 * t + 1:4 * t + 2);
end