function [] = two_link_draw(qs, h, filename)
%TWO_LINK_DRAW Draw the two link arm at a sequence of configurations
%   Detailed explanation goes here
fig = figure(1);
hold on;

for t = 1:size(qs, 2)
    [g_sl1, g_st] = two_link_kinematics(qs(1, t), qs(2, t));

    p0 = [0 0];
    p1 = [g_sl1(2, 4) g_sl1(3, 4)];
    p2 = [g_st(2, 4) g_st(3, 4)];

    plot([p0(1) p1(1) p2(1)], [p0(2) p1(2) p2(2)], 'r*');
    plot([p0(1) p1(1)], [p0(2) p1(2)], 'b');
    plot([p1(1) p2(1)], [p1(2) p2(2)], 'b');
    
    ptext = [0; 0; 0.05; 1];
    ptext = g_st * ptext;
    
%     text(p2(1), p2(2) + 0.5, strcat('t = ', num2str(t)));
    text(ptext(2), ptext(3), strcat('t = ', num2str(t - 1)));
    axis([-3 3 0 3]);
%     axis([-3 3 -3 3]);
    
    if ~isempty(filename)
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256); 
        if t == 1 
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf); 
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
        end
    end

    
    pause(h);
end

hold off;
end

