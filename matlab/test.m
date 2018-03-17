th1 = 0.0;
th2 = 0.0;

for i = 1:10
    draw_arm(th1, th2);
    th1 = th1 + 0.1;
    th2 = th2 + 0.05;
    pause(0.5);
end


function draw_arm(th1, th2)
[g_sl1, g_st] = two_link_kinematics(th1, th2);

p0 = [0 0];
p1 = [g_sl1(2, 4) g_sl1(3, 4)];
p2 = [g_st(2, 4) g_st(3, 4)];

figure(1);
hold on;
plot([p0(1) p1(1) p2(1)], [p0(2) p1(2) p2(2)], '*');
plot([p0(1) p1(1)], [p0(2) p1(2)]);
plot([p1(1) p2(1)], [p1(2) p2(2)]);
axis([-2 2 0 3]);
end
