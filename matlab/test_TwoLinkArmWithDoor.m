% Test the dynamics
arm = TwoLinkArmWithDoor(2, 1, 1, 1, 1, 1, 1, 1);

q_init = [pi / 2; pi / 2; -pi / 2];
dq_init = [0.01; 0.01; 0.01];
dt = 0.01;
T = 5000;

%path = '/home/abajcsy/hybrid_ws/src/contact_manip/matlab/';
filename = ''; %strcat(path,'doorsim.gif');

% arm.simulate(q_init, dq_init, dt, T, filename);

%% Unit test line segment interesection 
% assert(arm.intersect([0; 0], [1; 0], [0; 1], [1; 1]) == Inf);
% assert(arm.intersect([-1; 0], [1; 0], [0; -1], [0; 1]) == 0.5);
% assert(arm.intersect([0; -1], [0; 1], [-1; 0], [1; 0]) == 0.5);
% assert(arm.intersect([-1; 0], [1; 0], [1; -1], [1; 1]) == 1.0);
% assert(arm.intersect([1; -1], [1; 1], [-1; 0], [1; 0]) == 1.0);
% assert(arm.intersect([-1; 0], [1; 0], [-1; -1], [-1; 1]) == 0);
% assert(arm.intersect([-1; -1], [-1; 1], [-1; 0], [1; 0]) == 0);
arm.intersect([-1; -1], [-1; 1], [-1; 0], [1; 0])

figure(1);
hold on;

for i = 1:1
   l1_low = [rand(); rand()];
   l1_high = [rand(); rand()];
   l2_low = [rand(); rand()];
   l2_high = [rand(); rand()];
   
   plot([l1_low(1) l1_high(1)], [l1_low(2) l1_high(2)], 'b');
   plot([l2_low(1) l2_high(1)], [l2_low(2) l2_high(2)], 'g');
   
   fprintf("intersect((%f, %f), (%f, %f), (%f, %f), (%f, %f))?\n", l1_low(1), l1_low(2), l1_high(1), l1_high(2), l2_low(1), l2_low(2), l2_high(1), l2_high(2));
   s = arm.intersect(l1_low, l1_high, l2_low, l2_high)
   
   p = (1 - s) * l1_low + s * l1_high;
   scatter(p(1), p(2), 'rx');
end