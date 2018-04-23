arm = TwoLinkArmWithDoor(2, 1, 1, 1, 1, 1, 1, 2);

%% Unit test the arm and door dynamics
q_init = [pi / 2; pi / 2; -pi / 2];
dq_init = [0.01; 0.01; 0.01];
dt = 0.01;
T = 5000;

%path = '/home/abajcsy/hybrid_ws/src/contact_manip/matlab/';
filename = ''; %strcat(path,'doorsim.gif');

% arm.simulate(q_init, dq_init, dt, T, filename);

%% Unit test line segment interesection 
[s1, s2] = arm.intersection([0; 0], [1; 0], [0; 1], [1; 1]);
assert(s1 == Inf && s2 == Inf);

[s1, s2] = arm.intersection([-1; 0], [1; 0], [0; -1], [0; 1]);
assert(s1 == 0.5 && s2 == 0.5);

[s1, s2] = arm.intersection([-1; 0], [1; 0], [1; -1], [1; 1]);
assert(s1 == 1 && s2 == 0.5);

[s1, s2] = arm.intersection([-1; 0], [1; 0], [-1; -1], [-1; 1]);
assert(s1 == 0 && s2 == 0.5);

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
   [s1, s2] = arm.intersection(l1_low, l1_high, l2_low, l2_high);
   
   p1 = (1 - s1) * l1_low + s1 * l1_high;
   p2 = (1 - s2) * l2_low + s2 * l2_high;
   scatter(p1(1), p1(2), 'rx');
   scatter(p2(1), p2(2), 'rx');
end

%% Unit test for the signed distance function
q = [pi/2 + 0.1; -pi/8; -2*pi/5];
% arm.show_state_debug(q, 2);

q = [pi/2; 0; -pi/2];
% arm.show_state_debug(q, 3);

q_min = [0; 0; -pi];
q_max = [pi; pi; 0];
for i = 1:1
    q = rand_q(q_min, q_max);
    arm.show_state_debug(q, 2);
end

function q = rand_q(q_min, q_max)
    q = [q_min(1) + rand() * (q_max(1) - q_min(1));
         q_min(2) + rand() * (q_max(2) - q_min(2));
         q_min(3) + rand() * (q_max(3) - q_min(3))];
end