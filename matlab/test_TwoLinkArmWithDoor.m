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

%% Unit test for collision checks between states
q1 = [pi / 2; 0; -pi / 2; 0; 0; 0];
q2 = [3 * pi / 4; 0; -pi / 2; 0; 0; 0];
arm.show_state_debug(q1, 3);
arm.show_state_debug(q2, 3);
dist_btwn = arm.signed_dist_no_contact_between(q1, q2);
fprintf('dist btwn: %f\n', dist_btwn);

delta_q_min = [-pi / 4; -pi / 4; 0];
delta_q_max = [pi / 4; pi / 4; 0];
for i = 1:1
    q1 = rand_q(q_min, q_max);
    delta_q = rand_q(delta_q_min, delta_q_max);
    q2 = q1 + delta_q;
    arm.show_state_debug(q1, 4);
    arm.show_state_debug(q2, 4);
    dist_btwn = arm.signed_dist_no_contact_between(q1, q2);
    fprintf('dist btwn: %f\n', dist_btwn);
end

%% Unit tests for edge cases
qs = [0.785398, 0.288607, 0.517969, 0.766938, 0.707911, 0.814223, 0.822567, 0.896405, 0.936396, 0.907256, 0.896517, 1.471954, 1.858691, 1.876357, 1.858023, 2.135908, 2.133786, 2.596049, 3.098776, 3.108404, 3.141593; 
      0.000000, 0.747723, 0.648851, 0.543254, 0.600778, 0.294830, 0.298460, 0.340276, 0.290194, 0.242961, 0.397431, 0.357699, 0.602253, 0.683325, 0.690393, 0.141197, 0.529739, -0.053591, 0.045865, 0.024079, -0.000000; 
      -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997, -1.308997];
% q0 = [0.785398; 0.000000; -1.308997];
% q1 = [0.746128; 0.000000; -1.243547];
% q2 = [0.706858; 0.000000; -1.178097];
% q3 = [0.667588; 0.000000; -1.112647];
% q4 = [0.628319; 0.000000; -1.047198];
% q5 = [0.589049; 0.000000; -0.981748];
% q6 = [0.549779; 0.000000; -0.916298];
% q7 = [0.510509; 0.000000; -0.850848];
% q8 = [0.471239; 0.000000; -0.785398];
% q9 = [0.431969; 0.000000; -0.719948];
% q10 = [0.392699; 0.000000; -0.654498];
% q11 = [0.353429; 0.000000; -0.589049];
% q12 = [0.314159; 0.000000; -0.523599];
% q13 = [0.274889; 0.000000; -0.458149];
% q14 = [0.235619; 0.000000; -0.392699];
% q15 = [0.196350; 0.000000; -0.327249];
% q16 = [0.157080; 0.000000; -0.261799];
% q17 = [0.117810; 0.000000; -0.196350];
% q18 = [0.078540; 0.000000; -0.130900];
% q19 = [0.039270; 0.000000; -0.065450];
% q20 = [0.000000; 0.000000; 0.000000];
% qs = [q0, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17, q18, q19, q20];

for i = 0:20
    q = qs(:, i + 1);
    arm.show_state_debug(q, 5);
%     str = input('Press enter');
    pause(0.25);
end

function q = rand_q(q_min, q_max)
    q = [q_min(1) + rand() * (q_max(1) - q_min(1));
         q_min(2) + rand() * (q_max(2) - q_min(2));
         q_min(3) + rand() * (q_max(3) - q_min(3))];
end