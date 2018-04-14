% Test the dynamics
arm = TwoLinkArm(2, 1, 1, 1, 1, 1);

q_init = [2*pi/3; pi/4];
dq_init = [-0.1; 0];
dt = 0.01;
T = 20000;

arm.simulate_no_contact(q_init, dq_init, dt, T);
