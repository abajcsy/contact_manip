% Test the dynamics
arm = TwoLinkArmWithDoor(2, 1, 1, 1, 1, 1, 1, 1);

q_init = [pi / 2; pi / 2; -pi / 2];
dq_init = [0.01; 0.01; 0.01];
dt = 0.01;
T = 5000;

%path = '/home/abajcsy/hybrid_ws/src/contact_manip/matlab/';
filename = ''; %strcat(path,'doorsim.gif');

arm.simulate(q_init, dq_init, dt, T, filename);