function derive_two_link_dynamics()
%DERIVE_TWO_LINK_DYNAMICS Summary of this function goes here
%   Detailed explanation goes here
syms m1 m2 l1 l2 th1 th2 dth1 dth2 I1x I1y I1z I2x I2y I2z g real;

axis1 = [1; 0; 0];
q1 = [0; 0; l1];
axis2 = [1; 0; 0];
q2 = [0; 0; l1 + l2];

g_s1_init = [1 0 0 0; 0 1 0 0; 0 0 1 0.5 * l1; 0 0 0 1];
g_s2_init = [1 0 0 0; 0 1 0 0; 0 0 1 l1 + 0.5 * l2; 0 0 0 1];

twist1 = [-cross(axis1, q1); axis1];
twist2 = [-cross(axis2, q2); axis2];

J1 = [adjoint_inv(twist_angle_to_trans(twist1, th1) * g_s1_init) * twist1 zeros(6, 1)];
J2 = [adjoint_inv(twist_angle_to_trans(twist1, th1) * twist_angle_to_trans(twist2, th2) * g_s2_init) * twist1 adjoint_inv(twist_angle_to_trans(twist2, th2) * g_s2_init) * twist2];

%% Inertia matrix
M1 = sym(zeros(6, 6));
M1(1:3, 1:3) = m1 * eye(3);
M1(4, 4) = I1x;
M1(5, 5) = I1y;
M1(6, 6) = I1z;

M2 = sym(zeros(6, 6));
M2(1:3, 1:3) = m2 * eye(3);
M2(4, 4) = I2x;
M2(5, 5) = I2y;
M2(6, 6) = I2z;

disp(M1);
disp(M2);

disp(size(J1));
disp(size(J2));

M = J1' * M1 * J1 + J2' * M2 * J2;
disp(simplify(M));

%% Coriolis matrix
theta = [th1; th2];
dtheta = [dth1; dth2];
C = sym(zeros(2, 2));
for i = 1:2
    for j = 1:2
        for k = 1:2
            C(i, j) = C(i, j) + (diff(M(i, j), theta(k)) + diff(M(i, j), theta(j)) + diff(M(k, j), theta(i))) * dtheta(k);
        end
        
        C(i, j) = 0.5 * C(i, j);
    end
end

% disp(C);

%% Gravitational forces
% g = 9.81;
g1 = twist_angle_to_trans(twist1, th1) * g_s1_init;
g2 = twist_angle_to_trans(twist1, th1) * twist_angle_to_trans(twist2, th2) * g_s2_init;

V = m1 * g * g1(3, 4) + m2 * g * g2(3, 4);

N = sym(zeros(2, 1));
N(1, 1) = diff(V, th1);
N(2, 1) = diff(V, th2);

% disp(simplify(N));

%% Create functions to generate matrices
matlabFunction(M, 'File', 'two_link_inertia_matrix');
matlabFunction(C, 'File', 'two_link_coriolis_matrix');
matlabFunction(N, 'File', 'two_link_gravity_matrix');

end

