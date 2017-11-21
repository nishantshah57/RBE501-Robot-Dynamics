syms l1 l2 l3 q1 q2 q3 m1 m2 m3 dq1 dq2 dq3 ddq1 ddq2 ddq3 g real;
%% Position and Velocity
% Position is found using DH-Parameters for the 3-link massless arm
theta = [q1 q2 q3]';
alpha = [sym(pi/2) 0 0]';
a = [0 l2 l3]';
d = [l1 0 0]';

DH_parameters = [theta alpha a d];

T1 = [cos(q1)    0      sin(q1)    0;
      sin(q1)    0     -cos(q1)    0;
      0          1            0   l1;
      0          0            0   1];
  
T2 = [cos(q2)    -sin(q2)    0    l2*cos(q2);
      sin(q2)    cos(q2)     0    l2*sin(q2);
      0          0           1    0;
      0          0           0    1];
  
T3 = [cos(q3)    -sin(q3)    0    l3*cos(q3);
      sin(q3)    cos(q3)     0    l3*sin(q3);
      0          0           1    0;
      0          0           0    1];
 
T01 = T1;
T02 = T01*T2;
T03 = simplify(T02*T3);

v1 = jacobian(T01(1:3,4),[q1])* [dq1];
v2 = jacobian(T02(1:3,4),[q1, q2]) * [dq1; dq2];
v3 = jacobian(T03(1:3,4),[q1, q2, q3]) * [dq1; dq2; dq3];
%% Kinematic and Potential Energy
K1 = 0.5 * m1 * (v1.' * v1);
K2 = 0.5 * m2 * (v2.' * v2);
K3 = 0.5 * m3 * (v3.' * v3);

P1 = m1 * g * T01(3,4);
P2 = m2 * g * T02(3,4);
P3 = m3 * g * T03(3,4); % z-values taken as height

K = K1 + K2 + K3;
P = P1 + P2 + P3;
%% Lagrangian Equation
L = simplify(K - P);

syms t1(t) t2(t) t3(t) real;
%% Euler-Lagrange Equation
d_L_K1 = diff(L,dq1);
L_K1_t = subs(d_L_K1, [q1 q2 q3 dq1 dq2 dq3], [t1 t2 t3 diff(t1(t),t) diff(t2(t),t) diff(t3(t),t)]);
d_L_K1_t = diff(L_K1_t, t);
d_L_K1 = subs(d_L_K1_t, [t1 t2 t3 diff(t1(t),t) diff(t2(t),t) diff(t3(t),t)...
    diff(t1(t),t,t) diff(t2(t),t,t) diff(t3(t),t,t)], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

d_L_K2 = diff(L,dq2);
L_K2_t = subs(d_L_K2, [q1 q2 q3 dq1 dq2 dq3], [t1 t2 t3 diff(t1(t),t) diff(t2(t),t) diff(t3(t),t)]);
d_L_K2_t = diff(L_K2_t,t);
d_L_K2 = subs(d_L_K2_t, [t1 t2 t3 diff(t1(t),t) diff(t2(t),t) diff(t3(t),t)...
    diff(t1(t),t,t) diff(t2(t),t,t) diff(t3(t),t,t)], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

d_L_K3 = diff(L,dq3);
L_K3_t = subs(d_L_K3, [q1 q2 q3 dq1 dq2 dq3], [t1 t2 t3 diff(t1(t),t) diff(t2(t),t) diff(t3(t),t)]);
d_L_K3_t = diff(L_K3_t,t);
d_L_K3 = subs(d_L_K3_t, [t1 t2 t3 diff(t1(t),t) diff(t2(t),t) diff(t3(t),t)...
    diff(t1(t),t,t) diff(t2(t),t,t) diff(t3(t),t,t)], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

d_L_P1 = diff(L,q1);
d_L_P2 = diff(L,q2);
d_L_P3 = diff(L,q3);
 
Tau_1 = simplify(d_L_K1 - d_L_P1);
Tau_2 = simplify(d_L_K2 - d_L_P2);
Tau_3 = simplify(d_L_K3 - d_L_P3);

fprintf("\nThe dynamic equation of the 3 link arm robot is: \n\n");
Tau = [Tau_1; Tau_2; Tau_3];
pretty(Tau);

%% Dynamical Model
M11 = simplify(Tau_1 - subs(Tau_1, ddq1, 0)) / ddq1;
M12 = simplify(Tau_1 - subs(Tau_1, ddq2, 0)) / ddq2;
M13 = simplify(Tau_1 - subs(Tau_1, ddq3, 0)) / ddq3;

M21 = simplify(Tau_2 - subs(Tau_2, ddq1, 0)) / ddq1;
M22 = simplify(Tau_2 - subs(Tau_2, ddq2, 0)) / ddq2;
M23 = simplify(Tau_2 - subs(Tau_2, ddq3, 0)) / ddq3;

M31 = simplify(Tau_3 - subs(Tau_3, ddq1, 0)) / ddq1;
M32 = simplify(Tau_3 - subs(Tau_3, ddq2, 0)) / ddq2;
M33 = simplify(Tau_3 - subs(Tau_3, ddq3, 0)) / ddq3;

M = [M11 M12 M13;
    M21 M22 M23;
    M31 M32 M33];

G11 = subs(Tau_1, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);
G21 = subs(Tau_2, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);
G31 = subs(Tau_3, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);

G = [G11; 
    G21;
    G31];

C = ([Tau_1;Tau_2;Tau_3]) - (M*[ddq1;ddq2;ddq3] - G);

fprintf('\nThe inertia matrix M is: \n\n ');
pretty(M);
fprintf('\nThe centripetal vector C is: \n\n ');
pretty(C);
fprintf('\nThe gravity vector G is: \n\n ')
pretty(G);
