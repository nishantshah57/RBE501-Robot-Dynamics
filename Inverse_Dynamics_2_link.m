syms a1 a2 m1 m2 I1 I2 q1(t) q2(t) dq1 dq2 ddq1 ddq2 g real;

%% Position and Velocity 
% The position and velocity of the 2-link arm having masses at end of each
% link

x1 = 0.5*a1*cos(q1);
y1 = 0.5*a1*sin(q1);

x2 = a1*cos(q1) + 0.5*a2*cos(q1 + q2);
y2 = a1*sin(q1) + 0.5*a2*sin(q1 + q2);

v1_x = diff(x1,t);
v1_y = diff(y1,t);

v2_x = diff(x2,t);
v2_y = diff(y2,t);

v1 = [v1_x; v1_y];
v2 = [v2_x; v2_y];

%% Kinematic and Potential Energy

K1 = simplify((0.5*m1*(v1.' * v1)) + (0.5*I1*(diff(q1(t),t)^2)));
K2 = simplify((0.5*m2*(v2.' * v2)) + (0.5*I2*(diff(q1(t),t) + diff(q2(t),t))^2));

P1 = 0.5*m1*g*a1*sin(q1);
P2 = m2*g*a1*sin(q1) + 0.5*g*a2*m2*sin(q1 + q2);

%% Lagrangian Equation

L = simplify((K1 + K2) - (P1 + P2));

syms th1 th2 dth1 dth2 ddth1 ddth2 real;
%% Euler-Lagrange Equation

A1 = subs(L, diff(q1(t),t), dth1);
A1t = diff(A1,dth1);
temp_A1t = subs(A1t,[dth1],[diff(q1(t),t)]);
dA1t = diff(temp_A1t,t);
A1 = subs(dA1t, [q1 q2 diff(q1(t),t) diff(q2(t),t) diff(q1(t),t,t) diff(q2(t),t,t)],...
    [th1 th2 dth1 dth2 ddth1 ddth2]);

A2 = subs(L, diff(q2(t),t), dth2);
A2t = diff(A2,dth2);
temp_A2t = subs(A2t,[dth2],[diff(q2(t),t)]);
dA2t = diff(temp_A2t,t);
A2 = subs(dA2t, [q1 q2 diff(q1(t),t) diff(q2(t),t) diff(q1(t),t,t) diff(q2(t),t,t)],...
    [th1 th2 dth1 dth2 ddth1 ddth2]);

B1_temp = subs(L,[q1(t)],[th1]);
dB1 = diff(B1_temp, th1);
B1 = subs(dB1, [q1 q2], [th1 th2]);

B2_temp = subs(L, q2(t), th2);
dB2 = diff(B2_temp, th2);
B2 = subs(dB2, [q1 q2], [th1 th2]);

Tau_1 = simplify(A1 - B1);
Tau_2 = simplify(A2 - B2);

fprintf('\nThe dynamics of the robot is: T = \n\n');
Tau = [Tau_1; Tau_2];
pretty(Tau);
%% Dynamical Model 
M11 = simplify(Tau_1 - (subs(Tau_1, [ddth1], [0]))) / ddth1;
M12 = simplify(Tau_1 - (subs(Tau_1, [ddth2], [0]))) / ddth2;

M21 = simplify(Tau_2 - (subs(Tau_2, ddth1, 0))) / ddth1;
M22 = simplify(Tau_2 - (subs(Tau_2, ddth2, 0))) / ddth2;

M = [M11 M12;
    M21 M22];

G = simplify(subs(Tau, [dth1 dth2 ddth1 ddth2], [0 0 0 0]));

C = simplify(Tau - (M(t)*[ddth1;ddth2]) - G);

fprintf('\nThe inertia matrix M is: \n\n ');
pretty(M);
fprintf('\nThe centripetal vector C is: \n\n ');
pretty(C);
fprintf('\nThe gravity vector G is: \n\n ')
pretty(G);
