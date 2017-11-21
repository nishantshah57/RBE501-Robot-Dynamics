syms theta1 theta2 theta3 theta4 l1 l2 l3 

%% Question 3
% DH Parameters
DH = [[theta1 0  l1 pi/2]
      [theta2 0  l2 0]
      [theta3 0  0  pi/2]
      [theta4 l3 0 0]];

%% Question 4 
% Transformation Matrix
% Symbolic Frame Transformations
s_A1 = simplify(FrameTransform(DH(1,:)));
s_A2 = simplify(FrameTransform(DH(2,:)));
s_A3 = simplify(FrameTransform(DH(3,:)));
s_A4 = simplify(FrameTransform(DH(4,:)));

s_A01 = simplify(s_A1)
s_A02 = simplify(s_A01 * s_A2)
s_A03 = simplify(s_A02 * s_A3)
s_A04 = simplify(s_A03 * s_A4)

DH = subs(DH, [theta1 theta2 theta3 theta4 l1 l2 l3], [0 0 0 0 0 70 100]);
%% Question 5 
% Forward Kinematics when theta values are zero
% Frame Transformation at theta 0,0,0,0
A1 = FrameTransform(DH(1,:));
A2 = FrameTransform(DH(2,:));
A3 = FrameTransform(DH(3,:));
A4 = FrameTransform(DH(4,:));

A01 = A1;
A02 = A01 * A2;
A03 = A02 * A3;
A04 = vpa(A03 * A4)

Leg_Position = A04(1:3,4)

%% Question 6
v = [0; 0; 10; 1];
Frame_b = A04 * v

%% Question 7
% Inverse Kinematics
Px = 80;
Py = 0;
Pz = -100;
r = sqrt(Px^2 + Py^2);

theta_1 = atan2(Py,Px)
theta_1_1 = pi
% Here theta_1 value is zero,we can consider the value of theta_1 as "pi"
% as well and solve the further qequations

q = sqrt((r^2) + Pz^2);

beta = vpa(acos((70^2 + q^2 - 100^2)/(2*70*q)));
alpha = vpa(atan2(-100,80));

theta_2 = vpa(alpha + beta) % First value of theta_2
theta_2_2 = vpa(pi - alpha - beta) % Second value of theta_2


gamma = vpa(acos((70^2 + 100^2 - q^2)/(2*70*100)) );
theta_3 = gamma % Goes with first value of theta_2 
theta_3_3 = vpa(gamma - pi) %Goes with second value of theta_2_2

%% Question 8
%Symbolic Jacobian
Jv = jacobian([s_A04(1,4),s_A04(2,4),s_A04(3,4)],[theta1, theta2, theta3, theta4]);

k = [0; 0; 1];
Jw = [k, A01(1:3,1:3)*k, A02(1:3,1:3)*k, A03(1:3,1:3)*k];
J = [Jv; Jw]

%% Question 10
% Joint Velocities leading to foot linear velocity for all theta values equal to zero 
x_dot = [0; 0; 10; 0; 0; 0];
Jm = (subs(J, [theta1 theta2 theta3 theta4 l1 l2 l3], [0 0 0 0 0 70 100]));
theta_dot = pinv(Jm) * x_dot;

%% Question 9
% Singularities of leg 
Jb = (subs(J, [theta1 theta2 theta3 theta4 l1 l2 l3], [pi/2 pi/2 pi/2 0 0 70 100]));
Rank_Jb = rank(Jb); % Rank of J matrix when all theta values are zero
Rank_Jm = rank(Jm); % Rank of J Matrix when all theta values are pi/2
Rank_J = rank(J); % maximum rank of Jacobian Matrix

% For zero and pi/2 theta values, the rank is less than maximum rank of J
% matrix. Therefore singular configuration for zero and pi/2 theta values.
% If there are any theta values for which the rank goes less than the
% maximum rank of Jacobian Matrix, then for that set of theta value the
% configuration will be in Sigular condition.

%% Frame Transfromation Function
function A = FrameTransform(y)
u = y(1); d = y(2); a = y(3); v = y(4);

A = [cos(u), -sin(u)*cos(v), sin(u)*sin(v), a*cos(u);
     sin(u), cos(u)*cos(v), -cos(u)*sin(v), a*sin(v);
     0, sin(v), cos(v), d;
     0, 0, 0, 1];
end 
