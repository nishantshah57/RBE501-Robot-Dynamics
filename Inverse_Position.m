clear all;
clc;
close all;
syms t1 t2 t3 t4 t5 t6;
t2 = t2 + pi/2;
%% DH Parameters
DH = [t1 475 150   pi/2
      t2 0   600   0
      t3 720 120   pi/2
      t4 0   0     -pi/2
      t5 0   0     pi/2
      t6 85  0     0];

%% Transformation Matrix 
H = eye(4);  
for i = 1:6
    
    A{i} = [cos(DH(i,1)), -sin(DH(i,1))*cos(DH(i,4)), sin(DH(i,1))*sin(DH(i,4)), DH(i,3)*cos(DH(i,1));
                sin(DH(i,1)), cos(DH(i,1))*cos(DH(i,4)), -cos(DH(i,1))*sin(DH(i,4)), DH(i,3)*sin(DH(i,1));
                0, sin(DH(i,4)), cos(DH(i,4)), DH(i,2);
                0, 0, 0, 1];
            
    H = H*A{i};
end
EE_Position = H(1:3,4);
EE_Orientation = H(1:3,1:3);
F = A{1}*A{2}*A{3};
T01 = A{1};
T02 = A{1}*A{2};

%% Question A. 
%3D Workspace of the Robot 
plotMatrix = [];
index = 0;
for q1 = -pi/2:0.25:pi/2
    for q2 = -pi/3:0.25:pi/3
        for q3 = -pi/3:0.25:pi/3
            index = index + 1;
            x = subs(F(1,4),{t1,t2,t3},{q1,q2,q3});
            y = subs(F(2,4),{t1,t2,t3},{q1,q2,q3});
            z = subs(F(3,4),{t1,t2,t3},{q1,q2,q3});
            plotMatrix(index,:) = [x, y, z];
        end
    end
end

plot3(plotMatrix(:,1),plotMatrix(:,2),plotMatrix(:,3),'*');

%% Question A
%2D Side View/Workspace of the robot
% Workspace can also be derived by circle formulas such as 
%(x-h)^2 + (y-k)^2 = r^2
% where r is the radius or the link length can be taken.
% h,k is the centre or origin of the circle or the arc to be derived.

figure;
hold on
th1 = pi/2:-2*pi/300:-pi/6;
xunit1 = 1320*cos(th1) +0;
yunit1 = 1320*sin(th1) +475;

th2 = pi/2:-2*pi/300:-pi/6;
xunit2 = 720*cos(th2) + 0;
yunit2 = 720*sin(th2) + 1075;

th3 = 20*pi/180:-2*pi/300:-100*pi/180;
xunit3 = sqrt(600^2 + 720^2 - 2*600*720*cos(pi/3))*cos(th3) + 0;
yunit3 = sqrt(600^2 + 720^2 - 2*600*720*cos(pi/3))*sin(th3) + 475;

th4 = 210*pi/180:2*pi/300:330*pi/180;
xunit4 = 720*cos(th4) + 600*sin(pi/3);
yunit4 = 720*sin(th4) + 475-600*sin(pi/6);
plot(xunit1,yunit1,xunit2,yunit2,xunit3,yunit3,xunit4,yunit4);
hold off

%% Question D. and E.
% Jacobian Matrix derived and Joint Velocities obtained
syms t1 t2 t3 real;
Jv = simplify(jacobian([F(1,4),F(2,4),F(3,4)],[t1,t2,t3]))
k = [0 0 1]';
Jw3 = simplify(T02(1:3,1:3)*k);
Jw2 = simplify(T01(1:3,1:3)*k);
Jw1 = [0 0 1]';
 
J = simplify([Jv ; Jw1 Jw2 Jw3])
J1 = [ 720*cos(t1) - 150*sin(t1) + 600*sin(t1)*sin(t2) + 120*cos(t2)*sin(t1)*sin(t3) + 120*cos(t3)*sin(t1)*sin(t2), -120*cos(t1)*(cos(t2 + t3) + 5*cos(t2)), -120*cos(t2 + t3)*cos(t1);
150*cos(t1) + 720*sin(t1) - 600*cos(t1)*sin(t2) - 120*cos(t1)*cos(t2)*sin(t3) - 120*cos(t1)*cos(t3)*sin(t2), -120*sin(t1)*(cos(t2 + t3) + 5*cos(t2)), -120*cos(t2 + t3)*sin(t1);
                                                                                                           0,        - 120*sin(t2 + t3) - 600*sin(t2),         -120*sin(t2 + t3);
                                                                                                           0,                                 sin(t1),                   sin(t1);
                                                                                                           0,                                -cos(t1),                  -cos(t1);
                                                                                                           1,                                       0,                         0];
J_p = pinv(J1);
x_vel = [5 5 10 0 0 0]';
Joint_velocities = vpa(subs(J_p*x_vel,[t1,t2,t3],[0.1974,1.3353,1.4056]))
%% Question B. and C, 
%Transformation matrix G is given and Joint angles are obtained
% Theta values are obtained in Radians

syms r11 r12 r13 r14 r21 r22 r23 r24 r31 r32 r33 r34 r41 r42 r43 r44 real;
D = [r11 r12 r13 r14; r21 r22 r23 r24; r31 r32 r33 r34; r41 r42 r43 r44];
N = simplify(inv(F)*D);

G = [1 0 0 500;
     0 1 0 100;
     0 0 1 1500;
     0 0 0 1];

x1 = G(1,4);
y1 = G(2,4);

t1 = vpa(atan2(y1,x1));

r = sqrt(G(1,4)^2 + G(2,4)^2);
ax = sqrt(DH(3,2)^2 + DH(3,3)^2);
s3 = double((r-DH(1,3))^2 + G(3,4)^2 - DH(2,3)^2 - ax^2);
c3 = real(double(sqrt(1-s3^2)));
t3 = atan2(s3,c3) - double(atan(DH(3,3),DH(3,2)));

s2 = -double(((DH(2,3)+ax*s3)*(r-DH(1,3)) + ax*c3*G(3,4)));
c2 = double((DH(2,3)+ax*s3)*G(3,4) + ax*c3*(r-DH(1,3)));
t2 = atan2(s2,c2) + pi/2;

M  = subs(N,{r11,r12,r13,r14,r21,r22,r23,r24,r31,r32,r33,r34,r41,r42,r43,r44,t2,t3},{G(1,1),G(1,2),G(1,3),G(1,4),G(2,1),G(2,2),G(2,3),G(2,4),G(3,1),G(3,2),G(3,3),G(3,4),G(4,1),G(4,2),G(4,3),G(4,4),th2,th3});
t4 = 0;
% because it is aligned on the frame 3 and no rotation has happened

s6 = sin(t1 + t2 + t3)/2 - sin(t2 - t1 + t3)/2;
c6 = -(cos(t1 + t2 + t3)/2 + cos(t2 - t1 + t3)/2);
t6 = vpa(atan2(s6,c6));

c5 =  sin(t2 + t3);
s5 = sqrt(1-c5^2);
t5 = atan2(s5,c5);

theta1 = t1
theta2 = t2
theta3 = t3
theta4 = t4
theta5 = t5
theta6 = t6

% Question F. is derived and explained on paper and uploaded as a PDF
