clear all;
close all;
clc;

%% DH-Parameters Table
% Here the general values of theta are taken in the form of t.
syms t1 t2 t3 t4 t5 t6 real;

% Matlab considers all values in radian, so in order to convert the final
% answer I have multiplied it by pi/180, so that the sin and cos functions
% take in radian values equivalent to degree values.

t = [t1 (t2+90)*(pi/180) t3 t4 t5 t6]'; % Joint angle values
Theta = [0 90 0 0 0 0]'; %Just for the Table
d = [330 0 0 320 0 80]'; % Joint Distance values
a = [75 300 75 0 0 0]'; % Link Length values
ap = [pi/2 0 pi/2 -pi/2 pi/2 0]'; % Twist angle Values
alpha = rad2deg(ap);
DH_Parameter = table(Theta,d,a,alpha)

%% Homogenous Transformation Matrix
 for i = 1:6
    A(:,:,i) = Homogenous(t(i), d(i), a(i), alpha(i));
 end
 
T01 = A(:,:,1)  
T12 = A(:,:,2)
T23 = A(:,:,3)
T34 = A(:,:,4)
T45 = A(:,:,5)
T56 = A(:,:,6)

%% Composite Transformation as a series of Homogenous Transformation
T06 = simplify(T01*T12*T23*T34*T45*T56)

%% Substituting the values of theta 
T06_Given = double(simplify(subs(T06,[t1,t2,t3,t4,t5,t6],(pi/180)*[0,75,30,135,-45,60])));

% Position Forward Kinematics
T06_Given_Position_Kinematics = T06_Given(1:3,4)

%Orientation Forward Kinematics
T06_Given_Orientation_Kinematics = T06_Given(1:3,1:3)

%% Home Pose of the Robot
T06_Home = simplify(subs(T06,[t1,t2,t3,t4,t5,t6],[0,0,0,0,0,0]))

% Position Forward Kinematics For home Position
T06_Home_Postion_Kinematics = T06_Home(1:3,4)

% Orientation Forward Kinematics For home Position
T06_Home_Orientation_Kinematics = T06_Home(1:3,1:3)

%% Plotting of the Fanuc Robot 
TGiven = plotLink(0,75,30,135,-45,60);
THome = plotLink(0,0,0,0,0,0);
function Y = plotLink(th1, th2, th3, th4, th5, th6)
t = [th1*(pi/180), (th2+90)*(pi/180), th3*(pi/180), th4*(pi/180), th5*(pi/180), th6*(pi/180)]; %Input Theta values for plotting
d = [330 0 0 320 0 80]'; % Joint Distance values
a = [75 300 75 0 0 0]'; % Link Length values
ap = [pi/2 0 pi/2 -pi/2 pi/2 0]'; % Twist angle Values
alpha = rad2deg(ap); % Converted all alpha values from radian to degrees

 for i = 1:6
    A(:,:,i) = Homogenous(t(i), d(i), a(i), alpha(i));
 end
 
T(:,:,1) = A(:,:,1); 
for j=2:6   
    T(:,:,j) = T(:,:,j-1)*A(:,:,j); 
end

xPoint = zeros(1,6);
yPoint = zeros(1,6);
zPoint = zeros(1,6);
for j = 1:6
    xPoint(j+1) = T(1,4,j);
    yPoint(j+1) = T(2,4,j);
    zPoint(j+1) = T(3,4,j);
end

figure('Name','3D Plot of Fanuc Robot','NumberTitle','off') 
plot3(xPoint, yPoint, zPoint,'b','LineWidth',2); 
hold on 
scatter3(xPoint,yPoint,zPoint,100,'r','filled'); 
grid on 
title('3D Plot of Fanuc LR Mate 200 IC Robot Manipulator','FontSize',10) 
xlabel('x-axis') 
ylabel('y-axis') 
zlabel('z-axis') 
hold off
Y = T(:,:,6);
end
%% Function for calculating the Homogenous Matrix of the Robot
function Y = Homogenous(t,d,a,alpha)
Rot_z_theta = [cos(t) -sin(t) 0 0
     sin(t) cos(t) 0 0
     0 0 1 0
     0 0 0 1];
Trans_z_d = [1 0 0 0
     0 1 0 0
     0 0 1 d
     0 0 0 1];
Trans_x_a = [1 0 0 a
     0 1 0 0
     0 0 1 0
     0 0 0 1];
Rot_x_alpha = [1 0 0 0
     0 cosd(alpha) -sind(alpha) 0
     0 sind(alpha) cosd(alpha) 0
     0 0 0 1];
 
Y = Rot_z_theta * Trans_z_d * Trans_x_a * Rot_x_alpha;
end
 