%% Theoretical analysis of configuration



%
%% *** Robot (kinematic) model parameters *** 
clear all; 
close all; 
l(2) = 15.0; %% in cm   
l(4) = 4.0;
l(5) = 8.0;
h(1) = -10.0;
%% *** sampling period *** 
%% *** for the robot motion, kinematic simulation: 
dt = 0.001; %dt = 0.001; i.e. 1 msec)   

%% *** Create (or load from file) reference signals *** 
%% *** DESIRED MOTION PROFILE - TASK SPACE *** 
Tf=10.0; 	% 10sec duration of motion 
t=0:dt:Tf;  

%xd0,td0,yd1: initial/final end-point position --> desired task-space trajectory  
xd0 = 15;	
xd1 = 15; 
yd0 = 0; 
yd1 = 10;  
zd0 = h(1);
zd1 = h(1);
% Example of desired trajectory : linear segment (x0,y0)-->(x1,y1); Time duration: Tf; 
disp('Initialising Desired Task-Space Trajectory (Motion Profile) ...'); %% 
disp(' ');   
xd(1,1) = xd0; 
yd(1,1) = yd0; 
zd(1,1) = zd0;
lambda_x = (xd1-xd0)/Tf; 
lambda_y = (yd1-yd0)/Tf; 
lambda_z = (zd1-zd0)/Tf; 
kmax=Tf/dt + 1; 
for k=2:kmax;    
   xd(k,1) = xd(k-1,1) + lambda_x*dt;    
   yd(k,1) = yd(k-1,1) + lambda_y*dt; 
   zd(k,1) = zd(k-1,1) + lambda_z*dt;
end  
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% ****** KINEMATIC SIMULATION - Main loop ****** 
disp('Kinematic Simulation ...'); %% 
disp(' '); %%  

%% ***** INVESRE KINEMATICS  -->  DESIRED MOTION - JOINT SPACE ***** 
%% compute the reference joint-motion vectors: 
%% {qd(k,i), i=1,...,n (num of degrees of freedom), with k=1,..., kmax,} 
%% and reference joint (angular) velocities {qd_1(k,i)} 
qd(:,1) = 2*atan( (-h - sqrt(h^2 + xd(:).^2 - l(2)^2))./(l(2)+xd(:)) ); 
c1 = cos(real(qd(:,1)));
s1 = sin(real(qd(:,1)));
sub_x(:) = -zd(:).*c1-xd(:).*s1; %% subsystem of the last 2 joints
sub_y(:) = yd(:) - sub_x(:);     
rd2 = sub_x(:).^2 + sub_y(:).^2;
qd(:,3) = acos( (rd2(:)-l(4)^2-l(5)^2)./(2*l(4)*l(5)) );
s3 = sin(real(qd(:,3)));
c3 = cos(real(qd(:,3)));
qd(:,2) = atan( (sub_y(:))./sub_x(:)) - asin((l(5)*s3)./sqrt(rd2(:)));
%% ***** FORWARD KINEMATICS  JOINT MOTION -->  CARTESIAN POSITIONS ***** 

c2 = cos(real(qd(:,2)));
s2 = sin(real(qd(:,2)));
%%(xd1, yd1, zd1) : cartesian position of the 1st link's local reference frame 
xd1 = zeros(10001,1);   
yd1 = zeros(10001,1); 
zd1 = zeros(10001,1);
%%(xd2, yd2, zd2) : cartesian position of the 2nd link's local reference frame 
xd2 = l(2)*c1 ;   
yd2 = zeros(10001,1);  
zd2 = -l(2)*s1;
%%(xd3, yd3, zd3) : cartesian position of the 3rd link's local reference frame 
xd3 = l(2)*c1 - l(4).*c2.*s1;
yd3 = -2^(1/2)*l(4)*sin(qd(:,2) + pi/4);
zd3 = -l(2)*s1-l(4).*c1.*c2;
%%(xde, yde, zde) : cartesian position of the end effector's local reference frame 
xde = l(2)*c1 - l(4).*c2.*s1 - l(5)*c2.*c3.*s1 + l(5)*s1.*s2.*s3;
yde = - l(5)*(c2.*c3 - s2.*s3) - l(5)*(s2.*c3 + s3.*c2) - l(4)*c2 - l(4)*s2;
zde = l(5)*c1.*s2.*s3 - l(2)*s1 - l(4)*c1.*c2 - l(5)*c1.*c2.*c3 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%% *** SAVE and PLOT output data *** %%** use functions plot(...)  
save C:\Users\User\Desktop\Ρομποτική\results.mat;  %% --> save data to 'matlab.mat' file   
%% End effector position
fig1 = figure;  
subplot(4,4,1); 
plot(t,xd); 
ylabel('px (cm)'); 
xlabel('time t (sec)');  

subplot(4,4,2); 
plot(t,yd); 
ylabel('py (cm)'); 
xlabel('time t (sec)');  

subplot(4,4,3); 
plot(t,zd); 
ylabel('pz (cm)'); 
xlabel('time t (sec)');  
%% End effector velocity
subplot(4,4,4); 
plot(t,[0;diff(xde)]); 
ylabel('Vxe (m/s)'); 
xlabel('time t (sec)');    

subplot(4,4,5); 
plot(t,[0;diff(yde)]); 
ylabel('Vye (m/s)'); 
xlabel('time t (sec)');   

subplot(4,4,6); 
plot(t,[0;diff(zde)]); 
ylabel('Vze (m/s)'); 
xlabel('time t (sec)');   
%% Joint angles
subplot(4,4,7); 
plot(t,qd(:,1)); 
ylabel('q1 (rad)'); 
xlabel('time t (sec)');  

subplot(4,4,8); 
plot(t,qd(:,2)); 
ylabel('q2 (rad)'); 
xlabel('time t (sec)');    

subplot(4,4,9); 
plot(t,qd(:,3)); 
ylabel('q3 (rad)'); 
xlabel('time t (sec)');    

%% Velocities of joints
subplot(4,4,10); 
plot(t,[0;diff(qd(:,1))]); 
ylabel('v1 (rad/s)'); 
xlabel('time t (sec)');    

subplot(4,4,11); 
plot(t,[0;diff(qd(:,2))]); 
ylabel('v2 (rad/s)'); 
xlabel('time t (sec)');    

subplot(4,4,12); 
plot(t,[0;diff(qd(:,3))]); 
ylabel('v3 (rad/s)'); 
xlabel('time t (sec)');    





%%*** stick diagram --> animate robot motion ... (**optional**) 
%% within a for (or while) loop, use periodic plot(...) functions to draw the geometry (current pos)  
%% of the robot, and thus animate its motion ...  

fig2 = figure; 
axis([-15 15 -15 15 -15 15 ]) %%set xy plot axes (caution: square axes, i.e. dx=dy) 
axis on 
hold on 
xlabel('x (cm)'); 
ylabel('y (cm)'); 
zlabel ('z (cm)')
dtk=1000; %% plot robot position every dtk samples, to animate its motion 
plot3([0],[0],[0],'o'); 
for tk=1:dtk:kmax,    %%% 	
   pause(0.1);	%% pause motion to view successive robot configurations    
   plot3([0,xd1(tk)],[0,yd1(tk)],[0,zd1(tk)]);					
   plot3([xd1(tk)],[yd1(tk)],[zd1(tk)],'o');    
   plot3([xd1(tk),xd2(tk)],[yd1(tk),yd2(tk)],[zd1(tk),zd2(tk)]);	
   plot3([xd2(tk)],[yd2(tk)],[zd2(tk)],'.');    
   plot3([xd2(tk),xd3(tk)],[yd2(tk),yd3(tk)],[zd2(tk),zd3(tk)]);
   plot3([xd3(tk)],[yd3(tk)],[zd3(tk)],'+');    
   plot3([xd3(tk),xde(tk)],[yd3(tk),yde(tk)],[zd3(tk),zde(tk)]);
   plot3([xde(tk)],[yde(tk)],[zde(tk)],'*'); 
end       