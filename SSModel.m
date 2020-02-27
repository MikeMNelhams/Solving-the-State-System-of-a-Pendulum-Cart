% State Space model for the barrel-cart system
% Use:
% http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling
%
M = 30; % Mass of the Cart (kg)
m = 10; % Mass of the pendulum (kg)
b = 0.1; % Drag Coefficient (kgs^(-2))
l = 0.25; % Length of the Pendulum (m)
I = (1/3)*m*(l^2); % Moment of Inertia (kgm^4)
g = 9.8; % Gravitational Acceleration (ms^-2)

% Constants for state space Matrix
const1 = (-(I+m*l^2)*b)/(I*(M+m)+M*m*l^2);
const2 = (m^2*g*l^2)/(I*(M+m)+M*m*l^2);
const3 = (-m*l*b)/(I*(M+m)+M*m*l^2);
const4 = (m*g*l*(M+m))/(I*(M+m)+M*m*l^2);
const5 = (I+m*l^2)/(I*(M+m)+M*m*l^2);
const6 = (m*l)/(I*(M+m)+M*m*l^2);

% Matrix A: state-space matrix
A = [0 1 0 0;
     0 const1 const2 0;
     0 0 0 1;
     0 const3 const4 0];

% Input Matrix (Force)
B = [0;
     const5;
     0;
     const6];

% Output Matrix (x and theta)
C = [5 0 0 0;
     0 0 1 0];

% Feedthrough matrix (always 0)
D = [0;
     0];
 
% State Names
states = {'x' 'x_dot' 'theta' 'theta_dot'};
input = {'F'};
outputs = {'x'; 'theta'};

% Using the ss module to create a state-space
sys_ss = ss(A,B,C,D,'statename',states,'inputname',input,'outputname',outputs)

% Linear Quadratic Regular system
Q = C'*C; % Change C and Q to increase/decrease rate of stabilisation
R = 1; % R always 1
K = lqr(A,B,Q,R); % Linear Quadratic Regulator
A = (A-B*K) 

% Create a system class
sys_cl = ss(A,B,C,D,'statename',states,'inputname',input,'outputname',outputs); % state space - controlled system

% Simulate controlled system
t = 0:0.01:120; % Set time scale for simulation
r =-50*ones(size(t)); % Set distance objective for cart position (negative is positive x-direction)
[y,t,x]=lsim(sys_cl,r,t, [0 0 0 0]); % Simulate model over time t, with obj r and initial states []

% Plot x against t
% yyaxis left
plot(t,y(:,1), 'LineWidth',2)
xlabel('Time, t (s)')
ylabel('Displacement, x (m)')
title('Plotting the displacement of the cart against time (50m), Q_2')

% Plot theta against t
% yyaxis right
plot(t, y(:,2).*(180/pi), 'LineWidth',2)
xlabel('Time, t (s)')
ylabel('Angle, \theta (degrees)')
title('Plotting the angle of the rod against time (50m), Q_2')
