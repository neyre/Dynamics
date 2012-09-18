% Golf Ball Simulation
%
% Simulates the flight of a golf ball flying in wind using the ODE45 differential
% equation solver.  Incorporates air resistance with static wind blowing the ball.
%
% Function Inputs: None
% Function Outputs: None
% The function will plot the trajectory of the golf ball.
%
% Initial conditions and parameters are specified in the function for easy running.
%
% Nick Eyre - Sept 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 1 - Problem 2


function golfball
clear all

%System Parameters
W = [0 -15 0]; %Wind, mph
v = 150;      %Launch Speed, mph
a = 35;       %Angle from Horizontal, degrees
m = .0459;    %kg
Cd = .25;     %unitless
D = .04115;   %m
p = 1.02;     %kg/m^2
g = -9.8;     %m/s^2

%Simulation Parameters
Tmax = 100;   %seconds
dt = .02;     %seconds

%Calculated Parameters
c = pi*p*Cd*D^2/8/m;  %Combined Drag Constants
G = [0;0;g];  %Gravity (m/s^2)
W = W*.447;   %Convert to m/s
v = v*.447; %convert to m/s
V = [v*sind(a) 0 v*cosd(a)];

%State Variable Definitions
% Z=[R, V] (postition & its derivatives are vectors)

%Run Simulation
Z0 = [[0 0 0],V];  %Initial Conditions
T0 = [0:dt:Tmax];  %Max Time Vector
options = odeset('Events', @events);  %Setup Events
[T Z] = ode45(@equations, T0, Z0, options);  %Run Solver

%Plot Results
subplot(3,1,1)
  plot3(Z(:,1),Z(:,2),Z(:,3)) %3D Plot
  title('3D Trajectory')
  grid on
subplot(3,1,2)
  plot(Z(:,1),Z(:,3)) %2D Plot X vs Z
  title('Trajectory Right View'); xlabel('X Position (m)'); ylabel('Z Position (m)');
subplot(3,1,3)
  plot(Z(:,2),Z(:,3)) %2D Plot Y vs Z
  title('Trajectory Left View'); xlabel('Y Position (m)'); ylabel('Z Position (m)');

  %Function Definition
  function derivatives = equations(T, X)
    %Extract
    R = X(1:3);
    V = X(4:6);

    %Calculate Relative Velocity & Drag Force
    Vr = V-W';
    Fd = -c*norm(Vr)*Vr;

    %Calculate Acceleration
    A = G + Fd;

    %Return State Vector
    derivatives = [V;A];
  end

  %Events Function
  function [eventvalue,stopthecalc,eventdir] = events(T,X)
    %Stop when it hits ground
    eventvalue = X(3);
    stopthecalc = 1;
    eventdir = -1;
  end

end %projectile