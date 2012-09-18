% Exploding Projectile Simulation
%
% Simulates the flight of a projectile.  At the top of its flight, the projectile (mt)
% splits into two smaller projectiles (m1 & m2).  Each of the projectiles is given an
% additional instantaneous velocity at this point.  The simulation is implemented in two
% phases using the same set of differential equations with different stopping conditions.
% Phase one covers the flight of the combined projectile and ends when the projectile is
% at the top of its flight.  After phase one ends, the two smaller projectiles are given
% the additional velocities and simulated, ending when they hit the ground.
%
% Function Inputs: None
% Function Outputs: None
% The function will plot trajectories of the three projectiles over time.
%
% Initial conditions and parameters are specified in the function for easy running.
%
% Nick Eyre - Sept 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 1 - Problem 1


function projectile

%Define System Parameters
mt = 5;     %kg
v0 = [7 2]; %m/s
v1 = [0 0]; %m/s
v2 = [0 -2];%m/s
g = -9.8;   %m/s^2

%Simulation Parameters
Tmax = 100; %s
dt = .02;   %s

%Calculated Parameters
m1 = .8 * mt;
m2 = .2 * mt;
g = [0; g];  %Gravity (m/s^2)

%State Variable Definitions
% Z=[R, V] (postition & its derivatives are vectors)

%Phase 1: Run Simulation for combined projectile
Z0 = [[0 0],v0];  %Initial Conditions
T0 = [0:dt:Tmax];  %Max Time Vector
options = odeset('Events', @phase_one_events);  %Setup Events
[T1 Z1] = ode45(@equations, T0, Z0, options);   %Run Solver

%Phase 2: Run Simulations for each of the smaller projectiles
T0 = [T1(end):dt:Tmax];    %Max Time Vector
M1_Z0 = [Z1(end,1:2) Z1(end,3:4)+v1];   %Initial Conditions for M1 = Output of Phase 1 with Added Velocity
M2_Z0 = [Z1(end,1:2) Z1(end,3:4)+v2];   %Initial Conditions for M2 = Output of Phase 1 with Added Velocity
options = odeset('Events', @phase_two_events);
[M1_T2 M1_Z2] = ode45(@equations, T0, M1_Z0, options);   %Simulate M1
[M2_T2 M2_Z2] = ode45(@equations, T0, M2_Z0, options);   %Simulate M2

%Plot Results
clf,hold all
plot(Z1(:,1),Z1(:,2))             %Plot MT X vs Y
plot(M1_Z2(:,1),M1_Z2(:,2),':')   %Plot M1 X vs Y
plot(M2_Z2(:,1),M2_Z2(:,2),'-.')  %Plot M2 X vs Y
legend('Combined','Mass 1','Mass 2')
title('Projectile Trajectory')
xlabel('X Position (m)')
ylabel('Y Position (m)')

  % Function Defintions
  function derivatives = equations(T, X)
    %Extract
    R = X(1:2);
    V = X(3:4);

    %Calculate Derivatives
    A = g;

    %Return Derivatives
    derivatives = [V; A];
  end

  %Phase One (One Projectile) Events Function
  %Stop when at top of trajectory (Vy = 0)
  function [eventvalue,stopthecalc,eventdir] = phase_one_events(T,X)
    eventvalue = X(4);
    stopthecalc = 1;
    eventdir = -1;
  end

  %Phase Two (Two Projectiles) Events Function
  %Stop when hits ground
  function [eventvalue,stopthecalc,eventdir] = phase_two_events(T,X)
    eventvalue = X(2);
    stopthecalc = 1;
    eventdir = -1;
  end

end %projectile