% Hopper Simulation 2D
%
% Simulates the flight of a hopper built with two masses separated by a
% spring.  Uses two phases (hopper on ground, hopper off ground) with the
% ODE45 differential equation solver.  This simulation runs in two dimensions,
% using vector calculations in MATLAB.  The top mass of the hopper is 
% given an initial displacement downward and to the side before
% launching along a trajectory.
%
% Function Inputs: None
% Function Outputs: None
% The function will plot the trajectory of the hopper in the X-Y plane.
%
% Initial conditions and parameters are specified in the function for easy running.
%
% Nick Eyre - Sept 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 1 - Problem 4

function hopper_2d
clear all

%System Parameters
g = 9.81;      % gravitational acceleration in m/s^2
ml = .01;      % mass in kg
mu = .02;      % mass in kg
L0 = .05;       % rest length of spring
k = 9810;      % spring constant (N/m)
dx = .005;         %initial horiz displacement
dy = .01;       % initial vertical displacement

%Simulation Parameters
plotfriction = false;
Tmax = 25;

%State Variable Definitions (quantities are vectors of length 2)
% Z = [Ru Vu Rl Vl]

%Setup Friction Vector for recording friction over time
global Friction; Friction = [];

%PHASE ONE SIMULATION
Z0 = [[dx L0-dy],[0 0],[0 0],[0 0]];  %Hopper Initial Conditions
t_span = [0, Tmax];  %Max Time Span
options = odeset('Events', @phase_one_events); %Setup Events
[T1, Z1] = ode45(@phase_one_equations, t_span, Z0, options);  %Run Solver

%PHASE TWO SIMULATION
Z0 = Z1(end,1:8); %Initial Conditions are Phase 1 Output
t_span = [T1(end), Tmax];  %Max Time Span
options = odeset('Events', @phase_two_events); %Setup Events
[T2, Z2] = ode45(@phase_two_equations, t_span, Z0, options);  %Run Solver

%Combine Results & Calculate Center of Mass
Z = [Z1;Z2]; T=[T1;T2];
CoM = [(mu*Z(:,1)+ml*Z(:,5)),(mu*Z(:,2)+ml*Z(:,6))]./(mu+ml);

%Plot Results
clf,hold on
plot(Z(:,1),Z(:,2),'b')
plot(Z(:,5),Z(:,6),'r')
plot(CoM(:,1),CoM(:,2),'--k')
xlabel('X Distance (m)')
ylabel('Height (m)')
legend('Upper Mass','Lower Mass','Center of Mass')
title('Hopper Trajectory')

%Plot Friction
if plotfriction
    clf
    plot(Friction(:,1),Friction(:,2))
    xlabel('Time (s)')
    ylabel('Force (N)')
    title('Frictional Force during Launch')
end

    %Phase One (Hopper on Ground) Function Definition
    function derivatives = phase_one_equations(T, X)
        %Extract
        Ru = X(1:2);
        Vu = X(3:4);
        Rl = X(5:6);
        Vl = X(7:8);

        %Lower Mass Doesn't Move b/c it is on ground
        Al = [0;0];

        %Calculate Spring Focre
        Fs = k * (L0 - norm(Ru-Rl)) * (Ru-Rl)./norm(Ru-Rl);

        %Calculate Motion of Upper Mass
        Au = [0;-g] + Fs/mu;

        %Calculate Minimum Friction Force & add to Friction Vector
        Ff = Fs(1);
        Friction = [Friction;[T Ff]];

        %Return State Vector
        derivatives = [Vu;Au;Vl;Al];
    end

    %Phase One (Hopper on Ground) Events Function
    function [eventvalue,stopthecalc,eventdir] = phase_one_events(T,X)
        %Extract
        Ru = X(1:2);
        Vu = X(3:4);
        Rl = X(5:6);
        Vl = X(7:8);

        %Calculate Spring Focre
        Fs = k * (L0 - norm(Ru-Rl)) * (Ru-Rl)./norm(Ru-Rl);
        
        %Stop at Liftoff (Y Spring Force = Gravitational Force on Lower Mass)
        eventvalue = -Fs(2) - ml*g;  %  ‘Events’ are detected when eventvalue=0
        stopthecalc = 1;    %  stop if event occurs
        eventdir = 1;      %  Detect only events with dydt>0
    end

    %Phase Two (Hopper off Ground) Function Definition
    function derivatives = phase_two_equations(T, X)
        %Extract
        Ru = X(1:2);
        Vu = X(3:4);
        Rl = X(5:6);
        Vl = X(7:8);

        %Calculate Spring Focre
        Fs = k * (L0 - norm(Ru-Rl)) * (Ru-Rl)./norm(Ru-Rl);

        %Calculate Motion of Masses
        Au = [0;-g] + Fs/mu;
        Al = [0;-g] - Fs/ml;

        %Return State Vector
        derivatives = [Vu;Au;Vl;Al];
    end

    %Phase Two (Hopper off Ground) Events Function
    function [eventvalue,stopthecalc,eventdir] = phase_two_events(T,X)
        %Stop when lower mass hits ground (Yl=0)
        eventvalue = X(6);  %  ‘Events’ are detected when eventvalue=0
        stopthecalc = 1;    %  stop if event occurs
        eventdir = -1;      %  Detect only events with dydt<0
    end
end %function