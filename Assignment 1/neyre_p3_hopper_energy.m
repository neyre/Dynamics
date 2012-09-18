% Hopper Energy Simulation
%
% Simulates the flight of a hopper built with two masses separated by a
% spring.  Uses two phases (hopper on ground, hopper off ground) with the
% ODE45 differential equation solver.  The top mass of the hopper is 
% given an initial displacement downward before launching along 
% a trajectory.  This simulation plots the energy of the system as a
% function of time.
%
% Function Inputs: None
% Function Outputs: None
% The function will plot the energy of the hopper system over time.
%
% Initial conditions and parameters are specified in the function for easy running.
%
% Nick Eyre - Sept 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 1 - Problem 3

function hopper_energy
clear all

% Define system parameters
g = 9.81;      % gravitational acceleration in m/s^2
ml = .02;      % mass in kg
mu = .02;      % mass in kg
L0 = .05;       % rest length of spring
k = 2000;      % spring constant (N/m)
d = .005;       % initial height of upper

%Simulation Parameters
Tmax = 25;

%State Variable Definitions
% Z = [Yu, Vu, Yl, Vl]

%PHASE ONE SIMULATION
Z_0 = [L0-d, 0, 0, 0];
t_span = [0, Tmax];  %Max Time Span
options = odeset('Events', @phase_one_events); %Setup Events
[T1, Z1] = ode45(@phase_one_equations, t_span, Z_0, options);  %Run Solver

%PHASE TWO SIMULATION
Z_0 = Z1(end,:); %Initial Conditions are Phase 1 Output
t_span = [T1(end), Tmax];  %Max Time Span
options = odeset('Events', @phase_two_events); %Setup Events
[T2, Z2] = ode45(@phase_two_equations, t_span, Z_0, options);  %Run Solver

%Combine Data from Both Phases & Calculate Energy
%Potential Energy = PE Gravity Upper + PE Gravity Lower + PE Spring
%Kinetic Energy = KE Upper + KE Lower
Z = [Z1;Z2];   T = [T1;T2];
PE = mu*g*Z(:,1) + ml*g*Z(:,3) + .5*k*(Z(:,1)-Z(:,3)-L0).^2;
KE = .5*mu*Z(:,2).^2 + .5*ml*Z(:,4).^2;
E = PE+KE;

%Plot Energy
clf, hold on
plot(T,PE,'b');
plot(T,KE,'r');
plot(T,PE+KE,'--k');
legend('Potential Energy','Kinetic Energy','Total Energy','Location','SouthEast')
title('Hopper Simulation Energy')
ylabel('Energy (J)')
xlabel('Time (s)')
axis([0 T(end) 0 E(1)*1.1])

    %Phase One (Hopper on Ground) Function Definition
    function derivatives = phase_one_equations(T, X)
        %Extract
        Yu = X(1);
        Vu = X(2);
        Yl = X(3);
        Vl = X(4);

        %Lower Mass Doesn't Move b/c it is on ground
        Al = 0;

        %Calculate Motion of Upper Mass
        Fs = k * (L0 - (Yu-Yl));
        Au = (Fs - mu*g)/mu;

        %Return State Vector
        derivatives = [Vu;Au;Vl;Al];
    end

    %Phase One (Hopper on Ground) Events Function
    function [eventvalue,stopthecalc,eventdir] = phase_one_events(T,X)
        %Extract
        Yu = X(1);
        Yl = X(3);

        %Stop at Liftoff (Spring Force = Gravitational Force on Lower Mass)
        eventvalue = (Yu-Yl-L0)*k - ml*g;  %  ‘Events’ are detected when eventvalue=0
        stopthecalc = 1;    %  stop if event occurs
        eventdir = 1;      %  Detect only events with dydt>0
    end

    %Phase Two (Hopper off Ground) Function Definition
    function derivatives = phase_two_equations(T, X)
        %Extract
        Yu = X(1);
        Vu = X(2);
        Yl = X(3);
        Vl = X(4);

        %Calculate Motion of Upper Mass
        Fs = k * (L0 - (Yu-Yl));
        Au = (Fs - mu*g)/mu;
        Al = (-Fs - ml*g)/ml;

        %Return State Vector
        derivatives = [Vu;Au;Vl;Al];
    end

    %Phase Two (Hopper off Ground) Events Function
    function [eventvalue,stopthecalc,eventdir] = phase_two_events(T,X)
        %Extract
        Yl = X(3);

        %Stop when lower mass hits ground (Yl=0)
        eventvalue = Yl;  %  ‘Events’ are detected when eventvalue=0
        stopthecalc = 1;    %  stop if event occurs
        eventdir = -1;      %  Detect only events with dydt<0
    end
end
