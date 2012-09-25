% Coin Funnel
%
% Plot trajectory of coin in funnel
%
% Nick Eyre - Sept 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 2 - Problem 1
% Dynamics, Analysis and Design of Systems in Motion Problem 2.3.31

clear,clf

%Calculate Trajectory Given Equations from Problem
T = linspace(0,10,100);
Theta = 4*T;
Z = -T./(T+1);
R = exp(2*Z);

%Convert to Cartesian & Plot
[X,Y,Z] = pol2cart(Theta,R,Z);
plot3(X,Y,Z);

xlabel('X')
ylabel('Y')
zlabel('Z')
title('Coin in Funnel Trajectory')