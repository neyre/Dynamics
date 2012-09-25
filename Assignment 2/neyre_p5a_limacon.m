% Plot Limacon
%
% Plot a limacon defined by the following equation (polar)
% r = b - c*cos(theta)
%
% Nick Eyre - Sept 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 2 - Problem 5a

%Constants
b = 10;
c = -1;

%Calculate
theta = linspace(0,2*pi,100);
r = b - c*cos(theta);

%Plot
polar(theta,r)