% Marble in Bowl
%
% Simulates the trajectory of a marble sliding frictionless in a hemispherical bowl.
% The marble is released at the rim of the bowl with a tangential velocity.
% Gravity is the only force acting on the ball.  The simulation is implemented in
% spherical coordinates using the ODE45 differential equation solver.  After the
% simulation is conducted, the trajectory are converted to 3D cartesian and is plotted.
%
% Function Inputs: None
% Function Outputs: None
% The function will plot the trajectory of the marble in the X-Y-Z plane.
%
% Initial conditions and parameters are specified in the function for easy running.
%
% Nick Eyre - Sept 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 2 - Problem 6

function marble
clear all,clf

%System Parameters
g = 9.81;      % gravitational acceleration in m/s^2
r = .2;        % bowl radius in meters
v0 = 5;        % tangential initial velcotiy in m/s
phi_0 = pi/2;  % initial phi value in radians

%Simulation Parameters
Tmax = 10;
dt = .01;
animate = false;

% Theta and R have constant derivatives
v0 = v0/r;  % Initial Theta Dot
Radius_Dot = 0;

%State Variable Definitions
% Z = [Theta, d/dt Theta, Phi, d/dt Phi, R]

% Simulate & Visualize
Z0 = [0 v0 phi_0 0 r];  %Marble Initial Conditions
Ts = [0:dt:Tmax];  %Max Time Span
[T, Z] = ode45(@equations, Ts, Z0);  %Run Solver
visualize_trajectory(T,Z);

    % Function Definition
    % Input: Time, State; Output: Derivative
    function derivatives = equations(T, X)
        %Extract State Variables
        Theta     = X(1);
        Theta_Dot = X(2);
        Phi       = X(3);
        Phi_Dot   = X(4);
        Radius    = X(5);

        Phi_DDot = g*sin(Phi)/r + Theta_Dot^2 * sin(Phi)*cos(Phi);
        Theta_DDot = -2*Phi_Dot*Theta_Dot/tan(Phi);

        %Return Derivative Vector
        derivatives = [Theta_Dot; Theta_DDot; Phi_Dot; Phi_DDot; Radius_Dot];
    end

    % Plot the Trajectory of the Marble in 3D
    % Input: Output Vector from ODE45
    function visualize_trajectory(T,Z)
        %Extract
        Theta  = Z(:,1);
        Phi    = Z(:,3);
        Radius = Z(:,5);
        
        %Convert to Cartesian
        [Theta,Phi,Radius] = phi_convert(Theta,Phi,Radius);
        [x,y,z]=sph2cart(Theta,Phi,Radius);
        hold on
        
        %Animate It!
        if animate
            for i=1:length(Theta)
                %Plot
                plot3(x(i),y(i),z(i),'.');

                %Build Axes & Title
                xlabel('X (m)')
                ylabel('Y (m)')
                zlabel('Z (m)')
                title('Marble Trajectory')
                
                %Scale & Show
                axis([-r r -r r -r 0])
                drawnow
            end
        end

        %Plot Final Result
        plot3(x,y,z)
        xlabel('X (m)')
        ylabel('Y (m)')
        zlabel('Z (m)')
        title('Marble Trajectory')
        axis([-r r -r r -r 0])
    end

    % Converts Phi from Range of [0,Pi] to [-Pi/2,Pi/2] for sph2cart built-in function
    function [Theta,Phi,Radius] = phi_convert(Theta,Phi,Radius)
        %Deal with Overflow (Phi > Pi)
        for i=1:length(Theta)
            if Phi(i) > pi
                Phi(i) = 2*pi - Phi(i);
                Theta(i) = Theta(i) + pi;
            end
        end

        %Convert Range
        Phi = pi/2 - Phi;
    end

end %function