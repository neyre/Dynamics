% Roller Coaster
%
% The track of a roller coaster is given by a function y(x) where several
% parameters control the shape of the track.  This simulation simulates the
% trajectory of a car traveling on the roller coaster subject to gravity
% as well as aerodynamic drag.
%
% Nick Eyre - Oct 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 3 - Problem 2

function coaster

  %Parameters
  B = .23;    % 0 < B < 1 controls height of second peak
  h = 100;    % Height of Track (m)
  L = 500;    % Length of Track (m)
  m = 500;    % Mass of Car (kg)
  g = 9.81;   % Acceleration of Gravity (m/s^2)
  rho = 1.02; % Air Density (kg/m^3)
  Cd = .1;    % Drag Coefficient of Car
  Area = 4;   % Frontal Area of Car (m^2)
  
  %Initial Conditions of Car
  x = 0;
  y = h;
  vx = .1;
  vy = 0;

  %Simulation Parameters
  t = 100;    %Max Integration Time
  dt = .1;    %Time step
  tolerance = 1e-14;  %Relative Integration Tolerance
  animate = false;

  %Setup & Run ODE45 Simulation, Visualize Results
  X0 = [x y vx vy]; %Initial Condition Vector
  T = [0:dt:t]; %Time Vector
  options = odeset('Events', @events, 'RelTol', tolerance);  %Setup Events
  [T, Z] = ode45(@equations, T, X0, options);  %Run Solver
  visualize(T,Z);

  %Equations of Motion Function
  function res = equations(T,X)
    %Extract Variables
    x = X(1);
    y = X(2);
    vx = X(3);
    vy = X(4);

    %Get State Parameters
    a = AofX(x);  %Slope
    b = BofX(x);  %Bending
    c = rho * Cd * Area / 2 / m;  %Combined Air Drag Coefficient
    V2 = (vx^2 + vy^2)*sign(vx);  %Magnitude of v^2

    %Build A derivative matrix and V state vector
    A = [1 0 0 0 0;...
         0 1 0 0 0;...
         0 0 1 0 (a/(m*sqrt(1+a^2)));...
         0 0 0 1 (-1/(m*sqrt(1+a^2)));...
         0 0 a -1 0];
    V = [vx;...
         vy;...
         (-c*V2/sqrt(1+a^2));...
         (-a*c*V2/sqrt(1+a^2) - g);...
         (-b*vx^2)];

    %Calculate and Return Derivatives
    D = A\V;
    res = D(1:4); %Derivatives
    %N = D(5); %Normal Force
  end %equations

  %Events Function for ODE45 Solver
  function [eventvalue,stopthecalc,eventdir] = events(T,X)
    %Stop when x = L
    eventvalue = X(1)-L;
    stopthecalc = 1;
    eventdir = 0;
  end

  %Function to Visualize Results
  function visualize(T,Z)
    %Extract
    X = Z(:,1);
    Y = Z(:,2);

    %Calculate Y(x) path trajectory
    Xs = linspace(0,L,100);
    Ys = h*(1-B*Xs/L).*cos(3*pi*Xs/2/L).^2;

    %Animate Car position
    if animate
      for i=1:length(T)
        clf; hold all;
        plot(Xs,Ys)
        plot(X(i),Y(i),'k.','MarkerSize',50)
        drawnow
      end
    end

    %Plot X vs Y
    clf; hold all
    plot(Xs,Ys)     %Plot Simulation Result
    plot(X,Y,'--')  %Plot Y(x) exact path trajectory
    legend('Numerical Simulation','Exact Path Trajectory')
    xlabel('Horizontal Distance (m)')
    ylabel('Height (m)')
    title('Roller Coaster Trajectory')
    
  end %visualize

  %Calculate Slope of the path A(x)
  %Function for A(x) was calculated analytically
  function res = AofX(x)
    res = -h/L*cos(3*pi*x/2/L)*(3*pi*(1-B*x/L)*sin(3*pi*x/2/L)+B*cos(3*pi*x/2/L));
  end %AofX

  %Calculate Bending of the path B(x)
  %Function for B(x) was calculated analytically
  function res = BofX(x)
    res = 6*B*h*pi/L^2* cos(3*pi*x/(2*L))* sin(3*pi*x/(2*L)) + h*9*pi^2/(2*L^2)* (1-B*x/L)* ((sin(3*pi*x/(2*L)))^2 - ...
    (cos(3*pi*x/(2*L)))^2);
  end %BofX


end %coaster