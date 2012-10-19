% Spherical Spring Pendulum
%
% A mass is suspended on a massless linear spring to form a spherical
% pendulum.
%
% Note: Due to numerical integration errors, energy is not always conserved
% by the simulation.  This would be something to look into fixing in a v2.
%
% Nick Eyre - Oct 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 3 - Problem 4

function springpendulum

  %Parameters
  k  = 250;   %Spring Constant (N/m)
  m  = 1;     %Mass (kg)
  L0 = 1;     %Unstretched Spring Length (m)
  g = 9.8;    %Gravity (m/s2)
  
  %Initial Conditions
  r         = 1;  %Radius (m)
  theta     = 0;  %Azimuthal Angle
  phi       = pi/2;   %Polar Angle, 0>phi>pi to avoid gimbal lock
  r_dot     = 0;  %Radial Velocity
  theta_dot = 0;  %Aximuthal Angular Velocity
  phi_dot   = 0;  %Polar Angular Velocity

  %Simulation Parameters
  t = 4;  %Max Integration Time
  dt = .01; %Time Step
  animate = true; 
  plot2d = false;
  plotenergy = true;

  %Setup & Run Simulation, Visualize Results
  X0 = [r theta phi r_dot theta_dot phi_dot]; %Initial Condition Vector
  T = [0:dt:t]; %Time Vector
  options = odeset('RelTol',1e-15);
  [T, Z] = ode45(@equations, T, X0,options);  %Run Solver
  visualize(T,Z);

  %Equations of Motion
  function res = equations(T,X)
    %Extract Variables
    r         = X(1);
    theta     = X(2);
    phi       = X(3);
    r_dot     = X(4);
    theta_dot = X(5);
    phi_dot   = X(6);

    %Calculate Derivatives
    %Formulas for Derivatives were analytically derived
    r_ddot = k*(L0-r)/m - g*cos(phi) + r*phi_dot^2 + r*theta_dot^2*sin(phi)^2;
    theta_ddot = -2*r_dot*theta_dot/r - 2*theta_dot*phi_dot/tan(phi);
    phi_ddot = g*sin(phi)/r - 2*r_dot*phi/r + sin(phi)*cos(phi)*theta_dot^2;

    %Return Derivative Vector
    res = [r_dot; theta_dot; phi_dot; r_ddot; theta_ddot; phi_ddot];
  end %equations

  % Visualize Results
  % Input: Output Vector from ODE45
  function visualize(T,Z)
    %Extract
    Radius = Z(:,1);
    Theta  = Z(:,2);
    Phi    = Z(:,3);
    Radius_dot = Z(:,4);
    Theta_dot  = Z(:,5);
    Phi_dot    = Z(:,6);

    %Calculate Potential & Kinetic Energy of System
    potentialspring = .5*k*(L0-Radius).^2;
    potentialgrav = m*g*Radius.*cos(Phi);
    potential = potentialspring + potentialgrav;
    kinetic = .5*m*(Radius_dot.^2 + ...
            Radius.^2 .* Phi_dot.^2 + ...
            Radius.^2 .* sin(Phi).^2 .* Theta_dot.^2);
    
    %Convert to Cartesian
    [Theta,Phi,Radius] = phi_convert(Theta,Phi,Radius);
    [Theta_dot,Phi_dot,Radius_dot] = phi_convert(Theta_dot,Phi_dot,Radius_dot);
    [x,y,z]=sph2cart(Theta,Phi,Radius);
    [x_dot,y_dot,z_dot]=sph2cart(Theta_dot.*Radius.*sin(Phi), Phi_dot.*Radius, Radius_dot);
    
    clf
    
    %Animate
    if animate
      axislimits = [min(x) max(x) min(y) max(y) min(z) max(z)];
      for i=1:length(Theta)
        %Plot
        clf; hold on
        plot3(x,y,z)
        plot3(x(i),y(i),z(i),'k.','MarkerSize',20);

        %Build Axes & Title
        xlabel('X (m)')
        ylabel('Y (m)')
        zlabel('Z (m)')
        title('Pendulum Trajectory')
        
        %Scale & Show
        axis(axislimits)
        drawnow
      end
    end

    %Plot Final Result
    clf
    plot3(x,y,z)
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    title('Pendulum Trajectory')

    %Plot position in X-Z plane
    if plot2d
      clf
      plot(x,z)
      xlabel('X (m)')
      ylabel('Z (m)')
      title('Pendulum Trajectory')
    end %plot2d

    %Plot Energy as function of Time
    if plotenergy
      clf, hold all
      plot(T,kinetic);
      plot(T,potential);
      plot(T,kinetic+potential);
      legend('Kinetic','Potential','Total')
      title('Pendulum Energy')
      ylabel('Energy (J)')
      xlabel('Time (s)')
    end

  end

  % Converts Phi from Range of [-Inf,Inf] to [-Pi/2,Pi/2] for sph2cart built-in function
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

end %springpendulum