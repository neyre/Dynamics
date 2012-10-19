% Planar Pendulum
%
% A simple planar pendulum consists of a mass supported by a (rigid) string.
% Gravity, aerodynamic drag and viscous drag act on the pendulum.
%
% Nick Eyre - Oct 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 3 - Problem 3

function pendulum

  %Parameters
  m = 1;   %Mass (kg)
  g = 9.8; %Gravity (m/s^2)
  r = 1;   %Radius (m)
  dv = 0;  %Viscous Damping Coefficient
  da = 0;  %Aerodynamic Drag Coefficient
  
  %Initial Conditions
  theta = pi/2; %Initial Angular Position (rad)
  omega = 0;    %Initial Angular Velocity (rad/s)

  %Simulation Parameters
  t = 10;
  dt = .01;
  plotextras = false;  %Plot energy and tension?

  %Setup & Run ODE45 Solver, Visualize Results
  X0 = [theta omega]; %Initial Condition Vector
  T = [0:dt:t]; %Time Vector
  [T, Z] = ode45(@equations, T, X0);  %Run Solver
  visualize(T,Z);

  %Equations of Motion
  function res = equations(T,X)

    %Extract State Variables
    theta = X(1);
    omega = X(2);

    %Calculate Derivative
    alpha = -g*sin(theta)/r - dv*omega - da*omega^2*sign(omega);

    %Return Derivative Vector
    res = [omega; alpha];

  end %equations

  %Visualize Results
  function visualize(T,Z)
    %Extract
    theta = Z(:,1);
    omega = Z(:,2);

    alpha = -g * sin(theta) / r;

    %Plot Just Displacement, Velocity, Acceleration
    if plotextras == false
      clf
      hold all
      plot(T,theta);
      plot(T,omega);
      plot(T,alpha);
      legend('Displacement (rad)','Velocity (rad/s)','Acceleration (rad/s^2)')
      xlabel('Time (s)')
      title('Pendulum Motion under Aerodynamic Drag')
    end %plotextras

    %Plot Everything
    if plotextras
      clf

      %Plot Displacement, Velocity, Acceleration
      subplot(3,1,1)
      hold all
      plot(T,theta);
      plot(T,omega);
      plot(T,alpha);
      legend('Displacement (rad)','Velocity (rad/s)','Acceleration (rad/s^2)')
      xlabel('Time (s)')
      title('Pendulum Motion')

      %Calculate Tension and Energy
      tension = m*g*cos(theta) + m*r*omega.^2;
      tensionx = tension .* sin(theta);
      tensiony = tension .* cos(theta);
      kinetic = .5*m*r^2 * omega.^2;
      potential = (1-cos(theta))*r*m*g;

      %Plot Energy
      subplot(3,1,2)
      hold all
      plot(T,kinetic)
      plot(T,potential)
      plot(T,kinetic+potential)
      legend('Kinetic','Potential','Total')
      xlabel('Time (s)')
      ylabel('Energy (J)')
      title('Energy')

      %Plot Components of Tension
      subplot(3,1,3)
      hold all
      plot(T,tensionx)
      plot(T,tensiony)
      plot(T,tension)
      legend('X-Tension','Y-Tension','Total Tension')
      title('Tension')
      ylabel('Force (N)')
      xlabel('Time (s)')
    end %plotextras
  end %visualize

end %coaster