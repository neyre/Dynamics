% String Pull
%
% A mass is rotating around a fixed point, connected by a string.  The string
% is pulled at a constant acceleration through a hole, shortening the radius
% of rotation.  Plot the tension in the string as a function of time.
%
% Nick Eyre - Oct 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 3 - Problem 1
% Dynamics, Analysis and Design of Systems in Motion Problem 3.2.33

function stringpull
  clear

  %Parameters
  m = .4;

  %Initial Conditions
  % [Radius, Theta, Theta_dot, Tension]
  V = [2.5 0 8 0];

  %Simulation Setup
  t = 1.2;
  dt = .0002;
  T = [0:dt:t];

  %Run Forward Euler Simulation
  for i=2:length(T)
    %Calculate r, rdot and rddot (all are simple functions of time)
    r      = -T(i)^2 + 2.5;
    r_dot  = -T(i)*2;
    r_ddot = -2;

    %Upack State Variables
    theta      = V(i-1,2);
    theta_dot  = V(i-1,3);

    %Calculate Derivatives
    theta_ddot = -2*r_dot*theta_dot/r;
    
    %Apply Derivatives to State Variables
    theta_new     = theta     + theta_dot *dt;
    theta_dot_new = theta_dot + theta_ddot*dt;

    %Calculate Tension
    tension = m*r*theta_dot^2 - m*r_ddot;

    %Repack State Variables
    V = [V;[r theta_new theta_dot_new tension]];
  end %for

  %Plot Tension Results from Numerical Simulation Solution
  clf;hold all
  plot(T,V(:,4))
  xlabel('Time (s)')
  ylabel('Tension (N)')
  title('String Tension')
  axis([0 1.2 0 800])
  
  %Plot Analytical Solution
  T = '0.8 + 0.4*(2.5-x^2)*(10*sqrt(2)/((5-t^2)^2))^2';
  ezplot(T)

end %function