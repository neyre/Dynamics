% Double Pendulum
%
% Simulates the trajectory of a double pendulum with massless rods in a plane.
% The pendulum is released with given initial conditions and its trajectory is plotted.
% Through simulation parameters, the simulation can be configured to plot the tension in
% the rods, the energy of the system, the state variables (position & velocity) or
% animate the system.
%
% Function Inputs: None
% Function Outputs: None
%
% Initial conditions and parameters are specified in the function for easy running.
%
% Eyre, Holzgrafe, Kessler - October 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 4 - Problem 4

function pendulum

    %Parameters
    m1 = 1;  %Mass (kg)
    m2 = 1;  %Mass (kg)
    l1 = 1;  %Length (m)
    l2 = 1;  %Length (m)
    g = -9.8; %Gravity (m/s^2)
    
    %Initial Conditions
    theta1 = pi/2;
    theta2 = 0;
    omega1 = 0;
    omega2 = 0;

    %Simulation Parameters
    t = 10;
    dt = .01;
    skipframes = 14;
    animate = false;
    plotpath = true;
    plotstatevariables = false;
    plotenergy = true;
    plottension = true;

    %Setup & Run Simulation
    X0 = [theta1 theta2 omega1 omega2]; %Initial Condition Vector
    T = [0:dt:t]; %Time Vector
    [T, Z] = ode45(@equations, T, X0);  %Run Solver
    visualize(T,Z);

    % Equations for ODE Solver
    function res = equations(T,X)
        %Unpack State Variables
        theta1 = X(1);
        theta2 = X(2);
        omega1 = X(3);
        omega2 = X(4);

        %Calculate Derivatives
        alpha1 = (g*(2*m1+m2)*sin(theta1) + m2*(g*sin(theta1-2*theta2) - ...
            2*(l2*omega2^2 + l1*omega1^2*cos(theta1-theta2))*sin(theta1-theta2)))...
            /(l1*(2*m1+m2-m2*cos(2*(theta1-theta2))));
        alpha2 = 2*sin(theta1-theta2)*(l1*(m1+m2)*omega1^2 - ...
            g*(m1+m2)*cos(theta1) + l2*m2*omega2^2*cos(theta1-theta2))...
            /(l2*(2*m1 + m2 - m2*cos(2*(theta1-theta2))));

        %Return Derivatives of State Variables
        res = [omega1;omega2;alpha1;alpha2];

    end %equations

    % Visualize Simulation Results
    function visualize(T,Z)
        if animate
            runanimation(Z);
        end
        if plotpath
            plotthepath(Z);
        end
        if plotstatevariables
            plotstates(T,Z);
        end
    end

    function [X1 Y1 X2 Y2 X1d Y1d X2d Y2d] = converttocartesian(Z)
        %Extract
        theta1 = Z(:,1);
        theta2 = Z(:,2);
        omega1 = Z(:,3);
        omega2 = Z(:,4);

        %Position
        X1 = l1 * sin(theta1);
        Y1 =-l1 * cos(theta1);
        X2 = l2 * sin(theta2) + X1;
        Y2 =-l2 * cos(theta2) + Y1;

        %Velocity
        X1d = l1*omega1.*cos(theta1);
        Y1d =-l1*omega1.*sin(theta1);
        X2d = l2*omega2.*cos(theta2) + X1d;
        Y2d =-l2*omega2.*sin(theta2) + Y1d;
    end

    function runanimation(Z)
        %Convert to Cartesian
        [X1 Y1 X2 Y2 X1d Y1d X2d Y2d] = converttocartesian(Z);

        %Configure Axes
        lmax = l1 + l2;
        ax = [-lmax lmax -lmax lmax];

        for i=1:skipframes:length(theta1)
            clf, hold on

            %Plot Historical Position
            plot(X1(1:i),Y1(1:i),'b')
            plot(X2(1:i),Y2(1:i),'g')

            %Plot Current System Position
            plot([0 X1(i)],[0 Y1(i)],'k')
            plot([X1(i) X2(i)],[Y1(i) Y2(i)],'k')
            plot(0,0,'k.','MarkerSize',20)
            plot(X1(i),Y1(i),'k.','MarkerSize',20)
            plot(X2(i),Y2(i),'k.','MarkerSize',20)


            legend('M1 Path', 'M2 Path')
            axis(ax)
            axis square
            drawnow;
        end %for
    end %visualize

    function plotthepath(Z)
        %Convert to Cartesian
        [X1 Y1 X2 Y2 X1d Y1d X2d Y2d] = converttocartesian(Z);

        %Configure Axes
        lmax = l1 + l2;
        ax = [-lmax lmax -lmax lmax];

        %Plot
        clf; hold on
        plot(X1,Y1,'b')
        plot(X2,Y2,'g')
        plot([0 X1(end)],[0 Y1(end)],'k')
        plot([X1(end) X2(end)],[Y1(end) Y2(end)],'k')
        plot(0,0,'k.','MarkerSize',20)
        plot(X1(end),Y1(end),'k.','MarkerSize',20)
        plot(X2(end),Y2(end),'k.','MarkerSize',20)
        legend('M1 Path', 'M2 Path')
        axis(ax)   
        axis square
    end %plotthepath

    function plotstates(T,Z)
        %Extract
        theta1 = Z(:,1);
        theta2 = Z(:,2);
        omega1 = Z(:,3);
        omega2 = Z(:,4);

        %Calculate Accelerations
        alpha1 = [];
        alpha2 = [];
        for i=1:length(theta1)
            a1 = (g*(2*m1+m2)*sin(theta1(i)) + m2*(g*sin(theta1(i)-2*theta2(i)) - ...
                2*(l2*omega2(i)^2 + l1*omega1(i)^2*cos(theta1(i)-theta2(i)))*sin(theta1(i)-theta2(i))))...
                /(l1*(2*m1+m2-m2*cos(2*(theta1(i)-theta2(i)))));
            a2 = 2*sin(theta1(i)-theta2(i))*(l1*(m1+m2)*omega1(i)^2 - ...
                g*(m1+m2)*cos(theta1(i)) + l2*m2*omega2(i)^2*cos(theta1(i)-theta2(i)))...
                /(l2*(2*m1 + m2 - m2*cos(2*(theta1(i)-theta2(i)))));
            alpha1 = [alpha1;a1];
            alpha2 = [alpha2;a2];
        end
        
        clf

        %Determine Mode
        curplot = 3;
        numplots = 2;
        if plotenergy
            numplots = numplots + 1;
        end
        if plottension
            numplots = numplots + 2;
        end
        
        subplot(numplots,1,1)
        hold on
        plot(T,theta1,'r')
        plot(T,omega1,'b')
        plot(T,alpha1,'g')
        title('Mass One')
        legend('Angular Position (rad)','Angular Velocity (rad/s)','Angular Acceleration (rad/s^2')
        
        subplot(numplots,1,2)
        hold on
        plot(T,theta2,'r')
        plot(T,omega2,'b')
        plot(T,alpha2,'g')
        title('Mass Two')
        legend('Angular Position (rad)','Angular Velocity (rad/s)','Angular Acceleration (rad/s^2')
        
        if plotenergy
            %Potential Energy
            [X1 Y1 X2 Y2 X1d Y1d X2d Y2d] = converttocartesian(Z);
            potential = -g*(m1*Y1 + m2*Y2);

            %Kinetic Energy
            v1 = (X1d.^2 + Y1d.^2).^(1/2);
            v2 = (X2d.^2 + Y2d.^2).^(1/2);
            k1 = .5 * m1 * v1.^2;
            k2 = .5 * m2 * v2.^2; 
            kinetic = k1 + k2;

            %Total Energy
            energy = potential + kinetic;

            subplot(numplots,1,curplot)
            hold on
            plot(T,potential,'r')
            plot(T,kinetic,'b')
            plot(T,energy,'k')
            title('Energy')
            ylabel('Energy (J)')
            legend('Potential','Kinetic','Total')

            curplot = curplot + 1;
        end

        if plottension
            %Calculate Tensions
            T1 = 2*m1*(l1*(m1+m2)*omega1.^2 - g*(m1+m2)*cos(theta1) ...
                + l2*m2*omega2.^2.*cos(theta1-theta2))./(2*m1+m2-m2*cos(theta1-theta2));
            T2 = -m1*m2*(-2*l2*omega2.^2 - 2*l1*omega1.^2.*cos(theta1-theta2)...
                +g*cos(2*theta1-theta2) + g*cos(theta2))./(2*m1+m2-m2*cos(theta1-theta2));

            %Calculate Reaction Force Components
            Rx = T1 .* sin(theta1);
            Ry = -T1 .* cos(theta1);

            subplot(numplots,1,curplot)
            hold on
            plot(T,T1,'r')
            plot(T,T2,'b')
            title('Tension')
            ylabel('Force (N)')
            legend('T1','T2')

            curplot = curplot + 1;

            subplot(numplots,1,curplot)
            hold on
            plot(T,Rx,'r')
            plot(T,Ry,'b')
            title('Reaction Force')
            ylabel('Force (N)')
            legend('X Component','Y Component')

        end

        xlabel('Time (s)')

    end %plotstates

end %coaster