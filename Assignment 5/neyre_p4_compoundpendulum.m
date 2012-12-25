% Compound Double Pendulum
%
% Simulates the trajectory of a compound double pendulum in a plane.
% The pendulum is released with given initial conditions and its trajectory is plotted.
% Through simulation parameters, the simulation can be configured to plot the state variables
% or animate the trajectory of the system.  The second derivatives for the system were computed
% in Wolfram Mathematica.
%
% Function Inputs: None
% Function Outputs: None
%
% Initial conditions and parameters are specified in the function for easy running.
%
% Eyre, Holzgrafe, Kessler - October 2012
% Olin College ENGR 2340 - Dynamics
% Assignment 5 - Problem 4

function pendulum

    %Parameters
    m1 = 1;  %Mass (kg)
    m2 = 1;  %Mass (kg)
    l1 = 1;  %Length (m)
    l2 = 1;  %Length (m)
    g = 9.8; %Gravity (m/s^2)
    
    %Moments of Inertia About Pivot Points
    I1G = m1*l1^2/12;
    I1O = m1*l1^2/3;
    I2G = m2*l2^2/12;
    
    %Initial Conditions
    theta1 = pi/2;
    theta2 = pi/2;
    omega1 = 0;
    omega2 = 0;

    %Simulation Parameters
    t = 8.7;
    dt = .01;
    tolerance = 1e-14;  %Relative Integration Tolerance
    skipframes = 14;
    animate = false;
    plotpath = true;
    plotstatevariables = false;

    %Setup & Run Simulation
    X0 = [theta1 theta2 omega1 omega2]; %Initial Condition Vector
    T = [0:dt:t]; %Time Vector
    options = odeset('RelTol', tolerance);  %Setup Events
    [T, Z] = ode45(@equations, T, X0,options);  %Run Solver
    visualize(T,Z);

    function res = equations(T,X)
        %Unpack
        O1 = X(1);
        O2 = X(2);
        O1d = X(3);
        O2d = X(4);

        %Calculate
        O1dd=-(l1*(g*(l2^2*m2*(m1+m2)+4*I2G*(m1+2*m2))*sin(O1)+ ...
            l2*m2*(g*l2*m2*sin(O1-2*O2)+((4*I2G+l2^2*m2)*O2d^2+ ...
            2*l1*l2*m2*O1d^2*cos(O1-O2))*sin(O1-O2))))/ ...
            (2*I1O*(4*I2G+l2^2*m2)+l1^2*m2*(8*I2G+l2^2*m2)- ...
            l1^2*l2^2*m2^2*cos(2*(O1-O2)));

        O2dd=-(l2*m2*(4*l1*(I1O+l1^2*m2)*O1d^2*sin(O1-O2)+l1^2*l2*m2*O2d^2 ...
            *sin(2*(O1-O2)) + g*(l1^2*(m1+2*m2)*sin(2*O1-O2)+(-4*I1O+l1^2*(m1-2*m2)) ...
                *sin(O2))))/(-2*I1O*(4*I2G+l2^2*m2)-l1^2*m2*(8*I2G+l2^2*m2)+...
                l1^2*l2^2*m2^2*cos(2*(O1-O2)));
            
        %Output
        res = [O1d;O2d;O1dd;O2dd];

    end %equations

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

        xlabel('Time (s)')

    end %plotstates

end %coaster