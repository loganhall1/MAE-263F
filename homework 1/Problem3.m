%%%%%%%%% Problem 3 %%%%%%%%%
%%%%%%%%% MAE 263F %%%%%%%%%%
%%%%%%%%% Homework 1 %%%%%%%%
%%%%%%%%% Logan Hall %%%%%%%%
%%%%%% October 18 2023 %%%%%%
% Inspired by code from UCLA MAE 263F Lecture 4 by Professor Khalid Jawed



%% Reset MATLAB Workspace
clear;
close all;
clc;




%% Problem 3.1 - Plot Maximum Vertical Displacement and Compare to Theoretical Prediction
N1 = 50; %Number of nodes
P1 = 2000; %[N] - Point load, as given in problem statement
dt1 = 10^(-2);% [second] - Time step size

%Call function to run the beam simulation
[timeArray1, y_max1, y_max_sim1, y_max_euler1] = Problem3_simulation(N1,P1,dt1);

%%%%%
%Plot the time vs maximum beam deflection
plot1 = figure(1);
plot(timeArray1, y_max1, 'k-');
%Format plot
xlabel('Time, t [sec]');
ylabel('Maximum Beam Deflection, y_{max} [m]');
grid on;
title("Time vs Maximum Beam Deflection, y_{max} at P = "+ P1 + "N");
saveas(plot1,'Problem3.1.png');
%%%%%




%% Problem 3.2 Part 1 - Plot Maximum Vertical Displacement Under a Higher Load
N2 = 50; %Number of nodes
P2 = 20000;%[N] - Point load, as given in problem statement
dt2 = 10^(-2);%[second] - Time step size

%Call function to run the beam simulation with larger P
[timeArray2, y_max2, y_max_sim2, y_max_euler2] = Problem3_simulation(N2,P2,dt2);

%%%%%
%Plot the time vs maximum beam deflection
plot2 = figure(2);
plot(timeArray2, y_max2, 'k-');
%Format plot
xlabel('Time, t [sec]');
ylabel('Maximum Beam Deflection, y_{max} [m]');
grid on;
title("Time vs Maximum Beam Deflection, y_{max} at P = "+ P2 + "N");
saveas(plot2,'Problem3.2_HigherP.png');
%%%%%


%% Problem 3.2 Part 2 - Plot P vs y_max to Evaluate Where Solutions Diverge
%Define variables
N3 = 50; %Number of nodes
P_array = 100:100:20000;%[N] - Range of point loads
dt3 = 10^(-2);%[second] - Time step size

%Initialize vectors
y_max_sim3 = zeros(length(P_array),1);%Initialize array of y_max for the simulation
y_max_euler3 = zeros(length(P_array),1);%Initialize array of y_max for the euler method
theoretical_error = zeros(length(P_array),1);%Initialize array of theoretical error between simulation and euler method



%For loop to cycle through the range of point loads
for i = 1:length(P_array)

    %Call function to evaluate y_max for both simulation and theoretical methods
    [timeArray3(i,:), y_max3(:,i), y_max_sim3(i), y_max_euler3(i)] = Problem3_simulation(N3,P_array(i),dt3);
    
    %Calulate the percent error between simulation and theoretical methods
    theoretical_error(i) = abs(((y_max_sim3(i) - y_max_euler3(i))/y_max_sim3(i))*100);

end

%Determining the first load where the percent error between simulated and theoretical is greater than 3%
divergent_load = P_array(find(theoretical_error > 3,1));




%%%% Plot the Percent Error Between Theoretical and Simulation
plot3 = figure(3);
plot(P_array, theoretical_error, 'r-');% Plot Percent Error
hold on
plot([P_array(1) P_array(end)], [3 3], 'k-')% Plot a constant 3% error
%Format plot
xlabel('Load, P [N]');
ylabel('Percent Error [%]');
lgd4 = legend('Percent Error Between Theoretical and Simulation', '3% Error');
lgd4.Location = 'Northwest';
grid on;
title('Percent Error Between Theoretical and Simulation');
saveas(plot3,'Problem3.2_PercentError.png');
%%%%





%%%% Plot the Comparison Between Theoretical and Simulation - Zoomed Out
plot4 = figure(4);
plot(P_array, y_max_sim3, 'k-');
hold on;
plot(P_array, y_max_euler3, 'r-');
%Format Plot
xlabel('Load, P [N]');
ylabel('Maximum Beam Deflection, y_{max} [m]');
legend('y_{max, euler}', 'y_{max, simulation}')
grid on;
title('Comparing Theoretical and Simulation - Zoomed Out');
saveas(plot4,'Problem3.2_Zoomout.png');
%%%%





%%%% Plot the Comparison Between Theoretical and Simulation - Zoomed In
plot5 = figure(5);
plot(P_array, y_max_sim3, 'k-');
hold on;
plot(P_array, y_max_euler3, 'r-');
%Format Plot
xlabel('Load, P [N]');
ylabel('Maximum Beam Deflection, y_{max} [m]');
legend('y_{max, euler}', 'y_{max, simulation}')
xlim([0,5000])
grid on;
title('Percent Error Between Theoretical and Simulation - Zoomed In');
saveas(plot5,'Problem3.2_Zoomin.png');
%%%%







function [timeArray, y_max, y_max_sim, y_max_euler] = Problem3_simulation(N,Load,dt)
% Number of nodes

ndof = N * 2; % number of degrees of freedom 
RodLength = 1; %[m] given bar length
deltaL = RodLength / (N-1); %[m]

loadDistance = .75;%Distance point load is applied from left edge
loadNode = round(N * (loadDistance/RodLength));%Determine which Node to apply load
dofLoadNode = loadNode*2;%Determine which dof to apply load to

% Inner and outer radii of beam cross section
R_o = 0.013; %[m] Outer radii of beam cross section
r_i = 0.011; %[m] Inner radii of beam cross section




% Density
rho_metal = 2700; % kg/m^3

Y = 70*10^9; % Young's modulus (Y instead of E for clarity)

g = 9.8; % m/s^2 - gravity

totalTime = 1; % [sec] - total simulation time

% Utility parameter
I_inertia = (pi/4)*(R_o^4 - r_i^4); % [in^4] Calculate moment of inertia
EI = Y * I_inertia; % [Nm^2] - bending stiffness
EA = Y * pi * (R_o^2 - r_i^2); % [N]

% Geometry - initial configuration
nodes = zeros(N,2);
for c=1:N % Loop over all the nodes
    nodes(c,1) = (c-1) * deltaL; % x coordinates
    nodes(c,2) = 0;
end


% Point Load Vector, P

P = zeros(ndof, 1);%Initalize with no loads applied
P(dofLoadNode) = -Load;%Apply point load P at the dof the load is applied

% Mass, M
M = zeros(ndof, ndof);
for k=1:N
    M(2*k-1, 2*k-1) = (pi*(R_o^2 - r_i^2)*RodLength*rho_metal)/(N-1); % Mass for x_k
    M(2*k, 2*k) = M(2*k-1, 2*k-1); % Mass for y_k
end


% Initial DOF
q0 = zeros(ndof, 1);
for c=1:N % loop over nodes
    q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
    q0( 2*c ) = nodes(c,2); % y1, y2, y3
end

u0 = zeros(ndof, 1); % old velocity (initial velocity)

% tolerance
tol = EI/RodLength^2 * 1e-3; % small enouch force that can be neglected

% Time marching scheme
Nsteps = round(totalTime/dt);

% 
y_max = zeros(Nsteps, 1); %store all dof vector at each time step
y_max(1) = min(q0(2:2:end)); % Taking min

for c = 2:Nsteps
    %Printing out the time at each loop iteration has been commented out.
    %uncomment to view the time at each iteration
    %fprintf('Time = %f\n', (c-1) * dt);

    fixedDOF = [1; 2; ndof]; %Define which dof are fixed from boundary contitions
    freeDOF = 3:ndof-1; %Define which dof are free from boundary contitions
    boundaryConditionVector = [0;0;0];%[m] Define boundary contitions

    % Guess
    q = q0; % New DOFs are initialized to be equal to old DOFs
    q(fixedDOF) = boundaryConditionVector; % Imposing boundary conditions
    
    % Newton Raphson
    err = 10 * tol;
    while err > tol

        q_free = q(freeDOF);%Define free degrees of freedom

        f = M / dt * ( (q-q0)/dt - u0 ); %Calculate first part of f, Inertia
        J = M / dt^2; %Calculate Jacobian, Inertia


        %
        % Elastic forces
        %

        % Linear spring
        for k=1:N-1
            %Define x_k, y_k, y_(k+1), and x_(k+1)
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            l_k = deltaL;
            %Calculate derivative of stretching energy E_k^s
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            %Calculate the hessian of the stretching energy E_k^s
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            ind = [2*k-1, 2*k, 2*k+1, 2*k+2];
            %Update Jacobian and f
            f(ind) = f(ind) + dF;
            J(ind,ind) = J(ind,ind) + dJ;
        end

        % Bending spring
        for k=2:N-1
            %Define x_k, y_k, y_(k+1), x_(k+1), y_(k-1), x_(k-1)
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            l_k = deltaL;
            %Calculate the derivative of bending energy E_k^b
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            %Calculate the hessian of the bending energy E_k^b
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            ind = [2*k-3, 2*k-2, 2*k-1, 2*k, 2*k+1, 2*k+2];
            %Update Jacobian and f
            f(ind) = f(ind) + dF;
            J(ind, ind) = J(ind, ind) + dJ;
        end

        %Point Load
        f = f - P;

        %Define which elements of F and J are free to help with
        %calculations in next step
        f_free = f(freeDOF);
        J_free = J(freeDOF, freeDOF);

        % Update for the free components of dq and q
        dq_free = (J_free \ f_free);
        q_free = q_free - dq_free;
        
        %Recalculate err to see if it can exit while loop
        err = sum ( abs(f_free) );
        q(freeDOF) = q_free;
    end

    % New velocity
    u = (q - q0) / dt;

    % Store some information
    y_max(c,:) = min(q(2:2:end));


    % The drawnow plot function has been commented out to increase the
    % speed the code can run. It may be uncommented to visualize the
    % simulation.

    % Plot
    % figure(1);
    % plot( q(1:2:end), q(2:2:end), 'ro-');
    % 
    % xlabel('x [meter]');
    % ylabel('y [meter]');
    % drawnow
    % Update (new becomes old)
    q0 = q;
    u0 = u;
end



%% Plot maximum beam deflection
% plot1 = figure(2);
timeArray = (1:Nsteps) * dt;
% plot(timeArray, y_max, 'k-');
% xlabel('Time, t [sec]');
% ylabel('Maximum Beam Deflection, y_{max} [m]');



%% Theoretical Prediction from Euler Beam Theory
d_euler = loadDistance;
l_euler = RodLength;
c_euler = min(d_euler,l_euler-d_euler);
%Implement theoretical prediction of Euler beam theory
y_max_euler = -Load*c_euler*((l_euler^2-c_euler^2)^(1.5))/(9*sqrt(3)*EI*l_euler);
%Calculate the steady state maximum vertical deflection
y_max_sim = (y_max(end));
end

