%%%%%%%%% Problem 2 %%%%%%%%%
%%%%%%%%% MAE 263F %%%%%%%%%%
%%%%%%%%% Homework 1 %%%%%%%%
%%%%%%%%% Logan Hall %%%%%%%%
%%%%%% October 18 2023 %%%%%%
% Inspired by code from UCLA MAE 263F Lecture 3 by Professor Khalid Jawed





%% Reset MATLAB Workspace
clear;
close all;
clc;






%% Problem 2.1 and 2.2 - Plotting Vertical Position and Vertical Velocity
%Define parameters
N_1 = 21;
dt_1 = 10^(-2);%[s]
[midNode1, timeArray1, all_mid_v1, stored_q, terminal_velocity1] = problem2_simulation(N_1,dt_1);

% Plot middle node downward velocity
fig1 = figure(1);
plot(timeArray1, all_mid_v1, 'k-');
%format plot
xlabel('Time, t [sec]');
ylabel('Velocity (vertical) of Middle Node, v [meter/second]');
grid on;
title('Velocity of Middle Node (Along y-axis)')
saveas(fig1,'Problem2.1_velocity.png');

% Plot position of R2 along y-axis
fig2 = figure(2);
plot(timeArray1, stored_q(:,midNode1*2), 'k-');
%format plot
xlabel('Time, t [sec]');
ylabel('Position (Along y-axis) of Middle Node, y_{middle node} [meter]');
grid on;
title('Position of R_2 (Along y-axis)')
saveas(fig2,'Problem2.1_position.png');

% Plot final deformed shape of the beam
fig3 = figure(3);
plot(stored_q(end,1:2:end), stored_q(end,2:2:end), 'r-o');
%format plot
axis equal
xlabel('x [meter]');
ylabel('y [meter]');
grid on;
title('Final Deformed Shape in 2D Space at t = 10 sec')
saveas(fig3,'Problem2.2_finalshape.png');



%% Problem 2.3 Part 1 - Spacial Discretization
%Define parameters
N_2 = 3:2:33;
dt_2 = 10^(-2);

%initialize midNode and terminal_velocity
midNode2 = zeros(1,length(dt_2));
terminal_velocity2 = zeros(1,length(dt_2));

%for loop to calculate mid node and terminal velocity at each value of N
for z = 1:length(N_2)
    disp('Calculating for Node: ')
    disp(N_2(z))
    [midNode2(z), ~, ~, ~, terminal_velocity2(z)] = problem2_simulation(N_2(z),dt_2);
end

%Plot N vs terminal velocity
fig4 = figure(4);
plot(N_2, terminal_velocity2, 'k-');

%Format plot
xlabel('Number of Nodes, N');
ylabel('Terminal Velocity [meter/second]');
grid on;
title('Terminal Velocity vs Number of Nodes')
saveas(fig4,'Problem2.3_Nodes.png');


%% Problem 2.3 Part 2 - Temporal Discretization
%Define parameters
N_3 = 21;%Number of nodes
dt_3 = [.25 .1 .075 .050 .025 .010 .005]; %[s], time step

%initialize midNode and terminal_velocity
midNode3 = zeros(1,length(dt_3));
terminal_velocity3 = zeros(1,length(dt_3));

%for loop to calculate mid node and terminal velocity at each value of dt
for z = 1:length(dt_3)
    disp(dt_3(z))
    [midNode3(z), ~, ~, ~, terminal_velocity3(z)] = problem2_simulation(N_3,dt_3(z));
end

%Plot Temporal Discretization
fig5 = figure(5);
plot(dt_3, terminal_velocity3, 'k-');
%format plot
xlabel('Time Step, \Delta t');
ylabel('Terminal Velocity [meter/second]');
grid on;
title('Terminal Velocity vs \Delta t')
saveas(fig5,'Problem2.3_Steps.png');

%% Function to simulate for problem 2
function [midNode, timeArray, all_mid_v, stored_q, terminal_velocity] = problem2_simulation(N, dt)


ndof = N * 2; % number of degrees of freedom

RodLength = 0.1; %[m] given bar length
deltaL = RodLength / (N-1);%[m] discrete length

% Radii of spheres
R = zeros(N,1); % Vector of size N - Radius of N nodes
R(:) = deltaL/10; %[m]
midNode = (N+1)/2;
R(midNode) = 0.025;

% Density
rho_metal = 7000; % kg/m^3
rho_f = 1000; % fluid % kg/m^3
rho = rho_metal - rho_f; % kg/m^3

r0 = 0.001; % meter - rod radius
Y = 1e9; % Young's modulus (Y instead of E for clarity)

g = 9.8; % m/s^2 - gravity

visc = 1000; % pa-s

totalTime = 50; % second - total simulation time

% Utility parameter
ne = N - 1; % number of edges
EI = Y * pi * r0^4 / 4; % Nm^2 - bending stiffness
EA = Y * pi * r0^2; % Newton

% Geometry - initial configuration
nodes = zeros(N,2);
for c=1:N % Loop over all the nodes
    nodes(c,1) = (c-1) * deltaL; % x coordinates
    nodes(c,2) = 0;
end

% Mass, M
M = zeros(ndof, ndof);
for k=1:N
    M(2*k-1, 2*k-1) = 4/3*pi*R(k)^3*rho_metal; % Mass for x_k
    M(2*k, 2*k) = M(2*k-1, 2*k-1); % Mass for y_k
end

% Viscous damping matrix, C
C = zeros(ndof,ndof);
for k=1:N
    C(2*k-1, 2*k-1) = 6 * pi * visc * R(k);
    C(2*k, 2*k) = C(2*k-1, 2*k-1);
end
% Weight vector, W
W = zeros(ndof, 1);
for k=1:N
    W(2*k-1) = 0; % weight along x is zero
    W(2*k) = -4/3*pi*R(k)^3*rho*g;
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

% Storage for y-velocity of the middle node
all_mid_v = zeros(Nsteps, 1);
stored_q = zeros(Nsteps, ndof); %store all dof vector at each time step
stored_q(1,:) = q0;

for c = 2:Nsteps

    %fprintf('Time = %f\n', (c-1) * dt);

    % Guess
    q = q0; % New DOFs are initialized to be equal to old DOFs
    % Newton Raphson
    err = 10 * tol;
    while err > tol
        f = M / dt * ( (q-q0)/dt - u0 );
        J = M / dt^2;

        %
        % Elastic forces
        %
        % Linear spring
        for k=1:N-1
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
            f(ind) = f(ind) + dF;
            J(ind,ind) = J(ind,ind) + dJ;
        end

        % Bending spring
        for k=2:N-1
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

        % Viscous force
        f = f + C * (q-q0) / dt;
        J = J + C / dt;

        % Weight
        f = f - W;

        % Update
        q = q - J \ f;
        err = sum ( abs(f) );
    end

    % New velocity
    u = (q - q0) / dt;

    % Store some information
    all_mid_v(c) = u(2*midNode);
    stored_q(c,:) = q;

    % The drawnow plot function has been commented out to increase the
    % speed the code can run. It may be uncommented to visualize the
    % simulation.
    
    %Plot
    % fig1 = figure(1);
    % plot( q(1:2:end), q(2:2:end), 'ro-');
    % axis equal
    % 
    % xlabel('x [meter]');
    % ylabel('y [meter]');
    % drawnow

    % Update (new becomes old)
    q0 = q;
    u0 = u;
end
%grid on;
%title('Final Deformed Beam Shape at: t = 50 seconds')
%saveas(fig1,'Problem2.2.png');


timeArray = (1:Nsteps) * dt;
terminal_velocity = all_mid_v(end);
end


