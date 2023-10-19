%%%%%%%%% Problem 1 %%%%%%%%%
%%%%%%%%% MAE 263F %%%%%%%%%%
%%%%%%%%% Homework 1 %%%%%%%%
%%%%%%%%% Logan Hall %%%%%%%%
%%%%%% October 18 2023 %%%%%%
% Inspired by code from UCLA MAE 263F Lecture 3 by Professor Khalid Jawed


%% Reset MATLAB Workspace
clear;
close all;
clc;


%% Assignment parameters
dt_implicit = 10^(-2); %dt for implicit situation [sec]
dt_explicit = 10^(-5); %dt for explicit situation [sec]
dt_explicit_larger = 10^(-2);%dt larger for explicit situation [sec]
dt_implicit_larger = 10^(-1); %dt larger for implicit situation [sec]

% Radius of spheres
R1 = 0.005; %[m]
R2 = 0.025; %[m]
R3 = 0.005; %[m]



%% Simulate Implicitly
terminal_velocity_implicit = problem1_implicit_simulation(dt_implicit, R1, R2, R3);

%% Simulate Explicitly
terminal_velocity_explicit = problem1_explicit_simulation(dt_explicit, R1, R2, R3);

%% Simulate Implicitly with all Radii (R1, R2, R3) the same

terminal_velocity_implicit_Rconst = problem1_implicit_simulation(dt_implicit, R1, R1, R1);

%% Simulate Explicitly with different time step size
%Running this next line will result in an invalid simulation and create graphs that are not correct
terminal_velocity_explicit = problem1_explicit_simulation(dt_explicit_larger, R1, R2, R3);

%% Simulate Implicitly with different time step size
terminal_velocity_implicit = problem1_implicit_simulation(dt_implicit_larger, R1, R2, R3);




%% Implicit Simulation Function
function [terminal_velocity_implicit] = problem1_implicit_simulation(dt, R1, R2, R3)

N = 3;

ndof = N * 2; % Number of degrees of freedom

% Rod length
RodLength = .1; %[m] given bar length

%[m] discrete length
deltaL = RodLength / (N-1);

% Density
rho_metal = 7000; % kg/m^3
rho_f = 1000; % kg/m^3
rho = rho_metal - rho_f; % kg/m^3

% Rod radius
r0 = 0.001; %[m]

% Young's modulus
Y = 1e9; % Using Y instead of E to avoid ambiguity

% Gravity
g = 9.8; % m/s^2

% Viscosity
visc = 1000; % Pa-s

% Total time
totalTime = 10; % seconds

% Utility quantities
EI = Y * pi * r0^4 / 4; % [Nm^2] - bending stiffness
EA = Y * pi * r0^2; % [N]

% Geometry - initial configuration
nodes = zeros(N, 2);
for c = 1:N % Loop over all the nodes
    nodes(c,1) = (c-1) * deltaL; % x coordinates
    nodes(c,2) = 0;
end



% Mass matrix
M = zeros(ndof,ndof);
%Calculate mass at each node
M(1,1) = 4/3 * pi * R1^3 * rho_metal;
M(2,2) = 4/3 * pi * R1^3 * rho_metal;
M(3,3) = 4/3 * pi * R2^3 * rho_metal;
M(4,4) = 4/3 * pi * R2^3 * rho_metal;
M(5,5) = 4/3 * pi * R3^3 * rho_metal;
M(6,6) = 4/3 * pi * R3^3 * rho_metal;

% Viscous damping matrix
C = zeros(6,6);
%Calculate viscous damping terms
C1 = 6 * pi * visc * R1;
C2 = 6 * pi * visc * R2;
C3 = 6 * pi * visc * R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;

% Weight Vector, W
W = zeros(ndof,1);
%Calculate weight and only assign in y-direction
W(2) = -4/3 * pi * R1^3 * rho * g;
W(4) = -4/3 * pi * R2^3 * rho * g;
W(6) = -4/3 * pi * R3^3 * rho * g;

% Initial DOF vector
q0 = zeros(ndof,1);
for c=1:N % loop over nodes
    q0(2*c-1) = nodes(c,1); % x coordinate
    q0(2*c) = nodes(c,2); % y coordinate
end

% Initialize velocity vector
u0 = zeros(ndof,1); % Old Velocity (initial velocity)

%Initialize turning angle
turning_angle = zeros(ndof,1);

% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector

% Number of time steps
Nsteps = floor( totalTime / dt )+1;
all_mid_y = zeros( Nsteps, 1); % y-position of R2 all_mid_v = zeros( Nsteps, 1); % y-velocity of R2
all_mid_y(1) = q(4); %Store initial position
all_mid_v(1) = u(4);  %store initial velocity
stored_vals = zeros(Nsteps, ndof); %store all dof vector at each time step
stored_vals(1,:) = q; %Store initial dof

% Tolerance
tol = EI / RodLength^2 * 1e-3; %small enough force that can be neglected


% Time marching scheme
for c=2:Nsteps

    %Print time at each step
    fprintf('Time = %f\n', (c-1) * dt );

    q = q0; % Guess New DOFs
    u = (q-q0)/dt; % New Velocity

    % Newton Raphson
    err = 10 * tol;
    while err > tol
        % Inertia
        f = M / dt * ( (q-q0) / dt - u );
        J = M / dt^2;

        %
        % Elastic forces
        %
        % Linear spring 1 between nodes 1 and 2
        %Define x_k, y_k, y_(k+1), and x_(k+1)
        xk = q(1);
        yk = q(2);
        xkp1 = q(3);
        ykp1 = q(4);
        %Calculate derivative of stretching energy E_k^s
        dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
        %Calculate the hessian of the stretching energy E_k^s
        dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
        %Update Jacobian and f
        f(1:4) = f(1:4) + dF;
        J(1:4,1:4) = J(1:4,1:4) + dJ;

        % Linear spring 2 between nodes 2 and 3
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        %Calculate derivative of stretching energy E_k^s
        dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
        %Calculate the hessian of the stretching energy E_k^s
        dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
        %Update Jacobian and f
        f(3:6) = f(3:6) + dF;
        J(3:6,3:6) = J(3:6,3:6) + dJ;

        % Bending spring between nodes 1, 2, and 3
        %Define x_k, y_k, y_(k+1), x_(k+1), y_(k-1), x_(k-1)
        xkm1 = q(1);
        ykm1 = q(2);
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        curvature0 = 0;
        %Calculate the derivative of bending energy E_k^b
        dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
        %Calculate the hessian of the bending energy E_k^b
        dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
        %Update Jacobian and f
        f(1:6) = f(1:6) + dF;
        J(1:6,1:6) = J(1:6,1:6) + dJ;

        % Viscous force
        f = f + C * (q - q0)/dt;
        J = J + C / dt;

        % Weight
        f = f - W;

        % Update
        q = q - J \ f;

        err = sum( abs(f) );
    end

    % Update
    u = (q - q0) / dt; % Velocity

    %Define vectors and magnitudes between points to determine turning
    %angle
    v1 = [q0(3)-q0(1), q0(4)-q0(2)];
    v2 = [q0(5)-q0(3), q0(6)-q0(4)];
    mag_v1 = sqrt(v1(1)^2+v1(2)^2);
    mag_v2 = sqrt(v2(1)^2+v2(2)^2);

    %Determine turning angle
    turning_angle(c) = rad2deg(acos((dot(v1,v2))/(mag_v1*mag_v2)));


    q0 = q; % Old position



    


  

    %Live plot the simulation
    figure(1);
    plot( q(1:2:end), q(2:2:end), 'ro-'); axis equal
    axis equal
    xlabel('x [meter]')
    ylabel('y [meter]')
    drawnow

    % Store mid position, mid velocity, and all dof
    all_mid_y(c) = q(4);
    all_mid_v(c) = u(4);
    stored_vals(c,:) = q;
    
end

% Calculate terminal velocity
terminal_velocity_implicit = all_mid_v(end);

%% Plot velocity of R2 along y-axis
fig2 = figure(2);
timeArray = (1:Nsteps) * dt-dt;
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity (Along y-axis) of R_2, v [meter/sec]');
grid on;
title('Velocity of R_2 (Along y-axis)')
saveas(fig2,'Problem1.1_velocity.png');


%% Plot position of R2 along y-axis
fig3 = figure(3);
timeArray = (1:Nsteps) * dt;
plot(timeArray, stored_vals(:,6), 'k-');
xlabel('Time, t [sec]');
ylabel('Position (Along y-axis) of R_2, y_2 [meter]');
grid on;
title('Position of R_2 (Along y-axis)')
saveas(fig3,'Problem1.1_position.png');


%% Plot x-y coordinates at specific times

plot_times = [0, 0.01, 0.05, 0.10, 1.0, 10.0];
plot_index = plot_times/dt+1;

for ii = 1:length(plot_index)
    plot_index = find(timeArray == plot_times(ii));
    h = figure(ii+3);
    plot(stored_vals(plot_index, 1:2:end), stored_vals(plot_index,2:2:end), 'ro-');
    axis equal
    xlabel('x [meter]')
    ylabel('y [meter]')
    title('Rigid Spheres and Elastic Beam Falling in Viscous Flow: t = ' + string(plot_times(ii)) + ' [seconds]')
    grid on
    xlim([0,0.1]);
    ylim([-0.07,0.01]);
    saveas(h,sprintf('Problem1.1_%d.png',ii));
end

%% Plot turning angle
fig10 = figure(10);
timeArray = (1:Nsteps) * dt;
plot(timeArray, turning_angle, 'k-');
xlabel('Time, t [sec]');
ylabel('Turning Angle [degrees]');
grid on;
title('Turning Angle')
saveas(fig10,'Problem1.1_theta.png');
end




















function [terminal_velocity_explicit] = problem1_explicit_simulation(dt, R1, R2, R3)
%% Physical parameters
% Number of vertices
N = 3;

ndof = N * 2; % Number of degrees of freedom



% Rod length
RodLength = .1; % meter

% Discrete length
deltaL = RodLength / (N-1);



% Density
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;

% Rod radius
r0 = 0.001;

% Young's modulus
Y = 1e9; % Using Y instead of E to avoid ambiguity

% Gravity
g = 9.8; % m/s^2

% Viscosity
visc = 1000; % Pa-s

% Total time
totalTime = 10; % seconds

% Utility quantities
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;

% Geometry
nodes = zeros(N, 2);
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
    nodes(c,2) = 0;
end



% Mass matrix
M = zeros(ndof,ndof);
M(1,1) = 4/3 * pi * R1^3 * rho_metal;
M(2,2) = 4/3 * pi * R1^3 * rho_metal;
M(3,3) = 4/3 * pi * R2^3 * rho_metal;
M(4,4) = 4/3 * pi * R2^3 * rho_metal;
M(5,5) = 4/3 * pi * R3^3 * rho_metal;
M(6,6) = 4/3 * pi * R3^3 * rho_metal;

% Viscous damping matrix
C = zeros(6,6);
C1 = 6 * pi * visc * R1;
C2 = 6 * pi * visc * R2;
C3 = 6 * pi * visc * R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;

% Weight Vector, W
W = zeros(ndof,1);
W(2) = -4/3 * pi * R1^3 * rho * g;
W(4) = -4/3 * pi * R2^3 * rho * g;
W(6) = -4/3 * pi * R3^3 * rho * g;

% Initial DOF vector
q0 = zeros(ndof,1);
for c=1:N
    q0(2*c-1) = nodes(c,1); % x coordinate
    q0(2*c) = nodes(c,2); % y coordinate
end

% Initialize velocity vector
u0 = zeros(ndof,1); % Old Velocity (initial velocity)

% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector

% Number of time steps
Nsteps = round( totalTime / dt );
all_mid_y = zeros( Nsteps, 1); % y-position of R2 all_mid_v = zeros( Nsteps, 1); % y-velocity of R2
all_mid_y(1) = q(4);
all_mid_v(1) = u(4);
stored_q = zeros(Nsteps, ndof); %store all dof vector at each time step
stored_q(1,:) = q;

% Tolerance
tol = EI / RodLength^2 * 1e-3; %small enough force that can be neglected


% Time marching scheme
for c=2:Nsteps
    fprintf('Time = %f\n', (c-1) * dt );
    
    
    % Elastic forces
    % Linear spring 1 between nodes 1 and 2
    xk_12 = q0(1);
    yk_12 = q0(2);
    xkp1_12 = q0(3);
    ykp1_12 = q0(4);
    dF_12 = gradEs(xk_12, yk_12, xkp1_12, ykp1_12, deltaL, EA);
    
    
    % Linear spring 2 between nodes 2 and 3
    xk_23 = q0(3);
    yk_23 = q0(4);
    xkp1_23 = q0(5);
    ykp1_23 = q0(6);
    dF_23 = gradEs(xk_23, yk_23, xkp1_23, ykp1_23, deltaL, EA);
    
    
    % Bending spring between nodes 1, 2, and 3
    xkm1 = q0(1);
    ykm1 = q0(2);
    xk = q0(3);
    yk = q0(4);
    xkp1 = q0(5);
    ykp1 = q0(6);
    curvature0 = 0;
    dF_bend = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);

    F_elastic = zeros(6,1);
    F_elastic(1:4) = F_elastic(1:4) + dF_12;
    F_elastic(3:6) = F_elastic(3:6) + dF_23;
    F_elastic(1:6) = F_elastic(1:6) + dF_bend;

    F_external = C*u0 - W;

    % Explicit Method
    % Calculating dof directly via algebraic manipulation
    for i = 1:2*N
        q(i) = q0(i) + u0(i)*dt - ((dt^2)/M(i,i))*(F_elastic(i)+F_external(i));
    end

    u = (q - q0)/dt;
    % Update
    u0 = u; % Velocity
    q0 = q; % Old position
  

    %The active plotting for the explicit method has been commented out in
    %order to allow the code to run in a reasonable timeframe. This is an
    %example of how inefficient the explicit method is. Uncomment this
    %section of code to see how the simulation breaks with a smaller
    %timestep.
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(1);
    % plot( q(1:2:end), q(2:2:end), 'ro-');
    % axis equal
    % xlabel('x [meter]')
    % ylabel('y [meter]')
    % drawnow
    %%%%%%%%%%%%%%%%%%%%%%%%

    % Store values for plotting
    all_mid_y(c) = q(4);
    all_mid_v(c) = u(4);
    stored_q(c,:) = q;
    
end


%%
%Define array of time
timeArray = (1:Nsteps) * dt;
%Calculate terminal velocity
terminal_velocity_explicit = all_mid_v(end);
plot_times = [0, 0.01, 0.05, 0.10, 1.0, 10.0];
plot_index = plot_times/dt+1;

%Plot at requested times
for ii = 1:length(plot_index)
    plot_index = find(timeArray == plot_times(ii));
    h = figure(ii+3);
    plot(stored_q(plot_index, 1:2:end), stored_q(plot_index,2:2:end), 'ro-');
    %Format plot
    axis equal
    xlabel('x [meter]')
    ylabel('y [meter]')
    title('Rigid Spheres and Elastic Beam Falling in Viscous Flow: t = ' + string(plot_times(ii)) + ' [seconds]')
    grid on
    xlim([0,0.1]);
    ylim([-0.07,0.01]);
    saveas(h,sprintf('Problem1.1_%d.png',ii));
end

%%
figure(2);
%Plot mid node vertical velocity vs time
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity (vertical) of mid-node, v [meter/sec]');
end