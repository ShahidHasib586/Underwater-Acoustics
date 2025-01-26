clear all;
clc;
close all;

%% Environment Parameters
Zseabed = 3500;      % 3.5 km : seabed depth
gamma = 0.0163;      % Gradient
C0 = 1450;           % velocity at the water surface
Z0 = 500;            % Transmitter depth
theta0 = 2.2*pi/180; % initial angle , Q1: 2 , Q2: 2.2 OR 6.4

%% Simulation Parameters
num_rays = 1;        % 
% Create some angles (IN OUR CASE it is a single beam for single ray)
thetas = theta0; 

dt = 0.01;           % Time step size
ts = 2000;%10000;           % Total time steps
% (simulation time = dt*ts = 100 sec if ts=10000)

%% Ray Tracing Simulation
figure;
hold on;

for j = 1:num_rays
    % Initialize arrays for each ray
    X = zeros(1, ts);
    Z = zeros(1, ts);
    theta = zeros(1, ts);
    C = zeros(1, ts);
    
    % Initial conditions
    X(1) = 0;                   % Initial X position
    Z(1) = Z0;                  % Initial depth
    C(1) = C0 + gamma*Z(1);     % Initial sound speed
    theta(1) = thetas(j);  % Initial angle for this ray
    
    for i = 2:ts
        % Updating position using previous step (i-1)
        X(i) = X(i-1) + C(i-1)*dt*cos(theta(i-1));
        Z(i) = Z(i-1) + C(i-1)*dt*sin(theta(i-1));
        
        % Updating sound speed
        C(i) = C0 + gamma*Z(i);
        
        % Calculate refraction parameter
        K = cos(theta(i-1)) * C(i)/C(i-1);
        
        % Check for reflection/boundary conditions
        if abs(K) > 1 || Z(i) < 0 || Z(i) >= Zseabed
            % Handle reflection
            theta(i) = -theta(i-1);  % Reverse direction
            Z(i) = Z(i-1);           % Keep previous depth
            C(i) = C(i-1);           % Keep previous sound speed
        else
            % Update angle using Snell's Law
            theta(i) = acos(K) * sign(theta(i-1));
        end
    end
    
    % Plot the ray's trajectory
    plot(X, Z, 'b');
end
% Plot the transmitter location (0 km, 500 m depth)
plot(0, 500, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'none');
% Plot the receiver location (20 km, 500 m depth)
plot(20000, 500, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'none');
%% Plot Formatting
xlabel('Horizontal Distance (m)');
ylabel('Depth (m)');
title('Acoustic Ray Tracing in Arctic Ocean');
set(gca, 'YDir', 'reverse');
axis tight;
grid on;
hold off;