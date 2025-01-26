clear all;
clc;
close all;

%% Environment Parameters
Zseabed = 3500;       % 3.5 km : seabed depth
Zinv = 700;           % Sound speed interface depth
gamma_0 = 0.0163;     % Gradient below Zinv
gamma_1 = -0.026;     % Gradient above Zinv
C0 = 1450;            % velocity at the water surface
Cinv = C0 + gamma_1*Zinv;  % Sound speed at Zinv
Z0 = 500;             % Transmitter depth
theta0 = 2*pi/180;    % Initial angle 

%% Simulation Parameters
num_rays = 20;         % Single ray
dt = 0.01;            % Time step size
ts = 2000;           % Total time steps
thetas = linspace(-theta0, theta0, num_rays);
%% Ray Tracing Simulation
figure;
hold on;

for j = 1:num_rays
    % Initialize containers
    X = zeros(1, ts);
    Z = zeros(1, ts);
    C = zeros(1, ts);
    theta = zeros(1, ts);
    
    % Initial values
    X(1) = 0;                    
    Z(1) = Z0;   

    if Z(1) <= Zinv
        C(1) = C0 + gamma_1*Z(1);   % Above Zinv
    else
        C(1) = Cinv + gamma_0*(Z(1) - Zinv); % Below Zinv
    end
    
    theta(1) = thetas(j);%;   theta0        % Initial angle (one angle or indexed from set of many angles for multi beam tracing...)
    
    for i = 2:ts
        X(i) = X(i-1) + C(i-1)*dt*cos(theta(i-1));
        Z(i) = Z(i-1) + C(i-1)*dt*sin(theta(i-1));
        
        % Update sound speed based on depth     
        if Z(i) <= Zinv
            C(i) = C0 + gamma_1*Z(i);              % Above Zinv
        elseif Z(i) > Zinv
            C(i) = Cinv + gamma_0*(Z(i) - Zinv);  % Below Zinv        
        end
        K = cos(theta(i-1)) * C(i)/C(i-1);
        % Check for boundary reflection
        if K >= 1 || Z(i) < 0 || Z(i) >= Zseabed
            % Reverse angle and clamp to boundary
            theta(i) = -theta(i-2);
            Z(i) = Z(i-1); 
            C(i) = C(i-2);  
        else
            theta(i) = acos(K) * sign(theta(i-1));
        end
    end
    
    % Plot the ray's trajectory
    plot(X, Z, 'b');
end
% Plot the receiver location (10 km, 1000 m depth)
plot(10000, 1000, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'none');
%% Plot Formatting
xlabel('Horizontal Distance (m)');
ylabel('Depth (m)');
title('Acoustic Ray Tracing in Mediterranean Sea');
set(gca, 'YDir', 'reverse');
axis tight;
grid on;
hold off;