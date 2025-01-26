clear all;
clc;
close all;

%% Environment Parameters
Zseabed = 3500;       % Seabed depth
Zinv = 700;           % Sound speed interface depth
gamma_0 = 0.0163;     % Gradient below Zinv
gamma_1 = -0.026;     % Gradient above Zinv
C0 = 1450;            % Sound speed at surface
Cinv = C0 + gamma_1*Zinv;  
theta0_deg = 5;      % Max angle in degrees
theta0 = theta0_deg * pi/180; 

%% Simulation Parameters
depths_to_test = [100,200,300,400,500,600,700, 800, 900,1000,1100,1150,1180,1200,1300,1400,1500,1600,1700,1800];  % Transmitter depths
num_rays = 5000;                   % Angular resolution
dt = 0.01;                        % Time step
ts = 2000;                        % Total time steps
tolerance = 0.4;                    % Stricter depth tolerance (meters)

%% Storage for Results
results = [];  % Format: [depth, theta_deg, arrival_time]

%% Main Simulation Loop
for depth_idx = 1:length(depths_to_test)
    Z0 = depths_to_test(depth_idx);
    thetas = linspace(-theta0, theta0, num_rays);
    
    for j = 1:num_rays
        [X, Z, C, theta] = deal(zeros(1, ts));
        X(1) = 0; Z(1) = Z0;
        
        % Initialize sound speed
        if Z0 <= Zinv
            C(1) = C0 + gamma_1*Z0;
        else
            C(1) = Cinv + gamma_0*(Z0 - Zinv);
        end
        theta(1) = thetas(j);
        
        min_distance = Inf;
        best_time = Inf;
        
        for i = 2:ts
            % Update position
            X(i) = X(i-1) + C(i-1)*dt*cos(theta(i-1));
            Z(i) = Z(i-1) + C(i-1)*dt*sin(theta(i-1));
            
            % Update sound speed
            if Z(i) <= Zinv
                C(i) = C0 + gamma_1*Z(i);
            else
                C(i) = Cinv + gamma_0*(Z(i) - Zinv);
            end
            
            % Snell's law update (corrected reflection handling)
            K = cos(theta(i-1)) * C(i)/C(i-1);
            if K >= 1 || Z(i) < 0 || Z(i) >= Zseabed
                theta(i) = -theta(i-1);  
                Z(i) = Z(i-1);
                C(i) = C(i-1);
            else
                theta(i) = acos(K) * sign(theta(i-1));
            end
            
            % Track closest approach to receiver
            distance = sqrt((X(i)-10000)^2 + (Z(i)-1000)^2);
            if distance < min_distance
                min_distance = distance;
                best_time = i*dt;
            end
        end
        
        % Store result if within tolerance
        if min_distance <= tolerance
            results = [results; Z0, thetas(j)*180/pi, best_time];
        end
    end
end

%% Find Optimal Solution
if ~isempty(results)
    [min_time, idx] = min(results(:,3));
    best_depth = results(idx,1);
    best_theta = results(idx,2);
    
    fprintf('Optimal Solution:\n');
    fprintf('Transmitter Depth = %d m\n', best_depth);
    fprintf('Emission Angle    = %.2f°\n', best_theta);
    fprintf('Arrival Time      = %.2f seconds\n', min_time);
    
    %% Plot Fastest Ray (Corrected Reflection Logic)
    Z0 = best_depth;
    theta_opt = best_theta * pi/180;
    
    % Reinitialize containers
    [X, Z, C, theta] = deal(zeros(1, ts));
    X(1) = 0; Z(1) = Z0;
    
    if Z0 <= Zinv
        C(1) = C0 + gamma_1*Z0;
    else
        C(1) = Cinv + gamma_0*(Z0 - Zinv);
    end
    theta(1) = theta_opt;
    
    for i = 2:ts
        X(i) = X(i-1) + C(i-1)*dt*cos(theta(i-1));
        Z(i) = Z(i-1) + C(i-1)*dt*sin(theta(i-1));
        
        if Z(i) <= Zinv
            C(i) = C0 + gamma_1*Z(i);
        else
            C(i) = Cinv + gamma_0*(Z(i) - Zinv);
        end
        
        K = cos(theta(i-1)) * C(i)/C(i-1);
        if K >= 1 || Z(i) < 0 || Z(i) >= Zseabed
            theta(i) = -theta(i-2);  % Corrected reflection
            Z(i) = Z(i-1);
            C(i) = C(i-2);
        else
            theta(i) = acos(K) * sign(theta(i-1));
        end
    end
    
    figure;
    plot(X, Z, 'b');
    hold on;
    plot(10000, 1000, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'none');
    xlabel('Horizontal Distance (m)');
    ylabel('Depth (m)');
    title(sprintf('Fastest Ray: Depth=%dm, Angle=%.2f°', best_depth, best_theta));
    set(gca, 'YDir', 'reverse');
    grid on;
    axis tight;
    
else
    fprintf('No valid rays reached the receiver!\n');
end