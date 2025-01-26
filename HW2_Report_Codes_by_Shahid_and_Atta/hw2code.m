clear all % Clear all variables, functions, and classes from the workspace
load('C:/Users/shahi/Downloads/Received.mat'); % Load the received signal data from the file

%% Given Constants in the question
total_time = 6.4; % Total duration of the signals in seconds
dt = 2e-04; % Sampling time interval in seconds
c = 1500; % Speed of sound in water in meters per second
initial_pos = 25; % Initial depth increment for each receiver in meters
Reflections = 10; % Number of reflections to consider in the simulation

[num_rec, signal_Length] = size(RecSig); % Get the number of receivers and the length of the signal
%%
% Plotting All the Received Signals
t = linspace(0, total_time, signal_Length); % Create a time vector for plotting

ray1 = RecSig(1,:); % Signal from receiver 1 at 25 m depth
ray5 = RecSig(5,:); % Signal from receiver 5 at 125 m depth
ray9 = RecSig(9,:); % Signal from receiver 9 at 225 m depth

figure(1) % Create a new figure for plotting all received signals
plot(t, RecSig); % Plot all signals from the 9 receivers
title('The Plot of all Received Signals'); % Set the title of the plot
xlabel('time (s)'); % Label the x-axis as time in seconds
ylabel('Amplitude (m)'); % Label the y-axis as amplitude in meters
xlim([0.40 6.18]); % Set the x-axis limits
ylim([-1.71 1.19]); % Set the y-axis limits

%% Computing the time for receiving of the signals
peaks_sig = []; % Initialize an empty array to store the time of the first peaks
for i = 1:9 % Loop through each receiver
    first_peak = max(RecSig(i,11000:11500)); % Find the maximum value (peak) in a specific range
    idx = find(RecSig(i, :) == first_peak); % Find the index of the peak in the signal
    peaks = [peaks_sig, idx * dt]; % Convert the index to time and store it
    fprintf('Time taken to reach first peak %f for Receiver %d\n', idx * dt, i); % Print the time of the first peak
end
%%
% Plotting around first peak
figure(2) % Create a new figure for plotting specific signals
plot(t, ray1); % Plot the signal from receiver 1
hold on; % Hold the current plot to overlay additional signals
plot(t, ray5); % Overlay the signal from receiver 5
plot(t, ray9); % Overlay the signal from receiver 9
title('The Plot of Signal 1, 5, and 9'); % Set the title of the plot
xlabel('time (s)'); % Label the x-axis as time in seconds
ylabel('Amplitude (m)'); % Label the y-axis as amplitude in meters
legend('Signal 1', 'Signal 5', 'Signal 9'); % Add a legend for the signals
xlim([2.1920 2.3203]); % Set the x-axis limits
ylim([-2 1.5]); % Set the y-axis limits

 % Plotting the 5th signal
figure(3) % Create a new figure for plotting signal 5
plot(t, ray5); % Plot the signal from receiver 5
title('The Plot of 5th Signal'); % Set the title of the plot
xlabel('time (s)'); % Label the x-axis as time in seconds
ylabel('Amplitude (m)'); % Label the y-axis as amplitude in meters
xlim([2.048 2.584]); % Set the x-axis limits
ylim([-1.05 1.15]); % Set the y-axis limits


%% Estimated Depth
H_est = 358.3; % Estimated depth of the water channel in meters

% Exploration area (x,z axis)
X_Pos = linspace(980, 1080, 20); % Create a horizontal position grid (x-axis)
Z_Pos = linspace(100, 150, 12.5); % Create a depth grid (z-axis)

plotArea = zeros(length(Z_Pos), length(X_Pos)); % Initialize an array to store the results of the simulation

for i = 1:length(X_Pos) % Loop through each x-position
    for j = 1:length(Z_Pos) % Loop through each z-position
        receivedSignal = zeros(size(RecSig)); % Initialize the reconstructed signal matrix
        for k = 1:size(RecSig, 1) % Loop through each receiver
            flipSig = flip(RecSig(k, :)); % Time-reverse the signal
            out = Green(X_Pos(i), Z_Pos(j), 0, initial_pos * k, H_est, c, flipSig', dt, Reflections); % Call the Green's function
            receivedSignal(k, :) = out; % Store the reconstructed signal
        end
        receivedSignal = sum(receivedSignal, 1).^2; % Compute the squared sum of the signals
        totalSignal = sum(receivedSignal, "all"); % Compute the total energy of the signals
        totalSignal = 10 * log10(totalSignal); % Convert the total energy to dB
        plotArea(j, i) = totalSignal; % Store the result in the plot area
    end
end

%% Green's function definition
function [g] = Green(xr, zr, xs, zs, h, c, sig, dt, n_ref)
    r_even_j = @(x, z, j) sqrt((x - xs).^2 + (z - zs + (j) * h).^2); % Compute distance for even reflections
    r_odd_j = @(x, z, j) sqrt((x - xs).^2 + (z + zs - (j + 1) * h).^2); % Compute distance for odd reflections

    fs = 1 / dt; % Sampling frequency
    g = 0; % Initialize the Green's function

    for j = -n_ref:n_ref % Loop through each reflection
        if mod(j, 2) == 1 % Check if the reflection is odd
            eps = -1; % Set the reflection coefficient for odd reflections
            rj = r_odd_j(xr, zr, j); % Compute the distance for odd reflections
        else % Even reflection
            eps = 1; % Set the reflection coefficient for even reflections
            rj = r_even_j(xr, zr, j); % Compute the distance for even reflections
        end
        delta = delayseq(sig, rj / c, fs); % Apply delay to the signal based on the distance
        gj = eps / (4 * pi * rj) * delta; % Compute the contribution for the current reflection
        g = g - gj; % Accumulate the contribution
    end
end

% Plot the results of the source localization
figure(4); % Create a new figure for the localization plot
image(X_Pos, Z_Pos, plotArea, 'CDataMapping', 'scaled'); % Plot the localization heatmap
colorbar; % Add a colorbar to the plot
title('Plot of the Source Localisation'); % Set the title of the plot
xlabel('Position (m)'); % Label the x-axis as position in meters
ylabel('Depth (m)'); % Label the y-axis as depth in meters

% Plot all received signals with vertical offsets
figure(5);
offsets = (0:num_rec-1) * 2; % Define vertical offsets for each receiver
hold on;
for i = 1:num_rec
    plot(t, RecSig(i, :) + offsets(i), 'DisplayName', ['Rec ', num2str(i)]); % Add offset to each signal
end
hold off;

% Customize the plot
title('Received Signals at Different Depths'); % Title of the plot
xlabel('Time (s)'); % Label the x-axis
ylabel('Signal Amplitude'); % Label the y-axis
legend('show'); % Display legend with receiver names
grid on; % Add grid lines
ylim([-1, num_rec * 2 + 1]); % Adjust y-axis limits to fit all signals
