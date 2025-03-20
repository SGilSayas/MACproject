% MATLAB Code: Cabin Temperature Simulation with AC

% Initial Parameters
T_outside = 35; % Outside temperature (°C)
T_AC = 20;      % Air conditioning output temperature (°C)
T_cabin_front = 30; % Initial front cabin temperature (°C)
T_cabin_back = 32;  % Initial back cabin temperature (°C)
heat_transfer_rate = 0.1; % Heat transfer rate with outside (°C per time step)
cooling_rate = 0.2; % Cooling rate provided by AC (°C per time step)
time_steps = 100; % Number of simulation steps
dt = 1; % Time step duration

% Initialize arrays to store temperatures over time
T_cabin_front_history = zeros(1, time_steps);
T_cabin_back_history = zeros(1, time_steps);

% Simulation loop
for t = 1:time_steps
    % Calculate heat exchange with the outside environment
    heat_gain_front = heat_transfer_rate * (T_outside - T_cabin_front) * dt;
    heat_gain_back = heat_transfer_rate * (T_outside - T_cabin_back) * dt;

    % Cooling effect of the AC
    cooling_front = cooling_rate * (T_cabin_front - T_AC) * dt;
    cooling_back = cooling_rate * (T_cabin_back - T_AC) * dt;

    % Update cabin temperatures
    T_cabin_front = T_cabin_front + heat_gain_front - cooling_front;
    T_cabin_back = T_cabin_back + heat_gain_back - cooling_back;

    % Store temperatures for plotting
    T_cabin_front_history(t) = T_cabin_front;
    T_cabin_back_history(t) = T_cabin_back;
end

% Plot results
figure;
plot(1:time_steps, T_cabin_front_history, '-r', 'DisplayName', 'Front Cabin Temperature');
hold on;
plot(1:time_steps, T_cabin_back_history, '-b', 'DisplayName', 'Back Cabin Temperature');
xlabel('Time Steps');
ylabel('Temperature (°C)');
title('Cabin Temperature Simulation');
legend;
grid on;
