clc
clear all 
close all

% Define constants for the function parameters
distance_driven = 'N/A'; % Example distance driven in km
N_Humans = 2; % Number of humans
COP = 3.5; % Coefficient of performance
TC_cell_constant = 22; 
Irr_constant = 0; 
category = 'E';
fuel = 'PETROL';

% Load data from Excel file
filename = 'trip_info_priusPV.xlsx';
data = readtable(filename);

% Extract relevant columns
trip_id = data.trip_id;
date_time_raw = data.date_time;
cabin_temp_actual = data.T_cabin_C; 
TC_cell = data.T_amb_C; 
ac_consumption_actual = data.AC_comp_kW * 1000; % Convert to W
Irr = data.Irradiance_W_m2;

% Ensure datetime column is properly formatted
data.date_time = datetime(date_time_raw, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Ensure correct data types
if iscell(data.trip_id)
    data.trip_id = string(data.trip_id); % Convert to string if needed
end

% Initialize storage for results
all_trip_ids = [];
all_date_times = [];
all_cabin_temp_actual = [];
all_cabin_temp_simulated = [];
all_ac_consumption_actual = [];
all_ac_consumption_simulated = [];
all_av_trip_MAC_W = [];
all_tot_trip_CO2_g_km = [];

% Get unique trip IDs
unique_trip_ids = unique(trip_id);

% Loop through each trip
for t = 1:length(unique_trip_ids)
    % Get the current trip ID
    current_trip_id = unique_trip_ids(t);
    
    % Filter data for the current trip
    trip_indices = (trip_id == current_trip_id);
    trip_data = data(trip_indices, :);
    
    % Calculate Total_time for this trip
    trip_start_time = trip_data.date_time(1);
    trip_end_time = trip_data.date_time(end);
    duration_trip_s = seconds(trip_end_time - trip_start_time); % Duration in seconds
    Total_time = height(trip_data); % Number of timesteps

    % Extract inputs for the cabin model function
    TC_cell_trip = trip_data.T_amb_C;
    Irr_trip = trip_data.Irradiance_W_m2;

    % Call the cabin model function
    [Tcabin, Qcompressor, compressor_av_W, comp_Wh_WLTP, CO2_total] = cabinmodelfunction4( ...
        Total_time, distance_driven, N_Humans, COP, TC_cell_constant, ...
        Irr_constant, TC_cell_trip, Irr_trip, category, fuel);

    % Process and store results
    cabin_temp_simulated = Tcabin - 273.15; % Convert from Kelvin to Celsius
    ac_consumption_simulated = Qcompressor; % Simulated AC consumption
    av_trip_MAC_W = repmat(compressor_av_W, Total_time, 1); % Average trip power
    tot_trip_CO2_g_km = repmat(CO2_total, Total_time, 1); % Total trip CO2 emissions

    % Ensure cabin_temp_simulated and other vectors are column vectors
    cabin_temp_simulated = cabin_temp_simulated(:); % Convert to column vector
    ac_consumption_simulated = ac_consumption_simulated(:); % Convert to column vector
    av_trip_MAC_W = av_trip_MAC_W(:); % Convert to column vector
    tot_trip_CO2_g_km = tot_trip_CO2_g_km(:); % Convert to column vector
    
    % Append results to the master arrays
    all_trip_ids = [all_trip_ids; trip_data.trip_id]; % Already columnar
    all_date_times = [all_date_times; trip_data.date_time]; % Already columnar
    all_cabin_temp_actual = [all_cabin_temp_actual; trip_data.T_cabin_C]; % Already columnar
    all_cabin_temp_simulated = [all_cabin_temp_simulated; cabin_temp_simulated];
    all_ac_consumption_actual = [all_ac_consumption_actual; trip_data.AC_comp_kW * 1000]; % Convert to W
    all_ac_consumption_simulated = [all_ac_consumption_simulated; ac_consumption_simulated];
    all_av_trip_MAC_W = [all_av_trip_MAC_W; av_trip_MAC_W];
    all_tot_trip_CO2_g_km = [all_tot_trip_CO2_g_km; tot_trip_CO2_g_km];
end

% Save results to a new Excel file
output_table = table(all_trip_ids, all_date_times, all_cabin_temp_actual, ...
                     all_cabin_temp_simulated, all_ac_consumption_actual, ...
                     all_ac_consumption_simulated, all_av_trip_MAC_W, ...
                     all_tot_trip_CO2_g_km, ...
                     'VariableNames', {'Trip_ID', 'Date_Time', 'Actual_Cabin_Temp', ...
                                       'Simulated_Cabin_Temp', 'Actual_AC_Consumption', ...
                                       'Simulated_AC_Consumption', 'Avg_Trip_MAC_W', ...
                                       'Total_Trip_CO2_g_km'});

output_filename = 'results_cabin_prius_trips.xlsx';
writetable(output_table, output_filename);
disp(['Simulation results saved to ', output_filename]);

%%
% Ensure all actual and simulated arrays are column vectors
all_cabin_temp_actual = all_cabin_temp_actual(:); % Convert to column vector
all_cabin_temp_simulated = all_cabin_temp_simulated(:); % Convert to column vector
all_ac_consumption_actual = all_ac_consumption_actual(:); % Convert to column vector
all_ac_consumption_simulated = all_ac_consumption_simulated(:); % Convert to column vector

% Check sizes
if length(all_cabin_temp_actual) ~= length(all_cabin_temp_simulated)
    error('Cabin temperature arrays (actual and simulated) must have the same length.');
end

if length(all_ac_consumption_actual) ~= length(all_ac_consumption_simulated)
    error('AC compressor arrays (actual and simulated) must have the same length.');
end

% Define error metrics
calculate_errors = @(actual, simulated) struct( ...
    'MAE', mean(abs(actual - simulated)), ...
    'RMSE', sqrt(mean((actual - simulated).^2)), ...
    'R_squared', 1 - sum((actual - simulated).^2) / sum((actual - mean(actual)).^2));

% Calculate error metrics
temp_errors = calculate_errors(all_cabin_temp_actual, all_cabin_temp_simulated);
ac_errors = calculate_errors(all_ac_consumption_actual, all_ac_consumption_simulated);

% Display error metrics
disp('Error Metrics for Cabin Temperature:');
disp(temp_errors);

disp('Error Metrics for AC Compressor Consumption:');
disp(ac_errors);

% Plot actual vs simulated cabin temperature
figure;
scatter(all_cabin_temp_actual, all_cabin_temp_simulated, 'b', 'filled');
hold on;
coeff_temp = polyfit(all_cabin_temp_actual, all_cabin_temp_simulated, 1);
x_fit_temp = linspace(min(all_cabin_temp_actual), max(all_cabin_temp_actual), 100);
y_fit_temp = polyval(coeff_temp, x_fit_temp);
plot(x_fit_temp, y_fit_temp, 'r-', 'LineWidth', 1.5);
text(mean(all_cabin_temp_actual), mean(all_cabin_temp_simulated), ...
    sprintf('R² = %.3f', temp_errors.R_squared), 'FontSize', 12, 'Color', 'k');
xlabel('Actual Cabin Temperature (°C)');
ylabel('Simulated Cabin Temperature (°C)');
title('Actual vs Simulated Cabin Temperature');
grid on;
hold off;

% Plot actual vs simulated AC compressor consumption
figure;
scatter(all_ac_consumption_actual, all_ac_consumption_simulated, 'b', 'filled');
hold on;
coeff_ac = polyfit(all_ac_consumption_actual, all_ac_consumption_simulated, 1);
x_fit_ac = linspace(min(all_ac_consumption_actual), max(all_ac_consumption_actual), 100);
y_fit_ac = polyval(coeff_ac, x_fit_ac);
plot(x_fit_ac, y_fit_ac, 'r-', 'LineWidth', 1.5);
text(mean(all_ac_consumption_actual), mean(all_ac_consumption_simulated), ...
    sprintf('R² = %.3f', ac_errors.R_squared), 'FontSize', 12, 'Color', 'k');
xlabel('Actual AC Compressor Consumption (W)');
ylabel('Simulated AC Compressor Consumption (W)');
title('Actual vs Simulated AC Compressor Consumption');
grid on;
hold off;
