% Vehicle cabin model, SGS, 2024
% Application case: PRIUS PV roof real-world driving
clc
clear all 
close all

% Define constants for the function parameters
distance_driven = 'N/A'; % Example distance driven in km
N_Humans = 2; % Number of humans
COP = 3.5; % Coefficient of performance
TC_cell_constant = 'N/A'; 
Irr_constant = 'N/A'; 
category = 'E';
fuel = 'PETROL';

% Load data from Excel file
filename = 'trip_info_priusPV_sample.xlsx';
data = readtable(filename);

% Extract relevant columns
trip_id = data.trip_id;
datetime_raw = data.date_time;
cabin_temp_actual = data.T_cabin_C; 
ambient_temp = data.T_amb_C; 
ac_consumption_actual = data.AC_comp_kW * 1000; % Convert to W
irradiance = data.Irradiance_W_m2;

% Ensure datetime column is properly formatted
data.datetime = datetime(datetime_raw, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Ensure correct data types
if iscell(data.trip_id)
    data.trip_id = string(data.trip_id); % Convert to string if needed
end

% Get unique trip IDs
unique_trip_ids = unique(trip_id);

% Initialize variables for simulation results
cabin_temp_simulated = zeros(height(data), 1);
ac_consumption_simulated = zeros(height(data), 1);

% Loop through each trip
for t = 1:length(unique_trip_ids)
    % Get the current trip ID
    current_trip_id = unique_trip_ids(t);
    
    % Filter data for the current trip
    trip_indices = (trip_id == current_trip_id);
    trip_data = data(trip_indices, :);
    
    % Calculate Total_time for this trip (based on trip duration in seconds)
    trip_start_time = trip_data.datetime(1);
    trip_end_time = trip_data.datetime(end);
    Total_time = seconds(trip_end_time - trip_start_time); % Duration in seconds
    
    % Round to the nearest integer greater than or equal to Total time
    Total_time = ceil(Total_time)*2;

    % Loop through each data point in the current trip
    for i = 1:height(trip_data)
        % Inputs to the cabinmodelfunction3
        TC_cell = trip_data.T_amb_C(i); % Ambient temperature
        Irr = trip_data.Irradiance_W_m2(i); % Irradiance
        
        % Call cabin model function
        [Tcabin, compressor_av_W, comp_Wh_WLTP, CO2_total] = cabinmodelfunction3( ...
            Total_time, distance_driven, N_Humans, COP, TC_cell_constant, ...
            Irr_constant, TC_cell, Irr, category, fuel);
        
        % Store results
        cabin_temp_simulated(trip_indices(i)) = Tcabin; % Simulated cabin temperature
        ac_consumption_simulated(trip_indices(i)) = compressor_av_W; % Simulated AC consumption
    end
end

% Compare actual vs simulated cabin temperatures and air conditioner consumption
figure;
subplot(2, 1, 1); % First subplot
plot(cabin_temp_actual, 'b-', 'DisplayName', 'Actual Cabin Temp');
hold on;
plot(cabin_temp_simulated, 'r--', 'DisplayName', 'Simulated Cabin Temp');
xlabel('Time Point');
ylabel('Cabin Temperature (°C)');
legend;
title('Comparison of Actual and Simulated Cabin Temperatures');
hold off;
subplot(2, 1, 2); % Second subplot
plot(ac_consumption_actual, 'b-', 'DisplayName', 'Actual AC Consumption');
hold on;
plot(ac_consumption_simulated, 'r--', 'DisplayName', 'Simulated AC Consumption');
xlabel('Time Point');
ylabel('AC Consumption (W)');
legend;
title('Comparison of Actual and Simulated AC Consumption');
hold off;

% Calculate error metrics for cabin temperature
mae_cabin_temp = mean(abs(cabin_temp_actual - cabin_temp_simulated)); % Mean Absolute Error
rmse_cabin_temp = sqrt(mean((cabin_temp_actual - cabin_temp_simulated).^2)); % Root Mean Squared Error
r_squared_cabin_temp = 1 - sum((cabin_temp_actual - cabin_temp_simulated).^2) / sum((cabin_temp_actual - mean(cabin_temp_actual)).^2);

% Calculate error metrics for air conditioner consumption
mae_ac = mean(abs(ac_consumption_actual - ac_consumption_simulated)); % Mean Absolute Error
rmse_ac = sqrt(mean((ac_consumption_actual - ac_consumption_simulated).^2)); % Root Mean Squared Error
r_squared_ac = 1 - sum((ac_consumption_actual - ac_consumption_simulated).^2) / sum((ac_consumption_actual - mean(ac_consumption_actual)).^2);

% Display error metrics
disp('Error Metrics for Cabin Temperature:');
disp(['Mean Absolute Error (MAE): ', num2str(mae_cabin_temp), ' °C']);
disp(['Root Mean Squared Error (RMSE): ', num2str(rmse_cabin_temp), ' °C']);
disp(['R-squared (R²): ', num2str(r_squared_cabin_temp)]);

disp('Error Metrics for Air Conditioner Consumption:');
disp(['Mean Absolute Error (MAE): ', num2str(mae_ac), ' W']);
disp(['Root Mean Squared Error (RMSE): ', num2str(rmse_ac), ' W']);
disp(['R-squared (R²): ', num2str(r_squared_ac)]);

% Save results to a new Excel file
output_table = table(trip_id, datetime_raw, cabin_temp_actual, cabin_temp_simulated, ...
                     ac_consumption_actual, ac_consumption_simulated);
output_filename = 'cabin_prius_results_trips.xlsx';
writetable(output_table, output_filename);
disp(['Simulation results saved to ', output_filename]);


%%
% % % Define constants for the function parameters
% Total_time = 1800; % Total time in seconds
% distance_driven = 'N/A'; % Example distance driven in km
% N_Humans = 2; % Number of humans
% COP = 3.5; % Coefficient of performance
% TC_cell_constant = 'N/A'; 
% Irr_constant = 'N/A'; 
% category = 'E';
% fuel = 'PETROL';
% 
% % Load data from Excel file
% filename = 'trip_info_priusPV.xlsx';
% data = readtable(filename);
% 
% % Extract relevant columns
% trip_id = data.trip_id;
% datetime = data.datetime;
% cabin_temp_actual = data.T_cabin_C; 
% ambient_temp = data.T_amb_C; 
% ac_consumption_actual = data.AC_comp_kW;
% ac_consumption_actual = ac_consumption_actual*1000; % W
% irradiance = data.Irradiance_W_m2;
% 
% % Initialize variables for simulation results
% num_data_points = height(data);
% cabin_temp_simulated = zeros(num_data_points, 1);
% ac_consumption_simulated = zeros(num_data_points, 1);
% 
% % Loop through each data point and run the model
% for i = 1:num_data_points
%     % Inputs to the cabinmodelfunction3
%     TC_cell = ambient_temp(i); % Vector
%     Irr = irradiance(i); % Vector
% 
%     % Call cabin model function
%     [Tcabin,compressor_av_W,comp_Wh_WLTP,CO2_total] = cabinmodelfunction3(Total_time,distance_driven,N_Humans,COP,TC_cell_constant,Irr_constant,TC_cell,Irr,category,fuel);
% 
%     % Store results
%     cabin_temp_simulated(i) = Tcabin; % Simulated cabin temperature
%     ac_consumption_simulated(i) = compressor_av_W; % Simulated AC consumption
% end
% 
% Compare actual vs simulated cabin temperatures
figure;
subplot(2, 1, 1); % First subplot
plot(cabin_temp_actual, 'b-', 'DisplayName', 'Actual Cabin Temp');
hold on;
plot(cabin_temp_simulated, 'r--', 'DisplayName', 'Simulated Cabin Temp');
xlabel('Time Point');
ylabel('Cabin Temperature (°C)');
legend;
title('Comparison of Actual and Simulated Cabin Temperatures');
hold off;

% Compare actual vs simulated air conditioner consumption
subplot(2, 1, 2); % Second subplot
plot(ac_consumption_actual, 'b-', 'DisplayName', 'Actual AC Consumption');
hold on;
plot(ac_consumption_simulated, 'r--', 'DisplayName', 'Simulated AC Consumption');
xlabel('Time Point');
ylabel('AC Consumption (kW)');
legend;
title('Comparison of Actual and Simulated AC Consumption');
hold off;

% Calculate error metrics for cabin temperature
mae_cabin_temp = mean(abs(cabin_temp_actual - cabin_temp_simulated)); % Mean Absolute Error
rmse_cabin_temp = sqrt(mean((cabin_temp_actual - cabin_temp_simulated).^2)); % Root Mean Squared Error
r_squared_cabin_temp = 1 - sum((cabin_temp_actual - cabin_temp_simulated).^2) / sum((cabin_temp_actual - mean(cabin_temp_actual)).^2);

% Calculate error metrics for air conditioner consumption
mae_ac = mean(abs(ac_consumption_actual - ac_consumption_simulated)); % Mean Absolute Error
rmse_ac = sqrt(mean((ac_consumption_actual - ac_consumption_simulated).^2)); % Root Mean Squared Error
r_squared_ac = 1 - sum((ac_consumption_actual - ac_consumption_simulated).^2) / sum((ac_consumption_actual - mean(ac_consumption_actual)).^2);

% Display error metrics
disp('Error Metrics for Cabin Temperature:');
disp(['Mean Absolute Error (MAE): ', num2str(mae_cabin_temp), ' °C']);
disp(['Root Mean Squared Error (RMSE): ', num2str(rmse_cabin_temp), ' °C']);
disp(['R-squared (R²): ', num2str(r_squared_cabin_temp)]);

disp('Error Metrics for Air Conditioner Consumption:');
disp(['Mean Absolute Error (MAE): ', num2str(mae_ac), ' kW']);
disp(['Root Mean Squared Error (RMSE): ', num2str(rmse_ac), ' kW']);
disp(['R-squared (R²): ', num2str(r_squared_ac)]);
% 
% % Save results to a new Excel file
% output_table = table(trip_id, datetime, cabin_temp_actual, cabin_temp_simulated, ...
%                      ac_consumption_actual, ac_consumption_simulated);
% output_filename = 'cabin_prius_results.xlsx';
% writetable(output_table, output_filename);
% disp(['Simulation results saved to ', output_filename]);
