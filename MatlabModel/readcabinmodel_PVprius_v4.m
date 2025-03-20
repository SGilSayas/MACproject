clc
clear all 
close all
%% Read first:
% 1) SELECT GROUP FILE
% 2) CHANGE FILE NAMES OF FIGURES TO SAVE THEM BEFORE RUNNING THE CODE
%
% Gil-Sayas, S., 2025
% Reading Thermal Model to Simulate Vehicle's Cabin
% Working data: Toyota Prius PHEV with VIPV
% Verison read cabin model: 4 --> WORKING ON ADDING REF CYCLE
% Version cabin model function: 4 (cabinmodelfunction4)
%
%%
% Define constants for the function parameters
distance_driven = 'N/A'; % Example distance driven in km
N_Humans = 2; % Number of humans
COP = 3.5; % Coefficient of performance
TC_amb_constant = 'N/A'; 
Irr_constant = 'N/A'; 
category = 'F';
fuel = 'PETROL';

% Load data from Excel file
%filename = 'Prius_Trips_50_to_60_Minutes.xlsx'; %group 7
%filename = 'Prius_Trips_40_to_50_Minutes.xlsx'; %group 6 % only 2/3 validated
%filename = 'Prius_Trips_30_to_40_Minutes.xlsx'; %group 5 %only 2/3 validated
%filename = 'Prius_Trips_20_to_30_Minutes.xlsx'; %group 4
%filename = 'Prius_Trips_10_to_20_Minutes.xlsx'; %group 3
%filename = 'Prius_Trips_05_to_10_Minutes.xlsx'; %group 2
filename = 'Prius_Trips_01_to_05_Minutes.xlsx'; %group 1
%filename = 'trip_info_priusPV.xlsx';
data = readtable(filename);

% Extract relevant columns
trip_id = data.trip_id;
absolute_time = data.absolute_time;
TC_cabin = data.TC_cabin; 
TC_amb = data.TC_amb; 
ac_consumption_actual = data.ac_consumption_actual * 1000; % Convert to W
Irr = data.irr;

% Ensure correct data types
if iscell(trip_id)
    trip_id = string(trip_id); % Convert to string if needed
end

% Initialize containers for filtered data
filtered_trip_ids = [];

% Get unique trip IDs
unique_trip_ids = unique(trip_id);

% Filter trips based on average AC consumption
for j = 1:length(unique_trip_ids)
    % Get the current trip ID
    current_trip_id = unique_trip_ids(j);

    % Filter data for the current trip
    trip_indices = (trip_id == current_trip_id);
    trip_data = data(trip_indices, :);

    % Calculate average AC consumption for the trip
    avg_ac_consumption = mean(trip_data.ac_consumption_actual * 1000); % Convert to W

    % Calculate difference actul cabin and ambient temperature for the trip
    diff_temp_ini = abs(trip_data.TC_cabin(10)-trip_data.TC_amb(10));

    % Check if the trip meets the criteria
    if avg_ac_consumption >= 10 % && diff_temp_ini <=1
        filtered_trip_ids = [filtered_trip_ids; current_trip_id]; % Keep the trip
    end
end

% Filter the main dataset to only include valid trips
data = data(ismember(data.trip_id, filtered_trip_ids), :);

% % Set Irr=0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data.Irradiance_W_m2 = ones(height(data),1).*0; %2000;

% Get filtered unique trip IDs
filtered_unique_trip_ids = unique(filtered_trip_ids);

% Initialize storage for results
all_trip_ids = [];
all_relative_time = [];
all_irradiance = [];
all_ambient_temp = [];
all_cabin_temp_actual = [];
all_cabin_temp_simulated = [];
avg_cabin_temp_actual = [];
avg_cabin_temp_simulated = [];
all_ac_actual = [];
all_ac_simulated = [];
avg_ac_simulated = [];
all_tot_trip_CO2_g_km = [];
avg_ac_actual = [];
avg_ac_simulated = [];
trip_data.irr = [];

timestep = 0.5;

% Loop through each valid trip
for t = 1:length(filtered_unique_trip_ids)
    % Get the current trip ID
    current_trip_id = filtered_unique_trip_ids(t);

    % Filter data for the current trip (using the filtered data)
    trip_data = data(data.trip_id == current_trip_id, :);
    
    % Get the initial temp in the trip
    TC_cabin=trip_data.TC_cabin;
    TC_cabin_ini_vector = (TC_cabin(find(trip_data.TC_cabin,10)));
    TC_cabin_ini = TC_cabin_ini_vector(end);

    % Calculate Total_time for this trip (based on trip duration in seconds)
    Total_measured_time_s = height(trip_data)*timestep - timestep; % Total rows in the trip data
    
    % Create relative time for each trip
    relative_time = 0:timestep:height(trip_data)*timestep-timestep;

    % Round time to be used in cabin model function
    Total_simulation_time = height(trip_data);
    Total_simulation_time = round(Total_simulation_time);

    % Call the cabin model function
    [Tcabin, Qcompressor, compressor_av_W, comp_Wh_WLTP, CO2_total] = cabinmodelfunction4(Total_simulation_time,...
        distance_driven,N_Humans,COP,TC_cabin_ini,TC_amb_constant,Irr_constant,trip_data.TC_amb,trip_data.irr,category,fuel);

    % Store results
    cabin_temp_simulated = Tcabin - 273.15; % Simulated cabin temperature
    ac_consumption_simulated = Qcompressor; % Simulated AC consumption
    av_trip_MAC_W = repmat(compressor_av_W, Total_simulation_time, 1); % Average trip power
    tot_trip_CO2_g_km = repmat(CO2_total, Total_simulation_time, 1); % Total trip CO2 emissions
    
    % Calculate average actual and simulated AC consumption for this trip
    avg_ac_actual_trip = mean(trip_data.ac_consumption_actual * 1000); % Convert to W
    avg_ac_actual_trip = ones(height(trip_data),1).*avg_ac_actual_trip;
    % avg_cabin_temp_simulated_trip = mean(cabin_temp_simulated);
    % avg_cabin_temp_actual_trip = mean(trip_data.TC_cabin);

    % Ensure cabin_temp_simulated and other vectors are column vectors
    cabin_temp_simulated = cabin_temp_simulated(:); % Convert to column vector
    ac_consumption_simulated = ac_consumption_simulated(:); % Convert to column vector
    relative_time = relative_time(:); % Convert to column vector
    % av_trip_MAC_W = av_trip_MAC_W(:); % Convert to column vector
    % tot_trip_CO2_g_km = tot_trip_CO2_g_km(:); % Convert to column vector

    % Append results to the master arrays
    all_trip_ids = [all_trip_ids; trip_data.trip_id]; % Already columnar
    all_relative_time = [all_relative_time; relative_time]; % Already columnar
    all_irradiance = [all_irradiance; trip_data.irr]; % Already columnar
    all_ambient_temp = [all_ambient_temp; trip_data.TC_amb]; % Already columnar
    all_cabin_temp_actual = [all_cabin_temp_actual; trip_data.TC_cabin]; % Already columnar
    all_cabin_temp_simulated = [all_cabin_temp_simulated; cabin_temp_simulated];
    % avg_cabin_temp_simulated = [avg_cabin_temp_simulated; avg_cabin_temp_simulated_trip];
    % avg_cabin_temp_actual = [avg_cabin_temp_actual; avg_cabin_temp_actual_trip];
    all_ac_actual = [all_ac_actual; trip_data.ac_consumption_actual * 1000]; % Convert to W
    all_ac_simulated = [all_ac_simulated; ac_consumption_simulated];
    avg_ac_actual = [avg_ac_actual; avg_ac_actual_trip];
    avg_ac_simulated = [avg_ac_simulated; av_trip_MAC_W];
    all_tot_trip_CO2_g_km = [all_tot_trip_CO2_g_km; tot_trip_CO2_g_km];
end

% Save results to a new Excel file
output_table = table(all_trip_ids, all_relative_time, all_irradiance, all_ambient_temp, ...
    all_cabin_temp_actual, all_cabin_temp_simulated, all_ac_actual, ...
    all_ac_simulated, avg_ac_simulated, all_tot_trip_CO2_g_km, ...
    'VariableNames', {'Trip_ID','Absolute time','Irradiance','Ambient Temperature',...
    'Actual_Cabin_Temp','Simulated_Cabin_Temp','Actual_AC_Consumption', ...
    'Simulated_AC_Consumption','Avg_Trip_MAC_W','Total_Trip_CO2_g_km'});

%%
% Ensure all actual and simulated arrays are column vectors
all_cabin_temp_actual = all_cabin_temp_actual(:); % Convert to column vector
all_cabin_temp_simulated = all_cabin_temp_simulated(:); % Convert to column vector
all_ac_actual = all_ac_actual(:); % Convert to column vector
all_ac_simulated = all_ac_simulated(:); % Convert to column vector

% Check sizes
if length(all_cabin_temp_actual) ~= length(all_cabin_temp_simulated)
    error('Cabin temperature arrays (actual and simulated) must have the same length.');
end

if length(all_ac_actual) ~= length(all_ac_simulated)
    error('AC compressor arrays (actual and simulated) must have the same length.');
end
%%
% Cabin temperature Regression line & R2 value
coeff_temp = polyfit(all_cabin_temp_actual, all_cabin_temp_simulated, 1); % Coeff for first degree polynomial equation
y_fit_temp = polyval(coeff_temp, all_cabin_temp_actual);
SST_temp = sum((all_cabin_temp_simulated - mean(all_cabin_temp_simulated)).^2); % Total sum of squares
SSR_temp = sum((y_fit_temp - all_cabin_temp_simulated).^2); % Residual sum of squares
R_squared_temp = 1 - (SSR_temp / SST_temp);
% Cabin temperature Errors
MAE_cabin_temp = mean(abs(all_cabin_temp_actual - all_cabin_temp_simulated));
RMSE_cabin_temp = sqrt(mean((all_cabin_temp_actual - all_cabin_temp_simulated).^2));

% AC consumption Regression line and R2 value
coeff_ac = polyfit(avg_ac_actual, avg_ac_simulated, 1); % Linear fit (y = mx + c)
y_fit_ac = polyval(coeff_ac, avg_ac_actual); % Fitted values
SST_ac = sum((avg_ac_simulated - mean(avg_ac_simulated)).^2); % Total sum of squares
SSR_ac = sum((y_fit_ac - avg_ac_simulated).^2); % Residual sum of squares
R_squared_ac = 1 - (SSR_ac / SST_ac);
% AC Conumption Errors
MAE_ac = mean(abs(avg_ac_actual - avg_ac_simulated));
RMSE_ac = sqrt(mean((avg_ac_actual - avg_ac_simulated).^2));

% Display error metrics
fprintf('CABIN TEMPERATURE ERROR METRICS (ºC):\nMean Absolute Error (MAE): %.f\nRoot Mean Squared Error (RMSE): %.f\nR-squared: %.4f\n',...
    MAE_cabin_temp,RMSE_cabin_temp,R_squared_temp);
fprintf('AC COMPRESSOR CONSUMPTION ERROR METRICS (W):\nMean Absolute Error (MAE): %.f\nRoot Mean Squared Error (RMSE): %.f\nR-squared: %.4f\n',...
    MAE_ac,RMSE_ac,R_squared_ac);

%% PLOTS
figure(1);
scatter(all_cabin_temp_actual, all_cabin_temp_simulated, 'o');
hold on;
% Regression line
plot(all_cabin_temp_actual, y_fit_temp, 'r-', 'LineWidth', 1.5);
% Add a reference line y = x for comparison
plot([min(all_cabin_temp_actual) max(all_cabin_temp_actual)], [min(all_cabin_temp_actual) max(all_cabin_temp_actual)], 'k--');
hold on;
legend('Trips', 'Regression Line', 'y = x Line', 'Location', 'Best');
xlabel('Actual Cabin Temperature (°C)');
ylabel('Simulated Cabin Temperature (°C)');
    %title('Actual vs Simulated Cabin Temperature');
% Add R² and regression formula text
text(mean(all_cabin_temp_actual), mean(all_cabin_temp_simulated), ...
    sprintf('R² = %.3f', R_squared_temp), 'FontSize', 10, 'FontSize', 10,'HorizontalAlignment','right') 
grid on;
hold off;

% % % Plot actual vs simulated AC compressor consumption
% % figure;
% % scatter(all_ac_consumption_actual, all_ac_consumption_simulated, 'b', 'filled');
% % hold on;
% % coeff_ac = polyfit(all_ac_consumption_actual, all_ac_consumption_simulated, 1);
% % x_fit_ac = linspace(min(all_ac_consumption_actual), max(all_ac_consumption_actual), 100);
% % y_fit_ac = polyval(coeff_ac, x_fit_ac);
% % plot(x_fit_ac, y_fit_ac, 'r-', 'LineWidth', 1.5);
% % text(mean(all_ac_consumption_actual), mean(all_ac_consumption_simulated), ...
% %     sprintf('R² = %.3f', ac_errors.R_squared), 'FontSize', 12, 'Color', 'k');
% % xlabel('Actual AC Compressor Consumption (W)');
% % ylabel('Simulated AC Compressor Consumption (W)');
% % title('Actual vs Simulated AC Compressor Consumption');
% % grid on;
% % hold off;

% AC consumption plot
figure(2);
scatter(avg_ac_actual, avg_ac_simulated, 'o');
hold on;
% Plot regression line
plot(avg_ac_actual, y_fit_ac, 'r-', 'LineWidth', 1.5);
% Add a reference line y = x for comparison
plot([min(avg_ac_actual) max(avg_ac_actual)], [min(avg_ac_actual) max(avg_ac_actual)], 'k--');
xlabel('Average Actual AC Consumption (W)');
ylabel('Average Simulated AC Consumption (W)');
    % title('Comparison of Actual and Simulated Average AC Consumption');
legend('Trips', 'Regression Line', 'y = x Line', 'Location', 'Best');
grid on;
% Add R2 and regression formula text
text_x = min(avg_ac_actual) + 0.1 * range(avg_ac_actual);
text_y = max(avg_ac_simulated) - 0.1 * range(avg_ac_simulated);
text(text_x, text_y, sprintf('y = %.2fx + %.2f\nR² = %.2f', coeff_ac(1), coeff_ac(2), R_squared_ac), ...
     'FontSize', 10, 'BackgroundColor', 'w', 'Margin', 3);
hold off

% Time series of Temps and Irradiance
figure(3)
    % title('Teperatures and Irradiance for All Trips')
xlabel('Time (s)')
yyaxis left
grayColor = [.7 .7 .7];
plot(all_irradiance,'Color', grayColor)
set(gca, 'XColor','k', 'YColor','k')
ylabel('Irradiance (W/m^2)')
hold on
yyaxis right
ylabel('Temperature (ºC)')
plot(all_cabin_temp_simulated,'r--', 'LineWidth', 1)
plot(all_cabin_temp_actual,'b-.', 'LineWidth', 1)
plot(all_ambient_temp,'k-', 'LineWidth', 1)
hold on
grid on
set(gca, 'XColor','k', 'YColor','k')
legend('Irradiance','T cabin simulated','T cabin actual','T ambient actual')
hold off

% output_filename = 'filteredmodelresultsprius_group1.xlsx';
% writetable(output_table, output_filename);
% disp(['Filtered and simulation results saved to ', output_filename]);

saveas(figure(1),'group1_temps.png')
saveas(figure(2),'group1_compr.png')
saveas(figure(3),'group1_timetemps.png')