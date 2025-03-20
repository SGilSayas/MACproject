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
filename = 'trip_info_priusPV_sample.xlsx';
data = readtable(filename);

% Extract relevant columns
trip_id = data.trip_id;
date_time_raw = data.date_time;
cabin_temp_actual = data.T_cabin_C; 
TC_cell = data.T_amb_C; 
ac_consumption_actual = data.AC_comp_kW * 1000; % Convert to W
Irr = data.Irradiance_W_m2;

% Get unique trip IDs
unique_trip_ids = unique(trip_id);

% Loop through each trip
for t = 1:length(unique_trip_ids)
    % Get the current trip ID
    current_trip_id = unique_trip_ids(t);
    
    % Filter data for the current trip
    trip_indices = (trip_id == current_trip_id);
    trip_data = data(trip_indices, :);
    
    % Calculate Total_time for this trip (based on trip duration in seconds)
    trip_start_time = trip_data.date_time(1);
    trip_end_time = trip_data.date_time(end);
    duration_trip_s = seconds(trip_end_time - trip_start_time); % Duration in seconds
    Total_time = length(trip_data.date_time); % in 1/2 s
 
    [Tcabin,Qcompressor,compressor_av_W,comp_Wh_WLTP,CO2_total] = cabinmodelfunction4( ...
            Total_time, distance_driven, N_Humans, COP, TC_cell_constant, ...
            Irr_constant, TC_cell, Irr, category, fuel);

    % Store results
    cabin_temp_simulated = Tcabin-273.15; % Simulated cabin temperature
    ac_consumption_simulated = Qcompressor; % Simulated AC consumption
    trip_MAC_consumption_W = compressor_av_W; % Av power in the trip
    trip_emissions_g_km = CO2_total; % total CO2 emissions in the trip
    
    % Create columns with scalars
    av_trip_MAC_W = ones(1,Total_time).*trip_MAC_consumption_W;
    tot_trip_CO2_g_km = ones(1,Total_time).*trip_emissions_g_km;

    % Transpose values
    cabin_temp_simulated = cabin_temp_simulated';
    ac_consumption_simulated = ac_consumption_simulated';
    av_trip_MAC_W = av_trip_MAC_W';
    tot_trip_CO2_g_km = tot_trip_CO2_g_km';
end

% Save results to a new Excel file
output_table = table(trip_id, date_time_raw, cabin_temp_actual, cabin_temp_simulated, ...
                     ac_consumption_actual, ac_consumption_simulated,av_trip_MAC_W,tot_trip_CO2_g_km);
output_filename = 'cabin_prius_results_trips_sample.xlsx';
writetable(output_table, output_filename);
disp(['Simulation results saved to ', output_filename]);