% Vehicle cabin model, SGS, 2024
% Application case: MAC consumption for 12 cars at min/max EU temp (Jaime)
clc
clear all 
close all

% Input data files
filename1 = 'models_for_Susana_mod_SML.xlsx';
filename2 = 'table_Country_stats_TTminTmmax.xlsx';

% Extract data
table1 = readtable(filename1, 'ReadVariableNames', true);
Ft = table1.MS_Ft; % Fuel type
category = table1.category;
table2 = readtable(filename2, 'ReadVariableNames', true);
MS = table2.MS; % Member state
temp_max = table2.max_temp;
temp_min = table2.min_temp;

% Define constants for the function parameters
Total_time = 1800; % Total time in seconds
distance_driven = 'N/A'; % Example distance driven in km
N_Humans = 2; % Number of humans
COP = 3.5; % Coefficient of performance
Irr_constant = 0; 
TC_cell = 'N/A'; 
Irr = 'N/A';

% Initialize results
results = table();

% Loop through each vehicle and member state
rowIndex = 1; % Index for the results table
for i = 1:height(table1)
    for j = 1:height(table2)
        % Retrieve specific data for the current vehicle and state
        vehicle_category = category{i};
        fuel_type = Ft{i}; % Use fuel type directly from the input file
        member_state = MS{j};
        Tmin = temp_min(j);
        Tmax = temp_max(j);
        
        % Run for T_cell_constant = Tmin
        T_cell_constant = Tmin;
        [compressor_av_W_Tmin, comp_Wh_WLTP_Tmin, CO2_total_Tmin] = cabinmodelfunction2(...
            Total_time, distance_driven, N_Humans, COP, ...
            T_cell_constant, Irr_constant, TC_cell, Irr, vehicle_category, fuel_type);
        
        % Run for T_cell_constant = Tmax
        T_cell_constant = Tmax;
        [compressor_av_W_Tmax, comp_Wh_WLTP_Tmax, CO2_total_Tmax] = cabinmodelfunction2(...
            Total_time, distance_driven, N_Humans, COP, ...
            T_cell_constant, Irr_constant, TC_cell, Irr, vehicle_category, fuel_type);
        
        % Store the results in the table
        results(rowIndex, :) = {vehicle_category, fuel_type, member_state, Tmin, Tmax, ...
                                compressor_av_W_Tmin, comp_Wh_WLTP_Tmin, CO2_total_Tmin, ...
                                compressor_av_W_Tmax, comp_Wh_WLTP_Tmax, CO2_total_Tmax};
        rowIndex = rowIndex + 1; % Increment the row index
    end
end

% Define column names for the results table
results.Properties.VariableNames = {'Vehicle_Category', 'Fuel', 'Member_State', 'Tmin', 'Tmax', ...
                                    'Compressor_Avg_W_Tmin', 'Comp_Wh_WLTP_Tmin', 'CO2_Total_Tmin', ...
                                    'Compressor_Avg_W_Tmax', 'Comp_Wh_WLTP_Tmax', 'CO2_Total_Tmax'};

% Write results to an Excel file
output_filename = 'cabin_model_results_Tmin_Tmax_new.xlsx';
writetable(results, output_filename);

disp(['Results written to ', output_filename]);