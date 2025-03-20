% Data preprocessing Prius RWD 2024
clc
clear all 
close all

% Load data from Excel file
filename = 'trip_info_prius_original.csv';
data = readtable(filename,'NumHeaderLines',1); %'PreserveVariableNames',true); hdrs = data.Properties.VariableNames;

% Extract relevant columns as arrays using {}, as table using ()
trip_id = data{:,1};
times_input = data{:,2}; %data.date_time;
cabin_temp_actual = data{:,3}; %data.T_cabin_C; 
TC_cell = data{:,4};%data.T_amb_C; 
ac_consumption_actual = data{:,5}; %data.AC_comp_kW * 1000; % Convert to W
Irr = data{:,6}; %data.Irradiance_W_m2;

% Ensure trip_id is a string
if iscell(trip_id)
    trip_id = string(trip_id); % Convert to string if needed
end

% Convert time to duration (since mm:ss.s doesn't include date information)
time_durations = duration(times_input, 'InputFormat', 'mm:ss.S');

% Initialize containers for filtered data
filtered_trip_ids = [];

% % Get unique trip IDs
unique_trip_ids = unique(trip_id);

% Initialize array to store trip durations
trip_durations_inseconds = zeros(length(unique_trip_ids), 1);

% Calculate time duration of each trip
for i = 1:length(unique_trip_ids)
    % Get data corresponding to the current Trip ID
    current_trip_id = unique_trip_ids(i);
    mask = trip_id == current_trip_id;
    trip_times = time_durations(mask);
    trip_duration = max(trip_times) - min(trip_times);
    % Convert to seconds and store
    trip_durations_inseconds(i) = seconds(trip_duration);
end

% % % Display the results
% % disp('Trip durations in seconds:');
% % for i = 1:length(unique_trip_ids)
% %     fprintf('Trip ID: %d, Duration (seconds): %.2f\n', unique_trip_ids(i), trip_durations_inseconds(i));
% % end



% Categorize trips based on their duration
group_data = {...
    "Trips_1_to_5_Minutes.xlsx", [], []; % File name, trip IDs, and indices
    "Trips_5_to_10_Minutes.xlsx", [], [];
    "Trips_10_to_20_Minutes.xlsx", [], [];
    "Trips_20_to_30_Minutes.xlsx", [], [];
    "Trips_30_to_40_Minutes.xlsx", [], [];
    "Trips_40_to_50_Minutes.xlsx", [], [];
    "Trips_50_to_60_Minutes.xlsx", [], [];
    "Trips_Over_60_Minutes.xlsx", [], []};

% Initialize counters
counters = zeros(size(group_data, 1), 1);

% Categorize trips
for i = 1:length(trip_durations_inseconds)
    duration_min = trip_durations_inseconds(i) / 60; % Convert to minutes
    
    if duration_min < 1
        continue; % Skip trips shorter than 1 minute
    elseif duration_min <= 5
        idx = 1;
    elseif duration_min <= 10
        idx = 2;
    elseif duration_min <= 20
        idx = 3;
    elseif duration_min <= 30
        idx = 4;
    elseif duration_min <= 40
        idx = 5;
    elseif duration_min <= 50
        idx = 6;
    elseif duration_min <= 60
        idx = 7;
    else
        idx = 8;
    end
    
    % Append trip data to the group
    group_data{idx, 2} = [group_data{idx, 2}; unique_trip_ids(i)];
    group_data{idx, 3} = [group_data{idx, 3}; i];
    counters(idx) = counters(idx) + 1;
end

% Save the input data of each trip group to separate Excel files
for g = 1:size(group_data, 1)
    if ~isempty(group_data{g, 3})
        % Get indices of trips in the group
        group_indices = group_data{g, 3};
        
        % Filter input data for the trips in this group
        filtered_data = data(ismember(trip_id, group_data{g, 2}), :);
        
        % Convert times_input to duration (mm:ss.s format) for each row
        times_in_duration = duration(filtered_data.times_input, 'InputFormat', 'mm:ss.S');
        
        % Add the new 'trip_time' column
        trip_time_column = seconds(times_in_duration - min(times_in_duration)); % Time in seconds from 0 for each trip
        
        % Create new table with the required columns
        output_table = table(...
            filtered_data.trip_id, ...
            filtered_data.times_input, ...
            filtered_data.TC_cabin, ...
            filtered_data.TC_cell, ...
            filtered_data.ac_consumption_actual, ...
            filtered_data.irr, ...
            trip_time_column, ...
            'VariableNames', {'trip_id', 'times', 'TC_cabin', 'TC_cell', 'ac_consumption_actual', 'irr', 'trip_time'});
        
        % Write filtered data with new 'trip_time' column to Excel
        writetable(output_table, group_data{g, 1});
    end
end

% Display the results
fprintf('Trips by Duration:\n');
duration_labels = {
    '1-5 minutes', '5-10 minutes', '10-20 minutes', ...
    '20-30 minutes', '30-40 minutes', '40-50 minutes', ...
    '50-60 minutes', 'Over 60 minutes'};
for g = 1:size(group_data, 1)
    fprintf('%s: %d trips\n', duration_labels{g}, counters(g));
end

disp('Trips saved to Excel files successfully.');









%%
% % Filter trips based on average AC consumption
% for j = 1:length(unique_trip_ids)
%     % Get the current trip ID
%     current_trip_id = unique_trip_ids(j);
% 
%     % Filter data for the current trip
%     trip_indices = (trip_id == current_trip_id);
%     trip_data = data(trip_indices, :);
% 
%     % Calculate average AC consumption for the trip
%     avg_ac_consumption = mean(trip_data.AC_comp_kW * 1000); % Convert to W
% 
%     % Calculate difference actul cabin and ambient temperature for the trip
%     diff_temp_ini = abs(trip_data.T_cabin_C(10)-trip_data.T_amb_C(10));
% 
%     % Check if the trip meets the criteria
%     if avg_ac_consumption >= 100 && diff_temp_ini <=1
%         filtered_trip_ids = [filtered_trip_ids; current_trip_id]; % Keep the trip
%     end
% end
% 
% % Filter the main dataset to only include valid trips
% data = data(ismember(data.trip_id, filtered_trip_ids), :);

% % Set Irr=0 data.Irradiance_W_m2 = ones(height(data),1).*0; %2000;

% Get filtered unique trip IDs
% filtered_unique_trip_ids = unique(filtered_trip_ids);

%% Save data
% save filtered_prius_trips data??