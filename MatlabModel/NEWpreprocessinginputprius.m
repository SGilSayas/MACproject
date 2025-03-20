clc;
clear all;
close all;

% Load data from Excel file
filename = 'trip_info_prius_original.csv';
data = readtable(filename, 'NumHeaderLines', 1);

% Extract relevant columns
trip_id = data{:,1}; % Extract Trip ID as an array
times_input = data{:,2}; % Extract time as an array
TC_cabin = data{:,3}; % Cabin temperature
TC_amb = data{:,4}; % Cell/Ambient temperature
ac_consumption_actual = data{:,5}; % AC consumption
irr = data{:,6}; % Irradiance

% Ensure trip_id is a numeric or string array
if iscell(trip_id)
    trip_id = string(trip_id); % Convert to string if needed
end

% Ensure times_input is a valid string array
if iscell(times_input)
    times_input = string(times_input); % Convert cell array to string array
end

% Convert time to duration (mm:ss.s format)
time_durations = duration(times_input, 'InputFormat', 'mm:ss.S');

% Get unique Trip IDs
unique_trip_ids = unique(trip_id);

% Initialize array to store trip durations
trip_durations_inseconds = zeros(length(unique_trip_ids), 1);

% Calculate time duration of each trip
for i = 1:length(unique_trip_ids)
    % Get data corresponding to the current Trip ID
    current_trip_id = unique_trip_ids(i);
    mask = trip_id == current_trip_id; % Logical mask for the current trip
    trip_times = time_durations(mask); % Filter times for the current trip
    
    % Calculate duration of the trip
    trip_duration = max(trip_times) - min(trip_times);
    
    % Convert to seconds and store
    trip_durations_inseconds(i) = seconds(trip_duration);
end

% % Display trip duration in seconds
% disp('Trip durations in seconds:');
% for i = 1:length(unique_trip_ids)
%     fprintf('Trip ID: %d, Duration (seconds): %.2f\n', unique_trip_ids(i), trip_durations_inseconds(i));
% end

% Categorize trips based on their duration
group_data = {...
    "Prius_Trips_01_to_05_Minutes.xlsx", [], []; % File name, trip IDs, and indices
    "Prius_Trips_05_to_10_Minutes.xlsx", [], [];
    "Prius_Trips_10_to_20_Minutes.xlsx", [], [];
    "Prius_Trips_20_to_30_Minutes.xlsx", [], [];
    "Prius_Trips_30_to_40_Minutes.xlsx", [], [];
    "Prius_Trips_40_to_50_Minutes.xlsx", [], [];
    "Prius_Trips_50_to_60_Minutes.xlsx", [], [];
    "Prius_Trips_Over_60_Minutes.xlsx", [], []};

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

% Initialize cell array to store output tables for each group
group_output_tables = cell(size(group_data, 1), 1);

% Save the input data of each trip group to separate Excel files
for g = 1:size(group_data, 1)
    if ~isempty(group_data{g, 3})
        % Get indices of trips in the group
        group_indices = group_data{g, 3};
        
        % Filter input data for the trips in this group
        filtered_data = data(ismember(trip_id, group_data{g, 2}), :);
        
        % Set the new headers (column names) for filtered_data
        new_headers = {'trip_id', 'times_input', 'TC_cabin', 'TC_amb', 'ac_consumption_actual', 'irr'};
        filtered_data.Properties.VariableNames = new_headers;
        
        % Convert times_input to duration (mm:ss.s format) for each row
        times_in_duration = duration(filtered_data.times_input, 'InputFormat', 'mm:ss.S');
        
        % Add the new 'trip_time' column as a time series starting at 0 for each trip
        trip_time_column = (0:length(filtered_data.times_input)-1)' * 0.5; % Time in seconds starting from 0
        
        % Create new table with the required columns
        output_table = table(...
            filtered_data.trip_id, ...
            filtered_data.times_input, ...
            trip_time_column, ...
            filtered_data.TC_cabin, ...
            filtered_data.TC_amb, ...
            filtered_data.ac_consumption_actual, ...
            filtered_data.irr, ...
            'VariableNames', {'trip_id', 'times', 'absolute_time', 'TC_cabin', 'TC_amb', 'ac_consumption_actual', 'irr'});
        
        % Filter output_table for rows where ac_consumption_actual > 0.05  
        filtered_output_table = output_table(output_table.ac_consumption_actual > 0.05, :);
        output_table = filtered_output_table;

        % Write the output table to the cell array
        group_output_tables{g} = output_table;

        % Write filtered data with new 'trip_time' column to Excel
        % writetable(output_table, group_data{g, 1});       
        
    end
end
disp('Trips saved to Excel files successfully.');

% Display numbers of trips per duration group
fprintf('Trips by Duration:\n');
duration_labels = {
    '1-5 minutes', '5-10 minutes', '10-20 minutes', ...
    '20-30 minutes', '30-40 minutes', '40-50 minutes', ...
    '50-60 minutes', 'Over 60 minutes'};
for g = 1:size(group_data, 1)
    fprintf('%s: %d trips\n', duration_labels{g}, counters(g));
end

