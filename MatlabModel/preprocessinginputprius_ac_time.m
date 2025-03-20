clc;
clear all;
close all;

% Load data from Excel file
filename = 'trip_info_prius_original.csv';
data = readtable(filename, 'NumHeaderLines', 1);

% Extract data
trip_id = data{:,1}; % Extract Trip ID as an array
times_input = data{:,2}; % Extract time as an array
TC_cabin = data{:,3}; % Cabin temperature
TC_amb = data{:,4}; % Cell/Ambient temperature
ac_consumption_actual = data{:,5}; % AC consumption
irr = data{:,6}; % Irradiance

% Filter output_table for rows where ac_consumption_actual > 0.65 kW (error of 1A with 650V, max EM V)
ac_consumption_filtered = ac_consumption_actual;
ac_consumption_filtered(ac_consumption_filtered < 0.65) = 0;

% Ensure trip_id is a numeric or string array
if iscell(trip_id)
    trip_id = string(trip_id); % Convert to string if needed
end

% Ensure times_input is a valid string array
if iscell(times_input)
    times_input = string(times_input); % Convert cell array to string array
end

% Get unique Trip IDs
unique_trip_ids = unique(trip_id);

% Initialize invalid tirps table
invalid_trip_ids = [];

% Convert time to duration (mm:ss.s format)
time_durations = duration(times_input, 'InputFormat', 'mm:ss.S');

% Initialize array to store trip durations
trip_durations_inseconds = zeros(length(unique_trip_ids), 1);

% Calculate trips time duration and filter trips 
for i = 1:length(unique_trip_ids)
    % Get data corresponding to the current Trip ID
    current_trip_id = unique_trip_ids(i);
    mask = trip_id == current_trip_id; % Logical mask for the current trip
    trip_times = time_durations(mask); % Filter times for the current trip
    trip_ac_filtered = ac_consumption_filtered(mask);
    trip_TC_cabin = TC_cabin(mask);

    % Calculate duration of the trip
    trip_duration = max(trip_times) - min(trip_times);
    
    % Convert to seconds and store
    trip_durations_inseconds(i) = seconds(trip_duration);

    % Skip trips with less than 1 second duration
    if trip_durations_inseconds <= 60
        invalid_trip_ids = [invalid_trip_ids; current_trip_id];
        continue;
    end
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
% group_output_tables = cell(size(group_data, 1), 1);

% Initialize struct array to store output tables for each group
group_tables = struct('group_id', [], 'output_table', []);

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
        
        % Calculate total time where ac_consumption_actual > 0.65 for each trip
        unique_trips = unique(output_table.trip_id);
        active_time_per_trip = zeros(length(unique_trips), 1);
        for t = 1:length(unique_trips)
            trip_mask = output_table.trip_id == unique_trips(t);
            trip_data = output_table(trip_mask, :);
            
            % Calculate total time based on absolute_time
            if ~isempty(trip_data)
                active_time_per_trip(t) = max(trip_data.absolute_time) - min(trip_data.absolute_time);
            else
                active_time_per_trip(t) = 0;
            end
        end
        
        % Add the calculated active times to the filtered output table
        output_table.active_time = zeros(height(output_table), 1);
        for t = 1:length(unique_trips)
            trip_mask = output_table.trip_id == unique_trips(t);
            output_table.active_time(trip_mask) = active_time_per_trip(t);
        end
        
        % Save the output table for the current group in the struct array
        group_tables(g).group_id = g;
        group_tables(g).output_table = output_table;

        % % % Write the output table to the cell array
        % % output_table = output_table;
        % % group_output_tables{g} = output_table;

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

%% Access data:
% group_tables(1).output_table.ac_consumption_actual to extract the ac comp
% from the output table (with all data) of the group 1 of trips

%% Ensure group_tables has the correct size
num_groups = size(group_tables, 2); % Adjust dynamically based on the struct array

%% PLOTS
%% Bar plot of mean ac consumption per group
% 
% % Initialize an array to store the mean of ac_consumption_actual for each group
% mean_ac_consumption_per_group = nan(num_groups, 1); % Use NaN as default
% 
% % Calculate the mean ac_consumption_actual for each group
% for g = 1:num_groups
%     if isfield(group_tables, 'output_table') && ~isempty(group_tables(g).output_table)
%         % Get the filtered output table for the group
%         current_table = group_tables(g).output_table;
% 
%         % Calculate the mean of ac_consumption_actual for the group
%         mean_ac_consumption_per_group(g) = mean(current_table.ac_consumption_actual);
%     end
% end
% 
% % Plot the distribution of the mean ac_consumption_actual for each group
% figure;
% bar(mean_ac_consumption_per_group, 'FaceColor', [0.2, 0.6, 0.8]);
% xticks(1:num_groups);
% xticklabels({
%     '1-5 min', '5-10 min', '10-20 min', '20-30 min', ...
%     '30-40 min', '40-50 min', '50-60 min', 'Over 60 min'
% });
% xlabel('Trip Duration Groups');
% ylabel('Mean AC Consumption Actual');
% title('Distribution of Mean AC Consumption Actual by Trip Duration Groups');
% grid on;

%% Bar plot of ac consumption per trip for each group

% % Define the custom figure size
% figure_width = 1600; % Double the default width (in pixels)
% figure_height = 600; % Keep the default height
% figure_units = 'pixels'; % Specify units for figure size
% 
% % Loop through each group to calculate and plot the mean AC consumption per trip
% for g = 1:num_groups
%     if isfield(group_tables, 'output_table') && ~isempty(group_tables(g).output_table)
%         % Get the filtered output table for the group
%         current_table = group_tables(g).output_table;
% 
%         % Get unique trips in the group
%         unique_trips = unique(current_table.trip_id);
% 
%         % Calculate the mean AC consumption for each trip
%         mean_ac_per_trip = zeros(length(unique_trips), 1);
%         for t = 1:length(unique_trips)
%             % Filter data for the current trip
%             trip_mask = current_table.trip_id == unique_trips(t);
%             trip_data = current_table(trip_mask, :);
% 
%             % Calculate the mean AC consumption for the trip
%             mean_ac_per_trip(t) = mean(trip_data.ac_consumption_actual);
%         end
% 
%         % Convert unique_trips to a cell array of strings for xticklabels
%         unique_trips_str = cellstr(string(unique_trips));
% 
%         % Plot the mean AC consumption for all trips in the group
%         figure('Units', figure_units, 'Position', [100, 100, figure_width, figure_height]);
%         bar(mean_ac_per_trip, 'FaceColor', [0.2, 0.6, 0.8]);
%         xticks(1:length(unique_trips));
%         xticklabels(unique_trips_str);
%         xlabel('Trip ID');
%         ylabel('Mean AC Consumption Actual');
%         title(sprintf('Group %d: Mean AC Consumption Actual Per Trip', g));
%         grid on;
% 
%         % Save the figure as a PNG file
%         filename = sprintf('Group_%d_Mean_AC_Consumption.png', g);
%         saveas(gcf, filename);
%         fprintf('Figure for Group %d saved as %s\n', g, filename);
% 
%         % Close the figure to prevent clutter
%         close(gcf);
%     else
%         fprintf('Group %d has no data to plot.\n', g);
%     end
% end
% 
%% Histogram ac consumption per group
% 
% % Define the custom figure size
% figure_width = 1600; % Double the default width (in pixels)
% figure_height = 600; % Keep the default height
% figure_units = 'pixels'; % Specify units for figure size
% 
% % Loop through each group to calculate and plot the distribution of mean AC consumption
% for g = 1:num_groups
%     if isfield(group_tables, 'output_table') && ~isempty(group_tables(g).output_table)
%         % Get the filtered output table for the group
%         current_table = group_tables(g).output_table;
% 
%         % Get unique trips in the group
%         unique_trips = unique(current_table.trip_id);
% 
%         % Calculate the mean AC consumption for each trip
%         mean_ac_per_trip = zeros(length(unique_trips), 1);
%         for t = 1:length(unique_trips)
%             % Filter data for the current trip
%             trip_mask = current_table.trip_id == unique_trips(t);
%             trip_data = current_table(trip_mask, :);
% 
%             % Calculate the mean AC consumption for the trip
%             mean_ac_per_trip(t) = mean(trip_data.ac_consumption_actual);
%         end
% 
%         % Plot the distribution of mean AC consumption
%         figure('Units', figure_units, 'Position', [100, 100, figure_width, figure_height]);
%         histogram(mean_ac_per_trip, 'FaceColor', [0.2, 0.8, 0.4], 'EdgeColor', 'k', 'BinWidth', 0.01);
%         xlabel('Mean AC Consumption Actual');
%         ylabel('Frequency');
%         title(sprintf('Group %d: Distribution of Mean AC Consumption Actual', g));
%         grid on;
% 
%         % Save the figure as a PNG file
%         filename = sprintf('Group_%d_Mean_AC_Consumption_Distribution.png', g);
%         saveas(gcf, filename);
%         fprintf('Distribution figure for Group %d (Mean AC Consumption) saved as %s\n', g, filename);
% 
%         % Close the figure to prevent clutter
%         close(gcf);
%     else
%         fprintf('Group %d has no data to plot.\n', g);
%     end
% end
% 
% 
%% Bar plot active time ac compressor per trip and group
% 
% % Define the custom figure size
% figure_width = 1600; % Double the default width (in pixels)
% figure_height = 600; % Keep the default height
% figure_units = 'pixels'; % Specify units for figure size
% 
% % Loop through each group to calculate and plot the active time per trip
% for g = 1:num_groups
%     if isfield(group_tables, 'output_table') && ~isempty(group_tables(g).output_table)
%         % Get the filtered output table for the group
%         current_table = group_tables(g).output_table;
% 
%         % Get unique trips in the group
%         unique_trips = unique(current_table.trip_id);
% 
%         % Calculate the active time for each trip
%         active_time_per_trip = zeros(length(unique_trips), 1);
%         for t = 1:length(unique_trips)
%             % Filter data for the current trip
%             trip_mask = current_table.trip_id == unique_trips(t);
%             trip_data = current_table(trip_mask, :);
% 
%             % Calculate the active time (sum of time differences where AC > 0.05)
%             active_mask = trip_data.ac_consumption_actual > 0.05; % active time =1 when AC>50 W
%             active_times = trip_data.absolute_time(active_mask);
%             if ~isempty(active_times)
%                 active_time_per_trip(t) = active_times(end) - active_times(1);
%             end
%         end
% 
%         % Convert unique_trips to a cell array of strings for xticklabels
%         unique_trips_str = cellstr(string(unique_trips));
% 
%         % Plot the active time for all trips in the group
%         figure('Units', figure_units, 'Position', [100, 100, figure_width, figure_height]);
%         bar(active_time_per_trip, 'FaceColor', [0.8, 0.4, 0.2]);
%         xticks(1:length(unique_trips));
%         xticklabels(unique_trips_str);
%         xlabel('Trip ID');
%         ylabel('Active Time (seconds)');
%         title(sprintf('Group %d: Active Time Per Trip', g));
%         grid on;
% 
%         % Save the figure as a PNG file
%         filename = sprintf('Group_%d_Active_Time.png', g);
%         saveas(gcf, filename);
%         fprintf('Figure for Group %d (Active Time) saved as %s\n', g, filename);
% 
%         % Close the figure to prevent clutter
%         close(gcf);
%     else
%         fprintf('Group %d has no data to plot.\n', g);
%     end
% end
% 
%% Histograms for active time ac compressor
% 
% % Define the custom figure size
% figure_width = 1600; % Double the default width (in pixels)
% figure_height = 600; % Keep the default height
% figure_units = 'pixels'; % Specify units for figure size
% 
% % Loop through each group to calculate and plot the distribution of active time
% for g = 1:num_groups
%     if isfield(group_tables, 'output_table') && ~isempty(group_tables(g).output_table)
%         % Get the filtered output table for the group
%         current_table = group_tables(g).output_table;
% 
%         % Get unique trips in the group
%         unique_trips = unique(current_table.trip_id);
% 
%         % Calculate the active time for each trip
%         active_time_per_trip = zeros(length(unique_trips), 1);
%         for t = 1:length(unique_trips)
%             % Filter data for the current trip
%             trip_mask = current_table.trip_id == unique_trips(t);
%             trip_data = current_table(trip_mask, :);
% 
%             % Calculate the active time (sum of time differences where AC > 0.05)
%             active_mask = trip_data.ac_consumption_actual > 0.05;
%             active_times = trip_data.absolute_time(active_mask);
%             if ~isempty(active_times)
%                 active_time_per_trip(t) = active_times(end) - active_times(1);
%             end
%         end
% 
%         % Filter out zero values (trips with no active time)
%         active_time_per_trip = active_time_per_trip(active_time_per_trip > 0);
% 
%         % Plot the distribution of active time
%         figure('Units', figure_units, 'Position', [100, 100, figure_width, figure_height]);
%         histogram(active_time_per_trip, 'FaceColor', [0.4, 0.6, 0.8], 'EdgeColor', 'k', 'BinWidth', 10);
%         xlabel('Active Time (seconds)');
%         ylabel('Frequency');
%         title(sprintf('Group %d: Distribution of Active Time', g));
%         grid on;
% 
%         % Save the figure as a PNG file
%         filename = sprintf('Group_%d_Active_Time_Distribution.png', g);
%         saveas(gcf, filename);
%         fprintf('Distribution figure for Group %d saved as %s\n', g, filename);
% 
%         % Close the figure to prevent clutter
%         close(gcf);
%     else
%         fprintf('Group %d has no data to plot.\n', g);
%     end
% end
% 
%% Histogram mean ac consumption by cluster of 0.5 kW
% 
% % Loop through each group to calculate and plot the clustered mean AC consumption
% for g = 1:num_groups
%     if isfield(group_tables, 'output_table') && ~isempty(group_tables(g).output_table)
%         % Get the filtered output table for the group
%         current_table = group_tables(g).output_table;
% 
%         % Get unique trips in the group
%         unique_trips = unique(current_table.trip_id);
% 
%         % Calculate the mean AC consumption for each trip
%         mean_ac_per_trip = zeros(length(unique_trips), 1);
%         for t = 1:length(unique_trips)
%             % Filter data for the current trip
%             trip_mask = current_table.trip_id == unique_trips(t);
%             trip_data = current_table(trip_mask, :);
% 
%             % Calculate the mean AC consumption for the trip
%             mean_ac_per_trip(t) = mean(trip_data.ac_consumption_actual);
%         end
% 
%         % Calculate the maximum AC consumption for this group
%         max_ac_consumption = max(mean_ac_per_trip);
% 
%         % Define bin edges for clustering (0 to 0.5, 0.5 to 1, etc.)
%         bin_edges = 0:0.5:max_ac_consumption;
% 
%         % Plot the clustered distribution of mean AC consumption
%         figure('Units', figure_units);
%         histogram(mean_ac_per_trip, 'BinEdges', bin_edges, 'FaceColor', [0.2, 0.8, 0.4], 'EdgeColor', 'k');
%         xlabel('Mean AC Consumption Actual (clustered)');
%         ylabel('Frequency');
%         title(sprintf('Group %d: Clustered Distribution of Mean AC Consumption Actual', g));
%         grid on;
% 
%         % Save the figure as a PNG file
%         filename = sprintf('Group_%d_Clustered_Mean_AC_Consumption_Distribution.png', g);
%         saveas(gcf, filename);
%         fprintf('Clustered distribution figure for Group %d (Mean AC Consumption) saved as %s\n', g, filename);
% 
%         % Close the figure to prevent clutter
%         close(gcf);
%     else
%         fprintf('Group %d has no data to plot.\n', g);
%     end
% end
