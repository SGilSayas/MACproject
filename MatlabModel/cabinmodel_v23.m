% Vehicle's Cabin Model, Susana Gil-Sayas. Version: next of 15.2 new Qtarget, but
% not heat transfer coef validation added yet
clc
clear all 
close all

%% Test conditions
Irr = 0; % W/m2 
N_Humans = 1; % old: 1 human = 120 W, 0.66 = 80 W, 0.25 = 30 W;
heat_human = 0; % W
heat_equipment = 0 ; % 40; % W
T_ambient = 22; % [C] % NOT USED WHEN IMPORTING DATA
duration_WLTC = 1800; % s
timestep = 1; % s
t = 1; 
COP = 3.5;
heatpump = 1; % 0 for no, 1 for yes

%% What calibration? 
off_G8_22 = 0;
off_G7_22 = 0;
off_G8_35 = 0;
off_G7_35 = 0;
off_I5PV_35 = 0;
off_mustang_35 = 0;
off_tucsonCD_35 = 1;

bev = 1;
phev_cd = 0;
phev_cs = 0;
icev = 0;

hatchback = 0;
suv = 1; coef=1.5;

% Deletion high peak in G8 22C TC cell data?
deletion_peak = 1; % 1 for yes, 0 for no

% Want to reduce the amplitud of T amb? % Needed for G7_22
reduced_amplitude = 0; % 1 for yes, 0 for no
reduction = 5; 

%% Load input temperature
if(off_tucsonCD_35 == 1)
    filename = 'tucsonphev_35C_MACoff_TC.mat';
    loaded_data = load(filename);
    field_names = fieldnames(loaded_data);
    max_rows = max(cellfun(@(f) size(loaded_data.(f), 1), field_names));
    collected_data = cell(max_rows, numel(field_names));
    for i = 1:numel(field_names)
        % Extract the field data (already a double array)
        field_data = loaded_data.(field_names{i});
        % Remove NaN values and firs row (time = 0s) from field_data
        field_data = field_data(~isnan(field_data));
        % field_data(1, :) = [];
        num_rows = size(field_data, 1);
        % Store cleaned data in collected_data
        collected_data(1:num_rows, i) = num2cell(field_data);
        if num_rows < max_rows
            collected_data(num_rows+1:max_rows, i) = {NaN};
        end
        % Save each column as a separate variable without NaN values
        assignin('base', field_names{i}, field_data);
    end
elseif(off_mustang_35 == 1)
    filename = 'mustangbev_35C_MACoff_TC.mat';
    loaded_data = load(filename);
    field_names = fieldnames(loaded_data);
    max_rows = max(cellfun(@(f) size(loaded_data.(f), 1), field_names));
    collected_data = cell(max_rows, numel(field_names));
    for i = 1:numel(field_names)
        field_table = loaded_data.(field_names{i});
        field_data = table2array(field_table); % Extract data from table
        % Remove NaN values and firs row (time = 0s) from field_data
        field_data = field_data(~isnan(field_data));
        field_data(1, :) = [];
        num_rows = size(field_data, 1);
        % Store cleaned data in collected_data
        collected_data(1:num_rows, i) = num2cell(field_data);
        if num_rows < max_rows
            collected_data(num_rows+1:max_rows, i) = {NaN};
        end
        % Save each column as a separate variable without NaN values
        assignin('base', field_names{i}, field_data);
    end
elseif(off_I5PV_35 == 1)
    filename = 'ioniq5PV_35C_MACoff_TC.mat';
    loaded_data = load(filename);
    field_names = fieldnames(loaded_data);
    % Determine the maximum number of rows
    max_rows = max(cellfun(@(f) size(loaded_data.(f), 1), field_names));
    % Initialize a cell array to store column-aligned data
    collected_data = cell(max_rows, numel(field_names));
    % Fill the collected_data array with each field's values, padding with NaN if necessary
    for i = 1:numel(field_names)
        field_table = loaded_data.(field_names{i});
        field_data = table2array(field_table); % Extract data from table
        % Remove NaN values and firs row (time = 0s) from field_data
        field_data = field_data(~isnan(field_data));
        field_data(1, :) = [];
        num_rows = size(field_data, 1);
        % Store cleaned data in collected_data
        collected_data(1:num_rows, i) = num2cell(field_data);
        if num_rows < max_rows
            collected_data(num_rows+1:max_rows, i) = {NaN};
        end
        % Save each column as a separate variable without NaN values
        assignin('base', field_names{i}, field_data);
    end
elseif(off_G7_22==1)
    reduction = 0;  
    lab_input='G7_22C_MACoff.mat';
    struct=load(lab_input);
    table=struct2table(struct);
    lab_table=table2array(table);
    compressorspeed_rpms=table2array(lab_table(:,1));
    compressortorque_Nm=table2array(lab_table(:,2));
    TC_vent=table2array(lab_table(:,3)); % C
    TC_front=table2array(lab_table(:,4)); % C
    TC_driver=table2array(lab_table(:,5)); % C
    TC_pass=table2array(lab_table(:,6)); % C
    TC_back=table2array(lab_table(:,7)); % C
    TC_backpass=table2array(lab_table(:,8)); % C
    TC_cell=table2array(lab_table(:,9)); % C
    T_cell_OBD=table2array(lab_table(:,10)); % C
    T_evap=table2array(lab_table(:,11)); % C
    vehiclespeed_kmh=table2array(lab_table(:,12));
    n = height(lab_table)/(duration_WLTC); 
    TC_cell=TC_cell(1 : n : end);
    TC_vent=TC_vent(1 : n : end);
    TC_front=TC_front(1 : n : end);
    TC_driver=TC_driver(1 : n : end);
    TC_pass=TC_pass(1 : n : end);
    TC_back=TC_back(1 : n : end);
    TC_backpass=TC_backpass(1 : n : end);
    TC_cell=TC_cell(:,:);
    TC_vent=TC_vent(:,:);
    TC_front=TC_front(:,:);
    TC_driver=TC_driver(:,:);
    TC_pass=TC_pass(:,:);
    TC_back=TC_back(:,:);
    TC_backpass=TC_backpass(:,:);
elseif(off_G7_35==1)
    lab_input='G7_35C_MACoff_TC'; 
    struct=load(lab_input);
    table=struct2table(struct);
    lab_array=table2array(table);
    i=13;
    cell_array=table2array(lab_array(:,7));
    vent_array=table2array(lab_array(:,1));
    front_array=table2array(lab_array(:,2));
    driver_array=table2array(lab_array(:,3));
    pass_array=table2array(lab_array(:,4));
    back_array=table2array(lab_array(:,5)); 
    backpass_array=table2array(lab_array(:,6));
    n = height(lab_array)/(duration_WLTC); %10;
    TC_cell=cell_array(1:n:end-i);
    TC_vent=vent_array(1:n:end-i);
    TC_front=front_array(1:n:end-i);
    TC_driver=driver_array(1:n:end-i);
    TC_pass=pass_array(1:n:end-i);
    TC_back=back_array(1:n:end-i);
    TC_backpass=backpass_array(1:n:end-i);    
    TC_cell(1802:end,:)=[];
    TC_vent(1802:end,:)=[];
    TC_front(1802:end,:)=[];
    TC_driver(1802:end,:)=[];
    TC_pass(1802:end,:)=[];
    TC_back(1802:end,:)=[];
    TC_backpass(1802:end,:)=[];
        % For G7 35C: 
    TC_front=TC_vent;
elseif(off_G8_22==1)
    lab_input='G8_22C_MACoff_TCcell1'; 
    struct=load(lab_input);
    table=struct2table(struct);
    lab_table=table2array(table);
    array=table2array(lab_table(:,1));
    n0=height(lab_table)/duration_WLTC;
    TC_cell=array(1 : n0 : end);
    TC_cell=TC_cell(:,:); 
    lab_input='G8_22C_MACoff_TCcabin1';
    struct=load(lab_input);
    table=struct2table(struct);
    lab_table1=table2array(table);
    lab_table1=rmmissing(lab_table1); %delete NaN rows
    TC_vent=table2array(lab_table1(:,1)); 
    TC_front=table2array(lab_table1(:,2)); 
    TC_driver=table2array(lab_table1(:,3));
    TC_pass=table2array(lab_table1(:,4)); 
    TC_back=table2array(lab_table1(:,5));
    TC_backpass=table2array(lab_table1(:,6));
    n=height(lab_table1)/duration_WLTC;
    TC_vent=TC_vent(1 : n : end);
    TC_front=TC_front(1 : n : end);
    TC_driver=TC_driver(1 : n : end);
    TC_pass=TC_pass(1 : n : end);
    TC_back=TC_back(1 : n : end);
    TC_backpass=TC_backpass(1 : n : end);
    TC_vent=TC_vent(:,:);
    TC_front=TC_front(:,:);
    TC_driver=TC_driver(:,:);
    TC_pass=TC_pass(:,:);
    TC_back=TC_back(:,:);
    TC_backpass=TC_backpass(:,:);
    % Deletion high peak in G8 22C TC cell data:
    Total_time=length(TC_cell); 
    for k=1:Total_time
        if(k > 1600)
            TC_cell(k)=TC_cell(k-200)+0.4;
        end
    end
elseif(off_G8_35==1)
    lab_input='G8_35C_MACoff_TC';
    struct=load(lab_input);
    table=struct2table(struct);
    lab_array=table2array(table);
    cell_array=table2array(lab_array(:,7));
    vent_array=table2array(lab_array(:,1));
    front_array=table2array(lab_array(:,2));
    driver_array=table2array(lab_array(:,3));
    pass_array=table2array(lab_array(:,4));
    back_array=table2array(lab_array(:,5)); 
    backpass_array=table2array(lab_array(:,6));
    n=height(lab_array)/(duration_WLTC); %10;%8.5;
    TC_cell=cell_array(1 : n : end);
    TC_vent=vent_array(1 : n : end);
    TC_front=front_array(1 : n : end);
    TC_driver=driver_array(1 : n : end);
    TC_pass=pass_array(1 : n : end);
    TC_back=back_array(1 : n : end);
    TC_backpass=backpass_array(1 : n : end);
    TC_cell=TC_cell(:,:);
    TC_vent=TC_vent(:,:); 
    TC_front=TC_front(:,:); 
    TC_driver=TC_driver(:,:); 
    TC_pass=TC_pass(:,:);
    TC_back=TC_back(:,:);
    TC_backpass=TC_backpass(:,:);
        % For G8 35C: (error in frontal TCs measurements)
    TC_front=TC_backpass;
end
% save G8_35_cell_front_back TC_cell TC_front TC_back

%% Interpolate TC vectors to fit WLTC duration
n = length(TC_cell); % Get the total number of elements
if(n > duration_WLTC) % Select 1800 evenly spaced indices
    indices = round(linspace(1, n, 1800));
    TC_cell = TC_cell(indices);
    TC_vent = TC_vent(indices);
    TC_front = TC_front(indices);
    TC_driver = TC_driver(indices);
    TC_pass = TC_pass(indices);
    TC_back = TC_back(indices);
    TC_backpass = TC_backpass(indices);
end

%% Signals correction for unusual sinusoidal shape due to measurement problems
if(off_G7_22==1 || off_G7_35==1 )
    l=[length(TC_cell),
        length(TC_vent),
        length(TC_front),
        length(TC_driver),
        length(TC_pass),
        length(TC_back),
        length(TC_backpass)];
    T=[fft(TC_cell(:,:)/l(1)),fft(TC_vent(:,:)/l(2)),fft(TC_front(:,:)/l(3)),...
        fft(TC_driver(:,:)/l(4)),fft(TC_pass(:,:)/l(5)),fft(TC_back(:,:)/l(6)),fft(TC_backpass(:,:)/l(7))];
    k=0.1;
    Tfil=T(:,:);
    Tfil(2:15,:)=T(2:15,:)*k; 
    Tfil(l-13:l,:)=T(l-13:l,:)*k;
    TC_fil=[ifft(Tfil(:,1).*l(1),l(1)),ifft(Tfil(:,2).*l(2),l(2)),ifft(Tfil(:,3).*l(3),l(3)),...
        ifft(Tfil(:,4).*l(4),l(4)),ifft(Tfil(:,5).*l(5),l(5)),ifft(Tfil(:,6).*l(6),l(6)),ifft(Tfil(:,7).*l(7),l(7))];
    TC_cell=TC_fil(:,1);
    TC_vent=TC_fil(:,2);
    TC_front=TC_fil(:,3);
    TC_driver=TC_fil(:,4);
    TC_pass=TC_fil(:,5);
    TC_back=TC_fil(:,6);
    TC_backpass=TC_fil(:,7);
elseif(off_mustang_35==1) % step shape
    window_size = 50; % 50 %Adjust
    TC_cell = medfilt1(TC_cell, window_size);
    % TC_vent = medfilt1(TC_vent, window_size);
    % TC_front = medfilt1(TC_front, window_size);
    % TC_driver = medfilt1(TC_driver, window_size);
    % TC_pass = medfilt1(TC_pass, window_size);
    % TC_back = medfilt1(TC_back, window_size);
    % TC_backpass = medfilt1(TC_backpass, window_size);
elseif(off_I5PV_35 == 1)
    l=[length(TC_cell)];
    T=[fft(TC_cell(:,:)/l(1))];
    k=0.1;
    Tfil=T(:,:);
    Tfil(2:15,:)=T(2:15,:)*k; 
    Tfil(l-13:l,:)=T(l-13:l,:)*k;
    TC_fil=[ifft(Tfil(:,1).*l(1),l(1))];
    TC_cell=TC_fil(:,1);
    Total_time=length(TC_cell);
    time = 1:timestep:Total_time;
    time = time'; % Ensure it's a column vector

    % Extract the trend (data from 61 to 1800)
    time_trend = time(61:1770);
    TC_trend = TC_cell(61:1770);
    
    % Estimate frequency using FFT (optional but helps with accuracy)
    Y = fft(TC_trend - mean(TC_trend)); % Remove DC component
    P2 = abs(Y / length(TC_trend));
    P1 = P2(1:floor(length(TC_trend) / 2)); % Single-sided spectrum
    [~, idx] = max(P1(2:end)); % Find dominant frequency, ignoring DC component
    dominant_freq = (idx + 1) / length(TC_trend) * (0.4 * pi); % (2 * pi) %Convert to angular frequency
    
    % Define a sinusoidal fit model explicitly
    sinModel = fittype('a*sin(b*x + c) + d', 'independent', 'x', 'coefficients', {'a', 'b', 'c', 'd'});
    
    % Fit the data (adjust StartPoint based on estimated frequency)
    sin_fit = fit(time_trend, TC_trend, sinModel, ...
        'StartPoint', [0.5, dominant_freq, 0, mean(TC_trend)]);
    
    % Generate the first 60  and last 30 seconds using the fitted function
    time_first_60 = time(1:60);
    TC_substitute_ini = sin_fit.a * sin(sin_fit.b * time_first_60 + sin_fit.c) + sin_fit.d;
    TC_cell(1:60) = TC_substitute_ini;
    
    time_last_30 = time(1770:1800);
    TC_substitute_end = sin_fit.a * sin(sin_fit.b * time_last_30 + sin_fit.c) + sin_fit.d;
    TC_cell(1770:1800) = TC_substitute_end;
    
    figure(13);
    hold on;
    plot(time, TC_cell, 'b', 'DisplayName', 'Modified TC');
    plot(time_first_60, TC_substitute_ini, 'r', 'LineWidth', 2, 'DisplayName', 'Substituted First 60s');
    hold on
    plot(time_last_30, TC_substitute_end, 'r', 'LineWidth', 2, 'DisplayName', 'Substituted Last 30s');
    legend;
    xlabel('Time (s)');
    ylabel('Temperature (°C)');
    title('Temperature Signal with Replaced First 60 Seconds');
    hold off;
elseif(off_tucsonCD_35==1)
    % l=[length(TC_cell)];
    % T=[fft(TC_cell(:,:)/l(1))];
    % k=0.1;
    % Tfil=T(:,:);
    % Tfil(2:15,:)=T(2:15,:)*k; 
    % Tfil(l-13:l,:)=T(l-13:l,:)*k;
    % TC_fil=[ifft(Tfil(:,1).*l(1),l(1))];
    % TC_cell=TC_fil(:,1);
end

%% Definition of total time for the calulations
Total_time=length(TC_cell); 
time = 1:timestep:Total_time;

%% Load input Engine RPM (and XCU data)
lab_input2='G7_35C_MACoff_XCU'; 
struct2=load(lab_input2);
table2=struct2table(struct2);
modalxcu_array=table2array(table2);
ACcomp_speed_rpm = table2array(modalxcu_array(:,2));
ACcomp_torque_Nm = table2array(modalxcu_array(:,3));
T_Cell_XCU = table2array(modalxcu_array(:,4));
T_evap_C = table2array(modalxcu_array(:,5));
ICETorque_Nm = table2array(modalxcu_array(:,6));
engineSpeed_rpm = table2array(modalxcu_array(:,7));
fuelconsumption_lph = table2array(modalxcu_array(:,8));
refrigewrantFluidPress_bar = table2array(modalxcu_array(:,9));
exaustGasTempBank1_C = table2array(modalxcu_array(:,10));
n_xcu = height(modalxcu_array)/duration_WLTC;
ACcomp_rpm = ACcomp_speed_rpm(1 : round(n_xcu) : end);
ACcomp_Nm = ACcomp_torque_Nm(1 : round(n_xcu) : end);
engine_Nm = ICETorque_Nm(1 : round(n_xcu) : end);
engine_rpm = engineSpeed_rpm(1 : round(n_xcu) : end);
FC_lperh = fuelconsumption_lph(1 : round(n_xcu) : end);

vehiclespeed_struct = load("G7_vehiclespeed_kmh.mat","VehicleSpeed"); 
timespeed_struct = load("G7_vehiclespeed_kmh.mat","time_ModalXCU");
vehiclespeed_table = struct2table(vehiclespeed_struct);
timespeed_table = struct2table(timespeed_struct);
vehiclespeed_array = table2array(vehiclespeed_table);
timespeed_array = table2array(timespeed_table);
vehiclespeed_kmh = table2array(vehiclespeed_array(:,1));
time_vehiclespeed_s = table2array(timespeed_array(:,1));
vehiclespeed_n = height(vehiclespeed_array)/duration_WLTC;
time_vehiclespeed_s = time_vehiclespeed_s(1 : round(vehiclespeed_n) : end-1);
vehicle_kmh = vehiclespeed_kmh(1 : round(vehiclespeed_n) : end); %%% G7@35: end-1)
Total_time2=length(vehicle_kmh)-1; 
time2 = 0:timestep:Total_time2;

h_front = ones(length(vehicle_kmh),1);
h_side = ones(length(vehicle_kmh),1);
h_rear = ones(length(vehicle_kmh),1);
U_front = ones(length(vehicle_kmh),1);
U_back = ones(length(vehicle_kmh),1);
U_side = ones(length(vehicle_kmh),1);

%% Load input cell humidity (RH) and pressure (P)
load('G7_35_off_RH_P_amb')
RH_amb_mean=mean(RH_amb);
P_amb_kPa_mean=mean(P_amb_kPa);
RH_amb = ones(1,Total_time+1)*RH_amb_mean;
P_amb_kPa = ones(1,Total_time+1)*P_amb_kPa_mean;

%% Noise deletion
TC_cell = movmean(TC_cell,10);
TC_front = movmean(TC_front,20);
TC_back = movmean(TC_back,20);
TC_vent = movmean(TC_vent,10);
TC_driver = movmean(TC_driver,20);
TC_pass = movmean(TC_pass,10); % or middle
TC_backpass = movmean(TC_backpass,10); % or equipment

%% Shift of TC signals to same bias as Tamb
if(off_G7_22==1)
    TC_back= TC_back-abs(TC_back(1)-TC_cell(1)) ;
    TC_front= TC_front-abs(TC_back(6)-TC_front(6)) ;
    TC_vent= TC_vent - (TC_vent(60)-TC_cell(60));
    TC_driver= TC_driver - (TC_driver(60)-TC_cell(60));
    TC_backpass= TC_backpass - (TC_backpass(60)-TC_cell(60)); % or equipment
elseif(off_G7_35==1)
    TC_front= TC_front + abs(TC_front(6)-TC_cell(6))+0.2;
    TC_back= TC_back + abs(TC_back(2)-TC_cell(2));
    TC_vent= TC_vent + abs(TC_vent(2)-TC_cell(2));
    TC_driver= TC_driver + abs(TC_driver(2)-TC_cell(2));
    TC_backpass= TC_backpass - abs(TC_backpass(2)-TC_cell(2))+0.3; % or equipment
elseif(off_G8_22==1)
    TC_front= TC_front - abs(TC_front(260)-TC_cell(160))+0.3;
    TC_back= TC_back - abs(TC_back(260)-TC_cell(160))+0.3;
    TC_vent= TC_vent - abs(TC_vent(140)-TC_cell(105));
    TC_driver= TC_driver - abs(TC_driver(130)-TC_cell(105));
    TC_backpass= TC_backpass - abs(TC_backpass(182)-TC_cell(105)); 
elseif(off_G8_35==1)
    TC_front= TC_front - (TC_front(1)-TC_cell(1))+0.5+0.3; %new front
    TC_back= TC_back - (TC_back(1)-TC_cell(1))+0.3;
    % For G8 35C: (error in frontal TCs measurements)
    TC_vent= TC_vent - (TC_vent(1)-TC_cell(1))+0.5;
    TC_driver= TC_driver - (TC_driver(1)-TC_cell(1))+0.5;
    TC_backpass= TC_backpass - (TC_backpass(1)-TC_cell(1)); % or equipment
elseif(off_I5PV_35==1)
    TC_back = TC_back + abs(TC_back(1)-TC_cell(1)) ;
    TC_front = TC_front + abs(TC_front(1)-TC_cell(1));
    TC_vent = TC_vent + abs(TC_vent(1)-TC_cell(1));
    TC_driver = TC_driver + abs(TC_driver(1)-TC_cell(1));
    TC_backpass = TC_backpass + abs(TC_backpass(1)-TC_cell(1)); % or equipment
elseif(off_mustang_35 == 1)
    TC_back = TC_back - abs(TC_back(1)-TC_front(1));
    % TC_back = TC_back - abs(TC_back(120)-TC_cell(120));
    % TC_front = TC_front + abs(TC_front(120)-TC_cell(120));
elseif(off_tucsonCD_35 == 1)
    TC_back = TC_back + (TC_cell(1)-TC_back(1));
    TC_front = TC_front + (TC_cell(1)-TC_front(1));
    TC_vent = TC_vent + (TC_cell(1)-TC_vent(1));
    TC_driver = TC_driver + (TC_cell(1)-TC_driver(1));
    TC_backpass = TC_backpass + (TC_cell(1)-TC_backpass(1));
end

%% Transpose TC amb vector and initialization
T_cell = transpose(TC_cell)+273.15; % T_amb in K. Both TC_cell or T_cell_OBD
T_ini = T_cell(1);

%% Deletion peaks due to filtering at beginning and end:
deletion_filter_peaks = 1;
i=0;
if(deletion_filter_peaks == 1 && (off_G7_22==1 || off_G7_35==1))
    for n=1:Total_time
        if(n < 6)
            TC_front(n)=TC_back(n);
        end
        if(n < 33)
            TC_vent(n)=TC_back(n);
        elseif(n > 1750)
            T_cell(n)=T_cell(1750-i);
            TC_front(n)=TC_front(1750-i);
            TC_back(n)=TC_back(1750-i);
            TC_vent(n)=TC_vent(1750-i);
            TC_driver(n)=TC_driver(1750-i);
            TC_backpass(n)=TC_backpass(1750-i);
            i=i+1;
        end
    end
end

Irr = ones(1,Total_time)*Irr;
% if (Irr==0)
%     T_amb = T_ambient + 273.15; % [K]
% else
%     T_amb = ones(1,Total_time+1)*(T_ambient+273.15); %[K], vector
%     Irr = ones(1,Total_time+1)*Irr;
% end
    
%% Materials' thickness, e[m]
e_ws=0.006;
e_rw=0.005;
e_ABS=0.002;
e_steel=0.0005;
e_PU = 0.025; % polyurethane 
e_ceiling=e_ABS+e_steel+e_PU;
e_sidewindows=0.003;
e_doors=e_ceiling;

% Conduction Properties: thermoconductivity, k[W/(m*K)]
k_ws=0.8; % 0.8-1.4
k_rw=1.4; % 1.4-1.5
k_ABS=0.1; % Acrylonitrile butadiene styrene
k_steel=14.9;
k_PU = 0.022; % 0.022-0.028 W/mK polyurethane
k_ceiling=k_ABS*k_steel*k_PU;
k_sidewindows=1.4;
k_doors=k_ceiling;

% Initial Heat transfer Convection coefficients, W/(m2*K)
h_ws=5; 
h_sidewindows=5;
h_doors=5;
h_rw=5;
h_ceiling=5;
h_roof=20;
h_ext=20; 
h_base = 5;
h_cabin=5; %10;

%% Vehicle shape properties: Area, A (m2) ; Cabin geometry (m)
% if(hatchback == 1) %Golfs measurements
%     A_ws = 0.63*1.3;
%     A_rw = 0.29*1;
%     A_roof = 1.8*1.1;
%     A_sidewindows = 2*1.45*0.29;
%     A_doors = 2*1.45*0.29;
%     A_base = 2; % Marcos et al: 6m2
%     cabin_height =  0.5456;
%     cabin_width = 1.3;
%     cabin_roof_lenght = 1.8;
%     cabin_base_lenght = 1.45;
% elseif(suv == 1)
    
    A_ws = 0.63*1.3*coef;
    A_rw = 0.29*1*coef;
    A_roof = 1.8*1.1*coef;
    A_sidewindows = 2*1.45*0.29*coef;
    A_doors = 2*1.45*0.29*coef;
    A_base = 3*coef;
    cabin_height =  0.5456*coef;
    cabin_width = 1.3*coef;
    cabin_roof_lenght = 1.8*coef;
    cabin_base_lenght = 1.45*coef;
% end

A_dashboard = cabin_width*0.5; %m2
A_front = A_ws + A_dashboard; %3
A_back = A_rw; % 1.5
A_side = A_sidewindows + A_doors + A_roof;

seatcm3 = 150; %115; % cm3
seats_volume = 5*seatcm3*10^(-6); % m3

A_vehicle = A_ws + A_rw + A_roof + A_sidewindows + A_doors + A_base;
CabinVolume = (cabin_width*cabin_roof_lenght*cabin_height/2) + ...
    (cabin_width*cabin_base_lenght*cabin_height/2) - seats_volume; 
%%
    % Air properies
density_air = 1.18; % [Kg/m3] (at 25C)
cp_air = 1006;% + 10*1500; % [J/ kg*K]

    % Seats properties
mass_seats = 30; % 20-30 kg
cp_seats = 2000; % 500-2000 [J/kgK]

    % Human and equipment heat % ISO 8996, Alberto Viti Corsi (IDAE)
Qhuman = ones(1,Total_time)*N_Humans*heat_human; % [W]
Qequipment = ones(1,Total_time)*heat_equipment; % [W]

    % Radiation Properties 
sigma=0.0000000567037321; % Stefan-Boltzmann constant(sigma)[W/(m^2*K^4)]
    % epsilon:emissivity[]
epsilon_ws=0.9;
epsilon_rw=0.9;
epsilon_sidewindows=0.9; 
epsilon_doors=0.9;
epsilon_roof=0.9;
    % rho=reflectivity[]
rho_ws=0.246;
rho_rw=0.1;
rho_sidewindows=0.2;
rho_doors=0.74;
rho_roof=0.74;
rho_base=0.3;
    % tau=transmissivity[]
tao_ws=0.452;
tao_rw=0.311;
tao_sidewindows=0.475;
tao_doors=0;
tao_roof=0;
tao_base=0;
    % alpha=absorptivity[]
alfa_ws = 1-rho_ws-tao_ws;
alfa_rw = 1-rho_rw-tao_rw;
alfa_sidewindows = 1-rho_sidewindows-tao_sidewindows;
alfa_doors = 1-rho_doors-tao_doors;
alfa_roof = 1-rho_roof-tao_roof;
alfa_base = 1-rho_base-tao_base;

    % Global Horizontal Irradiance (Hypothesis: same in all windows)
G_roof_inc=Irr;
G_ws_inc=Irr; 
G_rw_inc=Irr;
G_sidewindows_inc=Irr;
G_doors_inc=Irr;
G_ws_r=rho_ws*G_ws_inc;
G_ws_a=alfa_ws*G_ws_inc; %bc
G_ws_t=tao_ws*G_ws_inc;
G_rw_r=rho_rw*G_rw_inc;
G_rw_a=alfa_rw*G_rw_inc; 
G_rw_t=tao_rw*G_rw_inc;
G_sidewindows_r=rho_sidewindows*G_sidewindows_inc;
G_sidewindows_a=alfa_sidewindows*G_sidewindows_inc; 
G_sidewindows_t=tao_sidewindows*G_sidewindows_inc;
G_doors_r=rho_doors*G_doors_inc;
G_doors_a=alfa_doors*G_doors_inc; 
G_doors_t=tao_doors*G_doors_inc;
G_roof_a=alfa_roof*G_roof_inc;

Qirr = Irr.*(A_front+A_side+A_back); % W

    % Conductances Matrix 
    % Tamb --K1-> Troof --K2-> Tceiling --K3-> Tcabin_back
K1=h_ext*A_roof; %h_roof?
K2=A_roof*k_ceiling/e_ceiling;
K3=h_ceiling*A_roof;
    % Tbase --K4-> Tcabin_back
K4=h_base*A_base;
    % Tamb --K5-> Tw_ext_ws --K6-> Tw_int_ws --K7-> Tcabin_front
K5=h_ext*A_ws;
K6=A_ws*k_ws/e_ws;
K7=h_ws*A_ws;
    % Tamb --K8-> Tw_ext_rw --K9-> Tw_int_rw --K10-> Tcabin_back
K8=h_ext*A_rw;
K9=A_rw*k_rw/e_rw;
K10=h_rw*A_rw;
    % Tamb --K11-> Tw_ext_sidewindows --K12-> Tw_int_sidewindows --K13-> Tcabin_back
K11=h_ext*A_sidewindows;
K12=A_sidewindows*k_sidewindows/e_sidewindows;
K13=h_sidewindows*A_sidewindows;
    % Tamb --K14-> Tw_ext_doors --K15-> Tw_int_doors --K16-> Tcabin_back
K14=h_ext*A_doors;
K15=A_doors*k_doors/e_doors;
K16=h_doors*A_doors;
    % Tcabin_front --K17-> Tcabin_back
K17=h_cabin*A_ws; %TBD
K = [1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    K1 -(K1+K2) K2 0 0 0 0 0 0 0 0 0 0 0;
    0 K2 -(K2+K3/2+K3/2) 0 0 0 0 0 0 0 0 0 (K3/2) (K3/2);
    0 0 0 -(K4/2+K4/2) 0 0 0 0 0 0 0 0 (K4/2) (K4/2);
    K5 0 0 0 -(K5+K6) 0 0 0 K6 0 0 0 0 0;
    K8 0 0 0 0 -(K8+K9) 0 0 0 K9 0 0 0 0;
    K11 0 0 0 0 0 -(K11+K12) 0 0 0 K12 0 0 0;
    K14 0 0 0 0 0 0 -(K14+K15) 0 0 0 K15 0 0;
    0 0 0 0 K6 0 0 0 -(K6+K7) 0 0 0 K7 0;
    0 0 0 0 0 K9 0 0 0 -(K9+K10) 0 0 0 K10;
    0 0 0 0 0 0 K12 0 0 0 -(K12+K13/2+K13/2) 0 (K13/2) (K13/2);
    0 0 0 0 0 0 0 K15 0 0 0 -(K15+K16/2+K16/2) (K16/2) (K16/2);
    0 0 (K3/2) (K4/2) 0 0 0 0 K7 0 (K13/2) (K16/2) -(K3/2+K4/2+K7+K16/2+K17) K17;
    0 0 (K3/2) (K4/2) 0 0 0 0 0 (K10) (K13/2) (K16/2) (K17) -(K3/2+K4/2+K10+K13/2+K16/2+K17)];

% % CabinVolume=2;
% front_volume = CabinVolume*1.5;
% back_volume = CabinVolume*0.15;

    % Initial Capacitance Matrix
C_amb = 0; %b.c.
C_roof = 0; %for now
C_ceiling = 0; %for now
C_base = mass_seats*cp_seats + 144240; %Capacitance, W/K  
C_ws=0;
C_rw=0;
C_sidewindows=0;
C_doors=0;
C_cabin_front = density_air*cp_air/timestep*CabinVolume*0.5;
C_cabin_back = density_air*cp_air/timestep*CabinVolume*0.5;
C=zeros(3);
C(1,1)=C_amb;
C(2,2)=C_cabin_front;
C(3,3)=C_cabin_back;

% C=zeros(14);
% C(1,1)=C_amb;
% C(2,2)=C_roof;
% C(3,3)=C_ceiling;
% C(4,4)=C_base;
% C(5,5)=C_ws;
% C(6,6)=C_rw;
% C(7,7)=C_sidewindows;
% C(8,8)=C_doors;
% C(9,9)=C_ws;
% C(10,10)=C_rw;
% C(11,11)=C_sidewindows;
% C(12,12)=C_doors;
% C(13,13)=C_cabin_front;
% C(14,14)=C_cabin_back;

    % Engine power and temperature:
% engine_kW=engine_Nm.*engine_rpm*2*pi/(60*1000);

    % Humidity ratio in gram of water per gram of dry air, X
Ps_22_kPa = 2.626; % kPa, Water saturation pressure at 22C
Ps_35_kPa = 5.63; % kPa, Water saturation pressure at 35C
Ps_kPa = Ps_22_kPa+(T_ini-22)/(35-22)*(Ps_35_kPa-Ps_22_kPa);
meanP_kPa_amb = mean(P_amb_kPa); % in kPa
X = 0.0732; % X = 0.62198*Ps_kPa.*RH_amb./(100.*P_amb_kPa - Ps_kPa.*RH_amb); % Faya, 2013

    % Enthalpies calculation. Humidity considered same in amb and cabin
e_amb = 1006.*T_cell + (2501000 + 1770.*T_cell)*X; % J/kg, Faya 2013

    % Air vent volume and flow rate:
vent_volumerate=0.0003; % m3/s, Faya 2013: 0.02; 0.001
leakage_volumerate = 0.0001; %% m3/s, Faya2014

    % CO2 emission calculation parameters
Vpe_petrol = 0.264; % l/kWh
CF_petrol = 2330; % gCO2/l
distance_driven = 23.25; % km

    % GENERAL EQUATIONS
    % % m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 
    % T_air(t) = (Q_windows(t)+Q_ceiling(t)+Q_base(t)+Q_human(t))*(timestep/(CabinVolume*cp_air*density_air)) + T_air(t-1);
    % Q_ws(t) = h_ws*A_ws*(T_ws(t)- T_air(t-1));
    % Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
    % Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
    % Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));
    % Q_windows(t)=Q_ws(t)+Q_sidewindows(t)+Q_doors(t)+Q_rw(t);

    %% INI
Troof(t)=T_ini;
Tceiling(t)=T_ini;
Tbase_int(t)=T_ini;
Tw_ext_ws(t)=T_ini;
Tw_ext_rw(t)=T_ini;
Tw_ext_sidewindows(t)=T_ini;
Tw_ext_doors(t)=T_ini;
Tw_int_ws(t)=T_ini;
Tw_int_rw(t)=T_ini;
Tw_int_sidewindows(t)=T_ini;
Tw_int_doors(t)=T_ini;
Tcabin_front(t)=T_ini;
Tcabin_back(t)=T_ini;
Tcabin(t)=(Tcabin_front(t)+Tcabin_back(t))/2;

Tcabin_front = ones(1,Total_time)*Tcabin_front(t);
Tcabin_back =  ones(1,Total_time)*Tcabin_back(t);
Tcabin = ones(1,Total_time)*(Tcabin_front(t)+Tcabin_back(t))/2; 

temperature=[T_ini;T_ini;T_ini];
prev_temp=[T_ini;T_ini;T_ini];

delta_temp_engine_mean = -ones(1,Total_time)*0.001;
delta_temp_tailpipe_mean = -ones(1,Total_time)*0.001;
check=0; check_f=0; check_b=0; check_t=0;
Qirr = ones(1,Total_time);
Qengine = ones(1,Total_time);
Qexhaust = ones(1,Total_time);
Qleakage = ones(1,Total_time);
Qtotal = ones(1,Total_time);
Qcv_emitted = ones(1,Total_time);
Qcompressor = ones(1,Total_time).*0;
CO2 = ones(1,Total_time).*0;
Qcv_emitted(t)=0;
t_active(t)=0;
compressor_time(t)=0;
heat_human = ones(1,Total_time);
T_engine(t)=T_ini;
T_exhaust(t)=T_ini;
e_cabin = ones(1,Total_time)*e_amb(1);

    % MAC Thermal Load parameters:
Ttarget = 22 + 273.15; %K
T0cabin= (Tcabin_front(1)+Tcabin_back(1))/2; %K
ttarget = ones(1,Total_time+1)* 900; %s, time to reach Ttarget from T0cabin -> TBD in the loop
tcons= ttarget/log(abs(T0cabin-Ttarget)); % pull-down constant: overall pull-down time

Qmac = ones(1,Total_time)*(CabinVolume*density_air*cp_air + ...
    mass_seats*cp_seats+144240)*(Tcabin(t)-Ttarget)/tcons(t); 
% Formula: Qmac_needed(t)= (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*(Tcabin(t)-Ttarget)/tcons;

    % Initial Overall heat transfer coefficients, UA (W/K):
UA_front(t) = (k_ws*A_ws)/e_ws + 0.5*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;
UA_back(t) = (k_rw*A_rw)/e_rw + 0.5*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;
UA_front=ones(1,Total_time).*UA_front(t);
UA_back=ones(1,Total_time).*UA_back(t);
U_front=ones(1,Total_time);
U_side=ones(1,Total_time);
U_back=ones(1,Total_time);

A_skin = 1.5; %m2
T_skin = 32 + 273.15; %K prev: 24.648 C %human skin is 4C warmer than ambient around 28C

h_engine = h_ext; % h_engine = f( velocidad aire que pasa por el motor )
A_engine = 0.01*A_ws; %Faya 2013: *0.01
A_exhaust = 0.01 * cabin_width*cabin_base_lenght; %Faya 2013
e_enginesheet = 0.001; %engine sheet thickness, m
e_base = 0.1; %m
U_engine = 1/( 1/h_engine + 1/h_cabin + e_enginesheet/k_steel );
U_exhaust = 1/( 1/h_engine + 1/h_cabin + e_base/k_steel );
% U_engine = h_ext;
% U_exhaust = h_cabin;

compressor_on_lower_limit = 0.5; % W, 


for t=1:Total_time
        
        % Overall heat transfer coefficients, UA (W/K):
    if t>1800-323 % extra high
        h_front(t) = 200;
        U_front(t) = 20;
    elseif t>1800-323-455 % high
        h_front(t) = 190;
        U_front(t) = 19;
    elseif t>1800-323-455-433 % medium
        h_front(t) = 170;
        U_front(t) = 17;
    else % low
        h_front(t) = 160;
        U_front(t) = 16;
    end

    h_side(t) = 0.9*h_front(t); % 0.95
    h_rear(t) = 0.1*h_front(t);
    U_side(t) = 0.9*U_front(t);
    U_back(t) = U_front(t);

    UA_front(t) = A_front*U_front(t) + 1*(A_side*U_side(t));
    UA_back(t) = A_back*U_back(t) + 1*(A_side*U_side(t));

    % Conductances Matrix 
    % Tamb --K1-> Troof --K2-> Tceiling --K3-> Tcabin_back
    K1=h_side(t)*A_roof; %h_roof?
    K2=A_roof*k_ceiling/e_ceiling;
    K3=h_cabin*A_roof;
    % Tbase --K4-> Tcabin_back
    K4=h_cabin*A_base;
    % Tamb --K5-> Tw_ext_ws --K6-> Tw_int_ws --K7-> Tcabin_front
    K5=h_front(t)*A_ws;
    K6=A_ws*k_ws/e_ws;
    K7=h_cabin*A_ws;
    % Tamb --K8-> Tw_ext_rw --K9-> Tw_int_rw --K10-> Tcabin_back
    K8=h_rear(t)*A_rw;
    K9=A_rw*k_rw/e_rw;
    K10=h_cabin*A_rw;
    % Tamb --K11-> Tw_ext_sidewindows --K12-> Tw_int_sidewindows --K13-> Tcabin_back
    K11=h_side(t)*A_sidewindows;
    K12=A_sidewindows*k_sidewindows/e_sidewindows;
    K13=h_cabin*A_sidewindows;
    % Tamb --K14-> Tw_ext_doors --K15-> Tw_int_doors --K16-> Tcabin_back
    K14=h_side(t)*A_doors;
    K15=A_doors*k_doors/e_doors;
    K16=h_cabin*A_doors;
    % Tcabin_front --K17-> Tcabin_back
    K17=h_cabin*cabin_width*cabin_height;%A_ws; %TBD

    UAroof = 1/( 1/K1 + 1/K2 + 1/K3 );
    UAbase = 1/( 1/K4 );
    UAws = 1/( 1/K5 + 1/K6 + 1/K7 );
    UArs = 1/( 1/K8 + 1/K9 + 1/K10 );
    UAsw = 1/( 1/K11 + 1/K12 + 1/K13 );
    UAdoors = 1/( 1/K14 + 1/K15 + 1/K16 );
    UAcabin = 1/( 1/K17 );

    Kfront = UAws + UAroof/2 + UAbase/2 + UAsw/2 + UAdoors/2;
    Kback = UArs + UAroof/2 + UAbase/2 + UAsw/2 + UAdoors/2;
    Kcabin = UAcabin;
    
    K = [1 0 0;
        Kfront -(Kfront+Kcabin) Kcabin;
        Kback Kcabin -(Kback+Kcabin)];

        % Cabin air enthalpy
    e_cabin(t) = 1006*Tcabin(t) + (2501000 + 1770*Tcabin(t))*X; % J/kg, Faya 2013
        
        % Engine and exhaust temperatures
    if(icev == 1)
        T_engine(t) = 273.15 + (-2*10^(-6)*engine_rpm(t)^(2) + 0.0355*engine_rpm(t) + 77.5); %K
        T_exhaust(t) = 273.15 + 0.138*engine_rpm(t) - 17;
            % Trials:
        % T_engine(t) = 273.15 - 2*10^(-6)*engine_rpm(t)^(2) + 0.08*engine_rpm(t) + 27.5; %K
        % T_exhaust(t) = 273.15 - 2*10^(-6)*engine_rpm(t)^(2) + 0.138*engine_rpm(t) - 17;

    elseif(bev == 1)
        T_engine(t) = 273.15 + engine_rpm(t)^(1.5)*3.07*10^(-2.8); % BEV case: EM: 150C. % 2600/(150+273.15)=6.1444
        T_exhaust(t) = 273.15 + engine_rpm(t)*0.03; % BEV case: HV bat 40C

    elseif(phev_cd == 1)
        T_engine(t) = 273.15 + TC_cell(t); 
        T_exhaust(t) = 273.15 + TC_cell(t);
    end

        % Heat Flows
    Qleakage(t) = - leakage_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    Qvent(t) = vent_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    Qbase(t) = alfa_base*(A_ws*G_ws_t(t)+A_rw*G_rw_t(t)+A_sidewindows*G_sidewindows_t(t)+A_sidewindows*G_doors_t(t));
    Qirr(t) = -A_roof*G_roof_a(t)-A_roof*G_roof_a(t)-A_rw*G_rw_a(t)-A_sidewindows*G_sidewindows_a(t)-A_doors*G_doors_a(t);
    heat_human(t) = h_cabin*A_skin*(Tcabin_front(:,t)-T_skin) ;
    Qhuman(t) = N_Humans.*heat_human(t);
    Qequipment(t) = heat_equipment;
    Qengine_calc(t) = A_engine*U_engine*(T_engine(t)-Tcabin_front(:,t));% - 11 - 44;
    Qexhaust_calc(t) = A_exhaust*U_exhaust*(T_exhaust(t)-Tcabin_back(:,t));%  + 17.77;

    Qengine(t) = Qengine_calc(t) ;
    Qexhaust(t) = Qexhaust_calc(t) ;

        % Boundary conditions vector
    % if(size(T_cell)==[1 1])
    %     Tbc=[T_cell;-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);...
    %         -A_sidewindows*G_sidewindows_a(t);-A_doors*G_doors_a(t);0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)];
    % else 
    %     Tbc=[T_cell(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_roof*G_roof_a(t);-A_rw*G_rw_a(t);...
    %     -A_sidewindows*G_sidewindows_a(t);-A_doors*G_doors_a(t);0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)];
    % end
    
    Tbc=[T_cell(t);
        -Qhuman(t)/2-Qirr(t)/2-Qbase(t)/2-Qengine(t)-Qvent(t); %-Qhuman(t)-Qengine(t)-Qvent(t)-Qirr(t)/2-Qbase(t)/2;
        -Qhuman(t)/2-Qirr(t)/2-Qbase(t)/2-Qexhaust(t)-Qleakage(t)]; %-Qequipment(t)-Qexhaust(t)-Qleakage(t)-Qirr(t)/2-Qbase(t)/2];

        % Temperatures calculation
    temperature=(inv(K - C))*(Tbc - C*prev_temp);
    prev_temp=temperature;

    %     % Exchange from the front and back:
    % if ( temperature(13)>T_cell(t) )
    %     temperature(13) = (density_air*CabinVolume/2*cp_air*prev_temp(13)/timestep + UA_front(t)*T_cell(t))/( density_air*CabinVolume/2*cp_air/timestep + UA_front(t) );
    %     prev_temp(13)=temperature(13);
    %     check_f=check_f+1;
    %     Qcv_emitted(t) = (temperature(13)-T_cell(t))*UA_front(t);
    %     Qcv_received(t) = 0;
    % else
    %     Qcv_emitted(t) = 0;
    %     Qcv_received(t) = (T_cell(t)-temperature(13))*UA_front(t);
    % end
    % 
    % if ( temperature(14)>T_cell(t) )
    %     temperature(14) = (density_air*CabinVolume/2*cp_air*prev_temp(14)/timestep + UA_back(t)*T_cell(t))/(density_air*CabinVolume/2*cp_air/timestep + UA_back(t));
    %     prev_temp(14)=temperature(14);
    %     check_b=check_b+1; 
    %     % Qcv_emitted(t) = Qcv_emitted(t)-(temperature(14)-T_cell(t))*UA_back(t);
    %     Qcv_emitted(t) = (temperature(14)-T_cell(t))*UA_back(t);
    %     Qcv_received(t) = 0;
    % else
    %     Qcv_emitted(t) = 0;
    %     %Qcv_received(t) = Qcv_received(t)+(T_cell(t)-temperature(14))*UA_back(t);
    %     Qcv_received(t) = (T_cell(t)-temperature(14))*UA_back(t);
    % end

    Tamb(:,t)=temperature(1); %T_amb: boundary condition, Tamb: simulated
    Tcabin_front(:,t)=temperature(2);
    Tcabin_back(:,t)=temperature(3);
    Tcabin(:,t)= (Tcabin_front(:,t)+Tcabin_back(:,t))/2;

    Qcv_emitted(t) = (Tcabin(:,t)-T_cell(t))*UA_front(t);
    Qcv_received(t) = (T_cell(t)-Tcabin(:,t))*UA_back(t);
    
        %% Thermal Loads
    Qcabin_received(t) = Qhuman(t) + Qequipment(t) + Qengine(t) + Qexhaust(t) + Qcv_received(t) + Qvent(t) + Qleakage(t);
    Qtotal(t) = Qcabin_received(t) + Qcv_emitted(t);
    ttarget(t) = 600;
    tcons= ttarget/abs(log(Ttarget-T0cabin));

    % Qtarget(t) = (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*(Ttarget-Tcabin(t))/tcons(t); %negative in cooling mode
    % if(heatpump==0)
    %     if(Tcabin(t)>Ttarget)
    %         % MAC heat needed to reach the T_target:
    %         Qtarget(t) = (CabinVolume*density_air*cp_air + mass_seats*cp_seats)*(Ttarget-Tcabin(t))/(ttarget(t));
    %         % Heat in the cabin that will need to be provided by MAC:
    %         Qmac_needed(t) = Qcabin_received(t) + Qcv_emitted(t) + Qtarget(t);
    %     else
    %         Qtarget(t) = 0;
    %         Qmac_needed(t) = 0;
    %     end
    % elseif(heatpump==1)
        Qtarget(t) = (CabinVolume*density_air*cp_air + mass_seats*cp_seats + 144240)*(Ttarget-Tcabin(t))/(ttarget(t));
        Qmac_needed(t) = Qcabin_received(t) + Qcv_emitted(t) + Qtarget(t);
    % end
    Qcompressor(t) = Qmac_needed(t)/COP;
    comp_cum_Wh(t) = cumtrapz(timestep.*Qcompressor(t)./3600);
    % energy_cum_comp_needed_Wh = cumtrapz(energy_comp_needed_Wh); %
    if abs(Qcompressor(t))>compressor_on_lower_limit
        t_active(t) = 1;
    end
    %% CO2 Emissions
    compressor_time(t) = ttarget(t); % OR sum(t_active)
    CO2(t) = abs(Qcompressor(t))*compressor_time(t)*Vpe_petrol/1000/3600*CF_petrol/distance_driven;

end

tot_t_acitve=sum(t_active);
total_active_time=tot_t_acitve(end)

Qcabin_received_plot=Qcabin_received;
Qmac_needed_plot=Qmac_needed;
Qemitted=Qleakage+Qcv_emitted;

% Delete NaN values
Qcabin_received=(Qcabin_received(~isnan(Qcabin_received)));
Qmac_needed=(Qmac_needed(~isnan(Qmac_needed)));
Qcompressor=Qcompressor(2:end);
Qcv_emitted=Qcv_emitted(2:end);
Qcabin_received=Qcabin_received(2:end);

    %% Energy and power
cabin_in_W_av = sum(Qcabin_received)/Total_time;
cabin_out_W_av = sum(Qcv_emitted)/Total_time;
cabin_W_av = cabin_in_W_av - abs(cabin_out_W_av);
target_av_W = sum(Qtarget)/Total_time;
evaporator_av_W = sum(Qmac_needed)/Total_time; % MAC needed=evaporator
compressor_av_W = sum(abs(Qcompressor))/Total_time

cabin_received_Wh = timestep.*Qcabin_received(2:end)./3600;
cabin_received_cum_Wh = cumsum(cabin_received_Wh);
cabin_received_Wh_WLTP=cabin_received_cum_Wh(end);
mac_needed_Wh = timestep.*Qmac_needed(2:end)./3600; % MAC needed=evaporator
mac_needed_cum_Wh = cumsum(mac_needed_Wh);
mac_needed_Wh_WLTP=mac_needed_cum_Wh(end);
comp_Wh = timestep.*Qcompressor(2:end)./3600;
comp_cum_Wh = cumtrapz(comp_Wh); %cumsum
comp_Wh_WLTP = comp_cum_Wh(end);


%% CO2 Emissions
CO2_total = CO2(end);
CO2_perc = CO2_total/95*100;

%% Save engine load to be applied in other scenarios:
%save Engine_thermal_loads_W_35G7_v9 Qengine Qexhaust

%% PLOTS

% figure(7)
% plot(comp_Wh)
% hold on
% yyaxis right
% plot(comp_cum_Wh)
% 
% figure(8)
% plot(Qcompressor)
% hold on
% yline(abs(compressor_av_W))

figure(3)
plot(T_engine-273.15)
hold on
plot(T_exhaust-273.15)

% Deletion of first 0 value for plotting (inicialization value):
x=T_cell(2); T_cell(1)=x;
% The boundary conditon 'T_amb' needs to be a vector to be plotted:
if (size(T_cell)==[1 1])
    T_cell = ones(1,Total_time)*(T_ambient+273.15); %[K], vector
end
if off_G7_22==1
    test = ['G7_22'];
elseif off_G7_35==1
    test = ['G7_35'];
elseif off_G8_22==1
    test = ['G8_22'];
elseif off_G8_35==1
    test = ['G8_35'];
end

% figure(13)
% plot(time,Tcabin-273.15)
% % legend('T target','T cabin','Location','best')

figure(12)
plot(time(3:Total_time),Qtarget(3:Total_time),':k','LineWidth',1.5)
hold on
plot(time(3:Total_time),Qmac_needed(3:Total_time),'-k','LineWidth',1.5)
grid minor
ylabel('MAC Heat [W]')
xlabel('Time [s]')
if(off_G7_22==1)
    ylim([-1200 200])
elseif(off_G8_22==1)
    ylim([-400 200])
end
yyaxis right
ax = gca;
ax.YColor = [1 0 0];
hold on
plot(time(3:Total_time),Qtotal(3:Total_time),'-r','LineWidth',1.5)
ylabel('Total Cabin Heat [W]','Color',[1 0 0])
if(off_G8_22==1 || off_G7_22==1 || off_G7_35==1 || off_G8_35==1)
    ylim([-5 130])
end
legend('Target','MAC','Cabin Total','Location','best')
% saveas(figure(12),strcat(test,'_4'),'png')

figure(11)
plot(time(3:Total_time),Qtotal(3:Total_time),'LineWidth',1.5)
hold on
plot(time(3:Total_time),Qcabin_received_plot(3:Total_time),'LineWidth',1.5)
plot(time(3:Total_time),Qemitted(3:Total_time),'LineWidth',1.5)
ylabel('Cabin Heat Loads [W]')
xlabel('Time [s]')
legend('Cabin Total','Cabin Received','Cabin Emitted','Location','best')
grid minor
% saveas(figure(11),strcat(test,'_3'),'png')

figure(10)
plot(time,Qtotal,'--k','LineWidth',1.5)
hold on
plot(time,Qcv_received,':','LineWidth',1.5)
plot(time,Qhuman,':','LineWidth',1.5)
plot(time,Qequipment,':','LineWidth',1.5)
plot(time,Qvent,':','LineWidth',1.5)
plot(time,Qengine,':','LineWidth',1.5,'Color','g')
plot(time,Qexhaust,':','LineWidth',1.5)
plot(time,Qcabin_received_plot,'-','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
ylabel('Thermal Loads [W]')
xlabel('Time [s]')
yyaxis right
ax = gca;
ax.YColor = [0.4940 0.1840 0.5560];
%ylim([-1000 600])
hold on
plot(time,Qmac_needed_plot,'-','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560])
ylabel('Q cooling MAC [W]','Color',[0.4940 0.1840 0.5560])
%set(gca,'Ydir','reverse')
legend('Ambient','Human','Equipment','Vents','Engine',...
    'Exhaust','Cabin Received','MAC','Location','northeastoutside')
grid minor
    % xlim([0 1800])
    % axis([1 Total_time -1500 1700])
    % axis([1 Total_time 0 1100])
    
% figure(9)
% colororder({'k','k'})
% plot(time,Qengine,':b','LineWidth',1.5)
% ylabel('Engine Heat [W]')
% xlabel('Time [s]')
% ylim([-1 400])
% hold on
% yyaxis right
% plot(time,engine_rpm,'-m','LineWidth',1)
% ylabel('Engine Speed [rpm]')
% legend('Engine Heat','Engine Speed','Location', 'Best')
% grid minor
% axis([1 Total_time 0 5500]) %745

% figure(8)
% ylabel('Calculated [W]')
% xlabel('Time [s]')
% plot(time,Qengine_calc)
% hold on
% plot(time,Qexhaust_calc)
% yyaxis right
% ylabel('Observed [W]')
% plot(time,Qengine_observed)
% plot(time,Qexhaust_observed)
% legend('Qengine calc','Qexhaust calc','Qengine observed','Qexhaust observed')

% figure(7)
% plot(time,Qcabin_received_plot) 
% xlabel('Time [s]') 
% hold on
% plot(time,-Qcv_emitted)
% ylabel('Cabin Heat [W]')
% legend('Received','Emitted')

% figure(6)
% plot(time,Qcabin_received_plot-abs(Qcv_emitted)) 
% xlabel('Time [s]') 
% hold on
% % plot(time,-Qcabin_emitted)
% ylabel('Cabin Heat [W]')
% legend('Received','Emitted')

% figure(5)
% colororder({'k','k'})
% plot(time,Qcabin_received_plot,':b','LineWidth',1.5)
% grid minor
% ylabel('Heat received (ambiebt to cabin) [W]')
% xlabel('Time [s]') 
% hold on
% yyaxis right
% plot(time,Qcv_emitted,':m','LineWidth',1.5)
% set(gca,'Ydir','reverse')
% xlim([0 1800])
% ylim([-700 0])
% ylabel('Heat emitted (cabin to ambient) [W]')
% legend('Received','Emitted','Location','Best')

% figure(4)
% xlabel('Time [s]') 
% plot(time,Qmac_needed_plot)
% ylabel('Heat needed to get T target [W]')

% figure(3)
% plot(time,Qengine)
% hold on
% plot(time,Qexhaust)
% ylabel('Heat [W]')
% xlabel('Time [s]') 
% yyaxis right
% plot(time,Eeng_kWh)
% hold on
% plot(time,Etail_kWh)
% ylabel('Energy [kWh]')
% legend('Q engine W','Q tailpipe W','E engine kWh','E tailpipe kWh')

% plot(time,Q_cv_cabin,':','LineWidth',1.5)
% hold on
% plot(time,Q_cv_front,time,Q_cv_back,time,Q_cv_top)
% hold on
% plot(time,Q_cv_left)
% legend('Cabin','Front','Back','Top','Left/Right','Location','best')
% xt=[1000 800 1000 1000 1000];
% yt = [Q_cv_front(end) Q_cv_back(end) Q_cv_left(end) Q_cv_top(end) Q_cv_cabin(end)];
% str = [string(Q_cv_front(end)) string(Q_cv_back(end)) string(Q_cv_left(end)) string(Q_cv_top(end)) string(Q_cv_cabin(end))];
% text(xt,yt,str)

figure(2)
plot(time,TC_cell,'-k','LineWidth',1) %T amb 
hold on
plot(time,TC_front,'-','LineWidth',1.5) 
plot(time,TC_back,'-','LineWidth',1.5)
plot(time,TC_vent,'-','LineWidth',1.5)
plot(time,TC_driver,'-','LineWidth',1.5)
plot(time,TC_backpass,'-','LineWidth',1.5)
hold off
legend('TC cell','TC front','TC back','TC vent','TC driver','TC equipment','Location','best')
grid minor
xlabel('Time [s]') 
ylabel('Temperature [ºC]')
if(off_G7_22==1 || off_G8_22==1)
    axis([1 Total_time 20 25])
elseif(off_G7_35==1 || off_G8_35==1)
    axis([1 Total_time 32 38])
end
% saveas(figure(2),strcat(test,'_2'),'png')

figure(1)
plot(time,Tcabin_front-273.15,':b','LineWidth',1.5)
hold on
plot(time,Tcabin_back-273.15,':m','LineWidth',1.5)
plot(time,TC_cell,'-k','LineWidth',1) %T amb 
plot(time,TC_front,'-b','LineWidth',1.5) 
plot(time,TC_back,'-m','LineWidth',1.5)
hold off
legend('T front simulated','T back simulated','TC ambient','TC front','TC back','Location','best')
grid minor
xlabel('Time [s]') 
ylabel('Temperature [ºC]')
% if(G7_22==1 || G8_22==1)
%     axis([1 1800 20 25])
% elseif(G7_35==1 || G8_35==1)
%     axis([1 1800 33 40])
% end     
% saveas(figure(1),strcat(test,'_1'),'png')

