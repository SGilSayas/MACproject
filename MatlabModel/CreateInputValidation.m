% Script creates input for ValidationJRCcabinmodel14nodes and ValidationJRC_CabinEnergyLoads, 
% Susana Gil-Sayas 2023-2024-2025
clc, clear all, close all


%% Hyundai Tucson PHEV, 2022 campaign
file_tucphev_MACoff = '35C_macOFF_Vela2-19_05_2022-1-HI006-WLTC Class3b-COLD.xlsx'; % CD, cold start
    thermocouples_tucsonphev = readtable(file_tucphev_MACoff,'Sheet','Continuous','ReadVariableNames',true);
    time_TC = thermocouples_tucsonphev(:,2);        % timeseries
    TC_vent = thermocouples_tucsonphev(:,20);       % vent           % Beckhoff BK5120 4DI 4DO 28AI 2AO AI13
    TC_front = thermocouples_tucsonphev(:,21);      % front ceiling  % Beckhoff BK5120 4DI 4DO 28AI 2AO AI14
    TC_driver = thermocouples_tucsonphev(:,22);     % driver         % Beckhoff BK5120 4DI 4DO 28AI 2AO AI15
    TC_pass = thermocouples_tucsonphev(:,23);       % passanger      % Beckhoff BK5120 4DI 4DO 28AI 2AO AI16
    TC_back = thermocouples_tucsonphev(:,24);       % back ceiling   % Beckhoff BK5120 4DI 4DO 28AI 2AO AI17
    TC_backpass = thermocouples_tucsonphev(:,25);   % back seat      % Beckhoff BK5120 4DI 4DO 28AI 2AO AI18
    TC_cell = thermocouples_tucsonphev(:,16);       % cell           % CellTemperature
    P_cell_hPa = thermocouples_tucsonphev(:,17);    % cell pressure  % Barometer
    save tucsonphev_35C_MACoff_TC time_TC TC_vent TC_front TC_driver TC_pass TC_back TC_backpass TC_cell

%% Ford Mustang Mach-e BEV, 2021 camapign
% file_mustangbev_MACoff = 'FordMustangMachE_2021-09-07_35cold_macOFF.xlsx'; % on: _____.xlsx
    % thermocouples_mustangbev = readtable(file_mustangbev_MACoff,'Sheet','Data1','ReadVariableNames',true);
    % time_TC = thermocouples_mustangbev(:,2); % timeseries
    % TC_vent = thermocouples_mustangbev(:,10); % vent
    % TC_front = thermocouples_mustangbev(:,11); % front ceiling
    % TC_driver = thermocouples_mustangbev(:,12); % driver
    % TC_pass = thermocouples_mustangbev(:,13); % passanger
    % TC_back = thermocouples_mustangbev(:,14); % back ceiling
    % TC_backpass = thermocouples_mustangbev(:,15); % back seat
    % TC_cell = thermocouples_mustangbev(:,25); % cell:38 and 25
    % save mustangbev_35C_MACoff_TC time_TC TC_vent TC_front TC_driver TC_pass TC_back TC_backpass TC_cell

%% Hyundai Ioniq 5 with PV, 2024 camapign
% on: V8_240126_04_HI009_WLTC_COLD_T35_MACoN_NM.xlsx
% file_ioniq_MACoff = 'V8_240126_01_HI009_1A_WLTC_COLD_T35_MACoFF_NC.xlsx';
    % 
    % thermocouples_ioniq = readtable(file_ioniq_MACoff,'Sheet','Modal Analog','ReadVariableNames',true);
    % time_TC = thermocouples_ioniq(:,1); % timeseries
    % TC_vent = thermocouples_ioniq(:,4); % vent
    % TC_front = thermocouples_ioniq(:,8); % front ceiling
    % TC_driver = thermocouples_ioniq(:,7); % driver
    % TC_pass = thermocouples_ioniq(:,9); % passanger
    % TC_back = thermocouples_ioniq(:,11); % back ceiling
    % TC_backpass = thermocouples_ioniq(:,12); % back seat
    % TC_cell = thermocouples_ioniq(:,2); % cell
    % save ioniq5PV_35C_MACoff_TC time_TC TC_vent TC_front TC_driver TC_pass TC_back TC_backpass TC_cell
    % 
    % OBD_ioniq = readtable(file_ioniq_MACoff,'Sheet','Modal XCU','ReadVariableNames',true);
    % % NOTE: OBD's outdoor and indoor temperatures are switched -> fixed:
    % time_OBD = OBD_ioniq(:,1);
    % Tevap_C = OBD_ioniq(:,8);
    % Tbat1_C = OBD_ioniq(:,16);
    % Tbat2_C = OBD_ioniq(:,17);
    % Tbat3_C = OBD_ioniq(:,18);
    % Tbat4_C = OBD_ioniq(:,19);
    % Tbat5_C = OBD_ioniq(:,20);
    % Tinterior_C = OBD_ioniq(:,3);
    % Toutdoor_C = OBD_ioniq(:,4);
    % Tws_C = OBD_ioniq(:,6);
    % Tfrontfeet_C = OBD_ioniq(:,7);
    % SOC = OBD_ioniq(:,9);
    % save ioniq5PV_35C_MACoff_OBD time_OBD Tevap_C Tbat1_C Tbat2_C Tbat3_C Tbat4_C Tbat5_C Tinterior_C Toutdoor_C Tws_C Tfrontfeet_C SOC 
    % 
    % PA_ioniq = readtable(file_ioniq_MACoff,'Sheet','Modal Power Analyzer','ReadVariableNames',true);
    % time_PA = PA_ioniq(:,1); % timeseries of power analyser
    % HVbat_A = PA_ioniq(:,3); % HV battery
    % convLVside_A = PA_ioniq(:,12); % DC/DC to 12V side 
    % auxplus_A = PA_ioniq(:,21); % AUX HV (+) 
    % auxminus_A = PA_ioniq(:,31); % AUX HV (-) 
    % compressorplus_A = PA_ioniq(:,40); % AC Compressor (+) 
    % compressorminus_A = PA_ioniq(:,49); % AC Compressor (-) 
    % PTCcabinplus_A = PA_ioniq(:,58); % Cabin PTC Heater (+) 
    % PTCcabinminus_A = PA_ioniq(:,69); % Cabin PTC heater (-) 
    % PTCbatplus_A = PA_ioniq(:,79); % Battery heater (+) 
    % PTCbatminus_A = PA_ioniq(:,88); % Battery heater (-) 
    % save ioniq5PV_35C_MACoff_PA time_PA HVbat_A convLVside_A auxplus_A ...
    %     auxminus_A compressorplus_A compressorminus_A PTCcabinplus_A ...
    %     PTCcabinminus_A PTCbatplus_A PTCbatminus_A

%% Engine speed and torque, G7 35C MAC off
% Total_time = 1800;
% lab_input2 = 'G7_35C_MACoff_XCU'; 
    % struct2 = load(lab_input2);
    % table2 = struct2table(struct2);
    % modalxcu_array = table2array(table2);
    % VehicleSpeed = table2array(modalxcu_array(:,1));
    % ACcomp_speed_rpm = table2array(modalxcu_array(:,2));
    % ACcomp_torque_Nm = table2array(modalxcu_array(:,3));
    % T_Cell_XCU = table2array(modalxcu_array(:,4));
    % T_evap_C = table2array(modalxcu_array(:,5));
    % ICETorque_Nm = table2array(modalxcu_array(:,6));
    % engineSpeed_rpm = table2array(modalxcu_array(:,7));
    % fuelconsumption_lph = table2array(modalxcu_array(:,8));
    % refrigewrantFluidPress_bar = table2array(modalxcu_array(:,9));
    % exaustGasTempBank1_C = table2array(modalxcu_array(:,10));
    % n=height(modalxcu_array)/Total_time;
    % vehicle_kmh = VehicleSpeed(1 : n : end);
    % ACcomp_rpm = ACcomp_speed_rpm(1 : n : end);
    % ACcomp_Nm = ACcomp_torque_Nm(1 : n : end);
    % engine_Nm = ICETorque_Nm(1 : n : end);
    % engine_rpm = engineSpeed_rpm(1 : n : end);
    % FC_lperh = fuelconsumption_lph(1 : n : end);
    % % Total_time2 = length(vehicle_kmh)-1; 
    % % time2 = 0:timestep:Total_time2;
    % ACcomp_W = ACcomp_Nm.*pi.*ACcomp_rpm/60;
    % save G7_35_off_rpm_Nm_FC engine_rpm engine_Nm FC_lperh

%% RH and P air in G7 35C MAC off
% lab_input3 = 'G7_35C_MACoff_Modaldata';
    % Total_time = 1800;
    % struct3NaN = load(lab_input3);
    % struct3NoNaN = structfun( @rmmissing , struct3NaN , 'UniformOutput' , false);
    % table3 = struct2table(struct3NoNaN);
    % modaldata_array = table2array(table3);
    % RH_amb = table2array(modaldata_array(:,1));
    % P_amb_kPa = table2array(modaldata_array(:,2));
    %     % Reshape RH and P data
    % extra1800RH = length(RH_amb)-Total_time;
    % extra1800P = length(P_amb_kPa)-Total_time;
    % RH_amb = RH_amb(1:end-extra1800RH,:);
    % P_amb_kPa = P_amb_kPa(1:end-extra1800P,:);
    % save G7_35_off_RH_P_amb RH_amb P_amb_kPa

%% HYUNDAI IONIQ 5 BEV
% 14C, MAC off
% file_name_XCU= 'XCU1.xlsx';
    % Modal_XCU=readtable(file_name_XCU,'Sheet','XCU','ReadVariableNames',true);
    % TC1=Modal_XCU(:,21);
    % TC2=Modal_XCU(:,22);
    % TC3=Modal_XCU(:,23);
    % TC4=Modal_XCU(:,24);
    % TC5=Modal_XCU(:,19);
    % TC6=Modal_XCU(:,20);
    % TC_Cell=Modal_XCU(:,11); %T CABIN!!
    % save HI5_14C_MACoff_XCU TC1 TC2 TC3 TC4 TC5 TC6 TC_Cell
    
    % 14C, MAC off

%% Humidities Golf 8
% file_vela_G835='Vela2-20_10_2021-1-VW057_20210630-WLTC Class3b-COLD.xlsx';
    % Continuous=readtable(file_vela_G835,'Sheet','Continuous','ReadVariableNames',true);
    % RH_G835=Continuous(:,18);
    % 
    % file_vela_G822='Vela2-29_10_2021-1-VW057_20210630-WLTC Class3b-COLD.xlsx';
    % Continuous=readtable(file_vela_G822,'Sheet','Continuous','ReadVariableNames',true);
    % RH_G822=Continuous(:,18);
    % 
    % file_vela_G735='VW054_WLTC_35.xlsx';
    % Modal_data=readtable(file_vela_G735,'Sheet','Modal data','ReadVariableNames',true);
    % time_Modaldata=Modal_data(:,1);
    % RH_G722=Modal_data(:,104);
    % 
    % file_vela_G722='';
    % Continuous=readtable(,'Sheet','Continuous','ReadVariableNames',true);
    % =Continuous(:,18);
    % 
    % save G8_humidities RH_35off RH_22off

%% Golf 8, 35C, MAC off
% Vela2-20_10_2021-1-VW057_20210630-WLTC Class3b-COLD
    % file_name_TC='Vela2-20_10_2021-1-VW057_20210630-WLTC Class3b-COLD.xlsx';
    % Continuous=readtable(file_name_TC,'Sheet','Continuous','ReadVariableNames',true);
    % TC1=Continuous(:,21);
    % TC2=Continuous(:,22);
    % TC3=Continuous(:,23);
    % TC4=Continuous(:,24);
    % TC5=Continuous(:,25);
    % TC6=Continuous(:,26);
    % TC_Cell=Continuous(:,19);
    % save G8_35C_MACoff_TC TC1 TC2 TC3 TC4 TC5 TC6 TC_Cell
% Cabin TC: 'XCU1.xlsx'
    % file_name_XCU= 'XCU1.xlsx';
    % Modal_XCU=readtable(file_name_XCU,'Sheet','XCU1','ReadVariableNames',true);
    % TC1=Modal_XCU(:,35);
    % TC2=Modal_XCU(:,36);
    % TC3=Modal_XCU(:,37);
    % TC4=Modal_XCU(:,32);
    % TC5=Modal_XCU(:,33);
    % TC6=Modal_XCU(:,34);
    % TC_Cell=Modal_XCU(:,26); %T CABIN!!
    % save G8_35C_MACoff_XCU TC1 TC2 TC3 TC4 TC5 TC6 TC_Cell

%% Golf 8, 22C, MAC off 
% file 1:Vela2-29_10_2021-1-VW057_20210630-WLTC Class3b-COLD(no cabin TC)
    % file_name_TC='Vela2-29_10_2021-1-VW057_20210630-WLTC Class3b-COLD.xlsx';
    % Continuous=readtable(file_name_TC,'Sheet','Continuous','ReadVariableNames',true);
    % TC_Cell=Continuous(:,19);
    % save G8_22C_MACoff_TCcell1 TC_Cell
    % 
    % file_name_XCU= 'XCU_CS_hot_start_heating_off_AC_off_cabinmodelpaper22G8.xlsx' ; 
    % %old: 'XCU_29_10_2021-1.xlsx';
    % XCU=readtable(file_name_XCU,'Sheet','XCU','ReadVariableNames',true);
    % TC1 = XCU(:,34);
    % TC2 = XCU(:,35);
    % TC3 = XCU(:,36);
    % TC4 = XCU(:,31);
    % TC5 = XCU(:,32);
    % TC6 = XCU(:,33);
    % save G8_22C_MACoff_TCcabin1 TC1 TC2 TC3 TC4 TC5 TC6 
    
    % file_name_XCU= 'XCU_CS_hot_start_heating_off_AC_off_cabinmodelpaper22G8.xlsx';
    % XCU=readtable(file_name_XCU,'Sheet','XCU','ReadVariableNames',true);
    % engine_rpm = XCU(:,3);
    % save G8_22C_MACoff_engineRPM engine_rpm

% file 3: 'Vela2-29_10_2021-3-VW057_20210630-WLTC Class3b-HOT.xlsx'
    % file_name_TC='Vela2-29_10_2021-3-VW057_20210630-WLTC Class3b-HOT.xlsx';
    % Continuous=readtable(file_name_TC,'Sheet','Continuous','ReadVariableNames',true);
    % TC_Cell=Continuous(:,19);
    % save G8_22C_MACoff_TCcell3 TC_Cell
    % file_name_XCU= 'XCU_29_10_2021-3.xlsx';
    % XCU=readtable(file_name_XCU,'Sheet','XCU','ReadVariableNames',true);
    % TC1 = XCU(:,35);
    % TC2 = XCU(:,36);
    % TC3 = XCU(:,37);
    % TC4 = XCU(:,32);
    % TC5 = XCU(:,33);
    % TC6 = XCU(:,34);
    % save G8_22C_MACoff_TCcabin3 TC1 TC2 TC3 TC4 TC5 TC6 

% file_lab = 'x';
    % Modal_Analog=readtable(file_lab,'Sheet','Modal Analog','ReadVariableNames',true);
    % TC1=Modal_Analog(:,11);
    % TC2=Modal_Analog(:,12);
    % TC3=Modal_Analog(:,13);
    % TC4=Modal_Analog(:,14);
    % TC5=Modal_Analog(:,15);
    % TC6=Modal_Analog(:,16);
    % TC_Cell=Modal_Analog(:,8);
    % Modal_XCU=readtable(file_V8,'Sheet','Modal XCU','ReadVariableNames',true);
    % save x_xC_MACoff_TClab TC1 TC2 TC3 TC4 TC5 TC6 TC_Cell

%% Golf 7
    % Golf 7, 35C, MAC off: 'Vela2-28_04_2022-2-VW054_35-WLTC Class3b-HOT.xlsx';
    % Golf 7, 35C, MAC on: 'Vela2-28_04_2022-4-VW054_35-WLTC Class3b-HOT.xlsx';
    % Golf 7, 22C, MAC off: 'LAB_20220511_05_MAC_off.xlsx' = V8_220511_05_VW054_WLTP_WARM_MAC_OFF_22C.xlsx
    % Golf 7, 22C, MAC on: LAB_20220511_06_MAC_on

%% Golf 7, 35C, MAC on, hot: 'VW054_WLTC_35.xlsx'
    % file_name_XCU= 'VW054_WLTC_35.xlsx';
    % Modal_XCU=readtable(file_name_XCU,'Sheet','Modal XCU','ReadVariableNames',true);
    % time_ModalXCU=Modal_XCU(:,1);
    % %%not recorded: VehicleSpeed=Modal_XCU(:,15);
    % ACcomp_speed_rpm=Modal_XCU(:,5);
    % ACcomp_torque_Nm=Modal_XCU(:,11);
    % T_Cell_XCU=Modal_XCU(:,14);
    % T_evap_C=Modal_XCU(:,8);
    % ICETorque_Nm=Modal_XCU(:,27);
    % engineSpeed_rpm=Modal_XCU(:,34);
    % fuelconsumption_lph=Modal_XCU(:,33);
    % refrigewrantFluidPress_bar=Modal_XCU(:,31);
    % exaustGasTempBank1_C=Modal_XCU(:,32);
    % save G7_35C_MACoff_XCU VehicleSpeed ACcomp_speed_rpm ACcomp_torque_Nm T_Cell_XCU T_evap_C ICETorque_Nm engineSpeed_rpm fuelconsumption_lph refrigewrantFluidPress_bar exaustGasTempBank1_C time_ModalXCU
    % 
    % Modal_data=readtable(file_name_XCU,'Sheet','Modal data','ReadVariableNames',true);
    % time_Modaldata=Modal_data(:,1);
    % RH=Modal_data(:,104);
    % Pressure_kPa=Modal_data(:,103);
    % save G7_35C_MACoff_Modaldata RH Pressure_kPa

%% Golf 7, 22C, MAC off: V8_220511_05_VW054_WLTP_WARM_MAC_OFF_22C.xlsx
% Golf 7, 23-28C, MAC on: 'VW054_WLTC_28.xlxs' not valid test
    % file_V8 = 'V8_220511_05_VW054_WLTP_WARM_MAC_OFF_22C.xlsx';
    % Modal_Analog=readtable(file_V8,'Sheet','Modal Analog','ReadVariableNames',true);
    % TC1=Modal_Analog(:,11);
    % TC2=Modal_Analog(:,12);
    % TC3=Modal_Analog(:,13);
    % TC4=Modal_Analog(:,14);
    % TC5=Modal_Analog(:,15);
    % TC6=Modal_Analog(:,16);
    % TC_Cell=Modal_Analog(:,8);

    % Modal_XCU=readtable(file_V8,'Sheet','Modal XCU','ReadVariableNames',true);
    % ACcomp_speed_rpm=Modal_XCU(:,5);
    % T_evap_C=Modal_XCU(:,8);
    % ACcomp_torque_Nm=Modal_XCU(:,11);
    % T_Cell_XCU=Modal_XCU(:,14);
    % VehicleSpeed=Modal_XCU(:,15);
    
    % save G7_22C_MACoff TC1 TC2 TC3 TC4 TC5 TC6 TC_Cell VehicleSpeed ACcomp_speed_rpm ACcomp_torque_Nm T_Cell_XCU T_evap_C

%% Golf 7, 28C, MAC on, hot: 'VW054_WLTC_28.xlsx'
% file_name_XCU= 'VW054_WLTC_28.xlsx';
    % Modal_XCU=readtable(file_name_XCU,'Sheet','Modal XCU','ReadVariableNames',true);
    % time_ModalXCU=Modal_XCU(:,1);
    % VehicleSpeed=Modal_XCU(:,15);
    % ACcomp_speed_rpm=Modal_XCU(:,5);
    % ACcomp_torque_Nm=Modal_XCU(:,11);
    % T_Cell_XCU=Modal_XCU(:,14);
    % T_evap_C=Modal_XCU(:,8);
    % ICETorque_Nm=Modal_XCU(:,27);
    % engineSpeed_rpm=Modal_XCU(:,34);
    % fuelconsumption_lph=Modal_XCU(:,33);
    % refrigewrantFluidPress_bar=Modal_XCU(:,31);
    % exaustGasTempBank1_C=Modal_XCU(:,32);
    % save G7_vehiclespeed_kmh time_ModalXCU VehicleSpeed
