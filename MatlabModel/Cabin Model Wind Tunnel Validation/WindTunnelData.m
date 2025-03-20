% Wind Tunnel Data, S. Gil-Sayas, 2024
clc,clear all,close all
%
WT = 1;
%% Load Wind Tunnel Data
if(WT==1)
    Ps_kPa = 5.63; % kPa, Water saturation pressure at 35C %for now
    WT_table = readtable("result_test_1.xlsx");
    TC_pass = table2array(WT_table(:,"x301_C_"));
    TC_vent = table2array(WT_table(:,"x302_C_"));
    TC_backpass = table2array(WT_table(:,"x303_C_"));
    TC_driver = table2array(WT_table(:,"x304_C_"));
    TC_front = table2array(WT_table(:,"x305_C_"));
    TC_back = table2array(WT_table(:,"x306_C_"));
    TC_cell = table2array(WT_table(:,"x307_C_"));
    timeWT = table2array(WT_table(:,"time"));
end
wt_temps=[timeWT TC_cell TC_vent TC_front TC_driver TC_pass TC_back TC_backpass];
%% Zone division
% Version slides: x1=1;x2=(751);x3=(760);x4=(911);x5=(1121);x6=(1491);x7=(1651);x8=(1961);x9=(2261);x10=length(wt_temps);
% Version1: x1=1;x2=(750);x3=(770);x4=(900);x5=(1120);x6=(1490);x7=(1650);x8=(1980);x9=(2280);x10=length(wt_temps);
% Version2:
x1=1;x2=(750);x3=(785);x4=(915);x5=(1127);x6=(1490);x7=(1646);x8=(1980);x9=(2280);x10=length(wt_temps);
x=[x1;x2;x3;x4;x5;x6;x7;x8;x9;x10];
zones(1)={'Constant temperature 20.7 C, wind off, cabin resistance 8W'};
zones(2)={'Fast cooling, wind 20 m/s, cabin resistance 8W'};
zones(3)={'Heating, wind 20 m/s, cabin resistance 8 W'};
zones(4)={'Heating, wind 25 m/s, cabin resistance 8 W'};
zones(5)={'Heating, wind 30 m/s, cabin resistance 8 W'};
zones(6)={'Cooling, wind 20 m/s, cabin resistance 8 W'};
zones(7)={'Cooling, wind 10 m/s, cabin resistance 8 W'};
zones(8)={'Cooling, wind off, cabin resistance 8 W'};
zones(9)={'Cooling, wind off, cabin resistance off'};
%% Moving mean:
% = movmean(tempconst_wind0,10);
%% Extract all data to Excel:
header = {'Time', 'TC_cell', 'TC_vent', 'TC_front', 'TC_driver', 'TC_pass', 'TC_back', 'TC_backpass'};
for i=1:9
    writecell(header,'windtunneldata.xlsx','Sheet',i,'Range','A1')
    writematrix(wt_temps(x(i):x(i+1),:),'windtunneldata.xlsx','Sheet',i,'Range','A2')
end
writecell({'Sheet 1';'Sheet 2';'Sheet 3';'Sheet 4';'Sheet 5';'Sheet 6';'Sheet 7';'Sheet 8';'Sheet 9'},'windtunneldata.xlsx','Sheet',i+1,'Range','A1')
writecell(transpose(zones),'windtunneldata.xlsx','Sheet',i+1,'Range','B1')
%% General Plot
figure(1)
plot(timeWT,TC_cell,':k','LineWidth',1.5)
hold on
plot(timeWT,TC_front,'LineWidth',1.5)
plot(timeWT,TC_back,'LineWidth',1.5)
plot(timeWT,TC_vent,'LineWidth',1.5)
plot(timeWT,TC_driver,'LineWidth',1.5)
plot(timeWT,TC_pass,'LineWidth',1.5)
plot(timeWT,TC_backpass,'LineWidth',1.5)
legend('Test Cell','Front','Back','Vent','Driver','Passenger','Back Passenger','Location','best')
grid minor
xlabel('Time [s]')
ylabel('Temperature [ºC]')
%% Extract only zones 3, 4, and 5: heating at increasing wind speeds
