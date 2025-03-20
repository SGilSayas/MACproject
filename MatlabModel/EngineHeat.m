% Engine and exhaust Heat Power characterisation, Susana Gil-Sayas 2023
% Data from Golf 7 tests at 35C, MAC off, 
% 'Vela2-28_04_2022-2-VW054_35-WLTC Class3b-HOT.xlsx'
clc, clear all, close all

% Load calculated heat powers
engineheatinput='Qengine_Qtailpipe_W_35G7_v6';
struct=load(engineheatinput);
table=struct2table(struct);
arr=table2array(table);
Qengine=arr(:,1);
Qtailpipe=arr(:,2);

% Load lab cell thermocouples
file_name2='Vela2-28_04_2022-2-VW054_35-WLTC Class3b-HOT.xlsx';
Continuous=readtable(file_name2,'Sheet','Continuous','ReadVariableNames',true);
TC_Cell=Continuous(:,36);
TCa_t=Continuous(:,40);
TCb_t=Continuous(:,41);
TCc_t=Continuous(:,42);
TCd_t=Continuous(:,43);
TCe_t=Continuous(:,44);
TCf_t=Continuous(:,45);
TC_Cell=table2array(TC_Cell);
TCa=table2array(TCa_t);
TCb=table2array(TCb_t);
TCc=table2array(TCc_t);
TCd=table2array(TCd_t);
TCe=table2array(TCe_t);
TCf=table2array(TCf_t);
n=10;
TC_amb=TC_Cell(1:n:end);
TC_a=TCa(1:n:end);
TC_b=TCb(1:n:end);
TC_c=TCc(1:n:end);
TC_d=TCd(1:n:end);
TC_e=TCe(1:n:end);
TC_f=TCf(1:n:end);

% Load engine data
file_name1='XCU_28_04_2022-2-VW54_35.xlsx';
XCUt=readtable(file_name1,'Sheet','All signals','ReadVariableNames',true);
XCU=table2array(XCUt);
engineSpeed=XCU(:,5);
fuelConsumption=XCU(:,6);
ICETorque=XCU(:,12);
engCoolTemp=XCU(:,24);
testingtime=XCU(:,1);
% delete NaN
engineSpeed_=engineSpeed(~isnan(engineSpeed));
fuelConsumption_=fuelConsumption(~isnan(engineSpeed));
ICETorque_=ICETorque(~isnan(engineSpeed));
engCoolTemp_=engCoolTemp(~isnan(engineSpeed));
testingtime_=testingtime(~isnan(engineSpeed));
Qengine_=Qengine(~isnan(engineSpeed));
% interpolation
timestep=1;
time = 1:timestep:1800;

Qengine_mesh=interp1(testingtime_,Qengine_,time,'linear','extrap');
engineSpeed_mesh = interp1(testingtime_,engineSpeed_,time,'linear','extrap');
ICETorque_mesh = interp1(testingtime_,ICETorque_,time,'linear','extrap');

Ones = ones(1800,1800);
Qengine_matrix= Qengine_mesh.*Ones;
mesh(engineSpeed_mesh,ICETorque_mesh,Qengine_matrix)

% %% Time series
% 
% figure(1)
% plot(Qengine)
% hold on
% plot(Qtailpipe)
% ylabel('Heat [W]')
% xlabel('Time [s]')
% yyaxis right
% plot(engineSpeed)
% ylabel('Engine Speed rpm')
% legend('Q engine W','Q tailpipe W','Engine Speed')
% 
% figure(2)
% plot(Qengine)
% hold on
% plot(Qtailpipe)
% ylabel('Heat [W]')
% xlabel('Time [s]')
% yyaxis right
% plot(time,fuelConsumption)
% ylabel('Fuel consumption')
% legend('Q engine W','Q tailpipe W','Fuel consumption')
% 
% figure(3)
% plot(Qengine)
% hold on
% plot(Qtailpipe)
% ylabel('Heat [W]')
% xlabel('Time [s]')
% yyaxis right
% plot(ICETorque)
% ylabel('Engine torque')
% legend('Q engine W','Q tailpipe W','Engine torque')
% 
% figure(4)
% plot(Qengine)
% hold on
% plot(Qtailpipe)
% ylabel('Heat [W]')
% xlabel('Time [s]')
% yyaxis right
% plot(engCoolTemp)
% ylabel('Engine coolant temperature [C]')
% legend('Q engine W','Q tailpipe W','Engine Coolant Temp C')
