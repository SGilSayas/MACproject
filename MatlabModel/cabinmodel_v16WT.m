% Vehicle's Cabin Model, Susana Gil-Sayas. Version: 15+wind tunnel
clc
clear all 
close all

%% Test conditions
Irr = -2.55; % IN NEGATIVE as it enters the cabin [W/(m2)] 
N_Humans = 1; % old: 1 human = 120 W, 0.66 = 80 W, 0.25 = 30 W;
heat_human = 80; %W
heat_equipment = 40 ; %40; % W
T_ambient = 22; % [C] % NOT USED WHEN IMPORTING DATA
duration_WLTC = 1800; % s
timestep = 1; % s
t = 1; 
COP = 3.5;
heatpump = 1; % 0 for no, 1 for yes
% What calibration? 
G7_22 = 0;
G7_35 = 0;
G8_22 = 0;
G8_35 = 0;
WT = 1;
% Deletion high peak in G8 22C TC cell data?
deletion_peak = 1; % yes/no
% Want to reduce the amplitud of T amb? % Needed for G7_22
reduced_amplitude = 0; % yes/no reduction = 5;
%% Load calculated engine and exhaust thermal load of test G7 at 35C
engineheatinput = 'Engine_thermal_loads_W_G7_35_v10'; 
load(engineheatinput);
Qengine_observed = Qengine_observed_G7_35;
Qexhaust_observed = Qexhaust_observed_G7_35;
%% Load wind tunnel data of zones 3, 4 and 5
if(WT==1)
    zone3 = readtable('windtunneldata.xlsx','Sheet','Sheet3');
    zone4 = readtable('windtunneldata.xlsx','Sheet','Sheet4');
    zone5 = readtable('windtunneldata.xlsx','Sheet','Sheet5');
    x3=2*height(zone3);
    x4=x3+2*height(zone4);
    wt_heating = [zone3;zone4;zone5];
    wt_heating=table2array(wt_heating);
    wt_heating_rs=ones(length(wt_heating)*2,8);
    j=1;
    while j<length(wt_heating)
        for k=1:length(wt_heating)*2-1 %1800
            if rem(k,2) ~= 0 % for odd (impar) k
                wt_heating_rs(k,:) = wt_heating(j,:);
            else % for even k
                wt_heating_rs(k,:) = wt_heating(j,:)+((wt_heating(j+1,:)-wt_heating(j,:))./2);
                j=j+1;
            end
            if k<=x3 % zone 3: 20 m/s = 72 km/h = medium
                vehicle_kmh(k)=72;
            elseif k>x3 && k<=x4 % zone 4: 25 m/s = 90 km/h = high
                vehicle_kmh(k)=90;
            elseif k>x4 % zone 5: 30 m/s = 108 km/h = extra high
                vehicle_kmh(k)=108;
             end
        end
    end
    TimeWT = wt_heating_rs(1:end-6,1);
    TC_cell = wt_heating_rs(1:end-6,2);
    TC_vent = wt_heating_rs(1:end-6,3);
    TC_front = wt_heating_rs(1:end-6,4);
    TC_driver = wt_heating_rs(1:end-6,5);
    TC_pass = wt_heating_rs(1:end-6,6);
    TC_back = wt_heating_rs(1:end-6,7);
    TC_backpass = wt_heating_rs(1:end-6,8);
    
    mean_TC_front=movmean(TC_front,20);
    mean_TC_back=movmean(TC_back,20);
end

%% Load input temperature
if(G7_22==1)
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
    n=height(lab_table)/(duration_WLTC); 
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
end
if(G7_35==1)
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
    n=height(lab_array)/(duration_WLTC); %10;
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
end
if(G8_22==1)
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
end
if(G8_35==1)
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
    %TC_front=TC_driver;
end

%% Filtering of Thermocouple Signals
if(G7_22==1 || G7_35==1)
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
end

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
n=height(modalxcu_array)/duration_WLTC;
ACcomp_rpm = ACcomp_speed_rpm(1 : n : end);
ACcomp_Nm = ACcomp_torque_Nm(1 : n : end);
engine_Nm = ICETorque_Nm(1 : n : end);
engine_rpm = engineSpeed_rpm(1 : n : end);
FC_lperh = fuelconsumption_lph(1 : n : end);

if WT == 0
    vehiclespeed_struct = load("G7_vehiclespeed_kmh.mat","VehicleSpeed"); 
    timespeed_struct = load("G7_vehiclespeed_kmh.mat","time_ModalXCU");
    vehiclespeed_table = struct2table(vehiclespeed_struct);
    timespeed_table = struct2table(timespeed_struct);
    vehiclespeed_array = table2array(vehiclespeed_table);
    timespeed_array = table2array(timespeed_table);
    vehiclespeed_kmh = table2array(vehiclespeed_array(:,1));
    time_vehiclespeed_s = table2array(timespeed_array(:,1));
    vehiclespeed_n = height(vehiclespeed_array)/duration_WLTC;
    time_vehiclespeed_s = time_vehiclespeed_s(1 : vehiclespeed_n : end-1);
    vehicle_kmh = vehiclespeed_kmh(1 : vehiclespeed_n : end); %%% G7@35: end-1)
    Total_time2=length(vehicle_kmh)-1; 
    time2 = 0:timestep:Total_time2;
end

%% Load engine rpms of G8 (22C case)
% if(G8_22==1 || G8_35==1)
%     lab_input4='G8_22C_MACoff_engineRPM'; 
%     struct4=load(lab_input4);
%     table4=struct2table(struct4);
%     modalxcu_array4=table2array(table4);
%     engine_rpm = table2array(modalxcu_array4(:,1));
%     n=height(modalxcu_array4)/duration_WLTC;
%     engine_rpm = engine_rpm(1 : n : end);
% end

%% Load input cell humidity and pressure
lab_input3='G7_35C_MACoff_Modaldata'; 
struct3NaN=load(lab_input3);
struct3NoNaN = structfun( @rmmissing , struct3NaN , 'UniformOutput' , false);
table3=struct2table(struct3NoNaN);
modaldata_array=table2array(table3);
RH_amb = table2array(modaldata_array(:,1));
P_amb_kPa = table2array(modaldata_array(:,2));

extra1800RH = length(RH_amb)-1801;
RH_amb = RH_amb(1:end-extra1800RH,:);
extra1800P = length(P_amb_kPa)-1801;
P_amb_kPa = P_amb_kPa(1:end-extra1800P,:);

Total_time=length(TC_cell); 
time = 1:timestep:Total_time;

% Noise deletion
TC_cell = movmean(TC_cell,10);
TC_front = movmean(TC_front,5);
TC_back = movmean(TC_back,5);
TC_vent = movmean(TC_vent,10);
TC_driver = movmean(TC_driver,10);
TC_pass = movmean(TC_pass,10); % or middle
TC_backpass = movmean(TC_backpass,10); % or equipment

% Shift of TC signals to same bias as Tamb
if(G7_22==1)
    TC_front= TC_front - abs(TC_front(107)-TC_cell(124))+abs(TC_front(10)-TC_back(10));
    TC_back= TC_back - (TC_back(1)-TC_cell(1));
    TC_vent= TC_vent - (TC_vent(60)-TC_cell(60));
    TC_driver= TC_driver - (TC_driver(60)-TC_cell(60));
    TC_backpass= TC_backpass - (TC_backpass(60)-TC_cell(60)); % or equipment
elseif(G7_35==1)
    TC_front= TC_front + abs(TC_front(81)-TC_cell(81))+0.6;
    TC_back= TC_back + abs(TC_back(2)-TC_cell(2));
    TC_vent= TC_vent + abs(TC_vent(2)-TC_cell(2));
    TC_driver= TC_driver + abs(TC_driver(2)-TC_cell(2));
    TC_backpass= TC_backpass - abs(TC_backpass(2)-TC_cell(2))+0.3; % or equipment
elseif(G8_22==1)
    TC_front= TC_front - abs(TC_front(260)-TC_cell(160))+0.3;
    TC_back= TC_back - abs(TC_back(260)-TC_cell(160))+0.3;
    TC_vent= TC_vent - abs(TC_vent(140)-TC_cell(105));
    TC_driver= TC_driver - abs(TC_driver(130)-TC_cell(105));
    TC_backpass= TC_backpass - abs(TC_backpass(182)-TC_cell(105)); 
elseif(G8_35==1)
    TC_front= TC_front - (TC_front(1)-TC_cell(1));%+0.8; %new front
    TC_back= TC_back - (TC_back(1)-TC_cell(1));%+0.3;
    % For G8 35C: (error in frontal TCs measurements)
    TC_vent= TC_vent - (TC_vent(1)-TC_cell(1))+0.5;
    TC_driver= TC_driver - (TC_driver(1)-TC_cell(1))+0.5;
    TC_backpass= TC_backpass - (TC_backpass(1)-TC_cell(1)); % or equipment
end

% Transpose TC amb vector and initialization
TC_cell=transpose(TC_cell)+273.15; % T_amb in K. Both TC_cell or T_cell_OBD
T_ini=TC_cell(1);

deletion_filter_peaks = 1;
i=0;
% Deletion peaks due to filtering at beginning and end:
if(deletion_filter_peaks == 1 && (G7_22==1 || G7_35==1))
    for n=1:Total_time
        if(n < 6)
            TC_front(n)=TC_back(n);
        end
        if(n < 33)
            TC_vent(n)=TC_back(n);
        elseif(n > 1750)
            TC_cell(n)=TC_cell(1750-i);
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

%% Initial conditions
T_sky = T_ini - 6; % [K]
cell_length = 7; %m
cell_width = 5.5; %m
cell_heigh = 4; %m       
V_cell = cell_length * cell_width * cell_heigh ;
%% Radiation Properties 
sigma=0.0000000567037321; % Stefan-Boltzmann constant(sigma)[W/(m^2*K^4)]
% epsilon:emissivity[], rho=reflectivity[], tau=transmissivity[], alpha=absorptivity[]
%Windshield
epsilon_ws=0.9; 
rho_ws=0.246;
tao_ws=0.452;
alfa_ws = 1-(rho_ws + tao_ws);
%Rear Window
epsilon_rw=0.9; 
rho_rw=0.1;
tao_rw=0.311;
alfa_rw = 1-rho_rw-tao_rw;
%Side windows, left and right
epsilon_lsw=0.9; 
epsilon_rsw=0.9;
rho_lsw=0.2;
tao_lsw=0.475;
alfa_lsw = 1-rho_lsw - tao_lsw;
rho_rsw=rho_lsw;
tao_rsw=tao_lsw;
alfa_rsw = alfa_lsw;
%Roof
epsilon_roof=0.9; 
rho_roof=0.74;
tao_roof=0;
alfa_roof = 1-rho_roof-tao_roof;
%Base
rho_base=0.3;
alfa_base=0.7;
%% Conduction Properties 
% A=area[m2]
A_ws=0.63*1.3; 
A_rw=0.29*1; 
A_roof=1.8*1.1;
A_sidewindows=2*1.45*0.29; 
A_doors=2*1.45*0.29;
A_base=6;
% k=thermoconductivity[W/(m*K)]
k_ws=0.8; % 0.8-1.4
k_rw=1.4; % 1.4-1.5
k_sidewindows=1.4;
k_ABS=0.1;%Acrylonitrile butadiene styrene
k_steel=14.9;
k_PU = 0.022; % 0.022-0.028 W/mK polyurethane
k_ceiling=k_ABS*k_steel*k_PU;
k_doors=k_ceiling;
% e=thickness[m]
e_ws=0.006;
e_rw=0.005;
e_sidewindows=0.003;
e_ABS=0.002;
e_steel=0.0005;
e_PU = 0.025; %polyurethane 
e_ceiling=e_ABS+e_steel+e_PU;
e_doors=e_ceiling;

%% Heat transfer Convection coefficients, W/(m2*K)
% exterior:
h_ws=5; 
h_rw=5;
h_base=5;
h_sidewindows=5;
h_doors=h_sidewindows;
h_roof=h_sidewindows;
h_engine=20; % = f( fan speed passing through the engine )
% interior:
h_cabin=5; %10;

% Geometry
cabin_height =  0.5456;
cabin_width = 1.3;
cabin_roof_lenght = 1.8;
cabin_base_lenght = 1.45;
CabinVolume = (cabin_width*cabin_roof_lenght*cabin_height/2)+...
    (cabin_width*cabin_base_lenght*cabin_height/2); 
A_dashboard = cabin_width*0.5; %m2

% Air properies
density_air = 1.18; % [Kg/m3] (at 25C)
cp_air = 1006;% + 10*1500; % [J/ kg*K]

% Seats properties
mass_seats = 30; % 20-30 kg
cp_seats = 2000; % 500-2000 [J/kgK]
% density_PE = 860; % [kg/m3] cp_PE = 2302.7; %polyethylene [J/ kg*K]

% Human and equipment heat % ISO 8996, Alberto Viti Corsi (IDAE)
Qhuman = ones(1,Total_time)*N_Humans*heat_human; % [W]
Qequipment = ones(1,Total_time)*heat_equipment; % [W]

% Global Horizontal Irradiance (Hypothesis: same in all windows)
G_roof_inc=Irr;
G_ws_inc=Irr; 
G_rw_inc=Irr;
G_lsw_inc=Irr;
G_rsw_inc=Irr;

% Solar radiations
%ws
G_ws_r=rho_ws*G_ws_inc;
G_ws_a=alfa_ws*G_ws_inc; %bc
G_ws_t=tao_ws*G_ws_inc;
%rw
G_rw_r=rho_rw*G_rw_inc;
G_rw_a=alfa_rw*G_rw_inc; 
G_rw_t=tao_rw*G_rw_inc;
%lsw
G_lsw_r=rho_lsw*G_lsw_inc;
G_lsw_a=alfa_lsw*G_lsw_inc; 
G_lsw_t=tao_lsw*G_lsw_inc;
%rsw
G_rsw_r=rho_rsw*G_rsw_inc;
G_rsw_a=alfa_rsw*G_rsw_inc; 
G_rsw_t=tao_rsw*G_rsw_inc;
%roof
G_roof_a=alfa_roof*G_roof_inc;

% Boundary conditions vector
Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ ...
    (A_sidewindows*G_lsw_t(t))+ (A_sidewindows*G_rsw_t(t)) );

Tbc=[TC_cell(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);...
    -A_rw*G_rw_a(t);-A_sidewindows*G_lsw_a(t);-A_doors*G_rsw_a(t);0;0;0;0;-Qhuman(t);-Qequipment(t)]; 

% Conductances Matrix 
% Tamb --K1-> Troof --K2-> Tceiling --K3-> Tcabin_back
K1=h_roof*A_roof; %h_roof?
K2=A_roof*k_ceiling/e_ceiling;
K3=h_cabin*A_roof;
% Tbase --K4-> Tcabin_back
K4=h_base*A_base;
% Tamb --K5-> Tw_ext_ws --K6-> Tw_int_ws --K7-> Tcabin_front
K5=h_ws*A_ws;
K6=A_ws*k_ws/e_ws;
K7=h_cabin*A_ws;
% Tamb --K8-> Tw_ext_rw --K9-> Tw_int_rw --K10-> Tcabin_back
K8=h_rw*A_rw;
K9=A_rw*k_rw/e_rw;
K10=h_cabin*A_rw;
% Tamb --K11-> Tw_ext_lsw --K12-> Tw_int_lsw --K13-> Tcabin_back
K11=h_sidewindows*A_sidewindows;
K12=A_sidewindows*k_sidewindows/e_sidewindows;
K13=h_cabin*A_sidewindows;
% Tamb --K14-> Tw_ext_rsw --K15-> Tw_int_rsw --K16-> Tcabin_back
K14=h_doors*A_doors;
K15=A_doors*k_doors/e_doors;
K16=h_cabin*A_doors;
% Tcabin_front --K17-> Tcabin_back
K17=h_cabin*A_ws; %TBD
K=[1 0 0 0 0 0 0 0 0 0 0 0 0 0;
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

% Engine heat
% ini-end of the cycle: delta_temp_engine_mean = (TC_front(end)-TC_front(1)+TC_back(end)-TC_back(1))/2;
% Qengine_W = cp_air*CabinVolume*density_air*delta_temp_engine_mean/Total_time; 
% E_engine_Ws = Qengine_W * 1800; % Ws
% E_engine_Wh = E_engine_Ws / 3600; % Wh

% Capacitance Matrix
C_amb=0; %b.c.
C_roof=0; %for now
C_ceiling=0; %for now
C_base= 144240 + (mass_seats*cp_seats); %Capacitance, W/K  
C_ws=0;
C_rw=0;
C_lsw=0;
C_rsw=0;
C_cabin_front=(density_air*(CabinVolume/2)*cp_air/timestep) ; % + ((mass_seats/2)*cp_seats/timestep);
C_cabin_back =(density_air*(CabinVolume/2)*cp_air/timestep) ; % + ((mass_seats/2)*cp_seats/timestep);
C=zeros(14);
C(1,1)=C_amb;
C(2,2)=C_roof;
C(3,3)=C_ceiling;
C(4,4)=C_base;
C(5,5)=C_ws;
C(6,6)=C_rw;
C(7,7)=C_lsw;
C(8,8)=C_rsw;
C(9,9)=C_ws;
C(10,10)=C_rw;
C(11,11)=C_lsw;
C(12,12)=C_rsw;
C(13,13)=C_cabin_front;
C(14,14)=C_cabin_back;

% Engine power and temperature:
engine_kW=engine_Nm.*engine_rpm*2*pi/(60*1000);
    % OP1: alpha=0.1; %(efficiency from mech power to thermal)/(mass flow warm air from engine compartment to cabin, m_aireng kg/s)
    % T_engine(t) = T_amb(t) + engine_kW*alpha/cp_air;
    % OP2: T_engine(t)=-2*10^(-6)*engine_rpm(t)^(2) + 0.0355*engine_rpm(t) + 77.5; %Fayazbakhsh 2013

% Humidity ratio in gram of water per gram of dry air, X
Ps_kPa = 2.626;
if(G7_22==1 || G8_22==1)
    Ps_kPa = 2.626; % kPa, Water saturation pressure at 22C
elseif(G7_35==1 || G8_35==1)
    Ps_kPa = 5.63; % kPa, Water saturation pressure at 35C
end
meanP_kPa_amb = mean(P_amb_kPa); % in kPa
X = 0.62198*Ps_kPa.*RH_amb./(100.*P_amb_kPa - Ps_kPa.*RH_amb); % Faya, 2013

% Enthalpies calculation. Humidity considered same in amb and cabin
e_amb = 1006.*TC_cell + (2501000 + 1770.*TC_cell).*X; % J/kg, Faya 2013
% e_cabin needs to be calc in the loop

% Air vent volume and flow rate:
vent_volumerate=0.0003; % m3/s, Faya 2013: 0.02; 0.001
leakage_volumerate = 0.0001; %% m3/s, Faya2014

% %% GENERAL EQUATIONS
% % m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 
% T_air(t) = (Q_windows(t)+Q_ceiling(t)+Q_base(t)+Q_human(t))*(timestep/(CabinVolume*cp_air*density_air)) + T_air(t-1);
% Q_ws(t) = h_ws*A_ws*(T_ws(t)- T_air(t-1));
% Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
% Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));
% Q_windows(t)=Q_ws(t)+Q_lsw(t)+Q_rsw(t)+Q_rw(t);

%% INI
Troof(t)=T_ini;
Tceiling(t)=T_ini;
Tbase_int(t)=T_ini;
Tw_ext_ws(t)=T_ini;
Tw_ext_rw(t)=T_ini;
Tw_ext_lsw(t)=T_ini;
Tw_ext_rsw(t)=T_ini;
Tw_int_ws(t)=T_ini;
Tw_int_rw(t)=T_ini;
Tw_int_lsw(t)=T_ini;
Tw_int_rsw(t)=T_ini;
Tcabin_front(t)=T_ini;
Tcabin_back(t)=T_ini;
Tcabin(t)=(Tcabin_front(t)+Tcabin_back(t))/2;
Tcabin=ones(1,Total_time)*(Tcabin_front(t)+Tcabin_back(t))/2; 
temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];
delta_temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];
delta_temp_engine_mean = -ones(1,Total_time)*0.001;
delta_temp_tailpipe_mean = -ones(1,Total_time)*0.001;
check=0; check_f=0; check_b=0; check_t=0;
Qengine = ones(1,Total_time); % [W]
Qexhaust = ones(1,Total_time);
Qleakage = ones(1,Total_time);
Qtotal = ones(1,Total_time);
Qcv_emitted(t)=0;
heat_human = ones(1,Total_time);
T_engine(t)=T_ini;
T_exhaust(t)=T_ini;
e_cabin = ones(1,Total_time)*e_amb(1);

% MAC Thermal Load parameters:
Ttarget = 22 + 273.15; %K
T0cabin= (Tcabin_front(1)+Tcabin_back(1))/2; %K
ttarget = ones(1,Total_time+1)* 900; %s, time to reach Ttarget from T0cabin -> TBD in the loop
tcons= ttarget/log(abs(T0cabin-Ttarget)); % pull-down constant: overall pull-down time
%Formula: Qmac_needed(t)= (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*(Tcabin(t)-Ttarget)/tcons;
Qmac = ones(1,Total_time) * ...
        (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*... 
        ((Tcabin(t)-Ttarget))/tcons(t); 
% % Overall heat transfer coefficients, UA (W/K):
% UA_front = (k_ws*A_ws)/e_ws + 0.95*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;
% UA_back = (k_rw*A_rw)/e_rw + 0.05*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;

A_skin=1.5; %m2
T_skin=24.648+273.15; %K

A_engine = 0.01*A_ws; %Faya 2013: *0.01
A_exhaust = 0.01 * cabin_width*cabin_base_lenght; %Faya 2013
e_enginesheet = 0.001; %engine sheet thickness, m
e_base = 0.1; %m
U_engine = 1/( 1/h_engine + 1/h_cabin + e_enginesheet/k_steel );
U_exhaust = 1/( 1/h_engine + 1/h_cabin + e_base/k_steel );
% U_engine = h_ext;
% U_exhaust = h_cabin;

for t=2:Total_time
    %% Overall heat transfer coefficients, UA (W/K):
    high_speed_index=0;
    if vehicle_kmh(t)<=56.5 % LOW
        f=0.5;
    elseif vehicle_kmh(t)>56.5 && vehicle_kmh(t)<=76.6 % MEDIUM
        f=0.8;
    elseif vehicle_kmh(t)>76.6 && vehicle_kmh(t)<=97.4 % HIGH
        f=0.8;
    elseif vehicle_kmh(t)>97.4 && vehicle_kmh(t)<=131.3 % EXTRA HIGH
        f=0.99;
    elseif vehicle_kmh(t)>131.3
        f=0.95;
        high_speed_index=high_speed_index+1;
    end
    % f = 1; % what % is considered "front"
    A_front = A_ws+f*(A_sidewindows+A_doors+A_roof);
    A_back = A_rw+(1-f)*(A_sidewindows+A_doors+A_roof);
    U_front = 1/h_ws+(k_ws/e_ws);
    U_back = 1/h_rw+(k_rw/e_rw);
    U_side= 1/h_sidewindows+(k_sidewindows/e_sidewindows)...
            + 1/h_doors+(k_doors/e_doors)...
            + 1/h_roof+(k_ceiling/e_ceiling);
    
    UA_front = A_front*U_front;
    UA_back = A_back*U_back;

    % Cabin air enthalpy
    e_cabin(t) = 1006*Tcabin(t) + (2501000 + 1770*Tcabin(t))*X(t); % J/kg, Faya 2013
        
    % Engine and exhaust temperatures
    T_engine(t) = 273.15 + (-2*10^(-6)*engine_rpm(t)^(2) + 0.0355*engine_rpm(t) + 77.5); %K
    T_exhaust(t) = 273.15 + 0.138*engine_rpm(t) - 17; 

    % Heats
    Qleakage(t) = - leakage_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    Qvent(t) = vent_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    heat_human(t) = h_cabin*A_skin*abs(temperature(13)-T_skin) ;
    Qhuman(t) = N_Humans.*heat_human(t);
    Qequipment(t) = heat_equipment;
    Qbase(t) = 0;
    Qengine_calc(t) = A_engine*U_engine*abs(T_engine(t)-temperature(13));% - 11 - 44;
    Qexhaust_calc(t) = A_exhaust*U_exhaust*abs(T_exhaust(t)-temperature(14));%  + 17.77;

    if ( (TC_front(t)+273.15-temperature(13)) < 0 )
        delta_temp_engine_mean(t) = -0.0001;
    else
        delta_temp_engine_mean(t) = (TC_front(t)+273.15-temperature(13)); % K
    end
    if ( (TC_back(t)+273.15-temperature(14)) < 0)
        delta_temp_exhaust_mean(t) = -0.0001;
    else
        delta_temp_exhaust_mean(t) = (TC_back(t)+273.15-temperature(14)); % K
    end
    % Qengine_observed(t)=0;
    % Qexhaust_observed(t)=0;
    % Qengine(t) = Qengine_observed(t) + Qengine_calc(t);
    % Qexhaust(t) = Qexhaust_observed(t) + Qexhaust_calc(t);
    
    % For G7_22:
    Qengine(t) = Qengine_calc(t) ;
    Qexhaust(t) = Qexhaust_calc(t) ;

    %Check: 
    % Qengine(t)=0;% 
    % Qexhaust(t)=0;% 
    % Qhuman(t)=0;
    % Qequipment(t)=0;% 
    % Qvent(t)=0;
    % Qleakage(t)=0;

    % Boundary conditions vector:
    if(Irr==0)
        if(size(TC_cell)==[1 1]) % IF T amb is a constant value:
            Tbc=[TC_cell;0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)];
        else % IF T amb is a vector from input (size(T_amb)==[1 1801]):
            Tbc=[TC_cell(t);0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)];
        end
    else  % IF there is irradiance and T amb is a vector
        Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_sidewindows*G_lsw_t(t)) + (A_sidewindows*G_rsw_t(t)) );
        Tbc=[TC_cell(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);...
            -A_sidewindows*G_lsw_a(t);-A_doors*G_rsw_a(t);0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)]; 
    end

    % Temperatures calculation
    temperature=(inv(K - C))*(Tbc - C*delta_temperature);
    delta_temperature=temperature;

   % Exchange from the front and back:
    if ( temperature(13)>TC_cell(t) )
        temperature(13)= ( density_air*CabinVolume/4*cp_air*delta_temperature(13)/timestep + UA_front*TC_cell(t) )/...
            ( density_air*CabinVolume/4*cp_air/timestep + UA_front );
        delta_temperature(13)=temperature(13);
       % CabinVolume=CabinVolume*(1+vent_volumerate);
        check_f=check_f+1;
        Qcv_emitted(t) = -(temperature(13)-TC_cell(t))*UA_front;
        Qcv_received(t) = 0;
    else
        Qcv_emitted(t) = 0;
        Qcv_received(t) = (TC_cell(t)-temperature(13))*UA_front;
    end

    if ( temperature(14)>TC_cell(t) )
        temperature(14) = (UA_back*TC_cell(t) + density_air*CabinVolume/4*cp_air*delta_temperature(14)/timestep)...
            / (UA_back + density_air*CabinVolume/4*cp_air/(timestep));
        delta_temperature(14)=temperature(14);
        check_b=check_b+1; % v14: V_cell=V_cell*(1+vent_volumerate);
%        Qcv_emitted(t) = Qcv_emitted(t)-(temperature(14)-TC_cell(t))*UA_back;
        Qcv_emitted(t) = -(temperature(14)-TC_cell(t))*UA_back;
        Qcv_received(t) = 0;
    else
        Qcv_emitted(t) = 0;
        %Qcv_received(t) = Qcv_received(t)+(TC_cell(t)-temperature(14))*UA_back;
        Qcv_received(t) = (TC_cell(t)-temperature(14))*UA_back;
    end

    
    % Temperatures extraction
    Tamb(:,t)=temperature(1); %T_amb: boundary condition, Tamb: simulated
    Troof(:,t)=temperature(2);
    Tceiling(:,t)=temperature(3);
    Tbase_int(:,t)=temperature(4);
    Tw_ext_ws(:,t)=temperature(5);
    Tw_ext_rw(:,t)=temperature(6);
    Tw_ext_lsw(:,t)=temperature(7);
    Tw_ext_rsw(:,t)=temperature(8);
    Tw_int_ws(:,t)=temperature(9);
    Tw_int_rw(:,t)=temperature(10);
    Tw_int_lsw(:,t)=temperature(11);
    Tw_int_rsw(:,t)=temperature(12);
    Tcabin_front(:,t)=temperature(13);
    Tcabin_back(:,t)=temperature(14);
    Tcabin(:,t)= (Tcabin_front(:,t)+Tcabin_back(:,t))/2;
    
    % Heats ambient-cabin, to and from
    % if Qcabin_emitted(t)==0
    %     Qcv_front(t)=-h_ws*A_ws*(Tw_int_ws(:,t)-Tcabin_front(:,t));
    %     Qcv_back(t)=-h_rw*A_rw*(Tw_int_rw(:,t)-Tcabin_back(:,t));
    %     Qcv_left(t)=-h_lsw*A_sidewindows*(Tw_int_lsw(:,t)-Tcabin_back(:,t));
    %     Qcv_right(t)=-h_rsw*A_doors*(Tw_int_rsw(:,t)-Tcabin_back(:,t));
    %     Qcv_top(t)=-h_ceiling*A_roof*(Tceiling(:,t)-Tcabin_back(:,t));
    % else
        % Qcv_front(t)=-h_ws*A_ws*(Tcabin_front(:,t)-Tw_int_ws(:,t));
        % Qcv_back(t)=-h_rw*A_rw*(Tcabin_back(:,t)-Tw_int_rw(:,t));
        % Qcv_left(t)=-h_lsw*A_sidewindows*(Tcabin_back(:,t)-Tw_int_lsw(:,t));
        % Qcv_right(t)=-h_rsw*A_doors*(Tcabin_back(:,t)-Tw_int_rsw(:,t));
        % Qcv_top(t)=-h_ceiling*A_roof*(Tcabin_back(:,t)-Tceiling(:,t));
    % end
    % Qcv_cabin(t)=Qcv_front(:,t)+Qcv_back(:,t)+Qcv_left(:,t)+Qcv_right(:,t)+Qcv_top(:,t);
    
    %% Thermal Loads
    Qambient(t) = Qcv_received(t); % Qambient(t) = abs(T_amb(t)-temperature(13))*UA_front + abs(T_amb(t)-temperature(14))*UA_back ;
    Qcabin_received(t) = Qhuman(t) + Qequipment(t) + Qengine(t) + Qexhaust(t) + Qambient(t) + Qvent(t) + Qleakage(t);

    Qtotal(t) = Qcabin_received(t) + Qcv_emitted(t);
    
    % Time to reach the T target, ttarget
        % ttarget(t) = log(abs(T0cabin-Ttarget)) * ...
        %             (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*...
        %             (Tcabin(t)-Ttarget)/(Qmac_needed(t));
    ttarget(t) = 600;
    
    tcons= ttarget/abs(log(Ttarget-T0cabin));
    
    if(heatpump==0)
        if(Tcabin(t)>Ttarget)
            % MAC heat needed to reach the T_target:
            Qtarget(t) = (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*(Ttarget-Tcabin(t))/tcons(t); %negative in cooling mode
            % Heat in the cabin that will need to be provided by MAC:
            Qmac_needed(t) = Qcabin_received(t) + Qcv_emitted(t) + Qtarget(t);
        else
            Qtarget(t) = 0;
            Qmac_needed(t) = 0;
        end
    elseif(heatpump==1)
        Qtarget(t) = (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*(Ttarget-Tcabin(t))/tcons(t); % negative: cooling, positive: heating
        Qmac_needed(t) = Qcabin_received(t) + Qcv_emitted(t) + Qtarget(t);
    end

end
text = ['Vehicle reached higher speed than WLTC limits during: ',num2str(high_speed_index),' seconds.'];
disp(text)

Qcabin_received_plot=Qcabin_received;
Qmac_needed_plot=Qmac_needed;
Qemitted=Qleakage+Qcv_emitted;

%% Energies
% Delete NaN values
Qcabin_received=(Qcabin_received(~isnan(Qcabin_received)));
Qmac_needed=(Qmac_needed(~isnan(Qmac_needed)));

% Cabin thermal load
energy_cabin_Wh = timestep.*Qcabin_received(2:end)./3600;
energy_cum_cabin_Wh = cumsum(energy_cabin_Wh);
total_energy_cum_cabin_Wh=energy_cum_cabin_Wh(end);
total_cabin_load_W=total_energy_cum_cabin_Wh * 3600 / duration_WLTC

% MAC needed cooling load
energy_mac_needed_Wh = timestep.*Qmac_needed(2:end)./3600;
energy_cum_mac_needed_Wh = cumsum(energy_mac_needed_Wh);
total_energy_cum_mac_needed_Wh=energy_cum_mac_needed_Wh(end);
total_mac_cooling_load_W=total_energy_cum_mac_needed_Wh * 3600 / duration_WLTC

% AC compressor power demand
compressor_electric_power_W = -total_mac_cooling_load_W/COP

%% Save engine load to be applied in other scenarios:
%save Engine_thermal_loads_W_35G7_v9 Qengine Qexhaust

%% PLOTS
% Deletion of first 0 value for plotting (inicialization value):
x=TC_cell(2); TC_cell(1)=x;
% The boundary conditon 'T_amb' needs to be a vector to be plotted:
if (size(TC_cell)==[1 1])
    TC_cell = ones(1,Total_time)*(T_ambient+273.15); %[K], vector
end
%%
figure(12)
plot(time(3:Total_time),Qtarget(3:Total_time),':k','LineWidth',1.5)
hold on
plot(time(3:Total_time),Qmac_needed(3:Total_time),'-k','LineWidth',1.5)
grid minor
ylabel('MAC Heat [W]')
xlabel('Time [s]')
if(G7_22==1)
    ylim([-1200 200])
elseif(G8_22==1)
    ylim([-400 200])
elseif(G7_35==1 || G8_35==1)
    ylim([-20500 -16500]) % G7_35: ylim([-21000 -18500]) % G8_35: ylim([-20000 -13000])
end
yyaxis right
ax = gca;
ax.YColor = [1 0 0];
hold on
plot(time(3:Total_time),Qtotal(3:Total_time),'-r','LineWidth',1.5)
ylabel('Total Cabin Heat [W]','Color',[1 0 0])
% if(G8_22==1)
%     ylim([-400 200])
% elseif(G7_22==1 || G7_35==1 || G8_35==1)
if(G8_22==1 || G7_22==1 || G7_35==1 || G8_35==1)
    ylim([-5 130])
end
legend('Target','MAC','Cabin Total','Location','best')

figure(11)
plot(time(3:Total_time),Qtotal(3:Total_time),'LineWidth',1.5)
hold on
plot(time(3:Total_time),Qcabin_received_plot(3:Total_time),'LineWidth',1.5)
plot(time(3:Total_time),Qemitted(3:Total_time),'LineWidth',1.5)
ylabel('Cabin Heat Loads [W]')
xlabel('Time [s]')
legend('Cabin Total','Cabin Received','Cabin Emitted','Location','best')
grid minor
%%
% figure(10)
% plot(time,Qtotal,'--k','LineWidth',1.5)
% hold on
% plot(time,Qambient,':','LineWidth',1.5)
% plot(time,Qhuman,':','LineWidth',1.5)
% plot(time,Qequipment,':','LineWidth',1.5)
% plot(time,Qvent,':','LineWidth',1.5)
% plot(time,Qengine,':','LineWidth',1.5,'Color','g')
% plot(time,Qexhaust,':','LineWidth',1.5)
% plot(time,Qcabin_received_plot,'-','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
% ylabel('Thermal Loads [W]')
% xlabel('Time [s]')
% yyaxis right
% ax = gca;
% ax.YColor = [0.4940 0.1840 0.5560];
% %ylim([-1000 600])
% hold on
% plot(time,Qmac_needed_plot,'-','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560])
% ylabel('Q cooling MAC [W]','Color',[0.4940 0.1840 0.5560])
% %set(gca,'Ydir','reverse')
% legend('Ambient','Human','Equipment','Vents','Engine',...
%     'Exhaust','Cabin Received','MAC','Location','northeastoutside')
% grid minor
    % xlim([0 1800])
    % axis([1 Total_time -1500 1700])
    % axis([1 Total_time 0 1100])
%%
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
%%
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
%%
% figure(7)
% plot(time,Qcabin_received_plot) 
% xlabel('Time [s]') 
% hold on
% plot(time,-Qcv_emitted)
% ylabel('Cabin Heat [W]')
% legend('Received','Emitted')
%%
% figure(6)
% plot(time,Qcabin_received_plot-abs(Qcv_emitted)) 
% xlabel('Time [s]') 
% hold on
% % plot(time,-Qcabin_emitted)
% ylabel('Cabin Heat [W]')
% legend('Received','Emitted')
%% 
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
%% 
% figure(4)
% xlabel('Time [s]') 
% plot(time,Qmac_needed_plot)
% ylabel('Heat needed to get T target [W]')
%%
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
%%
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
%%
figure(2)
plot(time,TC_cell-273.15,'-k','LineWidth',1) %T amb 
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
if(G7_22==1 || G8_22==1)
    axis([1 Total_time 20 25])
elseif(G7_35==1 || G8_35==1)
    axis([1 Total_time 32 38])
end
%%
figure(1)
plot(time,Tcabin_front-273.15,':b','LineWidth',1.5)
hold on
plot(time,Tcabin_back-273.15,':m','LineWidth',1.5)
plot(time,TC_cell-273.15,'-k','LineWidth',1) %T amb 
plot(time,TC_front,'-b','LineWidth',1.5) 
plot(time,TC_back,'-m','LineWidth',1.5)
grid minor
xlabel('Time [s]') 
ylabel('Temperature [ºC]')
if WT==1
    plot(time,mean_TC_front,'--k','LineWidth',1) 
    plot(time,mean_TC_back,'--k','LineWidth',1)
    xline(x3)
    xline(x4)
    legend('T front simulated','T back simulated','TC ambient','TC front','TC back','TC front mean','TC back mean','20 to 25 m/s','25 to 30 m/s','Location','best')
else
    xline(589)
    xline(589+433)
    xline(589+433+455)
    xline(589+433+455+323)
    legend('T front simulated','T back simulated','TC ambient','TC front','TC back','20 to 25 m/s','25 to 30 m/s','Location','best')
end
hold off
if(G7_22==1 || G8_22==1)
    axis([1 1800 20 25])
elseif(G7_35==1 || G8_35==1)
    axis([1 1800 33 38])
end
%%
% figure(13)
% plot(time,vehicle_kmh,'-b','LineWidth',1.5)
% xline(589)
% xline(589+433)
% xline(589+433+455)
% xline(589+433+455+323)
% grid minor
% xlabel('Time [s]') 
% ylabel('Vehicle Speed [km/h]')