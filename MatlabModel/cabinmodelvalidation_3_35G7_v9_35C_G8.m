% Vehicle's Cabin Model, Susana Gil-Sayas
clc
clear all 
close all

%% INPUTS:

% Test conditions
Irr = -0; % IN NEGATIVE as it enters the cabin [W/(m2)] 
N_Humans = 1; % old: 1 human = 120 W, 0.66 = 80 W, 0.25 = 30 W;
heat_human = 80; %W
heat_equipment = 40 ; %40; % W
T_ambient = 22; % [C] % NOT USED WHEN IMPORTING DATA

% Initial conditions
T_ini = 273.15 + T_ambient; % [K], initial value for simulation
T_sky = T_ini - 6; % [K]
cell_length = 7; %m
cell_width = 5.5; %m
cell_heigh = 4; %m       
V_cell = cell_length * cell_width * cell_heigh ;

% Radiation Properties 
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

% Conduction Properties 
% A=area[m2], k=thermoconductivity[W/(m*K)], e=thickness[m]
%Windshield
A_ws=0.63*1.3; 
k_ws=0.8; % 0.8-1.4
e_ws=0.006;
%Rear Window
A_rw=0.29*1; 
k_rw=1.4; % 1.4-1.5
e_rw=0.005;
%Roof and ceiling ONLY STEEL.
A_roof=1.8*1.1;
k_ABS=0.1;%Acrylonitrile butadiene styrene
k_steel=14.9;
k_PU = 0.022; % 0.022-0.028 W/mK polyurethane
k_ceiling=k_ABS*k_steel*k_PU;
e_ABS=0.002;
e_steel=0.0005;
e_PU = 0.025; %polyurethane 
e_ceiling=e_ABS+e_steel+e_PU;
%Side Windows and doors
A_sidewindows=2*1.45*0.29; 
A_doors=2*1.45*0.29;
k_sidewindows=1.4;
k_doors=k_ceiling;
e_sidewindows=0.003;
e_doors=e_ceiling;
%Base
A_base=6;

A_front = A_ws + A_sidewindows/2 + A_doors/2 + A_roof/2;
e_av_front = (e_ws + e_rw + e_sidewindows + e_doors + e_ceiling)/5 ;
k_av_front = (k_ws + k_rw + k_sidewindows + k_doors + k_ceiling)/5;

A_back = A_rw + A_sidewindows/2 + A_doors/2 + A_roof/2;
e_av_back = (e_rw + e_sidewindows + e_doors + e_ceiling)/4 ;
k_av_back = (k_rw + k_sidewindows + k_doors + k_ceiling)/4;

% Heat transfer Convection coefficients, W/(m2*K)
h_ws=5; 
h_lsw=5;
h_rsw=5;
h_rw=5;
h_ceiling=5;
h_roof=20;
h_ext=20; 
h_base = 5;
h_cabin=5; %10;

% Geometry
cabin_height =  0.5456;
cabin_width = 1.1;
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

duration_WLTC=1800; %s

% Load input for validation: 
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

n=height(lab_array)/duration_WLTC;; %10;%8.5;
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

TC_front=TC_backpass;

% Load calculated engine and exhaust thermal load of test G7 at 35C
engineheatinput='Engine_thermal_loads_W_35G7_v9';
load(engineheatinput);
%struct=load(engineheatinput);
% table=struct2table(struct);
% arr=table2array(table);
% Qengine_observed=arr(1,:);
% Qexhaust_observed=arr(2,:);

t = 1; timestep =  1; %s
Total_time=length(TC_cell); % 1800; % [s], WLTC duration
%i = Total_time-1800; % DEFINE i BEFORE LAB INPUT
time = 1:timestep:Total_time;
if (Irr==0)
    T_amb = T_ambient + 273.15; % [K]
else
    T_amb = ones(1,Total_time+1)*(T_ambient+273.15); %[K], vector
    Irr = ones(1,Total_time+1)*Irr;
end

%%
lab_input2='G7_35C_MACoff_XCU'; 
struct2=load(lab_input2);
table2=struct2table(struct2);
modalxcu_array=table2array(table2);
VehicleSpeed = table2array(modalxcu_array(:,1));
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
vehicle_kmh = VehicleSpeed(1 : n : end);
ACcomp_rpm = ACcomp_speed_rpm(1 : n : end);
ACcomp_Nm = ACcomp_torque_Nm(1 : n : end);
engine_Nm = ICETorque_Nm(1 : n : end);
engine_rpm = engineSpeed_rpm(1 : n : end);
FC_lperh = fuelconsumption_lph(1 : n : end);
Total_time2=length(vehicle_kmh)-1; 
time2 = 0:timestep:Total_time2;
%%
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
%%

T_amb=transpose(TC_cell)+273.15; % T_amb in K. Both TC_cell or T_cell_OBD
T_ini=T_amb(1);

% Shift of TC signals to same bias of Tamb and noise deletion
T_amb = movmean(T_amb,10);
TC_front = movmean(TC_front,5);
TC_back = movmean(TC_back,5);
TC_front= TC_front - abs(TC_front(80)-T_amb(80)+273.15)+0.5;
TC_back= TC_back - abs(TC_back(80)-T_amb(80)+273.15)+0.5;

% %TC_front= TC_front - abs(T_amb(2)-273.15-TC_front(6)); %%%%%% 22C G7 case %%%%%%%%

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

Tbc=[T_amb(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);...
    -A_rw*G_rw_a(t);-A_sidewindows*G_lsw_a(t);-A_doors*G_rsw_a(t);0;0;0;0;-Qhuman(t);-Qequipment(t)]; 

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
% Tamb --K11-> Tw_ext_lsw --K12-> Tw_int_lsw --K13-> Tcabin_back
K11=h_ext*A_sidewindows;
K12=A_sidewindows*k_sidewindows/e_sidewindows;
K13=h_lsw*A_sidewindows;
% Tamb --K14-> Tw_ext_rsw --K15-> Tw_int_rsw --K16-> Tcabin_back
K14=h_ext*A_doors;
K15=A_doors*k_doors/e_doors;
K16=h_rsw*A_doors;
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
C_base= 144240 + (mass_seats*cp_seats/timestep); %Capacitance, W/K  
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

meanP_kPa_amb = mean(P_amb_kPa); % in kPa
Ps35_kPa = 5.63; % kPa, Water saturation pressure at 35C
Ps22_kPa = 2.626; % kPa, Water saturation pressure at 22C

% Humidity ratio in gram of water per gram of dry air, X
X = 0.62198*Ps35_kPa.*RH_amb./(100.*P_amb_kPa - Ps35_kPa.*RH_amb); % Faya, 2013
%X = 0.62198*Ps22_kPa.*RH_amb./(100.*P_amb_kPa - Ps22_kPa.*RH_amb);

% Enthalpies calculation. Humidity considered same in amb and cabin
e_amb = 1006.*T_amb + (2501000 + 1770.*T_amb).*X; % J/kg, Faya 2013
% e_cabin needs to be calc in the loop

% Air vent volume and flow rate:
vent_volumerate=0.0002; % m3/s, Faya 2013: 0.02; 0.001
vent_massflowrate=vent_volumerate*density_air ; % m3/s * kg/m3 = kg/s

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

% %ctes:
% Ttarget = 22 + 273.15; %K
% T0cabin= (Tcabin_front(1)+Tcabin_back(1))/2; %K
% ttarget = 900; %s, time to reach the Ttarget from T0 cabin
% tcons= ttarget/log(abs(T0cabin-Ttarget)); % pull-down constant: overall pull-down time
% % Qmac_needed(t)= (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*...
% %     (Tcabin(t)-Ttarget)/tcons;


% Overall heat transfer coefficients, UA (W/K):
UA_front= (k_ws*A_ws)/e_ws + 0.5*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;
UA_back = (k_rw*A_rw)/e_rw + 0*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;
A_skin=1.5; %m2
T_skin=24.648+273.15; %K

for t=2:Total_time
    % Cabin air enthalpy:
    e_cabin(t) = 1006*Tcabin(t) + (2501000 + 1770*Tcabin(t))*X(t); % J/kg, Faya 2013
    % Heats:
        % Qhuman(t)=0;% Qequipment(t)=0;% Qengine(t)=0;% Qexhaust(t)=0;% Qvent(t)=0;
    Qleakage(t)=-0.02*density_air*(e_amb(t)-e_cabin(t)); %vent_massflowrate % W, J/s
    Qvent(t)=-1*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    heat_human(t) = h_cabin*A_skin*abs(temperature(13)-T_skin) ;
    Qhuman(t) = N_Humans.*heat_human(t);
    Qequipment(t) = heat_equipment;
    Qbase(t)=0;

    A_engine = 0.05*A_ws; %Faya 2013: *0.01
    U_engine = h_ext;
    A_exhaust = 0.01 * cabin_width*cabin_base_lenght; %Faya 2013
    U_exhaust = h_cabin;
    T_engine(t) = (273.15 + (-2*10^(-6)*engine_rpm(t)^(2) + 0.0355*engine_rpm(t) + 77.5)); %K
    T_exhaust(t) = 0.138*engine_rpm(t) - 17;
    Qengine_calc(t)= A_engine*U_engine*(T_engine(t)-temperature(13)) - 11 - 44;
    Qexhaust_calc(t)= A_exhaust*U_exhaust*(T_exhaust(t)-temperature(14)) + 17.77;
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
    % Qengine_observed(t) = cp_air*CabinVolume/2*density_air*(delta_temp_engine_mean(t))/timestep;
    % Qexhaust_observed(t) = cp_air*CabinVolume/2*density_air*(delta_temp_exhaust_mean(t))/timestep;
    
    Qengine(t)=1*Qengine_observed(t) + Qengine_calc(t) + 0*(cp_air*CabinVolume/2*density_air*(delta_temp_engine_mean(t))/timestep);
    Qexhaust(t)=1*Qexhaust_observed(t) + Qexhaust_calc(t) + 0*(cp_air*CabinVolume/2*density_air*(delta_temp_exhaust_mean(t))/timestep);

    % Boundary conditions vector:
    if(Irr==0)
        if(size(T_amb)==[1 1]) % IF T amb is a constant value:
            Tbc=[T_amb;0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)];
        else % IF T amb is a vector from input (size(T_amb)==[1 1801]):
            Tbc=[T_amb(t);0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)];
        end
    else  % IF there is irradiance and T amb is a vector
        Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_sidewindows*G_lsw_t(t)) + (A_sidewindows*G_rsw_t(t)) );
        Tbc=[T_amb(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);...
            -A_sidewindows*G_lsw_a(t);-A_doors*G_rsw_a(t);0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)]; 
    end

    % Temperatures calculation:
    temperature=(inv(K - C))*(Tbc - C*delta_temperature);
    delta_temperature=temperature;

   % Exchange from the front and back:
    %Formula:(temperature(13)-delta_temperature(13))*density_air*V_cell/timestep = (T_amb(t)*0.998-temperature(13))*UA_front;
    
    %T_amb(t)=mean(T_amb);

    if ( temperature(13)> T_amb(t) ) % si T_amb(t)*0.998 all ok
        temperature(13)= (UA_front*T_amb(t)*1 + density_air*V_cell*delta_temperature(13)/timestep)...
            / (UA_front + density_air*V_cell/(timestep));
        delta_temperature(13)=temperature(13);
        V_cell=V_cell*(1+vent_volumerate);
        check_f=check_f+1;
        % if (Qcv_emitted(t-1)==0)
        %     Qcv_emitted(t)=0;
        % else
            Qcv_emitted(t)=-(temperature(13)-T_amb(t))*UA_front;
        % end
        Qcv_received(t)=0;
    else
        Qcv_emitted(t)=0;
        Qcv_received(t)=(T_amb(t)-temperature(13))*UA_front;
    end

    if ( temperature(14)>T_amb(t) )
        temperature(14) = (UA_back*T_amb(t) + density_air*V_cell*delta_temperature(14)/timestep)...
            / (UA_back + density_air*V_cell/(timestep));
        delta_temperature(14)=temperature(14);
        check_b=check_b+1;
        V_cell=V_cell*(1+vent_volumerate);
        % if (Qcv_emitted(t-1)==0)
        %     Qcv_emitted(t)=0;
        % else
            Qcv_emitted(t)=Qcv_emitted(t)-(temperature(14)-T_amb(t))*UA_back;
        % end
        Qcv_received(t)=0;
    else
        Qcv_emitted(t)=0;
        Qcv_received(t)=Qcv_received(t)+(T_amb(t)-temperature(14))*UA_back;
    end

    
    % Temperatures extraction:
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
    
    
    % Thermal Loads
    Qambient(t)=-(T_amb(t)-temperature(13))*UA_front+(T_amb(t)-temperature(14))*UA_back;
    Qcabin_received(t)=(Qhuman(t) + Qequipment(t) + Qengine(t) + Qexhaust(t) + Qcv_received(t) + Qvent(t) + Qleakage(t));

    % Time to reach the T target, ttarget:
    Qmac_needed(t) = -(abs(Qcabin_received(t))-abs(Qcv_emitted(t)));
    ttarget(t) = log(abs(T0cabin-Ttarget)) * (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*...
        (Tcabin(t)-Ttarget)/abs(Qmac_needed(t));
    tcons= ttarget/log(abs(T0cabin-Ttarget));
    Qmac(t)= -(CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*... 
    (((temperature(13)+temperature(14))/2)-Ttarget)/tcons(t);
    % Qmac_needed(t) = Qmac_needed(t) + Qmac(t);
    % Qmac_needed(t) = Qcabin_received(t) + Qcv_emitted(t) + Qmac(t);
    Qmac_needed(t) = Qmac(t);

    % Energies
    Eengine_kWh(t) = trapz(Qengine(:,1:t))*t / (3600*1000);
    Eexhaust_kWh(t)= trapz(Qexhaust(:,1:t))*t / (3600*1000);
    Ecabin_kWh(t)= trapz(Qcabin_received(:,1:t))*t / (3600*1000);
    Emac_needed_kWh(t)= trapz(Qmac_needed(:,1:t))*t / (3600*1000);
    %cumW=trapz(Qengine_W(:,1:t);

end
Qcabin_received_plot=Qcabin_received;
Qmac_needed_plot=Qmac_needed;

% Delete NaN values
Qcabin_received=(Qcabin_received(~isnan(Qcabin_received)));
Qmac_needed=(Qmac_needed(~isnan(Qmac_needed)));

% Cabin thermal load
energy_cabin_Wh = timestep.*Qcabin_received(2:end)./3600;
energy_cum_cabin_Wh = cumsum(energy_cabin_Wh);
energy_cum_cabin_Wh(end)
% MAC needed cooling load
energy_mac_needed_Wh = timestep.*Qmac_needed(2:end)./3600;
energy_cum_mac_needed_Wh = cumsum(energy_mac_needed_Wh);
energy_cum_mac_needed_Wh(end)
% MAC cooling load
energy_mac_Wh = timestep.*Qmac(2:end)./3600;
energy_cum_mac_Wh = cumsum(energy_mac_Wh);
energy_cum_mac_Wh(end)


%% Save engine load to be applied in other scenarios:
%save Engine_thermal_loads_W_35G7_v9 Qengine Qexhaust

%% PLOTS
% Deletion of first 0 value for plotting (inicialization value):
x=T_amb(2); T_amb(1)=x;
% The boundary conditon 'T_amb' needs to be a vector to be plotted:
if (Irr==0) & (size(T_amb)==[1 1])
    T_amb = ones(1,Total_time)*(T_ambient+273.15); %[K], vector
end
%%
figure(11)
plot(time,Qcabin_received_plot,':k','LineWidth',1.5)
hold on
plot(time,-Qcv_emitted,':','LineWidth',1.5)
hold on
plot(time,-Qmac_needed_plot,'-','LineWidth',1.5)
plot(time,-Qmac,'-','LineWidth',1.5)
legend('Received cabin','Emitted','MAC needed','MAC','Location','northeastoutside')
grid minor
xlim([0 1800])

figure(10)
plot(time,Qcabin_received_plot,':k','LineWidth',1.5)
hold on
plot(time,Qambient,'-','LineWidth',1.5)
plot(time,Qhuman,':','LineWidth',1.5)
plot(time,Qequipment,':','LineWidth',1.5)
plot(time,Qengine,':','LineWidth',1.5)
plot(time,Qexhaust,':','LineWidth',1.5)
plot(time,Qvent,'-','LineWidth',1.5)
plot(time,Qleakage,':','LineWidth',1.5)
ylabel('Thermal Loads [W]')
%ylabel('Thermal Loads to the cabin [W]')
xlabel('Time [s]')
%yyaxis right
plot(time,Qcv_emitted,':','LineWidth',1.5)
hold on
plot(time,Qmac_needed_plot,'-','LineWidth',1.5)
%ylabel('Thermal Loads from the cabin [W]')
%set(gca,'Ydir','reverse')
legend('Total cabin','Ambient','Human','Equipment','Engine','Exhaust','Vents','Leakage','Emitted','MAC needed','Location','northeastoutside')
grid minor
xlim([0 1800])

%axis([1 Total_time 0 1100])
    
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
% figure(2)
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
figure(1)
plot(time,Tcabin_front-273.15,':b','LineWidth',1.5)
hold on
plot(time,Tcabin_back-273.15,':m','LineWidth',1.5)
plot(time,T_amb-273.15,'-k','LineWidth',1) %T amb 
plot(time,TC_front,'-b','LineWidth',1.5) 
plot(time,TC_back,'-m','LineWidth',1.5)
hold off
legend('T front simulated','T back simulated','TC ambient','TC front','TC back','Location','best')
grid minor
xlabel('Time [s]') 
ylabel('Temperature [ºC]')
%axis([1 Total_time 20 27])
axis([1 Total_time 34 38])
% hold on
% yyaxis right
% plot(engine_kW)
%axis([1 Total_time 0 3200])

