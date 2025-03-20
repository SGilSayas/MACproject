% Vehicle's Cabin Model, Susana Gil-Sayas. Version: next of 15.2 new Qtarget, but
% not heat transfer coef validation added yet
% FUNCTION VERSION
clc
clear all 
close all

function [Qcompressor_W_av,total_energy_cum_compr_needed_Wh,total_active_time] = cabinmodel(Total_time,distance_driven,N_Humans,COP,TC_cell_constant,Irr_constant,TC_cell,Irr)
    % Initial data
duration_WLTC = 1800; % s
timestep = 1; % s
t = 1;
heat_human = 80; %W
heat_equipment = 0 ; %40; % W
Vpe_petrol = 0.264; % l/kWh
CF_petrol = 2330; % gCO2/l
    % Inputs stored in CELL ARRAY:
    % Total_time = duration_WLTC; % s
    % distance_driven = 23.25; % km
    % N_Humans = 1;
    % COP = 3.5;
    % TC_cell_constant = 35; % C
    % Irr_constant = -0; % IN NEGATIVE (enters the cabin), W/m2
    % inputs = {Total_time,distance_driven,N_Humans,COP,TC_cell_constant,Irr_constant,TC_cell,Irr};
if Total_time == 'N/A' Total_time = 1800; end
if distance_driven == 'N/A' distance_driven = 23.25; end
if N_Humans == 'N/A' N_Humans = 1; end
if COP == 'N/A' COP = 3.5; end
if TC_cell_constant ~= 'N/A' TC_cell = ones(1,Total_time)*TC_cell_constant; end
if Irr_constant ~= 'N/A' Irr = ones(1,Total_time)*Irr_constant; end
if TC_cell == 'N/A' TC_cell = ones(1,Total_time)*TC_cell_constant; end
if Irr == 'N/A' Irr = ones(1,Total_time)*Irr_constant; end

    % Load input Engine speed and torque (and extra XCU data)
lab_input2 = 'G7_35C_MACoff_XCU'; 
struct2 = load(lab_input2);
table2 = struct2table(struct2);
modalxcu_array = table2array(table2);
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
Total_time2 = length(vehicle_kmh)-1; 
time2 = 0:timestep:Total_time2;

    % Load input cell humidity (RH) and pressure (P)
lab_input3 = 'G7_35C_MACoff_Modaldata'; 
struct3NaN = load(lab_input3);
struct3NoNaN = structfun( @rmmissing , struct3NaN , 'UniformOutput' , false);
table3 = struct2table(struct3NoNaN);
modaldata_array = table2array(table3);
RH_amb = table2array(modaldata_array(:,1));
P_amb_kPa = table2array(modaldata_array(:,2));
    % Reshape RH and P data
extra1800RH = length(RH_amb)-Total_time;
extra1800P = length(P_amb_kPa)-Total_time;
RH_amb = RH_amb(1:end-extra1800RH,:);
P_amb_kPa = P_amb_kPa(1:end-extra1800P,:);

Total_time = length(TC_cell); 
time = 1:timestep:Total_time;

T_cell = transpose(TC_cell)+273.15; % K, TC_cell in C
T_ini = T_cell(1); % K
T_sky = T_ini - 6; % K

    % Radiation Properties 
sigma = 0.0000000567037321; % [W/(m^2*K^4)], Stefan-Boltzmann constant
    % epsilon: emissivity[]
epsilon_ws=0.9;
epsilon_rw=0.9; 
epsilon_sidewindows=0.9; 
epsilon_doors=0.9;
epsilon_roof=0.9;
    % rho: reflectivity[]
rho_ws = 0.246;
rho_rw=0.1;
rho_sidewindows=0.2;
rho_doors=0.74;
rho_roof=0.74;
rho_base=0.3;
    % tau: transmissivity[]
tao_ws=0.452;
tao_rw=0.311;
tao_sidewindows=0.475;
tao_doors=0;
tao_roof=0;
    % alpha: absorptivity[]
alfa_ws = 1-(rho_ws + tao_ws);
alfa_rw = 1-rho_rw-tao_rw;
alfa_sidewindows = 1-rho_sidewindows-tao_sidewindows;
alfa_doors = 1-rho_doors-tao_doors;
alfa_roof = 1-rho_roof-tao_roof;
alfa_base=0.7;

%% Geometry (m) and areas (m2)
cabin_height =  0.5456;
cabin_width = 1.3;
cabin_roof_lenght = 1.8;
cabin_base_lenght = 1.45;
CabinVolume = (cabin_width*cabin_roof_lenght*cabin_height/2)+...
    (cabin_width*cabin_base_lenght*cabin_height/2); 
A_dashboard = cabin_width*0.5; %m2
A_ws=0.63*1.3;
A_rw=0.29*1;
A_roof=1.8*1.1;
A_sidewindows=2*1.45*0.29; 
A_doors=2*1.45*0.29;
A_base=6;
A_front = A_ws + A_sidewindows/2 + A_doors/2 + A_roof/2;
A_back = A_rw + A_sidewindows/2 + A_doors/2 + A_roof/2;
    % Seats properties
mass_seats = 30; % 20-30 kg
cp_seats = 2000; % 500-2000 [J/kgK]
%%
    % Conduction Properties. k=thermoconductivity[W/(m*K)], e=thickness[m]
k_ws=0.8; % 0.8-1.4
e_ws=0.006;
k_rw=1.4; % 1.4-1.5
e_rw=0.005;
k_ABS=0.1; %Acrylonitrile butadiene styrene
k_steel=14.9;
k_PU = 0.022; % 0.022-0.028 W/mK polyurethane
k_ceiling=k_ABS*k_steel*k_PU;
e_ABS=0.002;
e_steel=0.0005;
e_PU = 0.025; %polyurethane 
e_ceiling=e_ABS+e_steel+e_PU;
k_sidewindows=1.4;
k_doors=k_ceiling;
e_sidewindows=0.003;
e_doors=e_ceiling;
e_av_front = (e_ws + e_rw + e_sidewindows + e_doors + e_ceiling)/5 ;
k_av_front = (k_ws + k_rw + k_sidewindows + k_doors + k_ceiling)/5;
e_av_back = (e_rw + e_sidewindows + e_doors + e_ceiling)/4 ;
k_av_back = (k_rw + k_sidewindows + k_doors + k_ceiling)/4;

    % Initial Heat transfer Convection coefficients, W/(m2*K)
h_ws=5; 
h_lsw=5;
h_rsw=5;
h_rw=5;
h_ceiling=5;
h_roof=20;
h_ext=20; 
h_base = 5;
h_cabin=5; %10;

    % Air properies
density_air = 1.18; % [Kg/m3] (at 25C)
cp_air = 1006;% + 10*1500; % [J/ kg*K]

    % Human and equipment heat % ISO 8996, Alberto Viti Corsi (IDAE)
Qhuman = ones(1,Total_time)*N_Humans*heat_human; % [W]
Qequipment = ones(1,Total_time)*heat_equipment; % [W]

    % Global Horizontal Irradiance (Hypothesis: same in all windows)
G_roof_inc=Irr;
G_ws_inc=Irr; 
G_rw_inc=Irr;
G_sidewindows_inc=Irr;
G_doors_inc=Irr;
    % Extraction of all the Solar Irradiance fluxes
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

    % Ratiation emitted by a body:
%Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4);
% Ter1 = h_ext*A_ws*(((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))-T_ext(t))+epsilon_ws*sigma*A_ws*(((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))^4-(T_sky)^4);
% Ter2 = (A_ws*k_ws/e_ws)*(T_ws(t) - ((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))) +  A_ws*G_ws_a(t);

    % Boundary conditions vector
Qbase(t) = alfa_base*(A_ws*G_ws_t(t)+A_rw*G_rw_t(t)+A_sidewindows*G_sidewindows_t(t)+A_doors*G_doors_t(t));

Tbc=[T_cell(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);-A_sidewindows*G_sidewindows_a(t);-A_doors*G_doors_a(t);...
    0;0;0;0;-Qhuman(t);-Qequipment(t)]; 

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

    % Capacitance Matrix
C_amb=0; %b.c.
C_roof=0; %for now
C_ceiling=0; %for now
C_base= 144240 + (mass_seats*cp_seats); %Capacitance, W/K  
C_ws=0;
C_rw=0;
C_lsw=0;
C_rsw=0;
C_cabin_front=(density_air*(CabinVolume/2)*cp_air/timestep);
C_cabin_back =(density_air*(CabinVolume/2)*cp_air/timestep);
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
Ps_22_kPa = 2.626; % kPa, Water saturation pressure at 22C
Ps_35_kPa = 5.63; % kPa, Water saturation pressure at 35C
Ps_kPa = Ps_22_kPa+(T_ini-22)/(35-22)*(Ps_35_kPa-Ps_22_kPa);
meanP_kPa_amb = mean(P_amb_kPa); % in kPa
X = 0.62198*Ps_kPa.*RH_amb./(100.*P_amb_kPa - Ps_kPa.*RH_amb); % Faya, 2013

    % Enthalpies calculation. Humidity considered same in amb and cabin
e_amb = 1006.*T_cell + (2501000 + 1770.*T_cell).*X; % J/kg, Faya 2013
% e_cabin at loop

    % Air vent volume and flow rate:
vent_volumerate = 0.0003; % m3/s, Faya 2013: 0.02; 0.001
leakage_volumerate = 0.0001; %% m3/s, Faya2014

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
Qrd = ones(1,Total_time);
Qcompressor = ones(1,Total_time).*0;
CO2 = ones(1,Total_time).*0;
Qcv_emitted(t)=0;
t_active(t)=0; %=ones(1,Total_time).*0;
% compressor_time = ones(1,Total_time).*0;
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

%Formula: Qmac_needed(t)= (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*(Tcabin(t)-Ttarget)/tcons;
Qmac = ones(1,Total_time) * ...
        (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*... 
        ((Tcabin(t)-Ttarget))/tcons(t); 

% Overall heat transfer coefficients, UA (W/K):
UA_front= (k_ws*A_ws)/e_ws + 0.95*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;
UA_back = (k_rw*A_rw)/e_rw + 0.05*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;

A_skin=1.5; %m2
T_skin=24.648+273.15; %K

% Engine and exhaust heat transfer parameters
h_engine = h_ext; % h_engine = f( velocidad aire que pasa por el motor )
A_engine = 0.01*A_ws; %Faya 2013: *0.01
A_exhaust = 0.01 * cabin_width*cabin_base_lenght; %Faya 2013
e_enginesheet = 0.001; %engine sheet thickness, m
e_base = 0.1; %m
U_engine = 1/( 1/h_engine + 1/h_cabin + e_enginesheet/k_steel );
U_exhaust = 1/( 1/h_engine + 1/h_cabin + e_base/k_steel );

compressor_on_lower_limit = 0.5; % W, 

for t=2:Total_time
    
    % Cabin air enthalpy
    e_cabin(t) = 1006*Tcabin(t) + (2501000 + 1770*Tcabin(t))*X(t); % J/kg, Faya 2013
        
    % Engine and exhaust temperatures
    T_engine(t) = 273.15 + (-2*10^(-6)*engine_rpm(t)^(2) + 0.0355*engine_rpm(t) + 77.5); %K
    T_exhaust(t) = 273.15 + 0.138*engine_rpm(t) - 17; 

    % Heat Flows
    Qleakage(t) = - leakage_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    Qvent(t) = vent_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    heat_human(t) = h_cabin*A_skin*abs(temperature(13)-T_skin) ;
    Qhuman(t) = N_Humans.*heat_human(t);
    Qequipment(t) = heat_equipment;
    Qbase(t) = 0;
    Qengine_calc(t) = A_engine*U_engine*abs(T_engine(t)-temperature(13));% - 11 - 44;
    Qexhaust_calc(t) = A_exhaust*U_exhaust*abs(T_exhaust(t)-temperature(14));%  + 17.77;
    
    % Engin Heat
    Qengine_observed(t)=0;
    Qexhaust_observed(t)=0;
    Qengine(t) = Qengine_observed(t) + Qengine_calc(t);
    Qexhaust(t) = Qexhaust_observed(t) + Qexhaust_calc(t);
    
    % Boundary conditions vector:
    if(Irr==0)
        if(size(T_cell)==[1 1]) % IF T amb is a constant value:
            Tbc=[T_cell;0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)];
        else % IF T amb is a vector from input (size(T_amb)==[1 1801]):
            Tbc=[T_cell(t);0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)];
        end
    else  % IF there is irradiance and T amb is a vector
        Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_sidewindows*G_sidewindows_t(t)) + (A_doors*G_doors_t(t)) );
        Tbc=[T_cell(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);...
            -A_sidewindows*G_sidewindows_a(t);-A_doors*G_doors_a(t);0;0;0;0;-Qhuman(t)-Qengine(t)-Qvent(t);-Qequipment(t)-Qexhaust(t)-Qleakage(t)]; 
    end

    % Temperatures calculation
    temperature=(inv(K - C))*(Tbc - C*delta_temperature);
    delta_temperature=temperature;

   % Exchange from the front and back:
    if ( temperature(13)>T_cell(t) )
        temperature(13) = ( density_air*CabinVolume/2*cp_air*delta_temperature(13)/timestep + UA_front*T_cell(t) )/...
            ( density_air*CabinVolume/2*cp_air/timestep + UA_front );
        % UA_front = 3000;
        % Qcv_emitted(t) = -(temperature(13)-T_cell(t))*UA_front;
        % temperature(13) = delta_temperature(13) + timestep/(cp_air*density_air*CabinVolume/2)*(Qhuman(t)+Qengine(t)+Qvent(t)+Qequipment(t)+Qexhaust(t)+Qleakage(t)+Qcv_emitted(t));
        delta_temperature(13)=temperature(13);
        check_f=check_f+1;
        Qcv_emitted(t) = -(temperature(13)-T_cell(t))*UA_front;
        Qcv_received(t) = 0;
    else
        Qcv_emitted(t) = 0;
        Qcv_received(t) = (T_cell(t)-temperature(13))*UA_front;
    end

    if ( temperature(14)>T_cell(t) )
        temperature(14) = (UA_back*T_cell(t) + density_air*CabinVolume/2*cp_air*delta_temperature(14)/timestep)...
            / (UA_back + density_air*CabinVolume/2*cp_air/(timestep));
        delta_temperature(14)=temperature(14);
        check_b=check_b+1; 
        % Qcv_emitted(t) = Qcv_emitted(t)-(temperature(14)-T_cell(t))*UA_back; % cumulative
        Qcv_emitted(t) = -(temperature(14)-T_cell(t))*UA_back;
        Qcv_received(t) = 0;
    else
        Qcv_emitted(t) = 0;
        %Qcv_received(t) =
        %Qcv_received(t)+(T_cell(t)-temperature(14))*UA_back; % cumulative
        Qcv_received(t) = (T_cell(t)-temperature(14))*UA_back;
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
    
    % Time to reach the T target, ttarget
        % ttarget(t) = log(abs(T0cabin-Ttarget)) * ...
        %             (CabinVolume*density_air*cp_air + mass_seats*cp_seats+144240)*...
        %             (Tcabin(t)-Ttarget)/(Qmac_needed(t));
    ttarget(t) = 600;
    tcons= ttarget/abs(log(Ttarget-T0cabin));
    
    % Thermal Loads
    Qrd(t) = A_ws*G_ws_t(t) + A_rw*G_rw_t(t) + A_sidewindows*G_sidewindows_t(t) + A_doors*G_doors_t(t) + Qbase(t);
    % Qcv_received(t) = abs(T_amb(t)-temperature(13))*UA_front + abs(T_amb(t)-temperature(14))*UA_back ;
    Qcabin_received(t) = Qhuman(t) + Qequipment(t) + Qengine(t) + Qexhaust(t) + Qcv_received(t) + Qvent(t) + Qleakage(t) + Qrd(t);
    Qtotal(t) = Qcabin_received(t) + Qcv_emitted(t);
    Qtarget(t) = (CabinVolume*density_air*cp_air + mass_seats*cp_seats + 144240)*(Ttarget-Tcabin(t))/(ttarget(t));
    Qmac_needed(t) = Qcabin_received(t) + Qcv_emitted(t) + Qtarget(t);
    Qcompressor(t) = Qmac_needed(t)/COP;
    energy_cum_comp_needed_Wh(t) = cumtrapz(timestep.*Qcompressor(t)./3600);
        % energy_cum_comp_needed_Wh = cumtrapz(energy_comp_needed_Wh); %

    % CO2 Emissions
    compressor_time(t) = ttarget(t); % OR sum(t_active)
    CO2(t) = abs(Qcompressor(t))*compressor_time(t)*Vpe_petrol/1000/3600*CF_petrol/distance_driven;

    % Compressor activation time
    if abs(Qcompressor(t))>compressor_on_lower_limit
        t_active(t) = 1;
    end
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

%% Energies
% Cabin thermal load
energy_cabin_Wh = timestep.*Qcabin_received(2:end)./3600;
energy_cum_cabin_Wh = cumsum(energy_cabin_Wh);
total_energy_cum_cabin_Wh=energy_cum_cabin_Wh(end);
% MAC needed cooling load
energy_mac_needed_Wh = timestep.*Qmac_needed(2:end)./3600;
energy_cum_mac_needed_Wh = cumsum(energy_mac_needed_Wh);
total_energy_cum_mac_needed_Wh=energy_cum_mac_needed_Wh(end);
% Compressor
%%tot_energy_compressor = total_energy_cum_mac_needed_Wh/COP
energy_comp_needed_Wh = timestep.*Qcompressor(2:end)./3600;
energy_cum_comp_needed_Wh = cumtrapz(energy_comp_needed_Wh); %cumsum
total_energy_cum_compr_needed_Wh=energy_cum_comp_needed_Wh(end)

%% Powers
cabin_in_W_av = sum(Qcabin_received)/Total_time
cabin_out_W_av = sum(Qcv_emitted)/Total_time
Qtarget_W_av = sum(Qtarget)/Total_time;
Qcabin_W_av = cabin_in_W_av - abs(cabin_out_W_av)
% Qmac_needed_W_av = sum(abs(Qmac_needed))/Total_time %cooling load
Qcompressor_W_av = sum(abs(Qcompressor))/Total_time

%% CO2 Emissions
CO2_total = CO2(end);
CO2_perc = CO2_total/95*100;

end %% END OF FUNCTION 'CABINMODEL'




%% Save engine load to be applied in other scenarios:
%save Engine_thermal_loads_W_35G7_v9 Qengine Qexhaust

%% PLOTS
% % Deletion of first 0 value for plotting (inicialization value):
% x=T_cell(2); T_cell(1)=x;
% % The boundary conditon 'T_amb' needs to be a vector to be plotted:
% if (size(T_cell)==[1 1])
%     T_cell = ones(1,Total_time)*(T_amb); %[K], vector
% end
% if G7_22==1
%     test = ['G7_22'];
% elseif G7_35==1
%     test = ['G7_35'];
% elseif G8_22==1
%     test = ['G8_22'];
% elseif G8_35==1
%     test = ['G8_35'];
% end
% 
% figure(7)
% plot(energy_comp_needed_Wh)
% hold on
% yyaxis right
% plot(energy_cum_comp_needed_Wh)
% 
% figure(8)
% plot(Qcompressor)
% hold on
% yline(-Qcompressor_W_av)
% 
% % figure(13)
% % plot(time,Tcabin-273.15)
% % % legend('T target','T cabin','Location','best')
% 
% figure(12)
% plot(time(3:Total_time),Qtarget(3:Total_time),':k','LineWidth',1.5)
% hold on
% plot(time(3:Total_time),Qmac_needed(3:Total_time),'-k','LineWidth',1.5)
% grid minor
% ylabel('MAC Heat [W]')
% xlabel('Time [s]')
% if(G7_22==1)
%     ylim([-1200 200])
% elseif(G8_22==1)
%     ylim([-400 200])
% end
% yyaxis right
% ax = gca;
% ax.YColor = [1 0 0];
% hold on
% plot(time(3:Total_time),Qtotal(3:Total_time),'-r','LineWidth',1.5)
% ylabel('Total Cabin Heat [W]','Color',[1 0 0])
% if(G8_22==1 || G7_22==1 || G7_35==1 || G8_35==1)
%     ylim([-5 130])
% end
% legend('Target','MAC','Cabin Total','Location','best')
% % saveas(figure(12),strcat(test,'_4'),'png')
% 
% figure(11)
% plot(time(3:Total_time),Qtotal(3:Total_time),'LineWidth',1.5)
% hold on
% plot(time(3:Total_time),Qcabin_received_plot(3:Total_time),'LineWidth',1.5)
% plot(time(3:Total_time),Qemitted(3:Total_time),'LineWidth',1.5)
% ylabel('Cabin Heat Loads [W]')
% xlabel('Time [s]')
% legend('Cabin Total','Cabin Received','Cabin Emitted','Location','best')
% grid minor
% % saveas(figure(11),strcat(test,'_3'),'png')
% 
% % figure(10)
% % plot(time,Qtotal,'--k','LineWidth',1.5)
% % hold on
% % plot(time,Qcv_received,':','LineWidth',1.5)
% % plot(time,Qhuman,':','LineWidth',1.5)
% % plot(time,Qequipment,':','LineWidth',1.5)
% % plot(time,Qvent,':','LineWidth',1.5)
% % plot(time,Qengine,':','LineWidth',1.5,'Color','g')
% % plot(time,Qexhaust,':','LineWidth',1.5)
% % plot(time,Qcabin_received_plot,'-','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
% % ylabel('Thermal Loads [W]')
% % xlabel('Time [s]')
% % yyaxis right
% % ax = gca;
% % ax.YColor = [0.4940 0.1840 0.5560];
% % %ylim([-1000 600])
% % hold on
% % plot(time,Qmac_needed_plot,'-','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560])
% % ylabel('Q cooling MAC [W]','Color',[0.4940 0.1840 0.5560])
% % %set(gca,'Ydir','reverse')
% % legend('Ambient','Human','Equipment','Vents','Engine',...
% %     'Exhaust','Cabin Received','MAC','Location','northeastoutside')
% % grid minor
%     % xlim([0 1800])
%     % axis([1 Total_time -1500 1700])
%     % axis([1 Total_time 0 1100])
% 
% % figure(9)
% % colororder({'k','k'})
% % plot(time,Qengine,':b','LineWidth',1.5)
% % ylabel('Engine Heat [W]')
% % xlabel('Time [s]')
% % ylim([-1 400])
% % hold on
% % yyaxis right
% % plot(time,engine_rpm,'-m','LineWidth',1)
% % ylabel('Engine Speed [rpm]')
% % legend('Engine Heat','Engine Speed','Location', 'Best')
% % grid minor
% % axis([1 Total_time 0 5500]) %745
% 
% % figure(8)
% % ylabel('Calculated [W]')
% % xlabel('Time [s]')
% % plot(time,Qengine_calc)
% % hold on
% % plot(time,Qexhaust_calc)
% % yyaxis right
% % ylabel('Observed [W]')
% % plot(time,Qengine_observed)
% % plot(time,Qexhaust_observed)
% % legend('Qengine calc','Qexhaust calc','Qengine observed','Qexhaust observed')
% 
% % figure(7)
% % plot(time,Qcabin_received_plot) 
% % xlabel('Time [s]') 
% % hold on
% % plot(time,-Qcv_emitted)
% % ylabel('Cabin Heat [W]')
% % legend('Received','Emitted')
% 
% % figure(6)
% % plot(time,Qcabin_received_plot-abs(Qcv_emitted)) 
% % xlabel('Time [s]') 
% % hold on
% % % plot(time,-Qcabin_emitted)
% % ylabel('Cabin Heat [W]')
% % legend('Received','Emitted')
% 
% % figure(5)
% % colororder({'k','k'})
% % plot(time,Qcabin_received_plot,':b','LineWidth',1.5)
% % grid minor
% % ylabel('Heat received (ambiebt to cabin) [W]')
% % xlabel('Time [s]') 
% % hold on
% % yyaxis right
% % plot(time,Qcv_emitted,':m','LineWidth',1.5)
% % set(gca,'Ydir','reverse')
% % xlim([0 1800])
% % ylim([-700 0])
% % ylabel('Heat emitted (cabin to ambient) [W]')
% % legend('Received','Emitted','Location','Best')
% 
% % figure(4)
% % xlabel('Time [s]') 
% % plot(time,Qmac_needed_plot)
% % ylabel('Heat needed to get T target [W]')
% 
% % figure(3)
% % plot(time,Qengine)
% % hold on
% % plot(time,Qexhaust)
% % ylabel('Heat [W]')
% % xlabel('Time [s]') 
% % yyaxis right
% % plot(time,Eeng_kWh)
% % hold on
% % plot(time,Etail_kWh)
% % ylabel('Energy [kWh]')
% % legend('Q engine W','Q tailpipe W','E engine kWh','E tailpipe kWh')
% 
% % plot(time,Q_cv_cabin,':','LineWidth',1.5)
% % hold on
% % plot(time,Q_cv_front,time,Q_cv_back,time,Q_cv_top)
% % hold on
% % plot(time,Q_cv_left)
% % legend('Cabin','Front','Back','Top','Left/Right','Location','best')
% % xt=[1000 800 1000 1000 1000];
% % yt = [Q_cv_front(end) Q_cv_back(end) Q_cv_left(end) Q_cv_top(end) Q_cv_cabin(end)];
% % str = [string(Q_cv_front(end)) string(Q_cv_back(end)) string(Q_cv_left(end)) string(Q_cv_top(end)) string(Q_cv_cabin(end))];
% % text(xt,yt,str)
% 
% figure(2)
% plot(time,TC_cell,'-k','LineWidth',1) %T amb 
% hold on
% plot(time,TC_front,'-','LineWidth',1.5) 
% plot(time,TC_back,'-','LineWidth',1.5)
% plot(time,TC_vent,'-','LineWidth',1.5)
% plot(time,TC_driver,'-','LineWidth',1.5)
% plot(time,TC_backpass,'-','LineWidth',1.5)
% hold off
% legend('TC cell','TC front','TC back','TC vent','TC driver','TC equipment','Location','best')
% grid minor
% xlabel('Time [s]') 
% ylabel('Temperature [ºC]')
% if(G7_22==1 || G8_22==1)
%     axis([1 Total_time 20 25])
% elseif(G7_35==1 || G8_35==1)
%     axis([1 Total_time 32 38])
% end
% % saveas(figure(2),strcat(test,'_2'),'png')
% 
% figure(1)
% plot(time,Tcabin_front-273.15,':b','LineWidth',1.5)
% hold on
% plot(time,Tcabin_back-273.15,':m','LineWidth',1.5)
% plot(time,TC_cell,'-k','LineWidth',1) %T amb 
% plot(time,TC_front,'-b','LineWidth',1.5) 
% plot(time,TC_back,'-m','LineWidth',1.5)
% hold off
% legend('T front simulated','T back simulated','TC ambient','TC front','TC back','Location','best')
% grid minor
% xlabel('Time [s]') 
% ylabel('Temperature [ºC]')
% if(G7_22==1 || G8_22==1)
%     axis([1 1800 20 25])
% elseif(G7_35==1 || G8_35==1)
%     axis([1 1800 33 38])
% end     
% % saveas(figure(1),strcat(test,'_1'),'png')
