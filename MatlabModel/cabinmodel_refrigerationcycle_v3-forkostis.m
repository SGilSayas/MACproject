clc;
clear;
close all;
clear py;
%%
% Working on V4: 
% Trying to understand/fix the following:
% Front volume shares in C_cabin_front and _back influence a lot where the
% model is droping the simulated temperature to very bad resutls. In that
% valley, the PID error is dropping from 2C until 90C, afterwards starting
% from the original error but never reaching 0 (so our simulated cabin
% temperature is never equal to the set target tmperature (set-point)). In
% that valley, the heat exchanged from the front -calculated apart**- is 0,
% meaning that the heat needs to leave the front to work, and when not 
% happening model doesnt work. But physically, that heat shouldnt go out 
% during those seconds (when heat exchange lag is zero). 
% **: code line around 532 "%% Exchange from the front and back".
%
% Giorgos' suggestion: check units. Delete the section "exchange from the 
% front and back" and merge that if statement in the definition of the 
% conductances matrix, allowing the heat to go in or out of the cabin 
% already in that calculation. How it's now, the conductances matrix 
% calculate the eat going in, and if the T cabin is higher that the T 
% ambient, another heat flow appears going out.
%
% - Expected compressor behaviour: Compressor power having a peak, and 
% after slowly going down and remaining on at a lower power level. Never
% negative.
% - Expected refigerant mass flow (mf): beween 20 and 80 g/s, never
% negative
% - Expected PID error: going towards zero.
% - Expected Temepratures: starting around 35 C and getting to 22 C.
%
% V3:
% - Initial heat exchanged cabin Air with the HVAC defined with T target
% - Addition of PID CONTROLLER for Compressor Speed. Input = error
% between current cabin temp and target cabin temperature. Output =
% compressor rpm, which later is used to calculate the refrigerant mass
% flow, which is used to obtain compressor's power and evaporator's heat.
%
% V2: 
% - better check of compression ratio with an average minimum, 
% - POINT (2) is not saturated liquid, but a saturated mix, so input
%   parameter is s2is = s1, not T2is=Tcond (COP=1.7)
% - DECIDE: mass flow fixed or dependant on Q_cabi_req? 
%   fixed doesnt work but makes more physical sense
%
% V1:
% Combination cabin model and the refrigertion cycle.
% Cabin model: cabinmodelfunction4 (based on cabinmodel_v22)
% Refrigeration cycle: standard cycle (evap-comp-cond-valve), 
% using Python prompt CoolProp to get refrigerant properties (p-h diagram)
% Refrigerant used: R1234yf
%
% Gil-Sayas, S., 2025
%
%% If it appears a message not finding CoolProp: Redo link python-matlab
% in matlab terminal: 
% pyenv('Clear');
% pyenv('Version', 'C:\Users\susan\AppData\Python\Python310\python.exe')
% in matlab terminal: pyversion
% if all ok: check this on cmd window: pip install coolprop
% if all ok: check this on matlab terminal py.CoolProp.CoolProp.PropsSI('T','P',101325,'Q',0,'Water')

%% Test conditions
t = 1;
timestep = 1;
duration_WLTC = 1800; % s
Total_time = duration_WLTC;
time = 1:timestep:Total_time;
time = transpose(time);

N_Humans = 1;

fuel = 'PETROL'; %  or 'DIESEL';

% Cabin target temperature
T_target = 22 + 273.15; % K
T_target = 22 + 273.15;

% V1: If fixing time to reach target temperature
    % T0cabin= (Tcabin_front(1)+Tcabin_back(1))/2; % K
    % time_target = ones(1,Total_time+1)* 900; %s, time to reach T_target from T0cabin -> TBD in the loop
    % tcons= time_target/log(abs(T0cabin-T_target)); % pull-down constant: overall pull-down time

%% Standard copressor's isoentropic efficiency, eta_is = (h2s-h1)/(h2-h1)
eta_is = 0.8; % av: 0.8; 
compressor_on_lower_limit = 600; % W

%% T between HX-sources: 
% heat exchangers (evap&cond) and cold&hot sources (cabin&external air)
delta_T_evap = 7;   % Automotive: 7-12  ; HVAC: 5-10
delta_T_cond = 10;  % Automotive: 10-15 ; HVAC: 5-15

%% Ambient conditions
Irr_scalar = 0; % W/m2 
Irr = ones(1,Total_time+1)*Irr_scalar;

% If not loading data:
    % T_amb_constant = 35; %C
    % T_amb = ones(1,Total_time+1)*T_amb_constant + 273.15; %K % Row vecto

% Load ambient temperature
load G8_35_cell_front_back.mat TC_cell TC_front TC_back
T_amb = TC_cell + 273.15 ; %K
T_amb(Total_time) = T_amb(Total_time-1);
T_amb(Total_time+1) = T_amb(Total_time-1);
TC_front(Total_time) = TC_front(Total_time-1);
TC_front(Total_time+1) = TC_front(Total_time-1);
TC_back(Total_time) = TC_back(Total_time-1);
TC_back(Total_time+1) = TC_back(Total_time-1);

T_ini = T_amb(1); % alternative: T_ini = TC_cabin_ini + 273.15;

%% Select and load working fluid operation pressures
fluid = {'R1234yf'}; %{'Water'}; % %{'R134a'};
R1234yf_op_pres = readtable('R1234yf_operating_pressures.xlsx','VariableNamesRow',1);
    % Verify content: disp(R1234yf_op_pres);

%% Constant input data

% Load input cell humidity (RH) and pressure (P)
load('G7_35_off_RH_P_amb')
RH_amb_mean=mean(RH_amb);
P_amb_kPa_mean=mean(P_amb_kPa);
RH_amb = ones(1,Total_time+1)*RH_amb_mean;
P_amb_kPa = ones(1,Total_time+1)*P_amb_kPa_mean;

% Load input Engine and exhaust heat flows
load('heatflow_engine_exhaust_W.mat')
Qengine_av = ones(1,Total_time+1)*mean(Qengine);
Qexhaust_av = ones(1,Total_time+1)*mean(Qexhaust);

% Cabin dimensions for category C
vehicle_height = 1.46;
vehicle_width = 1.8;
vehicle_lenght = 4.36;
cat_factor = 1.1;
cabin_height = vehicle_height - 0.06 - 0.6215/2; % minus roof&base thickness and half of wheel
cabin_width = vehicle_width - 0.06 - 0.241/2; % minus doors thickness % mirrors=0.241, we use half bc some database width included mirrors and others did not.
cabin_roof_lenght = vehicle_lenght*0.43;
cabin_base_lenght = vehicle_lenght*0.35;
seatcm3 = 150*cat_factor; %115; % cm3
seats_volume = 5*seatcm3*10^(-6); % m3
V_cabin = (cabin_width*cabin_roof_lenght*cabin_height/2)+(cabin_width*cabin_base_lenght*cabin_height/2) - seats_volume; 
V_cabin = V_cabin*1;

% A=area[m2] 
A_ws=0.63*1.3*cat_factor;
A_rw=0.29*1*cat_factor;
A_roof=1.8*1.1*cat_factor;
A_sidewindows=2*1.45*0.29*cat_factor;
A_doors=2*1.45*0.29*cat_factor;
A_base = cabin_base_lenght*cabin_width*3;
A_dashboard = cabin_width*0.5; %m2
A_front = A_ws + A_dashboard; %3
A_back = A_rw; %1.5
A_side = A_sidewindows + A_doors + A_roof;
% A_front = A_ws + A_sidewindows/2 + A_doors/2 + A_roof/2;
% A_back = A_rw + A_sidewindows/2 + A_doors/2 + A_roof/2;

    % e=thickness[m]
e_ws=0.006;
e_rw=0.005;
e_ABS=0.002;
e_steel=0.0005;
e_PU = 0.025; % polyurethane 
e_ceiling=e_ABS+e_steel+e_PU;
e_sidewindows=0.003;
e_doors=e_ceiling;

% Fuel (addition of diesel and electric option for next dataset validation)
Vpe_petrol = 0.264; % l/kWh
CF_petrol = 2330; % gCO2/l
Vpe_diesel = 0.22; % l/kWh
CF_diesel = 2640; % gCO2/l
if strcmp(fuel, 'PETROL') || strcmp(fuel, 'PETROL/ELECTRIC')
    Vpe=Vpe_petrol;
    CF=CF_petrol;
elseif strcmp(fuel, 'DIESEL') || strcmp(fuel, 'DIESEL/ELECTRIC')
    Vpe=Vpe_diesel;
    CF=CF_diesel;
else
     disp('Specify PETROL or DIESEL as the last argument');
end

% Conduction Properties: k = thermoconductivity[W/(m*K)]
k_ws=0.8; % 0.8-1.4
k_rw=1.4; % 1.4-1.5
k_ABS=0.1; % Acrylonitrile butadiene styrene
k_steel=14.9;
k_PU = 0.022; % 0.022-0.028 W/mK polyurethane
k_ceiling=k_ABS*k_steel*k_PU;
k_sidewindows=1.4;
k_doors=k_ceiling;

% Heat transfer coefficient of interior cabin air, W/(K*m2)
h_cabin = 5; % 5-10

% Air properies
density_air = 1.18; % [Kg/m3] (at 25C)
cp_air = 1006;      % + 10*1500; % [J/ kg*K]

% Seats properties
mass_seats = 20;    % 20-30 kg
cp_seats = 2000;    % 500-2000 [J/kgK]

% Human and equipment heat % ISO 8996, Alberto Viti Corsi (IDAE)
heat_equipment = 40; % W
Qequipment = ones(1,Total_time+1)*heat_equipment; % W
A_skin = 1.5; % m2
T_skin = 24.648 + 273.15; % K

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
rho_roof=0.04; %G7:0.74;
rho_base=0.3;
    % tau=transmissivity[]
tao_ws=0.452;
tao_rw=0.311;
tao_sidewindows=0.475;
tao_doors=0;
tao_roof=0.9;%0;
tao_base=0;
    % alpha=absorptivity[]  
alpha_coef=1;
alpha_ws = alpha_coef*(1-rho_ws-tao_ws);
alpha_rw = alpha_coef*(1-rho_rw-tao_rw);
alpha_sidewindows = alpha_coef*(1-rho_sidewindows-tao_sidewindows);
alpha_doors = alpha_coef*(1-rho_doors-tao_doors);
alpha_roof = alpha_coef*(1-rho_roof-tao_roof);
alpha_base = alpha_coef*(1-rho_base-tao_base);

% Global Horizontal Irradiance (Hypothesis: same in all windows)
G_roof_inc=Irr;
G_ws_inc=Irr; 
G_rw_inc=Irr;
G_sidewindows_inc=Irr;
G_doors_inc=Irr;
G_ws_r=rho_ws*G_ws_inc;
G_ws_a=alpha_ws*G_ws_inc; %bc
G_ws_t=tao_ws*G_ws_inc;
G_rw_r=rho_rw*G_rw_inc;
G_rw_a=alpha_rw*G_rw_inc; 
G_rw_t=tao_rw*G_rw_inc;
G_sidewindows_r=rho_sidewindows*G_sidewindows_inc;
G_sidewindows_a=alpha_sidewindows*G_sidewindows_inc; 
G_sidewindows_t=tao_sidewindows*G_sidewindows_inc;
G_doors_r=rho_doors*G_doors_inc;
G_doors_a=alpha_doors*G_doors_inc; 
G_doors_t=tao_doors*G_doors_inc;
G_roof_a=alpha_roof*G_roof_inc;
% Not needed bc rad is added as G*A --> Qirr = Irr.*(A_front+A_side+A_back); % W

% Air vent volume and flow rate
vent_volumerate=0.0003; % m3/s, Faya 2013: 0.02; 0.001
leakage_volumerate = 0.0001; %% m3/s, Faya2014

% Average humidity ratio in gram of water per gram of dry air, X
Ps_22_kPa = 2.626; % kPa, Water saturation pressure at 22C
Ps_35_kPa = 5.63; % kPa, Water saturation pressure at 35C
Ps_kPa = Ps_22_kPa+(T_ini-22)/(35-22)*(Ps_35_kPa-Ps_22_kPa);
meanP_kPa_amb = mean(P_amb_kPa); % in kPa
X = 0.0732; % X = 0.62198*Ps_kPa.*RH_amb./(100.*P_amb_kPa - Ps_kPa.*RH_amb); % Faya, 2013

%% Calculation of Properties

% Enthalpies calculation. Humidity considered same in cabin and exterior
e_amb = 1006.*T_amb + (2501000 + 1770.*T_amb)*X; % J/kg, Faya 2013, % before ".*X"

% Capacitance Matrix
C_amb = 0; % boundary condition
C_roof = 0; %for now
C_ceiling = 0; %for now
C_base = mass_seats*cp_seats + 144240; %Capacitance, W/K  
C_ws=0;
C_rw=0;
C_sidewindows=0;
C_doors=0;
C_cabin_front = density_air*cp_air/timestep*V_cabin*0.8; % front air volume = 80% of total volume % other hypothesis: 0.9;
C_cabin_back = density_air*cp_air/timestep*V_cabin*0.2; % back air volume = 20% of total cabin volume % other hypothesis: 0.001;
C=zeros(14);
C(1,1)=C_amb;
C(2,2)=C_roof;
C(3,3)=C_ceiling;
C(4,4)=C_base;
C(5,5)=C_ws;
C(6,6)=C_rw;
C(7,7)=C_sidewindows;
C(8,8)=C_doors;
C(9,9)=C_ws;
C(10,10)=C_rw;
C(11,11)=C_sidewindows;
C(12,12)=C_doors;
C(13,13)=C_cabin_front;
C(14,14)=C_cabin_back;

%% Initialization of values
% Initial heat transfer coefficients, UA (W/K):
h_front(t) = 160;
h_side(t) = 0.9*h_front(t);
h_rear(t) = 0.1*h_front(t);
h_front = ones(1,Total_time+1)*h_front(t);
h_side = ones(1,Total_time+1)*h_side(t);
h_rear = ones(1,Total_time+1)*h_rear(t);

UA_front(t) = A_front*h_front(t) + 1*(A_side*h_side(t));
UA_back(t) = A_back*h_rear(t) + 1*(A_side*h_side(t));
UA_front = ones(1,Total_time+1).*UA_front(t);
UA_back = ones(1,Total_time+1).*UA_back(t);
    % UA_front(t) = (k_ws*A_ws)/e_ws + ...
    %     0.5*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;
    % UA_back(t) = (k_rw*A_rw)/e_rw + ...
    %     0.5*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;

%  Conductances Matrix 
K1=h_side(t)*A_roof;
K2=A_roof*k_ceiling/e_ceiling;
K3=h_cabin*A_roof;
K4=h_cabin*A_base;
K5=h_front(t)*A_ws;
K6=A_ws*k_ws/e_ws;
K7=h_cabin*A_ws;
K8=h_rear(t)*A_rw;
K9=A_rw*k_rw/e_rw;
K10=h_cabin*A_rw;
K11=h_side(t)*A_sidewindows;
K12=A_sidewindows*k_sidewindows/e_sidewindows;
K13=h_cabin*A_sidewindows;
K14=h_side(t)*A_doors;
K15=A_doors*k_doors/e_doors;
K16=h_cabin*A_doors;
K17=h_cabin*cabin_width*cabin_height;
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
    0 0 (K3/2) (K4/2) 0 0 0 0 0 (K10) (K13/2) (K16/2) (K17) ...
    -(K3/2+K4/2+K10+K13/2+K16/2+K17)];

% initial cabin surfaces tmperatures
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
Tcabin=ones(1,Total_time+1)*(Tcabin_front(t)+Tcabin_back(t))/2; 
temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;...
    T_ini;T_ini;T_ini;T_ini;T_ini];
prev_temp=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;...
    T_ini;T_ini;T_ini;T_ini;T_ini];

% Refrigerant cycle:
    % point 1: evaporator oulet = compressor inlet
    % point 2: compressor outlet = condenser inlet
    % point 3: condenser outlet = valve inlet
    % point 4: valve outlet = evaporator inlet
% initial parameters of refrigerant cycle followng p-h diagram
T1(t) = 0; % Temperature 1 [K]
P1(t) = 0; % Pressure 1 [Pa]
T2(t) = 0; % Temperature 2
P2(t) = 0; % Pressure 2 [Pa]
P2is(t) = 0; % Pressure 2' (isoentropic point) [Pa]
T3(t) = 0; % Temperature 3
P3(t) = 0; % Pressure 3 [Pa]
T4(t) = 0; % Temperature 4
P4(t) = 0; % Pressure 4 [Pa]
COP(t) = 0; % Coefficient of performance [-]

% initialization of heat flows
Qleakage = ones(1,Total_time+1);
Qhuman = ones(1,Total_time+1);
Qcabin_req = ones(1,Total_time+1);
Qcabin_received= ones(1,Total_time+1);
Qcv_emitted(t) = 0;
Qcv_emitted = ones(1,Total_time+1);
Qcabin_tot = zeros(1,Total_time+1);
Q_out_front = ones(1,Total_time+1);
Q_out_back = ones(1,Total_time+1);
Qtarget = ones(1,Total_time+1);
comp_cum_Wh = ones(1,Total_time+1);
Qcompressor = ones(1,Total_time+1).*0;
heat_human = ones(1,Total_time+1);
Q_evap_req = zeros(1,Total_time+1); % Heat exchanged in the evaporator [W]
Q_cond_req = zeros(1,Total_time+1); % Heat echanged in the condenser   [W]
W_comp(t) = 0; % Work performed by the compresor [W]
Q_evap(t) = 0; % Heat exchanged in the evaporator [W]
Q_cond(t) = 0; % Heat echanged in the condenser   [W]

%% v3: Initial heat exchanged of cabin Air with the HVAC, W
Q_MAC = ones(1,Total_time+1)*cp_air*V_cabin*density_air*(T_target - T_amb(1))/timestep;
%%
Q_MAC_zone = Q_MAC/2;

% previous versions initializations
cabin_received_cum_Wh = ones(1,Total_time+1);
cabin_received_Wh_WLTP = ones(1,Total_time+1);
mac_needed_Wh = ones(1,Total_time+1);
mac_needed_cum_Wh = ones(1,Total_time+1);
mac_needed_Wh_WLTP = ones(1,Total_time+1);
comp_Wh = ones(1,Total_time+1);
comp_cum_Wh = ones(1,Total_time+1);
comp_Wh_WLTP = ones(1,Total_time+1);

active_time(t) = 0;
active_time_simulated(t) = 0;
cooling_time_simulated(t) = 0;
heating_time_simulated(t) = 0; 
compressor_time(t)=0;

% initializations of flags, PID variables, and others
check_f = zeros(1,Total_time+1);
check_b = zeros(1,Total_time+1);
temp_error_prev = 0;
error_sum = 0;
temp_error(t) = 0; % Temperature error [-]
temp_error = zeros(1,Total_time+1);
mf(t) = 20/1000; % Compressor refrigerant mass flow [kg/s]
comp_speed = ones(1, Total_time+1)*1000; % Compressor speed [rpm]
e_cabin = ones(1,Total_time+1)*e_amb(1);
CO2 = ones(1,Total_time+1).*0;


%% CALCULATION LOOP
for t=2:Total_time+1 % better for for the mf
    
    % Overall heat transfer coefficients, UA W/(Km2):
    if t>1800-323 % extra high
        h_front(t) = 200;
    elseif t>1800-323-455 % high
        h_front(t) = 190;
    elseif t>1800-323-455-433 % medium
        h_front(t) = 170;
    else % low
        h_front(t) = 160;
    end
    h_side(t) = 0.9*h_front(t);
    h_rear(t) = 0.1*h_front(t);
    UA_front(t) = A_front*h_front(t) + 1*(A_side*h_side(t));
    UA_back(t) = A_back*h_rear(t) + 1*(A_side*h_side(t));

    %% Conductances Matrix 
    % Tamb --K1-> Troof --K2-> Tceiling --K3-> Tcabin_back
    K1=h_side(t)*A_roof;
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

    %% Heat Flows
    e_cabin(t) = 1006*Tcabin(t) + (2501000 + 1770*Tcabin(t))*X; % J/kg, Faya 2013 % Cabin air enthalpy
    Qleakage(t) = leakage_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    Qvent(t) = vent_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    Qbase(t) = alpha_base*(A_ws*G_ws_t(t)+A_rw*G_rw_t(t)+A_sidewindows*G_sidewindows_t(t)+A_sidewindows*G_doors_t(t));
    Qhuman(t) = N_Humans.*h_cabin*A_skin*abs(temperature(13)-T_skin);
    Qequipment(t) = heat_equipment;

    %% Boundary conditions vector for scalar or vector T amb inputs
    if(size(T_amb)==[1 1])
        Tbc=[T_amb;...
            -A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);...
            -A_sidewindows*G_sidewindows_a(t);-A_doors*G_doors_a(t);0;0;0;0;...
            -Qhuman(t)-Qengine_av(t)-Qvent(t)+Q_MAC_zone(t)+Qcv_emitted(t)/2;...
            -Qequipment(t)-Qexhaust_av(t)+Qleakage(t)+Q_MAC_zone(t)+Qcv_emitted(t)/2];
    else
        Tbc=[T_amb(t);...
            -A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);...
            -A_sidewindows*G_sidewindows_a(t);-A_doors*G_doors_a(t);0;0;0;0;...
            -Qhuman(t)-Qengine_av(t)-Qvent(t)+Q_MAC_zone(t)+Qcv_emitted(t)/2;...
            -Qequipment(t)-Qexhaust_av(t)+Qleakage(t)+Q_MAC_zone(t)+Qcv_emitted(t)/2];
    end

    %% Temperatures calculation
    temperature=(inv(K - C))*(Tbc - C*prev_temp);
    
    %% Exchange from the front and back:
    if (temperature(13) > T_amb(t)) 
        temperature(13) = (density_air*V_cabin/2*cp_air*prev_temp(13)/timestep + UA_front(t)*T_amb(t))/( density_air*V_cabin/2*cp_air/timestep + UA_front(t) );
        prev_temp(13)=temperature(13);
        check_f(t) = 1;
        Qcv_emitted(t) = (temperature(13)-T_amb(t))*UA_front(t);
        Qcv_received(t) = 0;
    else
        Qcv_emitted(t) = 0;
        Qcv_received(t) = (T_amb(t)-temperature(13))*UA_front(t);
    end
    if (temperature(14) > T_amb(t))
        temperature(14) = (density_air*V_cabin/2*cp_air*prev_temp(14)/timestep + UA_back(t)*T_amb(t))/(density_air*V_cabin/2*cp_air/timestep + UA_back(t));
        prev_temp(14)=temperature(14);
        check_b(t) = 1; 
        % Qcv_emitted(t) = Qcv_emitted(t)-(temperature(14)-T_amb(t))*UA_back(t);
        Qcv_emitted(t) = (temperature(14)-T_amb(t))*UA_back(t);
        Qcv_received(t) = 0;
    else
        Qcv_emitted(t) = 0;
        %Qcv_received(t) = Qcv_received(t)+(T_amb(t)-temperature(14))*UA_back(t);
        Qcv_received(t) = (T_amb(t)-temperature(14))*UA_back(t);
    end
%CHECKS
    Qcv_emitted(t) = 0;
    % to 0 does not change anything Qcv_received(t) = 0;

    % Heat requested from the MAC system (evaporator's or condenser's work), Q cabin
    Qcabin_received(t) = Qhuman(t) + Qequipment(t) + Qengine_av(t) + Qexhaust_av(t) + Qvent(t) + Qleakage(t) + Qcv_received(t);
    Qcabin_tot(t) = Qcabin_received(t) + Qcv_emitted(t);
    
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% v2 %%%%%%% 
    % Tcabin_av(t) = (temperature(13)+temperature(14))/2;
    % Qcabin_req(t) = cp_air*V_cabin*density_air*(T_target - Tcabin_av(t))/180; % 3 min
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% v3 %%%%%%% 
    % Qcabin_req(t)=Q_evap(t-1);
    % Q_evap_req(t) = abs(Qcabin_req(t-1));
%%
    % Define cooling or heating mode & evap/condenser temps
    if (T_amb(1) >= T_target)
       cooling = true;
       Q_evap_req(t) = abs(Qcabin_req(t-1));
       T_evap(t) = T_target - delta_T_evap;
       T_cond(t) = T_amb(t) + delta_T_cond;
    else % Heating
       cooling = false; 
       Q_cond_req(t) = abs(Qcabin_req(t-1));
       T_evap(t) = T_amb(t) - delta_T_evap;
       T_cond(t) = T_target + delta_T_cond;
    end
    
    % Call diagram p-h on CoolProp 
    Ref = char(fluid);
    
    % Check if T_amb is in the table
    idx = find(R1234yf_op_pres.T_ambient_C == T_amb(t));
    
    if ~isempty(idx) % If T_amb is found in the table
        P_LowSide_max_Pa = 1000*R1234yf_op_pres.P_LowSide_max_kPa(idx);
        P_HighSide_max_Pa = 1000*R1234yf_op_pres.P_HighSide_max_kPa(idx);
    else % If T_amb is not found, interpolate
        T_ambient = R1234yf_op_pres.T_ambient_C;
        P_LowSide_max_Pa = 1000*R1234yf_op_pres.P_LowSide_max_kPa;
        P_HighSide_max_Pa = 1000*R1234yf_op_pres.P_HighSide_max_kPa;
        
        P_LowSide_max_Pa = interp1(T_ambient, P_LowSide_max_Pa, T_amb, 'linear', 'extrap');
        P_HighSide_max = interp1(T_ambient, P_HighSide_max_Pa, T_amb, 'linear', 'extrap');
    end
    
    compression_ratio_min = 5.7357; % Average

    %% Refrigeration Cycle Points:
    % (1): evaporator outlet, salutared vapor; 
    % (2): compressor outlet, superheated vapor;
    % (3): condenser outlet, saturated fluid/liquid; 
    % (4): expansion device outlet, vapor-fluid mix;
    
    % (1) Evaporator's outlet = Compressor's inlet, Saturated vapor Q=1
    T1(t)= T_evap(t); % T_evap_out is known
    P1(t)=py.CoolProp.CoolProp.PropsSI('P','T',T1(t),'Q',1,Ref);
    P1(t) = double(P1(t));

    % Ensure evaporator pressure does not exceed the maximum low-side pressure
    if (P1(t) >= P_LowSide_max_Pa(t))
        P1(t) = P_LowSide_max_Pa(t);
        T1(t) = py.CoolProp.CoolProp.PropsSI('T','P',P1(t),'Q',1,Ref);
    end
    
    h1(t)=py.CoolProp.CoolProp.PropsSI('H','T',T1(t),'Q',1,Ref);     
    s1(t)=py.CoolProp.CoolProp.PropsSI('S','T',T1(t),'Q',1,Ref);
    
    % (2') Isoentropic compression (theoretical)
    s2is(t)=s1(t); 
    %% v2
    %v3: delete this h2is(t)=py.CoolProp.CoolProp.PropsSI('H','T',T_cond(t),'S',s2is(t),Ref); % falls inside of saturated mix area, so we can use T_cond
    P2is(t)=py.CoolProp.CoolProp.PropsSI('P','T',T_cond(t),'S',s2is(t),Ref);
    
    % (2) Compressor's outlet = Condenser's inlet
    P2(t)=P2is(t);

    % Actual compression ratio
    compression_ratio_actual(t)=P2(t)/P1(t);

    % Compression ratio needs to meet/exceed the minimum one
    if compression_ratio_actual(t) < compression_ratio_min
        % Adjust P2 so compression_ratio_actual = compression_ratio_min
        P2(t) = compression_ratio_min * P1(t);
        P2is(t)=P2(t);

        % Recompute T_cond based on adjusted P2
        T_cond(t) = py.CoolProp.CoolProp.PropsSI('T','P',P2(t),'Q',1,Ref);
        
        % Update delta_T_cond to reflect the new condenser saturation temperature
        delta_T_cond = T_cond(t) - T_amb(t);
    end
%%
    h2is(t)=py.CoolProp.CoolProp.PropsSI('H','P',P2(t),'S',s2is(t),Ref); % falls inside of saturated mix area, so we can use T_cond
    % or: h2is(t)=py.CoolProp.CoolProp.PropsSI('H','T',T_cond(t),'P',P2(t),Ref);

    % Get enthalpy from isoenthalic efficiency
    h2(t) = (h2is(t)-h1(t))/eta_is + h1(t);
    T2(t)=py.CoolProp.CoolProp.PropsSI('T','H',h2(t),'P',P2(t),Ref);
    s2(t)=py.CoolProp.CoolProp.PropsSI('S','H',h2(t),'P',P2(t),Ref);

    % (3) Condenser's outlet = Expansion vale's inlet, Saturated liquid
    % T3(t)=T_cond(t);
    P3(t)=P2(t);
    T3(t)=py.CoolProp.CoolProp.PropsSI('T','P',P3(t),'Q',0,Ref);
    h3(t)=py.CoolProp.CoolProp.PropsSI('H','P',P3(t),'Q',0,Ref);
    s3(t)=py.CoolProp.CoolProp.PropsSI('S','P',P3(t),'Q',0,Ref);

    % (4) Valve's outlet = Evaporator's inlet
    s4(t)=s3(t);
    P4(t)= P1(t);
    h4(t)=py.CoolProp.CoolProp.PropsSI('H','S',s4(t),'P',P4(t),Ref);
    T4(t)=py.CoolProp.CoolProp.PropsSI('T','S',s4(t),'P',P4(t),Ref);
    
    % State of each operating point
    Q1(t)=py.CoolProp.CoolProp.PropsSI('Q','S',s1(t),'P',P1(t),Ref);
    Q4(t)=py.CoolProp.CoolProp.PropsSI('Q','S',s4(t),'P',P1(t),Ref);
    
    % Coefficient of performance
    COP(t)=(h1(t)-h4(t))/(h2(t)-h1(t)); 
    % heat absorbed by the refrigerant in the evaporator/work done by the compressor

    %% v3: PID CONTROLLER for Compressor Speed
    % u(t) = Kp * e(t) + Ki * integral(e(t)) + Kd * derivative(e(t));
   
    % Controlled variable
    temp_error(t) = Tcabin(:,t) - T_target; % there's an opposite effect between error&comp_speed, so they are inversly proportional

    % PID gains
    Kp = 0.05;  % 0.05 Large immediate reaction on the output to bring the process value close to the set point
    Ki = 0.005; % The longer it takes for the process value to reach the set point, the more effect the integral will have on the output
    Kd = 0;     % If the process value is approaching the set point to fast, the derivative will limit the output to prevent the process value from overshooting the set point

    % Proportional, integral and derivative terms
    P = Kp * temp_error(t);
    error_sum = error_sum + temp_error(t) * timestep; % Total integral (accumulate error)
    I = Ki * error_sum;
    D = Kd * (temp_error(t) - temp_error(t-1)) / timestep;
    
    % Total control output
    comp_speed(t) = (P + I + D);

    % Relating Speed to Mass Flow
    k = 0.005; % (kg/s)/rmp, %%% 0.01 ; 0.0028 %%%%%%%%%%%%%%%%%%% TBD
    mf(t) = k*comp_speed(t); % Proportional relationship

    % trial:vLimitation of the mass flow -> doesnt help the simulation
        % if mf(t) > 50/1000
        %     mf(t) = 50/1000;
        %     comp_speed(t) = mf(t)/k;
        % elseif mf(t) < 20/1000
        %     mf(t) = 20/1000;
        %     comp_speed(t) = mf(t)/k;
        % end

    %% trial: not working
        % mf(t) = (P + I + D);
        % if mf(t) >= 5/1000 && mf(t) < 12/1000
        %     comp_speed(t) = 1000;
        % elseif mf(t) >= 12/1000 && mf(t) < 20/1000
        %     comp_speed(t) = 1500;
        % elseif mf(t) >= 20/1000 && mf(t) < 25/1000
        %     comp_speed(t) = 2000;
        % elseif mf(t) >= 25/1000 && mf(t) < 40/1000
        %     comp_speed(t) = 2500;
        % elseif mf(t) >= 40/1000 && mf(t) < 55/1000
        %     comp_speed(t) = 3000;
        % end

    %% v2:
    % % Mass flow calculation
    % if (cooling == true)
    %     mf(t) = Q_evap_req(t)/(h1(t) - h4(t));
    %     ratio = mf(t)/mf(t-1);
    % else
    %     mf(t) = Q_cond_req(t)/(h2(t) - h3(t));  
    %     ratio = (mf(t)- mf(t-1))/mf(t-1);
    % end
    % if mf(t-1) > 0
    %     % Keep mf changes between -5% and +5%
    %     if (ratio > 1.05)
    %         mf(t) = 1.05 * mf(t-1);
    %     elseif (ratio < 0.95)
    %         mf(t) = 0.95 * mf(t-1);
    %     end
    % end
    % % FIXED value: mf(t)=20*1000; %kg/s, small: 20-50 g/s, SUVS: 100-200, HDV: 200-500
    
    %% Components' heat loads calculation
    % Cooling capacity, J/s
    Q_evap(t) = mf(t)*(h1(t)-h4(t));
    % Heating capacity, J/s
    Q_cond(t) = mf(t)*(h3(t)-h2(t));
    % Compressor power, J/s
    W_comp(t) = mf(t)*(h2(t)-h1(t));

    if (cooling == true) Q_MAC(t) = Q_evap(t);
    else Q_MAC(t) = Q_cond(t);  
    end

    % Temperature recalculation after the MAC load:
        %Tcabin(t) = (Qcabin_req(t) - Q_MAC(t) )*(timestep/(V_cabin*cp_air*density_air)) + Tcabin(t-1);
    Q_MAC_zone(t) = Q_MAC(t)/2;
    temperature(13) = (Qcabin_tot(t) / 2 - Q_MAC_zone(t)) * (timestep / (V_cabin / 2 * cp_air * density_air)) + prev_temp(13);
    temperature(14) = (Qcabin_tot(t) / 2 - Q_MAC_zone(t)) * (timestep / (V_cabin / 2 * cp_air * density_air)) + prev_temp(14);

    % Save this timestep temperatures to be used in the following timestep
    prev_temp=temperature;

    % Add a timestep
    t = t + timestep;

    % Extract values
    Names = {'P1';'T1';'T2';'P2';'T3';'P3';'T4';'P4';'COP';'mf';'W_comp';'Q_evap';'Q_cond'};
    MAC_calcs = table(P1',T1',T2',P2',T3',P3',T4',P4',COP',mf',W_comp',Q_evap',Q_cond','VariableNames',Names);

    % Temperatures extraction
    Tamb(:,t)=temperature(1); %T_amb: boundary condition, Tamb: simulated
    Troof(:,t)=temperature(2);
    Tceiling(:,t)=temperature(3);
    Tbase_int(:,t)=temperature(4);
    Tw_ext_ws(:,t)=temperature(5);
    Tw_ext_rw(:,t)=temperature(6);
    Tw_ext_sidewindows(:,t)=temperature(7);
    Tw_ext_doors(:,t)=temperature(8);
    Tw_int_ws(:,t)=temperature(9);
    Tw_int_rw(:,t)=temperature(10);
    Tw_int_sidewindows(:,t)=temperature(11);
    Tw_int_doors(:,t)=temperature(12);
    Tcabin_front(:,t)=temperature(13);
    Tcabin_back(:,t)=temperature(14);

    Tcabin(:,t)= (Tcabin_front(:,t)+Tcabin_back(:,t))/2;
end
MAC_calcs = MAC_calcs(2:end,:);

%% Plots MAC components
figure(1)
plot(time, MAC_calcs.COP, 'LineWidth', 1);
xlabel('Time (s)')
ylabel('COP')
grid on

figure(2)
plot(time, MAC_calcs.mf*1000, 'LineWidth', 1); % kg/s to g/s
xlabel('Time (s)')
ylabel('Refrigerant mass flow, g/s')
grid on

figure(3)
plot(time, MAC_calcs.W_comp, 'LineWidth', 1);
xlabel('Time (s)')
ylabel('Compressor Work (W)')
grid on

figure(4)
plot(time, MAC_calcs.Q_evap, 'LineWidth', 1);
xlabel('Time (s)')
ylabel('Evaporator Heat Load (W)')
grid on

figure(5)
plot(time, MAC_calcs.Q_cond, 'LineWidth', 1);
xlabel('Time (s)')
ylabel('Condenser Heat Load (W)')
grid on

figure(6)
plot(time, Qcabin_received(2:end), 'LineWidth', 1);
xlabel('Time (s)')
ylabel('Cabin Heat Load (W)')
grid on

%% Plots Cabin Heat Loads
figure(7)
plot(time(3:Total_time),Qcabin_tot(3:Total_time),'LineWidth',1.5)
hold on
plot(time(3:Total_time),Qcabin_received(3:Total_time),'LineWidth',1.5)
plot(time(3:Total_time),Qcv_emitted(3:Total_time),'LineWidth',1.5)
ylabel('Cabin Heat Loads [W]')
xlabel('Time [s]')
legend('Cabin Total Req','Cabin Received','Cabin Emitted','Location','best')
grid minor

figure(8)
plot(time(3:Total_time),Qcabin_req(3:Total_time),':k','LineWidth',1.5)
hold on
plot(time(3:Total_time),Q_MAC(3:Total_time),'-k','LineWidth',1.5)
grid minor
ylabel('MAC Heat [W]')
xlabel('Time [s]')
yyaxis right
ax = gca;
ax.YColor = [1 0 0];
hold on
plot(time(3:Total_time),Qcabin_tot(3:Total_time),'-r','LineWidth',1.5)
ylabel('Total Cabin Heat [W]','Color',[1 0 0])
legend('Target Heat: Cabin & evap requested heat)','MAC Heat: Evap given heat','Cabin Total: present heat','Location','best')
% saveas(figure(12),strcat(test,'_4'),'png')

figure(9)
plot(Q_evap)
hold on
plot(-Qcabin_req)
legend('Evap heat given (Q MAC)','Cabin heat requested (= Q evap req)')

figure(10)
plot(time,Tcabin_front(3:end)-273.15,':b','LineWidth',1.5)
hold on
plot(time,Tcabin_back(3:end)-273.15,':m','LineWidth',1.5)
plot(time,T_amb(2:end)-273.15,'-k','LineWidth',1) 
% plot(time,TC_front(2:end),'-b','LineWidth',1.5) 
% plot(time,TC_back(2:end),'-m','LineWidth',1.5)
plot(time,Tcabin(3:end)-273.15,':k','LineWidth',1.5)
hold off
legend('T front simulated','T back simulated','TC ambient', 'T cabin simulated') %,'TC front','TC back','Location','best')
grid minor
xlabel('Time [s]') 
ylabel('Temperature [ºC]')

figure(11)
plot(time(3:Total_time),comp_speed(3:Total_time),':k','LineWidth',1.5)
hold on
plot(time(3:Total_time),100*MAC_calcs.mf(3:Total_time),'-k','LineWidth',1.5)
yyaxis right
plot(time, MAC_calcs.W_comp,'LineWidth',1);
grid minor
xlabel('Time [s]')
legend('Compressor speed [rpm]','Mass flow[kg/s]*100','Compressor work (Y)[W]')

figure(12)
plot(temp_error)
hold on
plot(check_f)
plot(check_b)
legend('error','front heat exchange flag','back heat exchange flag')