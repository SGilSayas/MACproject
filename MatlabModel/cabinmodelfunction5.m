function [Tcabin,Qcompressor,compressor_av_W,comp_Wh_WLTP,CO2_total] = cabinmodelfunction5(Total_time,distance_driven,N_Humans,COP,TC_cabin_ini,TC_cell_constant,Irr_constant,TC_cell,Irr,category,fuel)
    %% Version 5 includes:
    % WORKING ON refrigeration cycle (version: cabinmodel_refrigerationcycle)
    % a dynamic MAC active time, 
    % compressor lower limit now is 600 W (limit ins 1A: 1A*600V)
%%
    % Constant data
t = 1; timestep = 1; % s
heat_human = 80; %W
heat_equipment = 0 ; %40; % W
Vpe_petrol = 0.264; % l/kWh
CF_petrol = 2330; % gCO2/l
Vpe_diesel = 0.22; % l/kWh
CF_diesel = 2640; % gCO2/l

    %     % Inputs stored in CELL ARRAY:
    % Total_time = 1800; % s
    % distance_driven = 23.25; % km
    % N_Humans = 1;
    % COP = 3.5;
    % TC_cell_constant = 35; % C
    % Irr_constant = 0; % IN NEGATIVE (enters the cabin), W/m2

    % Filling the inpiut blanks with default values
if Total_time == 'N/A' Total_time = 1800; end
if distance_driven == 'N/A' distance_driven = 23.25; end
if N_Humans == 'N/A' N_Humans = 1; end
if COP == 'N/A' COP = 3.5; end
if TC_cell == 'N/A' TC_cell = ones(1,Total_time)*TC_cell_constant; end
if TC_cabin_ini == 'N/A'; TC_cabin_ini = TC_cell(1); end
if Irr == 'N/A' Irr = -ones(1,Total_time)*Irr_constant; end
if category == 'A'
    vehicle_height = 1.49;
    vehicle_width = 1.64;
    vehicle_lenght = 3.56;
    cat_factor = 0.9;
elseif category == 'B'
    vehicle_height = 1.47;
    vehicle_width = 1.73;
    vehicle_lenght = 4;
    cat_factor = 0.95;
elseif category == 'C'
    vehicle_height = 1.46;
    vehicle_width = 1.8;
    vehicle_lenght = 4.36;
    cat_factor = 1.1;
elseif category == 'D'
    vehicle_height = 1.457;
    vehicle_width = 1.832;
    vehicle_lenght = 1.457;
    cat_factor = 1.05;
elseif category == 'E'
    vehicle_height = 1.468;
    vehicle_width = 1.87;
    vehicle_lenght = 4.936;
    cat_factor = 1.1;
elseif category == 'F'
    vehicle_height = 1.479;
    vehicle_width = 1.903;
    vehicle_lenght = 5.172;
    cat_factor = 1.2;
elseif category == 'I'
    vehicle_height = 1.521;
    vehicle_width = 1.95;
    vehicle_lenght = 5.46;
    cat_factor = 1.3;
else
    % Validated Cabin Geometry
    vehicle_height = 0.5456 + 0.06 + 0.6215/2;
    vehicle_width = 1.3 + 0.06;
    vehicle_lenght = 4.1; % cabin_roof_lenght=1.8; cabin_base_lenght=1.45;
    cat_factor = 1;
end
if strcmp(fuel, 'PETROL') || strcmp(fuel, 'PETROL/ELECTRIC')
    Vpe=Vpe_petrol;
    CF=CF_petrol;
elseif strcmp(fuel, 'DIESEL') || strcmp(fuel, 'DIESEL/ELECTRIC')
    Vpe=Vpe_diesel;
    CF=CF_diesel;
else
     disp(['Specify PETROL or DIESEL as the last argument']);
end

cabin_height = vehicle_height - 0.06 - 0.6215/2; % minus roof&base thickness and half of wheel
cabin_width = vehicle_width - 0.06 - 0.241/2; % minus doors thickness % mirrors=0.241, we use half bc some database width included mirrors and others did not.
cabin_roof_lenght = vehicle_lenght*0.43;
cabin_base_lenght = vehicle_lenght*0.35;

seatcm3 = 150*cat_factor; %115; % cm3
seats_volume = 5*seatcm3*10^(-6); % m3
CabinVolume = (cabin_width*cabin_roof_lenght*cabin_height/2)+(cabin_width*cabin_base_lenght*cabin_height/2) - seats_volume; 
CabinVolume = CabinVolume*0.7;

    % Load input cell humidity (RH) and pressure (P)
load('G7_35_off_RH_P_amb')
RH_amb_mean=mean(RH_amb);
P_amb_kPa_mean=mean(P_amb_kPa);
RH_amb = ones(1,Total_time)*RH_amb_mean;
P_amb_kPa = ones(1,Total_time)*P_amb_kPa_mean;

    % Load input Engine and exhaust heat flows
load('heatflow_engine_exhaust_W.mat')
Qengine_av = mean(Qengine);
Qexhaust_av = mean(Qexhaust);

% Total_time = length(TC_cell); 
time = 1:timestep:Total_time;

    % Transpose TC amb vector and initialization
T_cell = transpose(TC_cell)+273.15; % T_amb in K. Both TC_cell or T_cell_OBD
% T_ini = T_cell(1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_ini = TC_cabin_ini + 273.15;
   
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

    % Conduction Properties: k=thermoconductivity[W/(m*K)]
k_ws=0.8; % 0.8-1.4
k_rw=1.4; % 1.4-1.5
k_ABS=0.1; % Acrylonitrile butadiene styrene
k_steel=14.9;
k_PU = 0.022; % 0.022-0.028 W/mK polyurethane
k_ceiling=k_ABS*k_steel*k_PU;
k_sidewindows=1.4;
k_doors=k_ceiling;

% e_av_front = (e_ws + e_rw + e_sidewindows + e_doors + e_ceiling)/5 ;
% k_av_front = (k_ws + k_rw + k_sidewindows + k_doors + k_ceiling)/5;
% e_av_back = (e_rw + e_sidewindows + e_doors + e_ceiling)/4 ;
% k_av_back = (k_rw + k_sidewindows + k_doors + k_ceiling)/4;

    % Heat transfer Convection coefficients, W/(m2*K)
h_ws=5; 
h_sidewindows=5;
h_doors=5;
h_rw=5;
h_ceiling=5;
h_roof=20;
h_ext=20; 
h_base = 5;
h_cabin=5; %10;

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
rho_roof=0.04; %G7:0.74;
rho_base=0.3;
    % tau=transmissivity[]
tao_ws=0.452;
tao_rw=0.311;
tao_sidewindows=0.475;
tao_doors=0;
tao_roof=0.9;%0;
tao_base=0;
    % alpha=absorptivity[] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
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
C_cabin_front = density_air*cp_air/timestep*CabinVolume*0.9;
C_cabin_back = density_air*cp_air/timestep*CabinVolume*0.001;
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

    % Humidity ratio in gram of water per gram of dry air, X
Ps_22_kPa = 2.626; % kPa, Water saturation pressure at 22C
Ps_35_kPa = 5.63; % kPa, Water saturation pressure at 35C
Ps_kPa = Ps_22_kPa+(T_ini-22)/(35-22)*(Ps_35_kPa-Ps_22_kPa); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanP_kPa_amb = mean(P_amb_kPa); % in kPa
X = 0.0732; %0.62198*Ps_kPa.*RH_amb./(100.*P_amb_kPa - Ps_kPa.*RH_amb); % Faya, 2013

    % Enthalpies calculation. Humidity considered same in amb and cabin
e_amb = 1006.*T_cell + (2501000 + 1770.*T_cell)*X; % J/kg, Faya 2013, % before .*X

    % Air vent volume and flow rate:
vent_volumerate=0.0003; % m3/s, Faya 2013: 0.02; 0.001
leakage_volumerate = 0.0001; %% m3/s, Faya2014

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
Tcabin=ones(1,Total_time)*(Tcabin_front(t)+Tcabin_back(t))/2; 
temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];
prev_temp=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];
check_f=0; check_b=0;
Qleakage = ones(1,Total_time);
Qtotal = ones(1,Total_time);
Qcabin_received= ones(1,Total_time);
Qcv_emitted = ones(1,Total_time);
Q_out_front = ones(1,Total_time);
Q_out_back = ones(1,Total_time);
Qtarget = ones(1,Total_time);
comp_cum_Wh = ones(1,Total_time);
Qcompressor = ones(1,Total_time).*0;
CO2 = ones(1,Total_time).*0;
Qcv_emitted(t) = 0;
active_time(t) = 0;
active_time_simulated(t) = 0;
cooling_time_simulated(t) = 0;
heating_time_simulated(t) = 0; 
compressor_time(t)=0;
heat_human = ones(1,Total_time);
e_cabin = ones(1,Total_time)*e_amb(1);
cabin_received_cum_Wh = ones(1,Total_time);
cabin_received_Wh_WLTP = ones(1,Total_time);
mac_needed_Wh = ones(1,Total_time);
mac_needed_cum_Wh = ones(1,Total_time);
mac_needed_Wh_WLTP = ones(1,Total_time);
comp_Wh = ones(1,Total_time);
comp_cum_Wh = ones(1,Total_time);
comp_Wh_WLTP = ones(1,Total_time);

    % MAC Thermal Load parameters:
Ttarget = 22 + 273.15; %K
T0cabin= (Tcabin_front(1)+Tcabin_back(1))/2; %K
ttarget = ones(1,Total_time+1)* 900; %s, time to reach Ttarget from T0cabin -> TBD in the loop
tcons= ttarget/log(abs(T0cabin-Ttarget)); % pull-down constant: overall pull-down time
Qmac = ones(1,Total_time);

    % Initial Overall heat transfer coefficients, UA (W/K):
UA_front(t) = (k_ws*A_ws)/e_ws + 0.5*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;
UA_back(t) = (k_rw*A_rw)/e_rw + 0.5*(k_sidewindows*A_sidewindows)/(e_sidewindows) ;
UA_front=ones(1,Total_time).*UA_front(t);
UA_back=ones(1,Total_time).*UA_back(t);
U_front=ones(1,Total_time);
U_side=ones(1,Total_time);
U_back=ones(1,Total_time);

A_skin=1.5; %m2
T_skin=24.648+273.15; %K
compressor_on_lower_limit = 600; % W, 

% temperature(13) = TC_cabin_ini;
% temperature(14) = TC_cabin_ini;

for t=2:Total_time
        
        % Overall heat transfer coefficients, UA (W/K):
    if t>1800-323 % extra high
        h_front(t) = 200;
        U_front(t) = 200;
    elseif t>1800-323-455 % high
        h_front(t) = 190;
        U_front(t) = 190;
    elseif t>1800-323-455-433 % medium
        h_front(t) = 170;
        U_front(t) = 170;
    else % low
        h_front(t) = 160;
        U_front(t) = 160;
    end

    h_side(t) = 0.9*h_front(t); % 0.95
    h_rear(t) = 0.1*h_front(t);
    U_side(t) = 0.9*U_front(t);
    U_back(t) = U_front(t);

    % h_side(t) = 0.9*h_front(t); % 0.95
    % h_rear(t) = 0.2*h_front(t);
    % U_front(t) = 1 / (1/h_front(t)+(e_ws/k_ws)); %*50000
    % U_back(t) = 1 / (1/h_rear(t)+(e_rw/k_rw));
    % U_side(t) = 1 / ( (1/h_side(t)+(e_sidewindows/k_sidewindows))...
    %                 + (1/h_side(t)+(e_doors/k_doors))...
    %                 + (1/h_side(t)+(e_ceiling/k_ceiling)) ); 
    %     % three surfaces: windows+doors+roof

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

    % UA_front = 500; %v15:185;
    % UA_back = 96;

        % Cabin air enthalpy
    e_cabin(t) = 1006*Tcabin(t) + (2501000 + 1770*Tcabin(t))*X; % J/kg, Faya 2013

        % Heat Flows
    Qleakage(t) = - leakage_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    Qvent(t) = vent_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    Qbase(t)=alpha_base*(A_ws*G_ws_t(t)+A_rw*G_rw_t(t)+A_sidewindows*G_sidewindows_t(t)+A_sidewindows*G_doors_t(t));
    heat_human(t) = h_cabin*A_skin*abs(temperature(13)-T_skin) ;
    Qhuman(t) = N_Humans.*heat_human(t);
    Qequipment(t) = heat_equipment;
    % if ( temperature(13)>T_cell(t) )
    %     Q_out_front(t) = UA_front(t)*(temperature(13)-T_cell(t))*500;
    % end
    % if (temperature(14)>T_cell(t))
    %     Q_out_back(t) = UA_back(t)*(temperature(14)-T_cell(t));
    % end

        % Boundary conditions vector
    if(size(T_cell)==[1 1])
        Tbc=[T_cell;-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);...
            -A_sidewindows*G_sidewindows_a(t);-A_doors*G_doors_a(t);0;0;0;0;...
            -Qhuman(t)-Qengine_av-Qvent(t)-Qmac(t);-Qequipment(t)-Qexhaust_av-Qleakage(t)];
    else 
        Tbc=[T_cell(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);...
        -A_sidewindows*G_sidewindows_a(t);-A_doors*G_doors_a(t);0;0;0;0;...
        -Qhuman(t)-Qengine_av-Qvent(t)-Qmac(t);-Qequipment(t)-Qexhaust_av-Qleakage(t)];
    end

        % Temperatures calculation
    temperature=(inv(K - C))*(Tbc - C*prev_temp);
    prev_temp=temperature;


        % Exchange from the front and back:
    if ( temperature(13)>T_cell(t) )
        temperature(13) = (density_air*CabinVolume/2*cp_air*prev_temp(13)/timestep + UA_front(t)*T_cell(t))/( density_air*CabinVolume/2*cp_air/timestep + UA_front(t) );
        prev_temp(13)=temperature(13);
        check_f=check_f+1;
        Qcv_emitted(t) = (temperature(13)-T_cell(t))*UA_front(t);
        Qcv_received(t) = 0;
    else
        Qcv_emitted(t) = 0;
        Qcv_received(t) = (T_cell(t)-temperature(13))*UA_front(t);
    end
    if ( temperature(14)>T_cell(t) )
        temperature(14) = (density_air*CabinVolume/2*cp_air*prev_temp(14)/timestep + UA_back(t)*T_cell(t))/(density_air*CabinVolume/2*cp_air/timestep + UA_back(t));
        prev_temp(14)=temperature(14);
        check_b=check_b+1; 
        % Qcv_emitted(t) = Qcv_emitted(t)-(temperature(14)-T_cell(t))*UA_back(t);
        Qcv_emitted(t) = (temperature(14)-T_cell(t))*UA_back(t);
        Qcv_received(t) = 0;
    else
        Qcv_emitted(t) = 0;
        %Qcv_received(t) = Qcv_received(t)+(T_cell(t)-temperature(14))*UA_back(t);
        Qcv_received(t) = (T_cell(t)-temperature(14))*UA_back(t);
    end

    % Qcv_emitted(t) = (temperature(13)-T_cell(t))*UA_front(t);
    % Qcv_received(t) = (T_cell(t)-temperature(13))*UA_front(t);
    % Qcv_emitted(t) = Qcv_emitted(t)+(temperature(14)-T_cell(t))*UA_back(t);
    % Qcv_received(t) = Qcv_received(t)+(T_cell(t)-temperature(14))*UA_back(t);

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
   
        % Thermal Loads
    Qcabin_received(t) = Qhuman(t) + Qequipment(t) + Qengine_av + Qexhaust_av + Qcv_received(t) + Qvent(t) + Qleakage(t);
    Qtotal(t) = Qcabin_received(t) + Qcv_emitted(t);
    
    % Record the simulated time that the MAC is on for heating and cooling
    if Tcabin(:,t) > Ttarget + 1
        cooling_time_simulated(t) = cooling_time_simulated(t) + 1;
    else
        cooling_time_simulated(t) = 0;
    end
    if Tcabin(:,t) < Ttarget + 1
        heating_time_simulated(t) = heating_time_simulated(t) + 1;
    else
        cooling_time_simulated(t) = 0;
    end
    active_time_simulated(t) = cooling_time_simulated(t) + heating_time_simulated(t);
    
    % Set a fixed target time
    ttarget(t) = 600; % Fhaya 2014: tcons= ttarget/abs(log(Ttarget-T0cabin));

    % Calculaton of target MAC load needed to reach Ttarget from Tcabin
    Qtarget(t) = (CabinVolume*density_air*cp_air + mass_seats*cp_seats + 144240)*(Ttarget-Tcabin(t))/(ttarget(t));

    % Definition of MAC work (evaporator) and compressor work and energy 
    Qmac(t) = Qcabin_received(t) + Qcv_emitted(t) + Qtarget(t);
    Qcompressor(t) = Qmac(t)/COP;
    comp_cum_Wh(t) = cumtrapz(timestep.*Qcompressor(t)./3600); % comp_cum_Wh = cumtrapz(comp_Wh); %

    if abs(Qcompressor(t))>compressor_on_lower_limit
        active_time(t) = 1;
    end
    compressor_time(t) = ttarget(t); % OR sum(t_active)
    
    CO2(t) = abs(Qcompressor(t))*compressor_time(t)*Vpe/1000/3600*CF/distance_driven;
end

% Compressor thermal work and actie timme that Qcomp > 0.5 kW
Qcompressor=abs(Qcompressor);
tot_t_acitve=sum(active_time);
total_active_time=tot_t_acitve(end);

Qcabin_received_plot=Qcabin_received;
Qmac_plot=Qmac;
Qemitted=Qleakage+Qcv_emitted;

    % Delete NaN values
Qcabin_received=(Qcabin_received(~isnan(Qcabin_received)));
Qmac=(Qmac(~isnan(Qmac)));
% Qcompressor=Qcompressor(2:end);
% Qcv_emitted=Qcv_emitted(2:end);
% Qcabin_received=Qcabin_received(2:end);

    % Power and Energy
cabin_in_W_av = sum(Qcabin_received)/Total_time;
cabin_out_W_av = sum(Qcv_emitted)/Total_time;
cabin_W_av = cabin_in_W_av - abs(cabin_out_W_av);
target_av_W = sum(Qtarget)/Total_time;
evaporator_av_W = sum(Qmac)/Total_time; % MAC needed=evaporator
compressor_av_W = sum(abs(Qcompressor))/Total_time;

cabin_received_Wh = timestep.*Qcabin_received(2:end)./3600;
cabin_received_cum_Wh = cumsum(cabin_received_Wh);
if ~isempty(cabin_received_cum_Wh)
    cabin_received_Wh_WLTP = cabin_received_cum_Wh(end);
else
    cabin_received_Wh_WLTP = 0; % Or another default value
end
mac_needed_Wh = timestep.*Qmac(2:end)./3600; % MAC needed=evaporator
mac_needed_cum_Wh = cumsum(mac_needed_Wh);
if ~isempty(mac_needed_cum_Wh)
    mac_needed_Wh_WLTP = mac_needed_cum_Wh(end);
else
    mac_needed_Wh_WLTP = 0; % Or another default value
end
comp_Wh = timestep.*Qcompressor(2:end)./3600;
comp_cum_Wh = cumtrapz(comp_Wh); %cumsum
if ~isempty(comp_cum_Wh)
    comp_Wh_WLTP = comp_cum_Wh(end);
else
    comp_Wh_WLTP = 0; % Or another default value
end

    % Final CO2 Emissions
CO2_total = CO2(end);
CO2_perc = CO2_total/95*100; %for 20-24 fleet targets
end