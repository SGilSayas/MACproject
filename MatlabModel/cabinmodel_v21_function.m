function [Qcompressor_W_av,total_energy_cum_compr_needed_Wh,CO2_total] = cabinmodelfunction(Total_time,distance_driven,N_Humans,COP,TC_cell_constant,Irr_constant,TC_cell,Irr)
%[a,b,c] = cabinmodelfunction(1800,'N/A',2,'N/A',35,500,'N/A','N/A')
   % Constant data
t = 1; timestep = 1; % s
heat_human = 80; %W
heat_equipment = 0 ; %40; % W
Vpe_petrol = 0.264; % l/kWh
CF_petrol = 2330; % gCO2/l

    % Inputs stored in CELL ARRAY:
% Total_time = 2500; % s
% distance_driven = 23.25; % km
% N_Humans = 1;
% COP = 3.5;
% TC_cell_constant = 35; % C
% Irr_constant = -0; % IN NEGATIVE (enters the cabin), W/m2

    % Filling the inpiut blanks with default values
if Total_time == 'N/A' Total_time = 1800; end
if distance_driven == 'N/A' distance_driven = 23.25; end
if N_Humans == 'N/A' N_Humans = 1; end
if COP == 'N/A' COP = 3.5; end
if TC_cell == 'N/A' TC_cell = ones(1,Total_time)*TC_cell_constant; end
if Irr == 'N/A' Irr = -ones(1,Total_time)*Irr_constant; end
   
    % Load input cell humidity (RH) and pressure (P)
load('G7_35_off_RH_P_amb')
RH_amb = ones(1,Total_time)*mean(RH_amb);
P_amb_kPa = ones(1,Total_time)*mean(P_amb_kPa);

    % Load input Engine and exhaust heat flows
load('heatflow_engine_exhaust_W.mat')

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

    % Geometry (m) and areas (m2)
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

    % Humidity ratio in gram of water per gram of dry air, X
Ps_22_kPa = 2.626; % kPa, Water saturation pressure at 22C
Ps_35_kPa = 5.63; % kPa, Water saturation pressure at 35C
Ps_kPa = Ps_22_kPa+(T_ini-22)/(35-22)*(Ps_35_kPa-Ps_22_kPa);
meanP_kPa_amb = mean(P_amb_kPa); % in kPa
X = 0.62198*Ps_kPa.*RH_amb./(100.*P_amb_kPa - Ps_kPa.*RH_amb); % Faya, 2013

    % Enthalpies calculation. Humidity considered same in amb and cabin
e_amb = 1006.*T_cell + (2501000 + 1770.*T_cell).*X; % J/kg

    % Air vent volume and flow rate:
vent_volumerate = 0.0003; % m3/s, Faya 2013: 0.02; 0.001
leakage_volumerate = 0.0001; %% m3/s, Faya2014

    % Initial values
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

compressor_on_lower_limit = 0.5; % W, 

for t=2:Total_time

    % Cabin air enthalpy
    e_cabin(t) = 1006*Tcabin(t) + (2501000 + 1770*Tcabin(t))*X(t); % J/kg, Faya 2013

    % Heat Flows
    Qleakage(t) = - leakage_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    Qvent(t) = vent_volumerate*density_air*(e_amb(t)-e_cabin(t)); % W, J/s
    heat_human(t) = h_cabin*A_skin*abs(temperature(13)-T_skin) ;
    Qhuman(t) = N_Humans.*heat_human(t);
    Qequipment(t) = heat_equipment;
    
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

    if ( temperature(13)>T_cell(t) )
        temperature(13) = ( density_air*CabinVolume/2*cp_air*delta_temperature(13)/timestep + UA_front*T_cell(t) )/...
            ( density_air*CabinVolume/2*cp_air/timestep + UA_front );
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
        Qcv_emitted(t) = -(temperature(14)-T_cell(t))*UA_back;
        Qcv_received(t) = 0;
    else
        Qcv_emitted(t) = 0;
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

    % Energies
    % Cabin thermal load
energy_cabin_Wh = timestep.*Qcabin_received(2:end)./3600;
energy_cum_cabin_Wh = cumsum(energy_cabin_Wh);
total_energy_cum_cabin_Wh = energy_cum_cabin_Wh(end);
    % MAC needed cooling load
energy_mac_needed_Wh = timestep.*Qmac_needed(2:end)./3600;
energy_cum_mac_needed_Wh = cumsum(energy_mac_needed_Wh);
total_energy_cum_mac_needed_Wh = energy_cum_mac_needed_Wh(end);
    % Compressor
%%tot_energy_compressor = total_energy_cum_mac_needed_Wh/COP
energy_comp_needed_Wh = timestep.*Qcompressor(2:end)./3600;
energy_cum_comp_needed_Wh = cumtrapz(energy_comp_needed_Wh); %cumsum
total_energy_cum_compr_needed_Wh = energy_cum_comp_needed_Wh(end)

    % Powers
cabin_in_W_av = sum(Qcabin_received)/Total_time
cabin_out_W_av = sum(Qcv_emitted)/Total_time
Qtarget_W_av = sum(Qtarget)/Total_time;
Qcabin_W_av = cabin_in_W_av - abs(cabin_out_W_av)
% Qmac_needed_W_av = sum(abs(Qmac_needed))/Total_time %cooling load
Qcompressor_W_av = sum(abs(Qcompressor))/Total_time

    % CO2 Emissions
CO2_total = CO2(end);
CO2_perc = CO2_total/95*100;

end