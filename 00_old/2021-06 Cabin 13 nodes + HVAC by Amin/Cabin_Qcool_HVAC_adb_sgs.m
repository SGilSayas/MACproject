% Vehicle Cabin Model
% April 2021, One node model. Created by Alejandor Vega and Amin Dreif 
% May 2021, upgrade to lumped model (13 and 14 nodes) by Susana Gil-Sayas
% June 2021, addition HVAC and Cool Prop 
clc
clear all 
close all
%% INPUTS
timestep =  1; %s

%Ambient
T_ini = 273.15 + 25.52; 
T_sky=T_ini-6; %T_sky = 293;

%Radiation Properties 
sigma=0.0000000567; %5.67037321*10^-8 W/(m^2*K^4) %Stefan-Boltzmann constant
% epsilon=emissivity[-], rho=reflectivity[-], tau=transmissivity[-], alpha=absorptivity[-]
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

%Conduction Properties
% A=area[m2], k=thermoconductivity[W/(m*K)], e=thickness[m]
    %Windshield
A_ws=0.63*1.3; 
k_ws=1.4;
e_ws=0.006;
    %Rear Window
A_rw=0.29*1; 
k_rw=1.4;
e_rw=0.005;
    %Side Windows
A_lsw=1.45*0.29; 
A_rsw=1.45*0.29;
k_lsw=1.4;
k_rsw=1.4;
e_lsw=0.003;
e_rsw=0.003;
    %Roof and ceiling ONLY STEEL.
A_roof=1.8*1.1; 
k_ceiling=14.9;
e_ceiling=0.0005;
    %Base
A_base=6;
C_base=144240; %Capacitance, W

%Introducir correlaciones
h_ws=5; 
h_lsw=5;
h_rsw=5;
h_rw=5;
h_ceiling=5;
h_roof=20;
h_ext=20; 
h_base = 5;

%Geometry aprox
cabin_h =  0.5456;
cabin_w = 1.1;
cabin_r_l = 1.8;
cabin_b_l = 1.45; 
V_cabin = cabin_w*(cabin_r_l + cabin_b_l)*cabin_h; % = 1.45*1.3*1;
Ttarget =23+273.15; %K

% Air properies
rho_air = 1.18; %Kg/m3 (25ºC)
cp_air = 1006; %J/kg*K  "

% Load input IRRADIANCE 6/21/2013
t = 1;
load('Inputs.mat');
load('T_Air_mod_Ref.mat');
    % Variable time
Total_time = Irradiancia(end,1);
time = 0:timestep:Total_time;
    % Interpolacion para el vector tiempo
Irr = interp1(Irradiancia(:,1),Irradiancia(:,2),time);
T_ext = 273.15 + interp1(Temp_Ext(:,1),Temp_Ext(:,2),time);
Tamb = T_ext;
% Vel = interp1(Veh_speed(:,1),Veh_speed(:,2),time);
Vel = time *0 + 5;

% Human heat
N_Humans=1;
Qhuman(t)=N_Humans * 120; 

    %Global Horizontal Irradiance (Hipótesis: same all windows)
G_roof_inc=Irr;
G_ws_inc=Irr; 
G_rw_inc=Irr;
G_lsw_inc=Irr;
G_rsw_inc=Irr;

    %Different radiations:
    %ws
G_ws_r=rho_ws*G_ws_inc;
G_ws_a=alfa_ws*G_ws_inc; 
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

% Boundary conditions
Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_lsw*G_lsw_t(t)) + (A_lsw*G_rsw_t(t)) );

% Conductances Matrix 
% Tamb --K1-> Troof --K2-> Tceiling --K3-> Tcabin
K1=h_ext*A_roof; %h_roof por qué existe?
K2=A_roof*k_ceiling/e_ceiling;
K3=h_ceiling*A_roof;
% Tbase --K4-> Tcabin_back
K4=h_base*A_base;
% Tamb --K5-> Tw_ext_ws --K6-> Tw_int_ws --K7-> Tcabin
K5=h_ext*A_ws;
K6=A_ws*k_ws/e_ws;
K7=h_ws*A_ws;
% Tamb --K8-> Tw_ext_rw --K9-> Tw_int_rw --K10-> Tcabin
K8=h_ext*A_rw;
K9=A_rw*k_rw/e_rw;
K10=h_rw*A_rw;
% Tamb --K11-> Tw_ext_lsw --K12-> Tw_int_lsw --K13-> Tcabin
K11=h_ext*A_lsw;
K12=A_lsw*k_lsw/e_lsw;
K13=h_lsw*A_lsw;
% Tamb --K14-> Tw_ext_rsw --K15-> Tw_int_rsw --K16-> Tcabin
K14=h_ext*A_rsw;
K15=A_rsw*k_rsw/e_rsw;
K16=h_rsw*A_rsw;
K=[1 0 0 0 0 0 0 0 0 0 0 0 0;
    K1 -(K1+K2) K2 0 0 0 0 0 0 0 0 0 0;
    0 K2 -(K2+K3) 0 0 0 0 0 0 0 0 0 K3;
    0 0 0 -(K4) 0 0 0 0 0 0 0 0 K4;
    K5 0 0 0 -(K5+K6) 0 0 0 K6 0 0 0 0;
    K8 0 0 0 0 -(K8+K9) 0 0 0 K9 0 0 0;
    K11 0 0 0 0 0 -(K11+K12) 0 0 0 K12 0 0;
    K14 0 0 0 0 0 0 -(K14+K15) 0 0 0 K15 0;
    0 0 0 0 K6 0 0 0 -(K6+K7) 0 0 0 K7;
    0 0 0 0 0 K9 0 0 0 -(K9+K10) 0 0 K10;
    0 0 0 0 0 0 K12 0 0 0 -(K12+K13) 0 K13;
    0 0 0 0 0 0 0 K15 0 0 0 -(K15+K16) K16;
    0 0 K3 K4 0 0 0 0 K7 K10 K13 K16 -(K3+K4+K7+K10+K13+K16)];

% Capacitance Matrix
C_amb=0; %b.c.
C_roof=0; %de momento
C_ceiling=0; %de momento
C_base=144240; %Capacitance, W
C_ws=0;
C_rw=0;
C_lsw=0;
C_rsw=0;
C_air=rho_air*(V_cabin)*cp_air/timestep;
C_cabin=C_air;
C=zeros(13);
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
C(13,13)=C_cabin;

%INI
G_ws_r(t)=0;
G_ws_a(t)=0;
G_ws_t(t)=0;

h_ext_v(t) = h_ext; % Convection coeficient vector
h_in_v(t) = h_ws;
h_in_h(t) = h_ceiling;

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
Tcabin(t)=T_ini;
temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];
delta_temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];

Q_windows(t) = 0;
trip_duration=Total_time;
i=1;

for t=2:Total_time+1

    Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_lsw*G_lsw_t(t)) + (A_lsw*G_rsw_t(t)) );
    Qhuman(t)=N_Humans * 120;

    Qcool_req(t)=rho_air*V_cabin*cp_air*(temperature(13)-Ttarget)/timestep; %trip_duration; %Tcabin=324.7876
    Tbc = [Tamb(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);-A_lsw*G_lsw_a(t);-A_rsw*G_rsw_a(t);0;0;0;0;-Qhuman(t)+Qcool_req(t)]; 
    
    temperature=(inv(K - C))*(Tbc - C*delta_temperature);
    delta_temperature=temperature;
    
    Tamb(:,t)=temperature(1); %T_amb es condición de contorno, Tamb la calculada
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
    Tcabin(:,t)=temperature(13);

Qcool_actual(t)=rho_air*V_cabin*cp_air*(Tcabin(:,t)-Ttarget)/timestep;%trip_duration;
Qcool_actual_cum=Qcool_actual(t-1)+Qcool_actual(t);
end
error=Qcool_actual-Qcool_req; %fila
% %if loop: 
% Qcool_req=Qcool_actual;
% E(i,:)=error;
% TCABIN(i,:)=Tcabin;
% i=i+1; %contador de cada run

T_cabin_final=Tcabin(:,trip_duration); %RUN#1: 52ºC; %RUN#2: 324.8446K=51.6946ºC ; %RUN#3: 324.8446

% figure(1)
% plot(time,Tcabin-273.15,'b',T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2),'k',T_Air_Exp(:,1),T_Air_Exp(:,2),'--')
% legend('T cabin','T cabin model (Paper)','T cabin thermocouple (Paper)','Location','best')
% grid on
% grid minor
% xlabel('Time [s]') 
% ylabel('Temperature [ºC]')

figure(2)
yyaxis left
plot(time,Tcabin-273.15,'b')
hold on
plot(time,Troof-273.15,time,Tamb-273.15)
grid on
grid minor
ylabel('Temperature [ºC]')
xlabel('Time [s]') 
yline(23)
yyaxis right
plot(time,Qcool_req)
hold on
plot(time,Qcool_actual,'k')
ylabel('Cooling heat [W]')
legend('T cabin','T roof','T ambient (b.c.)','T target','Qcool req','Qcool actual','Location','NorthWest')

% figure(3)
% plot(time,error)

%% HVAC
Q_air_cabin=Qcool_actual; %sgs
T_air(t) = Tcabin(t);
Q_windows(t) = 0;

%% HVAC Model
% This model calculates the HVAC system of an xEV vehicle
% Sign criteria ----->   [+] Positive heat == Heat absorbed by the refrigerant
%               ----->   [-] Negative heat == Heat rejected by the refrigerant
%% Set PATH FOR PYTHON and CoolProp Installation
% pcPythonExe = 'C:\Users\susan\AppData\Local\Programs\Python\Python36\python.exe';
% [ver, exec, loaded]	= pyversion(pcPythonExe);
% pyversion;
% 
% % CoolProp Instalation library instalation -python
% [v,e] = pyversion; system([e,' -m pip install --user -U CoolProp'])

%% Boundary conditions
% Load Inputs from cabin model
% load Q_Air_Cabin_T1.mat %adb
% load Inputs06212013.mat %adb

% GENERAL EQUATION
% m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 

N_human = 1;
Q_human = 120 * N_human;
Q_air_cabin = Q_air_cabin + Q_human;
heat = max(Q_air_cabin) ;

% % Variable time
% timestep = 1; % 1s
% Total_time = Irradiancia(end,1);
% time = 0:timestep:Total_time;
% t = 1;
% 
% % Interpolacion para el vector tiempo
% Irr = interp1(Irradiancia(:,1),Irradiancia(:,2),time);
% % Vel = interp1(Veh_speed(:,1),Veh_speed(:,2),time);

% % Sources (Temp_conf && T_ext)
% T_ext = 273 + interp1(Temp_Ext(:,1),Temp_Ext(:,2),time);
% Temp_conf = 23 + 273;
% T_ext_const = mean(T_ext);

T_ext = Tamb;
Temp_conf = Ttarget;
T_ext_const = mean(Tamb);


% % Load Inputs from cabin model
Delta_T = 5; % Delta T between the heat exchangers (Evap. + Cond.) and the sources (Cold and Hot)


% % Refrigerant fluids % R123a 1 % Water 2
Fluid_Number = 1;
fluids = {'R134a','Water'}; 
% Tamb => [18 29 43] ÂºC  //  Pevap = [172 310 345] kPa
% Tcond => [37 56 74] ÂºC // Pcond = [931 1550 2344] kPa   
Low_press_max = 345e3; % For R134a
Comp_ratio_min = 5.7357; % " Average"

% Mode True =  SUMMER (Cooling the cabin) /// False = WINTER (Heating the
% cabin)
Mode = true;
Temp_conf = Ttarget; %23 + 273;

% % Initialization HVAC
Q_evap_req(t) = 0; % Heat exchanged in the evaporator [W]
Q_cond_req(t) = 0; % Heat echanged in the condenser   [W]
W_comp(t) = 0; % Work performed by the compresor [W]
COP(t) = 0; % Coefficient of performance of the HP [-]
mf(t) = 0; % Compressor refrigerant mass flow [kg/s]
eta_is(t) = 0.8; % Compressor efficiency [-]
T1(t) = 0; % Temperature 1 [K]
P1(t) = 0; % Pressure 1 [Pa]
T2(t) = 0; % Temperature 2
P2(t) = 0; % Pressure 2 [Pa]
T3(t) = 0; % Temperature 3
P3(t) = 0; % Pressure 3 [Pa]
T4(t) = 0; % Temperature 4
P4(t) = 0; % Pressure 4 [Pa]
COP(t) = 0; % Coefficient of performance [-]
mf(t) = 0; % Mass flow [Kg/s]
Q_evap(t) = 0; % Heat exchanged in the evaporator [W]
Q_cond(t) = 0; % Heat echanged in the condenser   [W]
Q_HVAC(t) = 0; % Heat exchanged with cabin Air with the HVAC [W]

t = 1 + timestep;
while (t <= Total_time)
    %Different radiations:
    %ws
G_ws_r(t)=rho_ws*G_ws_inc(t);
G_ws_a(t)=alfa_ws*G_ws_inc(t); 
G_ws_t(t)=tao_ws*G_ws_inc(t);
G_ws_inc(t)=G_ws_r(t)+G_ws_a(t)+G_ws_t(t);
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

%    HVAC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_ceiling(t) = h_ceiling*A_roof*(Tceiling(t)- Tcabin(t-1));
% Q_base(t)=h_base*A_base*(T_base(t)-T_air(t-1));

 %HVAC function
Q_air_cabin(t) = Q_windows(t)+Q_ceiling(t)+Qbase(t)+Qhuman(t);
 %Summer Mode
% if (Mode == true)
 %Heat exchangers
% T_evap(t) = Temp_conf;
% % T_cond(t) = T_ext(t);
% T_cond(t) = T_ext_const;
% else
 %Winter Mode
% T_evap(t) = T_ext_const;
% % T_evap(t) = T_ext(t);
% T_cond(t) = Temp_conf;
% end

%% Calculating Loop
% t = 1 + timestep;
% while (t <= Total_time)
% Low and high pressure accordind to temperature - R132a
%
 if (mean(Q_air_cabin)>= 0)
   Mode = true;
   Q_evap_req(t) = abs(Q_air_cabin(t));
%    Q_evap_req(t) = heat;
else
   Mode = false;
   Q_cond_req(t) = abs(Q_air_cabin(t));
%    Q_cond_req(t) = heat;
end

% Summer Mode
if (Mode == true)
% Heat exchangers
T_evap(t) = Temp_conf - Delta_T ;
T_cond(t) = T_ext(t) + Delta_T;
% T_cond(t) = T_ext_const + Delta_T;
else
% Winter Mode
% T_evap(t) = T_ext_const - Delta_T;
T_evap(t) = T_ext(t) - Delta_T ;
T_cond(t) = Temp_conf + Delta_T;
end   
 
% py.CoolProp.py.CoolProp.CoolProp.PropsSI('T','P',101325,'Q',0,char(fluids(1)));

% Code Diagram p-h 
Ref = char(fluids(1)); 

%% Controlando P2
% T1(t)= T_evap(t); 
% P1(t)=py.CoolProp.CoolProp.PropsSI('P','T',T1(t),'Q',1,Ref);      
% h1(t)=py.CoolProp.CoolProp.PropsSI('H','T',T1(t),'Q',1,Ref);     
% s1(t)=py.CoolProp.CoolProp.PropsSI('S','T',T1(t),'Q',1,Ref);       
% P2(t)=10e6;    
% beta(t)=P2(t)/P1(t);    
% s2is(t)=s1(t); 
% h2is(t)=py.CoolProp.CoolProp.PropsSI('H','P',P2(t),'S',s2is(t),Ref);
% eta_is(t)=0.8;
% h2(t)=(h2is(t)-h1(t))/eta_is(t) + h1(t);
% T2(t)=py.CoolProp.CoolProp.PropsSI('T','H',h2(t),'P',P2(t),Ref);
% s2(t)=py.CoolProp.CoolProp.PropsSI('S','H',h2(t),'P',P2(t),Ref);
% T3(t)=T_cond(t);
% P3(t)=P2(t);
% h3(t)=py.CoolProp.CoolProp.PropsSI('H','T',T3(t),'P',P3(t),Ref);
% s3(t)=py.CoolProp.CoolProp.PropsSI('S','T',T3(t),'P',P3(t),Ref);
% s4(t)=s3(t);
% h4(t)=py.CoolProp.CoolProp.PropsSI('H','S',s4(t),'P',P1(t),Ref);
% T4(t)=py.CoolProp.CoolProp.PropsSI('T','S',s4(t),'P',P1(t),Ref);
% Q4(t)=py.CoolProp.CoolProp.PropsSI('Q','S',s4(t),'P',P1(t),Ref);
%% Tomand T1 y T2 como conocidas (Punto 1, vapor saturado) (Punto 3, liquido saturad0)

T1(t)= T_evap(t); 
P1(t)=py.CoolProp.CoolProp.PropsSI('P','T',T1(t),'Q',1,Ref);

if (P1(t) >= Low_press_max)
P1(t) = Low_press_max;
T1(t) = py.CoolProp.CoolProp.PropsSI('T','P',P1(t),'Q',1,Ref);
end
h1(t)=py.CoolProp.CoolProp.PropsSI('H','T',T1(t),'Q',1,Ref);     
s1(t)=py.CoolProp.CoolProp.PropsSI('S','T',T1(t),'Q',1,Ref);       

T2sat(t) = T_cond(t);
P2(t)=py.CoolProp.CoolProp.PropsSI('P','T',T2sat(t),'Q',1,Ref);

if (P2(t) <= Comp_ratio_min*P1(t))
P2(t) = Comp_ratio_min*P1(t);
end
% P2(t)=10e6;    
beta(t)=P2(t)/P1(t);    
s2is(t)=s1(t); 
h2is(t)=py.CoolProp.CoolProp.PropsSI('H','P',P2(t),'S',s2is(t),Ref);
eta_is(t)=0.8; 
h2(t)=(h2is(t)-h1(t))/eta_is(t) + h1(t);
T2(t)=py.CoolProp.CoolProp.PropsSI('T','H',h2(t),'P',P2(t),Ref);
s2(t)=py.CoolProp.CoolProp.PropsSI('S','H',h2(t),'P',P2(t),Ref);
% T3(t)=T_cond(t);
P3(t)=P2(t); %% Liquido saturado (Punto 3)
T3(t)=py.CoolProp.CoolProp.PropsSI('T','P',P3(t),'Q',0,Ref);
h3(t)=py.CoolProp.CoolProp.PropsSI('H','P',P3(t),'Q',0,Ref);
s3(t)=py.CoolProp.CoolProp.PropsSI('S','P',P3(t),'Q',0,Ref);
s4(t)=s3(t);
h4(t)=py.CoolProp.CoolProp.PropsSI('H','S',s4(t),'P',P1(t),Ref);
T4(t)=py.CoolProp.CoolProp.PropsSI('T','S',s4(t),'P',P1(t),Ref);
Q4(t)=py.CoolProp.CoolProp.PropsSI('Q','S',s4(t),'P',P1(t),Ref);
P4(t)= P1(t);
%%
COP(t)=(h2(t)-h3(t))/(h2(t)-h1(t));

if (Mode == true)
    mf(t) = Q_evap_req(t)/(h1(t) - h4(t));
    ratio = mf(t)/mf(t-1);
    if (ratio > 1.05) && (ratio~= inf)
     mf(t) = 1.05*mf(t-1);
    elseif (ratio < 0.95)
     mf(t) = 0.95*mf(t-1);
    end   
else
    mf(t) = Q_cond_req(t)/(h2(t) - h3(t));  
     ratio = (mf(t)- mf(t-1))/mf(t-1);
    if (ratio > 1.05) && (ratio~= inf)
     mf(t) = 1.05*mf(t-1);
    elseif (ratio < 0.95)
     mf(t) = 0.95*mf(t-1);
    end   
end
Q_evap(t) = mf(t)*(h1(t) - h4(t));
Q_cond(t) = mf(t)*(h3(t) - h2(t));
W_comp(t) = mf(t)*(h2(t) - h1(t));

if (Mode == true)
    Q_HVAC(t) = Q_evap(t);
else
    Q_HVAC(t) = Q_cond(t);  
end

Names = {'P1';'T1';'T2';'P2';'T3';'P3';'T4';'P4';'COP';'mf';'W_comp';'Q_evap';'Q_cond'};
HVAC_Out = table(P1',T1',T2',P2',T3',P3',T4',P4',COP',mf',W_comp',Q_evap',Q_cond','VariableNames',Names);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GENERAL EQUATION
% m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 
Tcabin(t) = (Q_windows(t)+Q_ceiling(t)+Qbase(t)+Qhuman(t) - Q_HVAC(t) )*(timestep/(V_cabin*cp_air*rho_air)) + Tcabin(t-1);

t = t + timestep;
end

%PLOTER HVAC
plottemplate1p(time(2:end-1), HVAC_Out.COP(2:end-1),'COP [-]','HP COP');
tufteAxesOrig;
plottemplate1p(time(2:end-1), HVAC_Out.mf(2:end-1)*1000,'Mass flow [g/s]','Mass flow');
tufteAxesOrig;
plottemplate1p(time(2:end-1), HVAC_Out.W_comp(2:end-1),'Power [W]','Compressor Power');
tufteAxesOrig;
plottemplate1p(time(2:end-1), HVAC_Out.Q_evap(2:end-1),'Heat [W]','Evaporator heat');
tufteAxesOrig;
plottemplate1p(time(2:end-1), HVAC_Out.Q_cond(2:end-1),'Heat [W]','Condenser heat');
tufteAxesOrig;
plottemplate1p(time(2:end-1), Q_air_cabin(2:end-1),'Heat [W]','Cabin thermal load');
tufteAxesOrig;
