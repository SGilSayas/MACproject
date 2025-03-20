%% June 2021
clc
clear all 
close all
%% INPUTS
% CARGAMOS FICHERO INPUTS CON LA IRRADIANCIA DEL DIA DEL TEST

load('Inputs06212013.mat');
% load('Inputs05072013.mat');
veh_mov = false;


%Ambient
T_ini = 273 + T_Air_Exp(1,2);

%Radiation Properties
epsilon_ws=0.9; %Windshield
ro_ws=0.246;
tao_ws=0.452;

epsilon_rw=0.9; %Rear Window
ro_rw=0.1;
tao_rw=0.311;

epsilon_lsw=0.9; %Side windows
epsilon_rsw=0.9;
ro_lsw=0.2;
tao_lsw=0.475;


epsilon_roof=0.9; %Roof
ro_roof=0.74;
tao_roof=0;

%Conduction Properties
A_ws=0.63*1.3; %Windshield
k_ws=1.4;
e_ws=0.006;

A_rw=0.29*1; %Rear Window
k_rw=1.4;
e_rw=0.005;

A_lsw=1.45*0.29; %Side Windows
A_rsw=1.45*0.29;
k_lsw=1.4;
k_rsw=1.4;
e_lsw=0.003;
e_rsw=0.003;

A_roof=1.8*1.1; %Roof and ceiling ONLY STEEL.
k_ceiling=14.9;
e_ceiling=0.0005;

%Lower base properties
A_base=6;
ro_base=0.3;
alfa_base=0.7;
C_base=144240;

%CONVECTION INPUTS

%Introducir correlaciones---
%CONVECTION INPUTS


h_ws=5; %Aquí supongo todas las h interiores iguales a 5
h_lsw=5;
h_rsw=5;
h_rw=5;
h_ceiling=5;
h_roof=20;
h_ext=20; %Aquí supongo, de momento, una h exterior contsante de 20
h_base = 5;
%RADIATION INPUTS
sigma=0.0000000567; %Constante de stefan boltzmann
%Geometry aprox

cabin_h =  0.5456;
cabin_w = 1.1;
cabin_r_l = 1.8;
cabin_b_l = 1.45;
CabinVolume = 2*cabin_w*(cabin_r_l + cabin_b_l)*cabin_h/2;

% CabinVolume = 1.45*1.3*1;


density_air = 1.18; %Kg/m3 (25ºC)
cp_air = 1006; %J/kg*K  "

%Humans

N_Humans = 0;


% ASIGNAMOS VARIABLES


% Determinacion de la duracion del ciclo e interpolacion de las variabels de entrada

% Variable time
Total_time = Irradiancia(end,1);
timestep =  1; %s
time = 0:timestep:Total_time;

% Interpolacion para el vector tiempo
Irr = interp1(Irradiancia(:,1),Irradiancia(:,2),time);
% Vel = interp1(Veh_speed(:,1),Veh_speed(:,2),time);
Vel = time *0 + 5;
T_ext = 273 + interp1(Temp_Ext(:,1),Temp_Ext(:,2),time);

%Radiación incidente en todas las ventanas (vector igual para todas, sacado
%del artículo, del point plotter de una gráfica suya, Global Horizontal
%Irradiance)
G_ws_inc = Irr; % Irradiancia es la misma para todas las ventanas (hipotesis)
G_lsw_inc = Irr;
G_rsw_inc = Irr;
G_rw_inc= Irr;
G_roof_inc= Irr;


% INICIALIZACION (en t=1 que valen los vectores? T_air(t=1)= T_ext(1,1))
T_sky = 293; % 0ºC // = T_amb_ext
alfa_ws = 1 - (ro_ws + tao_ws);
alfa_lsw = 1 - ro_lsw - tao_lsw;
alfa_rsw = alfa_lsw;
alfa_rw = 1 - ro_rw - tao_rw;
alfa_roof = 1 - ro_roof - tao_roof;
t = 1;
G_ws_r(t)=0;
G_ws_a(t)=0;
G_ws_t(t)=0;
T_ws(t) = T_ini;
T_lsw(t) = T_ini;
T_rsw(t) = T_ini;
T_rw(t) = T_ini;
T_ceiling(t)= T_ini;
T_base(t)  = T_ini;

h_ext_v(t) = h_ext; % Convection coeficient vector
h_in_v(t) = h_ws;
h_in_h(t) = h_ceiling;

T_air(t) = T_ini;
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

Delta_T = 5; % Delta T between the heat exchangers (Evap. + Cond.) and the sources (Cold and Hot)

% Refrigerant fluids
% R123a 1    
% Water 2

Fluid_Number = 1;

fluids = {'R134a','Water'}; 
% Tamb => [18 29 43] ºC  //  Pevap = [172 310 345] kPa
% Tcond => [37 56 74] ºC // Pcond = [931 1550 2344] kPa   

Low_press_max = 345e3; % For R134a
Comp_ratio_min = 5.7357; % " Average"

% Mode True =  SUMMER (Cooling the cabin) /// False = WINTER (Heating the
% cabin)
Mode = true;

Temp_conf = 23 + 273;

%% Initialization HVAC
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
%% CALCULO DE COEFICIENTES DE CONVECCION
Pr_air=0.7;
Kf_air=0.026;
vis_air = 15.52e-6;
L_roof=1.45;

%% Exterior convection coeficient
% % Roof
% Re_roof_ex = Re(Vel(t),vis_air,L_roof);
% if (Re_roof_ex < 5e5) % Flujo Laminar
%   
%     Nu_roof_ex = 0.664*Re_roof_ex^(0.5)*Pr_air^(1.5);
%     h_roof_ex=(Nu_roof_ex*Kf_air)/L_roof;
% else                  % Flujo Turbulento
%   Nu_roof_ex = (0.037*Re_roof_ex^(4/5) - 871)*Pr_air^(1/3);
%   h_roof_ex=(Nu_roof_ex*Kf_air)/L_roof;
% end
% 
% if(h_roof_ex<5)
%     h_roof_ex = 20;
% end
% h_ext_v(t) = h_roof_ex;
% 
% h_roof = h_roof_ex;
% h_ext = h_roof_ex;

% %% Interior convection coeficient
% %ROOF
% % L_roof=1.45;
% % Gr_roof=(9.81*0.0034*(T_roof(t-1)-T_ext(t-1))*(L_roof^3))/((0.000016^2));
% % Nusselt_roof=0.15*((Gr_roof*Pr_air)^0.333); %Consideramos flujo turbulento
% % h_roof=(Nusselt_roof*Kf_air)/L_roof;
% % h_ext=(Nusselt_roof*Kf_air)/L_roof;
% %CEILING
% L_ceiling=1.45;
% Gr_ceiling=(9.81*0.0034*(T_ceiling(t-1)-T_air(t-1))*(L_ceiling^3))/((0.000016^2));
% Nusselt_ceiling=0.15*((abs(Gr_ceiling)*Pr_air)^0.333);
% h_ceiling=(Nusselt_ceiling*Kf_air)/L_ceiling;
% Gr_ceiling_vec(t) = Gr_ceiling;
% h_in_h(t) = h_ceiling;
% if(h_ceiling<5)
% h_ceiling = 5;
% end
% %WINDSHIELD
% L_ws=1.3;
% Gr_ws=(9.81*cosd(60)*0.0034*(T_ws(t-1)-T_air(t-1))*(L_ws^3))/((0.000016^2));
% Nusselt_ws=(0.825+((0.837*((abs(Gr_ws*Pr_air))^0.1667))/((1+((0.492*Pr_air)^0.5625))^0.4706))^2);
% h_ws=(Nusselt_ws*Kf_air)/L_ws;
% h_in_v(t) = h_ws;
% if(h_ws<5)
% h_ws = 5;
% end
% %REAR WINDOW
% L_rw=1;
% Gr_rw=(9.81*cosd(55)*0.0034*(T_rw(t-1)-T_air(t-1))*(L_rw^3))/((0.000016^2));
% Nusselt_rw=(0.825+((0.837*(abs(Gr_rw*Pr_air)^0.1667))/((1+((0.492*Pr_air)^0.5625))^0.4706))^2);
% h_rw=(Nusselt_rw*Kf_air)/L_rw;
% if(h_rw<5)
% h_rw = 5;
% end
% %LEFT SIDE WINDOW
% L_lsw=1.45;
% Gr_lsw=(9.81*cosd(20)*0.0034*(T_lsw(t-1)-T_air(t-1))*(L_lsw^3))/((0.000016^2));
% Nusselt_lsw=(0.825+((0.837*(abs((Gr_lsw*Pr_air))^0.1667))/((1+((0.492*Pr_air)^0.5625))^0.4706))^2);
% h_lsw=(Nusselt_lsw*Kf_air)/L_lsw;
% if(h_lsw<5)
% h_lsw = 5;
% end
% %RIGHT SIDE WINDOW
% L_rsw=1.45;
% Gr_rsw=(9.81*cosd(20)*0.0034*(T_rsw(t-1)-T_air(t-1))*(L_rsw^3))/((0.000016^2));
% Nusselt_rsw=(0.825+((0.837*(abs((Gr_rsw*Pr_air))^0.1667))/((1+((0.492*Pr_air)^0.5625))^0.4706))^2);
% h_rsw=(Nusselt_rsw*Kf_air)/L_rsw;
% if(h_rsw<5)
% h_rsw = 5;
% end
% %BASE
% h_base=h_ceiling;

% T_Sky set to Ambient Temperature
T_sky = T_ext(t);
%% WINDSHIELD

% Q_cd(t)=Q_ws(t); %2nd Eq

G_ws_r(t)=ro_ws*G_ws_inc(t);
G_ws_a(t)=alfa_ws*G_ws_inc(t);
G_ws_t(t)=tao_ws*G_ws_inc(t);

G_ws_inc(t)=G_ws_r(t)+G_ws_a(t)+G_ws_t(t);

% Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
% Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));

% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_ws*G_ws_a(t); %1er Eq

% 1er Eq (Sustitucion)
% h_ext*A_ws*(T_ws_ext(t)-T_ext(t))+epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4) = (A_ws*k_ws/e_ws)*(T_ws_ext(t) - T_ws(t)) +  A_ws*G_ws_a(t);

% 2nd Eq 
% T_ws_ext(t) = (h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t);

% Iteracion
 It_max = 1000;
 err_prev = 0;
 err = 0;
 T_ws(t) = T_ws(t-timestep);
 rest  = false;
 DT = 0.01;
 for i = 0:0.02:It_max
     
     Ter1 = h_ext*A_ws*(((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))-T_ext(t))+epsilon_ws*sigma*A_ws*(((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))^4-(T_sky)^4);
     Ter2 = (A_ws*k_ws/e_ws)*(T_ws(t) - ((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))) +  A_ws*G_ws_a(t);
     err = abs(Ter2 - Ter1);
     
     if (err < err_prev)&&(rest == true)
         T_ws(t) = T_ws(t) - DT;
         rest = true; 
     elseif (err < err_prev)&&(rest == false)
         T_ws(t) = T_ws(t) + DT;
         rest = false;
     elseif (err >= err_prev)&&(rest == true)
         T_ws(t) = T_ws(t) + DT;
         rest = false; 
     elseif (err >= err_prev)&&(rest == false)
         T_ws(t) = T_ws(t) - DT;
         rest = true;
     end
     err_prev = err;
     if err < 0.01
         i = It_max;
     end
 end
 
%Calcula T_ws_ext(t)
T_ws_ext(t) = (h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t);
% Comprobacion del balance
% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_ws*G_ws_a(t); %1er Eq
% Q_cd(t)=Q_ws(t); %2nd Eq
Q_ws(t) = h_ws*A_ws*(T_ws(t)- T_air(t-1));
Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));

Conv_Rad_ws(t) = Q_cv_ext(t) + Q_rd_ext(t);
Cond_RadAbs_ws(t) = - Q_cd(t) + A_ws*G_ws_a(t);
Cond_ws(t) = Q_cd(t);
Convin_ws(t) = Q_cd(t);
Error_WS(t) = err;
%% LEFT SIDE WINDOW

G_lsw_r(t)=ro_lsw*G_lsw_inc(t);
G_lsw_a(t)=alfa_lsw*G_lsw_inc(t);
G_lsw_t(t)=tao_lsw*G_lsw_inc(t);

G_lsw_inc(t)=G_lsw_r(t)+G_lsw_a(t)+G_lsw_t(t);

% Iteracion
 It_max = 1000;
 err_prev = 0;
 err = 0;
 T_lsw(t) = T_lsw(t-timestep);
 rest  = false;
 DT = 0.01;
 for i = 0:0.02:It_max
     
     Ter1 = h_ext*A_lsw*(((h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t))-T_ext(t))+epsilon_lsw*sigma*A_lsw*(((h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t))^4-(T_sky)^4);
     Ter2 = (A_lsw*k_lsw/e_lsw)*(T_lsw(t) - ((h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t))) +  A_lsw*G_lsw_a(t);
     err = abs(Ter2 - Ter1);
     
     if (err < err_prev)&&(rest == true)
         T_lsw(t) = T_lsw(t) - DT;
         rest = true; 
     elseif (err < err_prev)&&(rest == false)
         T_lsw(t) = T_lsw(t) + DT;
         rest = false;
     elseif (err >= err_prev)&&(rest == true)
         T_lsw(t) = T_lsw(t) + DT;
         rest = false; 
     elseif (err >= err_prev)&&(rest == false)
         T_lsw(t) = T_lsw(t) - DT;
         rest = true;
     end
     err_prev = err;
     if err < 0.01
         i = It_max;
     end
 end
 
%Calcula T_lsw_ext(t)
 T_lsw_ext(t) = (h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t);
% Comprobacion del balance
% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_lsw*G_lsw_a(t); %1er Eq
% Q_cd(t)=Q_lsw(t); %2nd Eq
Q_lsw(t) = h_lsw*A_lsw*(T_lsw(t)- T_air(t-1));
Q_cv_ext_lsw(t)=h_ext*A_lsw*(T_lsw_ext(t)-T_ext(t)); %Exterior convection
Q_rd_ext_lsw(t)=epsilon_lsw*sigma*A_lsw*((T_lsw_ext(t))^4-(T_sky)^4); %Exterior radiation
Q_cd_lsw(t)=((k_lsw*A_lsw)/e_lsw)*(T_lsw_ext(t)-T_lsw(t));

Conv_Rad_lsw(t) = Q_cv_ext_lsw(t) + Q_rd_ext_lsw(t);
Cond_RadAbs_lsw(t) = - Q_cd_lsw(t) + A_lsw*G_lsw_a(t);
Cond_lsw(t) = Q_cd_lsw(t);
Convin_lsw(t) = Q_cd_lsw(t);

Error_LSW(t) = err;
%% RIGHT SIDE WINDOW

G_rsw_r(t)=ro_lsw*G_rsw_inc(t);
G_rsw_a(t)=alfa_rsw*G_rsw_inc(t);
G_rsw_t(t)=tao_lsw*G_rsw_inc(t);

G_rsw_inc(t)=G_rsw_r(t)+G_rsw_a(t)+G_rsw_t(t);

% Iteracion
 It_max = 1000;
 err_prev = 0;
 err = 0;
 T_rsw(t) = T_rsw(t-timestep);
 rest  = false;
 DT = 0.01;
 for i = 0:0.02:It_max
     
     Ter1 = h_ext*A_rsw*(((h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t))-T_ext(t))+epsilon_rsw*sigma*A_rsw*(((h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t))^4-(T_sky)^4);
     Ter2 = (A_rsw*k_rsw/e_rsw)*(T_rsw(t) - ((h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t))) +  A_rsw*G_rsw_a(t);
     err = abs(Ter2 - Ter1);
     
     if (err < err_prev)&&(rest == true)
         T_rsw(t) = T_rsw(t) - DT;
         rest = true; 
     elseif (err < err_prev)&&(rest == false)
         T_rsw(t) = T_rsw(t) + DT;
         rest = false;
     elseif (err >= err_prev)&&(rest == true)
         T_rsw(t) = T_rsw(t) + DT;
         rest = false; 
     elseif (err >= err_prev)&&(rest == false)
         T_rsw(t) = T_rsw(t) - DT;
         rest = true;
     end
     err_prev = err;
     if err < 0.01
         i = It_max;
     end
 end
 
%Calcula T_rsw_ext(t)
 T_rsw_ext(t) = (h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t);
% Comprobacion del balance
% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_rsw*G_rsw_a(t); %1er Eq
% Q_cd(t)=Q_rsw(t); %2nd Eq
Q_rsw(t) = h_rsw*A_rsw*(T_rsw(t)- T_air(t-1));
Q_cv_ext_rsw(t)=h_ext*A_rsw*(T_rsw_ext(t)-T_ext(t)); %Exterior convection
Q_rd_ext_rsw(t)=epsilon_rsw*sigma*A_rsw*((T_rsw_ext(t))^4-(T_sky)^4); %Exterior radiation
Q_cd_rsw(t)=((k_rsw*A_rsw)/e_rsw)*(T_rsw_ext(t)-T_rsw(t));

Conv_Rad_rsw(t) = Q_cv_ext_rsw(t) + Q_rd_ext_rsw(t);
Cond_RadAbs_rsw(t) = - Q_cd_rsw(t) + A_rsw*G_rsw_a(t);
Cond_rsw(t) = Q_cd_rsw(t);
Convin_rsw(t) = Q_cd_rsw(t);

Error_RSW(t) = err;

%% REAR WINDOW

G_rw_r(t)=ro_rw*G_rw_inc(t);
G_rw_a(t)=alfa_rw*G_rw_inc(t);
G_rw_t(t)=tao_rw*G_rw_inc(t);

G_rw_inc(t)=G_rw_r(t)+G_rw_a(t)+G_rw_t(t);
% Iteracion
 It_max = 1000;
 err_prev = 0;
 err = 0;
 T_rw(t) = T_rw(t-timestep);
 rest  = false;
 DT = 0.01;
 for i = 0:0.02:It_max
     
     Ter1 = h_ext*A_rw*(((h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t))-T_ext(t))+epsilon_rw*sigma*A_rw*(((h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t))^4-(T_sky)^4);
     Ter2 = (A_rw*k_rw/e_rw)*(T_rw(t) - ((h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t))) +  A_rw*G_rw_a(t);
     err = abs(Ter2 - Ter1);
     
     if (err < err_prev)&&(rest == true)
         T_rw(t) = T_rw(t) - DT;
         rest = true; 
     elseif (err < err_prev)&&(rest == false)
         T_rw(t) = T_rw(t) + DT;
         rest = false;
     elseif (err >= err_prev)&&(rest == true)
         T_rw(t) = T_rw(t) + DT;
         rest = false; 
     elseif (err >= err_prev)&&(rest == false)
         T_rw(t) = T_rw(t) - DT;
         rest = true;
     end
     err_prev = err;
     if err < 0.01
         i = It_max;
     end
 end
 
%Calcula T_rw_ext(t)
 T_rw_ext(t) = (h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t);
% Comprobacion del balance
% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_rw*G_rw_a(t); %1er Eq
% Q_cd(t)=Q_rw(t); %2nd Eq
Q_rw(t) = h_rw*A_rw*(T_rw(t)- T_air(t-1));
Q_cv_ext_rw(t)=h_ext*A_rw*(T_rw_ext(t)-T_ext(t)); %Exterior convection
Q_rd_ext_rw(t)=epsilon_rw*sigma*A_rw*((T_rw_ext(t))^4-(T_sky)^4); %Exterior radiation
Q_cd_rw(t)=((k_rw*A_rw)/e_rw)*(T_rw_ext(t)-T_rw(t));

Conv_Rad_rw(t) = Q_cv_ext_rw(t) + Q_rd_ext_rw(t);
Cond_RadAbs_rw(t) = - Q_cd_rw(t) + A_rw*G_rw_a(t);
Cond_rw(t) = Q_cd_rw(t);
Convin_rw(t) = Q_cd_rw(t);

Error_RW(t) = err;


%% GENERAL WINDOW EQUATION

Q_windows(t)=Q_ws(t)+Q_lsw(t)+Q_rsw(t)+Q_rw(t); %% Listo

%% ROOF AND CEILING

G_roof_a(t)=alfa_roof*G_roof_inc(t);

% Iteracion
 It_max = 1000;
 err_prev = 0;
 err = 0;
 T_ceiling(t) = T_ceiling(t-timestep);
 rest  = false;
 DT = 0.01;
 for i = 0:0.02:It_max
     
     Ter1 = h_roof*A_roof*(((h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t))-T_ext(t))+epsilon_roof*sigma*A_roof*(((h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t))^4-(T_sky)^4);
     Ter2 = (A_roof*k_ceiling/e_ceiling)*(T_ceiling(t) - ((h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t))) +  A_roof*G_roof_a(t);
     err = abs(Ter2 - Ter1);
     
     if (err < err_prev)&&(rest == true)
         T_ceiling(t) = T_ceiling(t) - DT;
         rest = true; 
     elseif (err < err_prev)&&(rest == false)
         T_ceiling(t) = T_ceiling(t) + DT;
         rest = false;
     elseif (err >= err_prev)&&(rest == true)
         T_ceiling(t) = T_ceiling(t) + DT;
         rest = false; 
     elseif (err >= err_prev)&&(rest == false)
         T_ceiling(t) = T_ceiling(t) - DT;
         rest = true;
     end
     err_prev = err;
     if err < 0.01
         i = It_max;
     end
 end
 
%Calcula T_rw_ext(t)
 T_roof(t) = (h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t);
% Comprobacion del balance
% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_roof*G_roof_a(t); %1er Eq
% Q_cd(t)=Q_ceiling(t); %2nd Eq
Q_ceiling(t) = h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1));
Q_cv_ext_roof(t)=h_ext*A_roof*(T_roof(t)-T_ext(t)); %Exterior convection
Q_rd_ext_roof(t)=epsilon_roof*sigma*A_roof*((T_roof(t))^4-(T_sky)^4); %Exterior radiation
Q_cd_ceiling(t)=((k_ceiling*A_roof)/e_ceiling)*(T_roof(t)-T_ceiling(t));

Conv_Rad_ceiling(t) = Q_cv_ext_roof(t) + Q_rd_ext_roof(t);
Cond_RadAbs_ceiling(t)= - Q_cd_ceiling(t) + A_roof*G_roof_a(t);
Cond_ceiling(t) = Q_cd_ceiling(t);
Convin_ceiling(t) = Q_ceiling(t);

Error_CR(t) = err;



%% BASE
% 1a Eq: (C_base*(T_base(t)-T_base(t-1))/(timestep)) = alfa_base*((A_ws*G_ws_t)+(A_rw*G_rw_t)+(A_lsw*G_lsw_t)+(A_lsw*G_rsw_t))-Q_base;

% 2a Eq:   Q_base=h_base*A_base*(T_base(t)-T_air(t-1));

% Sustituyendo 2a en 1a
T_base(t) = (timestep/(C_base+h_base*A_base)) * (alfa_base*((A_ws*G_ws_t(t))+(A_rw*G_rw_t(t))+(A_lsw*G_lsw_t(t))+(A_lsw*G_rsw_t(t))) + C_base*T_base(t-1) + h_base*A_base*T_air(t-1));
Q_base(t)=h_base*A_base*(T_base(t)-T_air(t-1));

%% BEINGS
Q_human(t)=N_Humans * 120;


%% HVAC function
Q_air_cabin(t) = Q_windows(t)+Q_ceiling(t)+Q_base(t)+Q_human(t);

%% HVAC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Summer Mode
% if (Mode == true)
% % Heat exchangers
% T_evap(t) = Temp_conf;
% % T_cond(t) = T_ext(t);
% T_cond(t) = T_ext_const;
% else
% % Winter Mode
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
%%



%% GENERAL EQUATION
% m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 

T_air(t) = (Q_windows(t)+Q_ceiling(t)+Q_base(t)+Q_human(t) - Q_HVAC(t) )*(timestep/(CabinVolume*cp_air*density_air))  + T_air(t-1);



t = t + timestep
end


%% PLOTER CABIN
% Temp Air
% plot(time(2:end),T_air - 273,'b',T_Air_Exp(:,1),T_Air_Exp(:,2),'r',T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2),'g');
plottemplate3p(time(2:end),T_air - 273,T_Air_Exp(:,1),T_Air_Exp(:,2),T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2));
tufteAxesOrig;
plottemplate4p(time(2:end),Q_windows,time(2:end),Q_ceiling,time(2:end),Q_base,time(2:end),Q_human);
tufteAxesOrig;
plottemplate1p(Irradiancia(:,1), Irradiancia(:,2),'Irradiance [W/m2]','Irradiance')
tufteAxesOrig;
% plottemplate1p(Veh_speed(:,1), Veh_speed(:,2),'Speed [m/s]','Vehicle Speed')
% tufteAxesOrig;

%% PLOTER HVAC
plottemplate1p(time(2:end-1), HVAC_Out.COP(2:end),'COP [-]','HP COP');
tufteAxesOrig;
plottemplate1p(time(2:end-1), HVAC_Out.mf(2:end)*1000,'Mass flow [g/s]','Mass flow');
tufteAxesOrig;
plottemplate1p(time(2:end-1), HVAC_Out.W_comp(2:end),'Power [W]','Compressor Power');
tufteAxesOrig;
plottemplate1p(time(2:end-1), HVAC_Out.Q_evap(2:end),'Heat [W]','Evaporator heat');
tufteAxesOrig;
plottemplate1p(time(2:end-1), HVAC_Out.Q_cond(2:end),'Heat [W]','Condenser heat');
tufteAxesOrig;
