
clc
clear all 
%% INPUTS
%Ambient

T_ini = 273 + 25.52;
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
CabinVolume = cabin_w*(cabin_r_l + cabin_b_l)*cabin_h/2;

% CabinVolume = 1.45*1.3*1;


density_air = 1.18; %Kg/m3 (25ºC)
cp_air = 1006; %J/kg*K  "

%Humans

N_Humans = 0;

% CARGAMOS FICHERO INPUTS CON LA IRRADIANCIA DEL DIA 6/21/2013
load('Inputs.mat');
load('T_Air_mod_Ref.mat');
% ASIGNAMOS VARIABLES


% Determinacion de la duracion del ciclo e interpolacion de las variabels de entrada

% Variable time
Total_time = Irradiancia(end,1);
timestep =  1; %s
time = 0:timestep:Total_time;

% Interpolacion para el vector tiempo
Irr = interp1(Irradiancia(:,1),Irradiancia(:,2),time);
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

T_air(t) = T_ini;

Q_windows(t) = 0;

t = 1 + timestep;
while (t <= Total_time)
%% CALCULO DE COEFICIENTES DE CONVECCION

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
%% GENERAL EQUATION
% m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 

T_air(t) = (Q_windows(t)+Q_ceiling(t)+Q_base(t)+Q_human(t))*(timestep/(CabinVolume*cp_air*density_air))  + T_air(t-1);

t = t +timestep;
end


%% PLOTER
% Temp Air
plot(time(2:end),T_air - 273,'b',T_Air_Exp(:,1),T_Air_Exp(:,2),'r',T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2),'g');

