% Vehicle's Cabin Model
% April 2021, One node model. Created by Alejandor Vega and Amin Dreif 
% May 2021, upgrade to lumped model (13 and 14 nodes) by Susana Gil-Sayas
clc
clear all 
close all
%% INPUTS
timestep =  1; %s

%Ambient
T_ini = 273.15 + 25.52; 
T_sky=T_ini-6; %T_sky = 293;

%Radiation Properties 
sigma=0.0000000567037321; % W/(m^2*K^4) %Stefan-Boltzmann constant(sigma)
% epsilon:emissivity[], rho=reflectivity[], tau=transmissivity[], alpha=absorptivity[]
    %Windshield
epsilon_ws=0.9; 
ro_ws=0.246;
tao_ws=0.452;
alfa_ws = 1-(ro_ws + tao_ws);
    %Rear Window
epsilon_rw=0.9; 
ro_rw=0.1;
tao_rw=0.311;
alfa_rw = 1-ro_rw-tao_rw;
    %Side windows, left and right
epsilon_lsw=0.9; 
epsilon_rsw=0.9;
ro_lsw=0.2;
tao_lsw=0.475;
alfa_lsw = 1-ro_lsw - tao_lsw;
ro_rsw=ro_lsw;
tao_rsw=tao_lsw;
alfa_rsw = alfa_lsw;
    %Roof
epsilon_roof=0.9; 
ro_roof=0.74;
tao_roof=0;
alfa_roof = 1-ro_roof-tao_roof;
    %Base
ro_base=0.3;
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

%Introducir correlaciones---
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
CabinVolume = cabin_w*(cabin_r_l + cabin_b_l)*cabin_h; % = 1.45*1.3*1;

% Air properies
density_air = 1.18; %Kg/m3 (25ÂºC)
cp_air = 1006; %J/kg*K  "

% Load input IRRADIANCE 6/21/2013
t = 1;
load('Inputs.mat');
load('T_Air_mod_Ref.mat');
    % Variable time
Total_time = Irradiancia(end,1);
time = 0:timestep:Total_time;
    % Interpolacion para el vector tiempo
% Irr = interp1(Irradiancia(:,1),Irradiancia(:,2),time);
% T_ext = 273.15 + interp1(Temp_Ext(:,1),Temp_Ext(:,2),time);
T_amb = 50 + 273.15;

% Human heat
N_Humans=1;
Qhuman(t)=N_Humans * 120; 

    %Global Horizontal Irradiance (Hipótesis: same all windows)
Irr=0;
G_roof_inc=Irr;
G_ws_inc=Irr; 
G_rw_inc=Irr;
G_lsw_inc=Irr;
G_rsw_inc=Irr;

    %Different radiations:
    %ws
G_ws_r=ro_ws*G_ws_inc;
G_ws_a=alfa_ws*G_ws_inc; 
G_ws_t=tao_ws*G_ws_inc;
    %rw
G_rw_r=ro_rw*G_rw_inc;
G_rw_a=alfa_rw*G_rw_inc; 
G_rw_t=tao_rw*G_rw_inc;
    %lsw
G_lsw_r=ro_lsw*G_lsw_inc;
G_lsw_a=alfa_lsw*G_lsw_inc; 
G_lsw_t=tao_lsw*G_lsw_inc;
    %rsw
G_rsw_r=ro_rsw*G_rsw_inc;
G_rsw_a=alfa_rsw*G_rsw_inc; 
G_rsw_t=tao_rsw*G_rsw_inc;
    %roof
G_roof_a=alfa_roof*G_roof_inc;

% Boundary conditions vector
Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_lsw*G_lsw_t(t)) + (A_lsw*G_rsw_t(t)) );
Tbc=[T_amb;-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);-A_lsw*G_lsw_a(t);-A_rsw*G_rsw_a(t);0;0;0;0;-Qhuman(t)]; 

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

%14 nodes: Tcabin_front --K17-> Tcabin_back
% K17=h_ws*A_ws; %TBD
% K=[1 0 0 0 0 0 0 0 0 0 0 0 0 0;
%     K1 -(K1+K2) K2 0 0 0 0 0 0 0 0 0 0 0;
%     0 K2 -(K2+K3) 0 0 0 0 0 0 0 0 0 0 K3;
%     0 0 0 -(K4) 0 0 0 0 0 0 0 0 0 K4;
%     K5 0 0 0 -(K5+K6) 0 0 0 K6 0 0 0 0 0;
%     K8 0 0 0 0 -(K8+K9) 0 0 0 K9 0 0 0 0;
%     K11 0 0 0 0 0 -(K11+K12) 0 0 0 K12 0 0 0;
%     K14 0 0 0 0 0 0 -(K14+K15) 0 0 0 K15 0 0;
%     0 0 0 0 K6 0 0 0 -(K6+K7) 0 0 0 K7 0;
%     0 0 0 0 0 K9 0 0 0 -(K9+K10) 0 0 0 K10;
%     0 0 0 0 0 0 K12 0 0 0 -(K12+K13) 0 0 K13;
%     0 0 0 0 0 0 0 K15 0 0 0 -(K15+K16) 0 K16;
%     0 0 0 0 0 0 0 0 K7 0 0 0 -(K7+K17) K17;
%     0 0 K3 K4 0 0 0 0 0 K10 K13 K16 K17 -(K3+K4+K10+K13+K16+K17)];

% Capacitance Matrix
C_amb=0; %b.c.
C_roof=0; %de momento
C_ceiling=0; %de momento
C_base=144240; %Capacitance, W
C_ws=0;
C_rw=0;
C_lsw=0;
C_rsw=0;

C_cabin=density_air*(CabinVolume)*cp_air/timestep;
C=zeros(13);
% C_cabin_front=density_air*(CabinVolume/2)*cp_air/timestep;
% C_cabin_back=density_air*(CabinVolume/2)*cp_air/timestep;
% C=zeros(14);

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
% C(13,13)=C_cabin_front;
% C(14,14)=C_cabin_back;

% %% GENERAL EQUATIONS
% % m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 
% T_air(t) = (Q_windows(t)+Q_ceiling(t)+Q_base(t)+Q_human(t))*(timestep/(CabinVolume*cp_air*density_air))  + T_air(t-1);
% Q_ws(t) = h_ws*A_ws*(T_ws(t)- T_air(t-1));
% Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
% Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));
% Q_windows(t)=Q_ws(t)+Q_lsw(t)+Q_rsw(t)+Q_rw(t); %% Listo

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
Tcabin(t)=T_ini;
temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];
delta_temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];

%%14 nodes:
% Tcabin_front(t)=T_ini;
% Tcabin_back(t)=T_ini;
% temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];
% delta_temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];

Q_windows(t) = 0;

for t=2:Total_time+1
        
    % Boundary conditions vector
    Qbase(t)=0; %alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_lsw*G_lsw_t(t)) + (A_lsw*G_rsw_t(t)) );
    Qhuman(t)=N_Humans * 120;
       
    Tbc=[T_amb;0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t)]; 
    %Tbc=[T_amb(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);-A_lsw*G_lsw_a(t);-A_rsw*G_rsw_a(t);0;0;0;0;-Qhuman(t);0]; 
    %Tbc=[T_amb;-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);-A_lsw*G_lsw_a(t);-A_rsw*G_rsw_a(t);0;0;0;0;-Qhuman(t)]; 
    
    %calculates different temperatures
    temperature=(inv(K - C))*(Tbc - C*delta_temperature);
    delta_temperature=temperature;
    
    %obtaining different temperatures of different elements
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
%     Tcabin_front(:,t)=temperature(13);
%     Tcabin_back(:,t)=temperature(14);
end

%% PLOTS 13 nodes
figure(1)
plot(time,Tcabin-273.15,'b',T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2),'k',T_Air_Exp(:,1),T_Air_Exp(:,2),'--')
legend('T cabin','T cabin model (Paper)','T cabin thermocouple (Paper)','Location','best')
grid on
grid minor
xlabel('Time [s]') 
ylabel('Temperature [ºC]')

figure(2)
plot(time,Tcabin-273.15,'b',T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2),'k',T_Air_Exp(:,1),T_Air_Exp(:,2),'--')
hold on
plot(time,Troof-273.15,time,T_amb-273.15)
legend('T cabin','T cabin model (Paper)','T cabin thermocouple (Paper)','T roof','T ambient (b.c.)','Location','best')
grid on
grid minor
xlabel('Time [s]') 
ylabel('Temperature [ºC]')

figure(3)
plot(time,Tcabin)
legend('T cabin','Location','best')
grid on
grid minor
xlabel('Time [s]') 
ylabel('Temperature [K]')

% %% PLOTS 14 nodes
% figure(1)
% plot(time,Tcabin_front-273.15,'b',time,Tcabin_back-273.15,'r',T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2),'k',T_Air_Exp(:,1),T_Air_Exp(:,2),'--')
% legend('T cabin front','T cabin back','T cabin model (Paper)','T cabin thermocouple (Paper)','Location','best')
% grid on
% grid minor
% xlabel('Time [s]') 
% ylabel('Temperature [ºC]')
% %axis([0 500 -350 100]) %xx yy
% 
% figure(2)
% plot(time,Tcabin_front-273.15,'b',time,Tcabin_back-273.15,'r',T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2),'k',T_Air_Exp(:,1),T_Air_Exp(:,2),'--')
% hold on
% plot(time,Troof-273.15,time,T_amb-273.15)
% legend('T cabin front','T cabin back','T cabin model (Paper)','T cabin thermocouple (Paper)','T roof','T ambient (b.c.)','Location','best')
% grid on
% grid minor
% xlabel('Time [s]') 
% ylabel('Temperature [ºC]')
% %axis([0 100 -350 100]) %xx yy



%% FIRST VERSION: iteration-substitution of Energy Balances Equations: running 8 min! 
% while (t <= Total_time)
% %% CALCULO DE COEFICIENTES DE CONVECCION
% % T_Sky set to Ambient Temperature
% T_sky = T_ext(t);
% %% WINDSHIELD
% % Q_cd(t)=Q_ws(t); %2nd Eq
% G_ws_r(t)=ro_ws*G_ws_inc(t);
% G_ws_a(t)=alfa_ws*G_ws_inc(t); %añadir como bc
% G_ws_t(t)=tao_ws*G_ws_inc(t);
% G_ws_inc(t)=G_ws_r(t)+G_ws_a(t)+G_ws_t(t);
% % Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
% % Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
% % Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));
% 
% % Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_ws*G_ws_a(t); %1er Eq
% 
% % 1st Eq (Sustitucion)
% % h_ext*A_ws*(T_ws_ext(t)-T_ext(t))+epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4) = (A_ws*k_ws/e_ws)*(T_ws_ext(t) - T_ws(t)) +  A_ws*G_ws_a(t);
% 
% % 2nd Eq 
% % T_ws_ext(t) = (h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t);
% 
% % Iteracion
%  It_max = 1000;
%  err_prev = 0;
%  err = 0;
%  T_ws(t) = T_ws(t-timestep);
%  rest  = false;
%  DT = 0.01;
%  for i = 0:0.02:It_max
%      
%      Ter1 = h_ext*A_ws*(((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))-T_ext(t))+epsilon_ws*sigma*A_ws*(((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))^4-(T_sky)^4);
%      Ter2 = (A_ws*k_ws/e_ws)*(T_ws(t) - ((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))) +  A_ws*G_ws_a(t);
%      err = abs(Ter2 - Ter1);
%      
%      if (err < err_prev)&&(rest == true)
%          T_ws(t) = T_ws(t) - DT;
%          rest = true; 
%      elseif (err < err_prev)&&(rest == false)
%          T_ws(t) = T_ws(t) + DT;
%          rest = false;
%      elseif (err >= err_prev)&&(rest == true)
%          T_ws(t) = T_ws(t) + DT;
%          rest = false; 
%      elseif (err >= err_prev)&&(rest == false)
%          T_ws(t) = T_ws(t) - DT;
%          rest = true;
%      end
%      err_prev = err;
%      if err < 0.01
%          i = It_max;
%      end
%  end
%  
% %Calcula T_ws_ext(t)
% T_ws_ext(t) = (h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t);
% % Comprobacion del balance
% % Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_ws*G_ws_a(t); %1er Eq
% % Q_cd(t)=Q_ws(t); %2nd Eq
% Q_ws(t) = h_ws*A_ws*(T_ws(t)- T_air(t-1));
% Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
% Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));
% 
% Conv_Rad_ws(t) = Q_cv_ext(t) + Q_rd_ext(t);
% Cond_RadAbs_ws(t) = - Q_cd(t) + A_ws*G_ws_a(t);
% Cond_ws(t) = Q_cd(t);
% Convin_ws(t) = Q_cd(t);
% Error_WS(t) = err;
% %% LEFT SIDE WINDOW
% G_lsw_r(t)=ro_lsw*G_lsw_inc(t);
% G_lsw_a(t)=alfa_lsw*G_lsw_inc(t);
% G_lsw_t(t)=tao_lsw*G_lsw_inc(t);
% G_lsw_inc(t)=G_lsw_r(t)+G_lsw_a(t)+G_lsw_t(t);
% 
% % Iteracion
%  It_max = 1000;
%  err_prev = 0;
%  err = 0;
%  T_lsw(t) = T_lsw(t-timestep);
%  rest  = false;
%  DT = 0.01;
%  for i = 0:0.02:It_max
%      
%      Ter1 = h_ext*A_lsw*(((h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t))-T_ext(t))+epsilon_lsw*sigma*A_lsw*(((h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t))^4-(T_sky)^4);
%      Ter2 = (A_lsw*k_lsw/e_lsw)*(T_lsw(t) - ((h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t))) +  A_lsw*G_lsw_a(t);
%      err = abs(Ter2 - Ter1);
%      
%      if (err < err_prev)&&(rest == true)
%          T_lsw(t) = T_lsw(t) - DT;
%          rest = true; 
%      elseif (err < err_prev)&&(rest == false)
%          T_lsw(t) = T_lsw(t) + DT;
%          rest = false;
%      elseif (err >= err_prev)&&(rest == true)
%          T_lsw(t) = T_lsw(t) + DT;
%          rest = false; 
%      elseif (err >= err_prev)&&(rest == false)
%          T_lsw(t) = T_lsw(t) - DT;
%          rest = true;
%      end
%      err_prev = err;
%      if err < 0.01
%          i = It_max;
%      end
%  end
%  
% %Calcula T_lsw_ext(t)
%  T_lsw_ext(t) = (h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t);
% % Comprobacion del balance
% % Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_lsw*G_lsw_a(t); %1er Eq
% % Q_cd(t)=Q_lsw(t); %2nd Eq
% Q_lsw(t) = h_lsw*A_lsw*(T_lsw(t)- T_air(t-1));
% Q_cv_ext_lsw(t)=h_ext*A_lsw*(T_lsw_ext(t)-T_ext(t)); %Exterior convection
% Q_rd_ext_lsw(t)=epsilon_lsw*sigma*A_lsw*((T_lsw_ext(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd_lsw(t)=((k_lsw*A_lsw)/e_lsw)*(T_lsw_ext(t)-T_lsw(t));
% 
% Conv_Rad_lsw(t) = Q_cv_ext_lsw(t) + Q_rd_ext_lsw(t);
% Cond_RadAbs_lsw(t) = - Q_cd_lsw(t) + A_lsw*G_lsw_a(t);
% Cond_lsw(t) = Q_cd_lsw(t);
% Convin_lsw(t) = Q_cd_lsw(t);
% 
% Error_LSW(t) = err;
% %% RIGHT SIDE WINDOW
% 
% G_rsw_r(t)=ro_lsw*G_rsw_inc(t);
% G_rsw_a(t)=alfa_rsw*G_rsw_inc(t);
% G_rsw_t(t)=tao_lsw*G_rsw_inc(t);
% 
% G_rsw_inc(t)=G_rsw_r(t)+G_rsw_a(t)+G_rsw_t(t);
% 
% % Iteracion
%  It_max = 1000;
%  err_prev = 0;
%  err = 0;
%  T_rsw(t) = T_rsw(t-timestep);
%  rest  = false;
%  DT = 0.01;
%  for i = 0:0.02:It_max
%      
%      Ter1 = h_ext*A_rsw*(((h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t))-T_ext(t))+epsilon_rsw*sigma*A_rsw*(((h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t))^4-(T_sky)^4);
%      Ter2 = (A_rsw*k_rsw/e_rsw)*(T_rsw(t) - ((h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t))) +  A_rsw*G_rsw_a(t);
%      err = abs(Ter2 - Ter1);
%      
%      if (err < err_prev)&&(rest == true)
%          T_rsw(t) = T_rsw(t) - DT;
%          rest = true; 
%      elseif (err < err_prev)&&(rest == false)
%          T_rsw(t) = T_rsw(t) + DT;
%          rest = false;
%      elseif (err >= err_prev)&&(rest == true)
%          T_rsw(t) = T_rsw(t) + DT;
%          rest = false; 
%      elseif (err >= err_prev)&&(rest == false)
%          T_rsw(t) = T_rsw(t) - DT;
%          rest = true;
%      end
%      err_prev = err;
%      if err < 0.01
%          i = It_max;
%      end
%  end
%  
% %Calcula T_rsw_ext(t)
%  T_rsw_ext(t) = (h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t);
% % Comprobacion del balance
% % Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_rsw*G_rsw_a(t); %1er Eq
% % Q_cd(t)=Q_rsw(t); %2nd Eq
% Q_rsw(t) = h_rsw*A_rsw*(T_rsw(t)- T_air(t-1));
% Q_cv_ext_rsw(t)=h_ext*A_rsw*(T_rsw_ext(t)-T_ext(t)); %Exterior convection
% Q_rd_ext_rsw(t)=epsilon_rsw*sigma*A_rsw*((T_rsw_ext(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd_rsw(t)=((k_rsw*A_rsw)/e_rsw)*(T_rsw_ext(t)-T_rsw(t));
% 
% Conv_Rad_rsw(t) = Q_cv_ext_rsw(t) + Q_rd_ext_rsw(t);
% Cond_RadAbs_rsw(t) = - Q_cd_rsw(t) + A_rsw*G_rsw_a(t);
% Cond_rsw(t) = Q_cd_rsw(t);
% Convin_rsw(t) = Q_cd_rsw(t);
% 
% Error_RSW(t) = err;
% 
% %% REAR WINDOW
% 
% G_rw_r(t)=ro_rw*G_rw_inc(t);
% G_rw_a(t)=alfa_rw*G_rw_inc(t);
% G_rw_t(t)=tao_rw*G_rw_inc(t);
% 
% G_rw_inc(t)=G_rw_r(t)+G_rw_a(t)+G_rw_t(t);
% % Iteracion
%  It_max = 1000;
%  err_prev = 0;
%  err = 0;
%  T_rw(t) = T_rw(t-timestep);
%  rest  = false;
%  DT = 0.01;
%  for i = 0:0.02:It_max
%      
%      Ter1 = h_ext*A_rw*(((h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t))-T_ext(t))+epsilon_rw*sigma*A_rw*(((h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t))^4-(T_sky)^4);
%      Ter2 = (A_rw*k_rw/e_rw)*(T_rw(t) - ((h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t))) +  A_rw*G_rw_a(t);
%      err = abs(Ter2 - Ter1);
%      
%      if (err < err_prev)&&(rest == true)
%          T_rw(t) = T_rw(t) - DT;
%          rest = true; 
%      elseif (err < err_prev)&&(rest == false)
%          T_rw(t) = T_rw(t) + DT;
%          rest = false;
%      elseif (err >= err_prev)&&(rest == true)
%          T_rw(t) = T_rw(t) + DT;
%          rest = false; 
%      elseif (err >= err_prev)&&(rest == false)
%          T_rw(t) = T_rw(t) - DT;
%          rest = true;
%      end
%      err_prev = err;
%      if err < 0.01
%          i = It_max;
%      end
%  end
%  
% %Calcula T_rw_ext(t)
%  T_rw_ext(t) = (h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t);
% % Comprobacion del balance
% % Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_rw*G_rw_a(t); %1er Eq
% % Q_cd(t)=Q_rw(t); %2nd Eq
% Q_rw(t) = h_rw*A_rw*(T_rw(t)- T_air(t-1));
% Q_cv_ext_rw(t)=h_ext*A_rw*(T_rw_ext(t)-T_ext(t)); %Exterior convection
% Q_rd_ext_rw(t)=epsilon_rw*sigma*A_rw*((T_rw_ext(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd_rw(t)=((k_rw*A_rw)/e_rw)*(T_rw_ext(t)-T_rw(t));
% 
% Conv_Rad_rw(t) = Q_cv_ext_rw(t) + Q_rd_ext_rw(t);
% Cond_RadAbs_rw(t) = - Q_cd_rw(t) + A_rw*G_rw_a(t);
% Cond_rw(t) = Q_cd_rw(t);
% Convin_rw(t) = Q_cd_rw(t);
% 
% Error_RW(t) = err;
% 
% 
% %% GENERAL WINDOW EQUATION
% 
% Q_windows(t)=Q_ws(t)+Q_lsw(t)+Q_rsw(t)+Q_rw(t); %% Listo
% 
% %% ROOF AND CEILING
% 
% G_roof_a(t)=alfa_roof*G_roof_inc(t);
% 
% % Iteracion
%  It_max = 1000;
%  err_prev = 0;
%  err = 0;
%  T_ceiling(t) = T_ceiling(t-timestep);
%  rest  = false;
%  DT = 0.01;
%  for i = 0:0.02:It_max
%      
%      Ter1 = h_roof*A_roof*(((h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t))-T_ext(t))+epsilon_roof*sigma*A_roof*(((h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t))^4-(T_sky)^4);
%      Ter2 = (A_roof*k_ceiling/e_ceiling)*(T_ceiling(t) - ((h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t))) +  A_roof*G_roof_a(t);
%      err = abs(Ter2 - Ter1);
%      
%      if (err < err_prev)&&(rest == true)
%          T_ceiling(t) = T_ceiling(t) - DT;
%          rest = true; 
%      elseif (err < err_prev)&&(rest == false)
%          T_ceiling(t) = T_ceiling(t) + DT;
%          rest = false;
%      elseif (err >= err_prev)&&(rest == true)
%          T_ceiling(t) = T_ceiling(t) + DT;
%          rest = false; 
%      elseif (err >= err_prev)&&(rest == false)
%          T_ceiling(t) = T_ceiling(t) - DT;
%          rest = true;
%      end
%      err_prev = err;
%      if err < 0.01
%          i = It_max;
%      end
%  end
%  
% %Calcula T_rw_ext(t)
%  T_roof(t) = (h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t);
% % Comprobacion del balance
% % Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_roof*G_roof_a(t); %1er Eq
% % Q_cd(t)=Q_ceiling(t); %2nd Eq
% Q_ceiling(t) = h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1));
% Q_cv_ext_roof(t)=h_ext*A_roof*(T_roof(t)-T_ext(t)); %Exterior convection
% Q_rd_ext_roof(t)=epsilon_roof*sigma*A_roof*((T_roof(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd_ceiling(t)=((k_ceiling*A_roof)/e_ceiling)*(T_roof(t)-T_ceiling(t));
% 
% Conv_Rad_ceiling(t) = Q_cv_ext_roof(t) + Q_rd_ext_roof(t);
% Cond_RadAbs_ceiling(t)= - Q_cd_ceiling(t) + A_roof*G_roof_a(t);
% Cond_ceiling(t) = Q_cd_ceiling(t);
% Convin_ceiling(t) = Q_ceiling(t);
% 
% Error_CR(t) = err;
% 
% 
% 
% %% BASE
% % 1a Eq: (C_base*(T_base(t)-T_base(t-1))/(timestep)) = alfa_base*((A_ws*G_ws_t)+(A_rw*G_rw_t)+(A_lsw*G_lsw_t)+(A_lsw*G_rsw_t))-Q_base;
% 
% % 2a Eq:   Q_base=h_base*A_base*(T_base(t)-T_air(t-1));
% 
% % Sustituyendo 2a en 1a
% T_base(t) = (timestep/(C_base+h_base*A_base)) * (alfa_base*((A_ws*G_ws_t(t))+(A_rw*G_rw_t(t))+(A_lsw*G_lsw_t(t))+(A_lsw*G_rsw_t(t))) + C_base*T_base(t-1) + h_base*A_base*T_air(t-1));
% Q_base(t)=h_base*A_base*(T_base(t)-T_air(t-1));
% 
% %% BEINGS
% Q_human(t)=N_Humans * 120;
% %% GENERAL EQUATION
% % m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 
% 
% T_air(t) = (Q_windows(t)+Q_ceiling(t)+Q_base(t)+Q_human(t))*(timestep/(CabinVolume*cp_air*density_air))  + T_air(t-1);
% 
% t = t +timestep;
% end

% PLOT

% % Temp Air
% plot(time(2:end),T_air - 273,'b',T_Air_Exp(:,1),T_Air_Exp(:,2),'r',T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2),'g');

