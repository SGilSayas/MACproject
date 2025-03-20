% Vehicle Cabin Model
% April 2021, One node model. Created by Alejandor Vega and Amin Dreif 
% May 2021, upgrade to lumped model (13 and 14 nodes) by Susana Gil-Sayas
clc
clear all 
close all
%% INPUTS
timestep =  1; %s

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
%S_cabin=cabin_h*

% Air properies
rho_air = 1.18; %Kg/m3 (25ºC)
cp_air = 1006; %J/kg*K  "

% Load input IRRADIANCE 6/21/2013
t = 1;
load('Inputs.mat');
load('T_Air_mod_Ref.mat');
    % Variable time
Total_time = Irradiancia(end,1);
time = 0:timestep:Total_time;%-1;
    % Interpolacion: Irr = interp1(Irradiancia(:,1),Irradiancia(:,2),time);
    
Total_time=Total_time*2;
time = 0:timestep:Total_time;
Irr=zeros(Total_time);

% Ambient temperature
T_ext = 273.15 + interp1(Temp_Ext(:,1),Temp_Ext(:,2),time);
Tamb=T_ext; T_ini=T_ext(1);

% Tamb_C = -30; %T ambient in celsius
% for j=2:Total_time+1
% Tamb(j) = Tamb_C + 273.15;
% %Irr(j)=0;
% end
% T_ini=Tamb(2);%Tamb_C + 273.15;

%T_sky=T_ini-6; %
T_sky = 293;
N_Humans=0;

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
%Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_lsw*G_lsw_t(t)) + (A_lsw*G_rsw_t(t)) );

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
delta_temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini+1];

% Troof(t)=T_ini(t);
% Tceiling(t)=T_ini(t);
% Tbase_int(t)=T_ini(t);
% Tw_ext_ws(t)=T_ini(t);
% Tw_ext_rw(t)=T_ini(t);
% Tw_ext_lsw(t)=T_ini(t);
% Tw_ext_rsw(t)=T_ini(t);
% Tw_int_ws(t)=T_ini(t);
% Tw_int_rw(t)=T_ini(t);
% Tw_int_lsw(t)=T_ini(t);
% Tw_int_rsw(t)=T_ini(t);
% Tcabin(t)=T_ini(t);
% temperature=[T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1)];
% delta_temperature=[T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1);T_ini(1)];

%Q_windows(t) = 0;

for t=2:Total_time+1
%
Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_lsw*G_lsw_t(t)) + (A_lsw*G_rsw_t(t)) );
Qhuman(t)=N_Humans * 120;

Q_req(t)=rho_air*V_cabin*cp_air*(temperature(13)-Ttarget)/timestep; %Tcabin(:,t-1)=324.7876
%Tbc = [Tamb(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);-A_lsw*G_lsw_a(t);-A_rsw*G_rsw_a(t);0;0;0;0;-Qhuman(t)+Qcool_req(t)]; 
%No IRR:
Tbc = [Tamb(t);-A_roof*G_roof_a(t);0;0;-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);-A_lsw*G_lsw_a(t);-A_rsw*G_rsw_a(t);0;0;0;0;-Qhuman(t)+Q_req(t)]; %if 50:+Qcool

    temperature=(inv(K - C))*(Tbc - C*delta_temperature);
    delta_temperature=temperature;
    
    Tamb(:,t)=temperature(1); %T_amb es condición de contorno, T la calculada
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

Q_used(:,t)=rho_air*V_cabin*cp_air*(Tcabin(:,t)-Ttarget)/timestep;
Q_used_cum=Q_used(t-1)+Q_used(t);
end
%Tamb(1)=Tamb(2);
Q_req(:,1)=Q_req(:,2);
Q_used(:,1)=Q_used(:,2);

%error=Qcool_actual-Qcool_req; %fila
% %if loop: 
% Qcool_req=Qcool_actual;
% E(i,:)=error;
% TCABIN(i,:)=Tcabin;
% i=i+1; %contador de cada run
%T_cabin_final=Tcabin(:,trip_duration); %RUN#1: 52ºC; %RUN#2: 324.8446K=51.6946ºC ; %RUN#3: 324.8446

% figure(1)
% plot(time,Tcabin-273.15,'b',T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2),'k',T_Air_Exp(:,1),T_Air_Exp(:,2),'--')
% legend('T cabin','T cabin model (Paper)','T cabin thermocouple (Paper)','Location','best')
% grid on
% grid minor
% xlabel('Time [s]') 
% ylabel('Temperature [ºC]')

% figure(2)
% yyaxis left
% plot(time,Tcabin-273.15,'b')
% hold on
% plot(time,Troof-273.15,time,Tamb-273.15)
% grid on
% grid minor
% ylabel('Temperature [ºC]')
% xlabel('Time [s]') 
% %yline(23)
% yyaxis right
% plot(time,Q_req)
% hold on
% plot(time,Q_used,'k')
% ylabel('Cooling heat [W]')
% legend('T cabin','T roof','T ambient (b.c.)','T target','Q req','Q used','Location','NorthWest')


ymax=max(Tamb)-273.15+0.1;
ymin=min(Tamb)-273.15-0.1;
figure(4)
title('T target = 23ºC; Q for heating <0; Q for cooling >0')
yyaxis left
plot(time,Tcabin-273.15,'ro',time,Tamb-273.15,'.') %antes: time
ylim([ymin ymax])
xlim([-100 13300])
grid on
grid minor
disp(Q_used_cum)
ylabel('Temperature [ºC]')
xlabel('Time [s]') 
%yline(Ttarget-273.15)
yyaxis right
plot(time,Q_req)
hold on
plot(time,Q_used,'k')
ylabel('Heat [W] = f(dT)')
legend('T cabin','T ambient (b.c.)','Qcool req','Qcool actual','Location','Best')
%legend('T cabin','T ambient (b.c.)','T target','Qcool req','Qcool actual','Location','NorthWest')
