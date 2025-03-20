clear all 
close all 
clc

% Material data
% density (kg/m3)
density_copper = 8966;
density_iron = 7877;
density_aluminium = 2700.787;
density_air = 1.28;
density_glass=15; %sgs (check: 24.6;)
% thermoconductivity (mW/(m*K))
thermconductivity_copper=413;thermconductivity_iron=94;thermconductivity_aluminium=237;
thermconductivity_air=26; %this data varied according to the temperature
% thermal conductivities %W/mK PAPER
k_glass=1.4;k_steel=14.9;k_air=26.3*10^-3;k_cotton=0.06; 
% specific heat capacity (J/KG.k)
heat_capacity_copper=385;heat_capacity_iron=450;heat_capacity_aluminium=902;
heat_capacity_air=1010;heat_capacity_glass=840; %sgs
% Heat capacities (Cp,Cv): simple model: constant; real model: Cp=f(T)
% Cp_glass=0.19; % Cp_air=0.24; %kcal/(kg*ºC)
% Cv_glass=2.85; Cv_air=0.29; %kcal/(m3*ºC)

%% INPUT data
t=1;
area_windows=1.5295; %m2, D. Marcos et al.
volume_cabin=1.8*1.1*1.1; %m3, approximation from roof size in D. Marcos et al.
h_interiorair=5;  %W(m2*K), natural convection, no interior ventilation activated
h_exteriorair=5; %W(m2*K), natural convection, when vehicle 0-8 km/h
temp_ext(:,t)=273.15+25; %K

%Convective conductances:
Kinw=area_windows*h_interiorair;
Kexw=area_windows*h_exteriorair;

% Values for inicialization
temp_prev=temp_ext;
temp_int(:,t)=temp_prev;
temp_window(:,t)=temp_prev;
delta_temp=[temp_prev;temp_prev;temp_prev];

% Heat and rediation 
people=1; %number of people in the vehicle
G_hum=-(65*people); %W, ref: UNE-EN ISO 8996 and Alberto Viti Corsi (IDAE)
G_rad=0; %calor por radiación, metido como generación del cristal

Cw=heat_capacity_glass; % specific heat capacity (J/KG.k)
Cin=heat_capacity_air; % specific heat capacity (J/KG.k)

% Matrix definition
K=[Kinw -Kinw 0;Kinw -(Kinw+Kexw) Kexw;0 0 1]; %Convective conductances
C=[Cin 0 0;0 Cw 0;0 0 0]; %Capacitances
temp_bc=[G_hum;G_rad;temp_prev]; %Boundary conditions

%% Heat transfer
%Paper %%%% We look for: temp_air (temp air inside cabin)
Aws=0.63*1.3; Arw=0.29*1; Asw=1.45*0.29; % length,m * width,m
Aw=[Aws Arw Asw];
thermconductivity_air=0.026; %W/(m*K)
L= [0.63 0.29 1.45] %lenght, m
Pr_air=0.71; %Prandtl
kinem_visc_air=0.00001349; %m^2/s (0ºC, 1 bar), kinematic viscosity
vol_thermal_exp_coef_air=0.00369; %1/K, volumetric thermal expansion coefficient air (beta)

%initial values to calculate h:
temp_w=[temp_ext temp_ext temp_ext]; temp_w_ext=[temp_ext temp_ext temp_ext]; temp_air=[temp_ext temp_ext temp_ext];
% -> with h -> balances
% -> from balances: new temps to be used next loop AND temp_air (cabin) 
% -> write temp_air(t)
g=9.81; %m/s2
theta=0; %degrees of vertical windows, base case: totally vertical
Gr=g*cos(theta)*vol_thermal_exp_coef_air*(temp_w-temp_air).*L^3/kinem_visc_air^2;
Ra=Pr_air.*Gr;
Nu_vertical_flat=( 0.825+ 0.837.*Ra^(1/6)/( 1+(0.492*Pr_air)^(9/16) )^(8/17) )^2;
Nu=[Nu_vertical_flat Nu_vertical_flat Nu_vertical_flat];
hw=thermconductivity_air.*Nu./L; 
hws=hw(1);hrw=hw(2);hlsw=hw(3);hrsw=hw(3);hsw=hlsw+hrsw;

%external convective heat transfer coefficient, vehicle speed=0-8 km/h
he=5;

%heat transfer from the windows (vectors)
Qw=hw.*Aw.*(temp_w-temp_air);
Qwindows=Qw(1)+Qw(2)+2*Qw(3)=Qws+Qrw+2*Qsw %tot heat transfer w->cabin

%external convective heat flow:
Qcv_ext=he.*Aw.*(temp_w_ext-temp_ext)
%conductive heat flow interior windows:
Qcd=k_glass*Aw.*(temp_w_ext-temp_w)./thickness_w
%longwave radiative heat flow:
Qrd_ext=emissivity_glass*constant_SB*Aw.*((temp_w_ext)^4-(temp_sky)^4)
%irradiance reflected by the windows:
Gw_r=reflectivity_w.*Gw_inc
%irradiance absorbed by the windows:
Gw_a=absorptivity_w.*Gw_inc
%irradiance transmitted by the windows:
Gw_t=transmissivity_w.*Gw_inc

% Window heat balances to obtain temp_w_ext
%  Qcv_ext + Qrd_ext = Qcd + Aw.*Gw_a
%  Qcd = Qw
%  Gw_inc = Gw_r + Gw_a + Gw_t 
% m_air*cp_air*dtemp_air/dt=Qw+Qceiling+Qbase+Qhuman;

% Developed balances:
%he.*Aw.*(temp_w_ext-temp_ext)+emissivity_glass*constant_SB*Aw.*((temp_w_ext)^4-(temp_sky)^4)...
%    =k_glass*Aw.*(temp_w_ext-temp_w)./thickness_w + Aw.*(absorptivity_w.*Gw_inc)
%temp_air=temp_w-k_glass*Aw.*(temp_w_ext-temp_w)./(thickness_w.*hw.*Aw);
%1 = reflectivity_w + absorptivity_w + transmissivity_w 
%m_air*cp_air*dtemp_air/dt=Qw %probar primero con solo windows

%%
%initial values to calculate h:
temp_w=[temp_ext temp_ext temp_ext]; 
temp_w_ext=[temp_ext temp_ext temp_ext]; 
temp_air=[temp_ext temp_ext temp_ext];
%
thickness_ws =6*10^-3; thickness_rw =5*10^-3; thickness_sw =3*10^-3; 
thickness_w=[thickness_ws thickness_rw thickness_sw]; %m
constant_SB=5.67037321*10^-8; %W/(m^2*K^4) %Stefan-Boltzmann constant(sigma):
emissivity_glass=0.9;
absorptivity_glass=0.2;
absorptivity_ws=absorptivity_glass; absorptivity_rw=absorptivity_glass;absorptivity_sw=absorptivity_glass;
absorptivity_w=[absorptivity_ws absorptivity_rw absorptivity_sw];
temp_ext=298; %25ºC
temp_sky=temp_ext-6;
temp_w_ext=temp_ext+1; %inicialización
Aws=0.63*1.3; Arw=0.29*1; Asw=1.45*0.29; % length,m * width,m
Aw=[Aws Arw Asw];
thermconductivity_air=0.026; %W/(m*K)
L=[0.63 0.29 1.45]; %lenght, m
Pr_air=0.71; %Prandtl
kinem_visc_air=0.00001349; %m^2/s (0ºC, 1 bar), kinematic viscosity
vol_thermal_exp_coef_air=0.00369; %1/K, volumetric thermal expansion coefficient air (beta)
g=9.81; %m/s2
theta=0; %degrees of vertical windows, base case: totally vertical
he=5;
k_glass=1.4; % thermal conductivities %W/mK PAPER

j=1; % 1=windshield, 2=rear window, 3=side indows
% Gr wth (temp_w-temp_air) bc we look for the h inide the cabin
Gr(j)=g*cos(theta)*vol_thermal_exp_coef_air*(temp_w(j)-temp_air(j))*L(j)^3/((kinem_visc_air)^2);
Ra(j)=Pr_air.*Gr(j);
Nu_vertical_flat(j)=( 0.825 + 0.837*Ra(j)^(1/6) / ( 1+(0.492*Pr_air)^(9/16) )^(8/17) )^2;
Nu=[Nu_vertical_flat Nu_vertical_flat Nu_vertical_flat];
hw=thermconductivity_air*Nu(j)/L(j); 
density_air = 1.28; %(kg/m3)
volume_cabin=1.8*1.1*1.1; %m3,
m_air=volume_cabin*density_air; %kg, mass of indoor air
t=1; %min
heat_capacity_air=1010;
Gw_inc=120; %W/m^2
temp_air_previous=temp_ext; %inicialización
%ADD temp_air_previous=temp_air; %before loop

a1(j)=he*Aw(j);
a2(j)=k_glass*Aw(1)/thickness_w(j); 
a3(j)=emissivity_glass*constant_SB*Aw(j);
a4(j)=Aw(j)*absorptivity_w(j)*Gw_inc;
a5(j)=emissivity_glass*constant_SB*Aw(j)*temp_sky^4;
a6(j)=Aw(j)*hw(j);
a7=m_air*heat_capacity_air/t;
a8=a7*temp_air_previous;

% syms x1 x2 x3
% eq1 = a6(j)*x2 -(a7+a6(j))*x3 +a8;
% eq2 = (a1(j)+a2(j))*x1 +a3(j)*x1^4 + (a2(j)-a1(j))*x2 -(a4(j)+a5(j));
% eq3 = (a6(j)+a2(j))*x2 -a6(j)*x3 -a2(j)*x1;
% sol = solve(eq1,eq2,eq3)
%sol.x1;

A=[a1 a2 a3 0;-a7 0 a5 -a6;0 0 a6 -a9];
b=[-a8;a5;0];
n=length(b);
d=det(A);
x=zeros(n,1);
for i=1:n
    Ab=[A(:,1:i-1),b,A(:,i+1:n)];
    x(i)=det(Ab)/d;
end
disp('x')
disp(x)
%%

% Matrix definition
K=[Kinw -Kinw 0;Kinw -(Kinw+Kexw) Kexw;0 0 1]; %Convective conductances
C=[Cin 0 0;0 Cw 0;0 0 0]; %Capacitances
temp_bc=[G_hum;G_rad;temp_prev]; %Boundary conditions
%%

temp_w=[temp_ws temp_rw temp_sw]; %unknown %interior surface window temp 
temp_w_ext=[temp_ws_ext temp_rw_ext temp_sw_ext]; %unknown %exterior surface window temp   

% temp_sky: SKY TEMPERATURE MODELISATION AND APPLICATIONS IN BUILDING SIMULATION. L. ADELARD, F. PIGNOLET-TARDAN, T. MARA, P. LAURET, F. GARDE, H. BOYER
%sky equivalent temperatue (method in ref 17)

reflectivity_ws=0.246; reflectivity_rw=0.1; reflectivity_sw=0.2;
reflectivity_w=[reflectivity_ws reflectivity_rw reflectivity_sw];
transmissivity_ws=0.452; transmissivity_rw=0.311; transmissivity_sw=0.475;
transmissivity_w=[transmissivity_ws transmissivity_rw transmissivity_sw]

% abroptivity between 0.2 and 0.5, 
%https://www.researchgate.net/post/Absorption-coefficient-and-emissivity-of-glass

Qw=[Qws Qrw Qsw]; %unknown
Gw_inc=[Gws_inc Grw_inc Gsw_inc]; %incident irradiance (horizontal global);[W/m2]
Gw_r = [Gws_r Grw_r Gsw_r];
Gw_a = [Gws_a Grw_a Gsw_a];
Gw_t = [Gws_t Grw_t Gsw_t];

%%
for t=2:1000 %seconds
    
    % Calculation different temperatures
    temp=(inv(K - C))*(temp_bc - C*delta_temp);

    delta_temp=temp;
    results_temp(:,t)=temp;
    results_bc(:,t)=temp_bc;

    temp_int(:,t)=temp(1);
    temp_window(:,t)=temp(2);
    temp_ext(:,t)=temp(3);
end

figure(1)
plot(temp_int)

figure(2)
plot(temp_window)

figure(3)
plot(temp_ext)
