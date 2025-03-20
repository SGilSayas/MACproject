% Vehicle's Cabin Model, April 2021, Susana Gil-Sayas
clc
clear all 
close all

%% INPUTS:

% Test conditions
Irr = 0; % [W/(m2*K)] 
N_Humans = 0.25;
T_ambient = 22; % [C] % NOT USED WHEN IMPORTING DATA

% Initial conditions
T_ini = 273.15 + T_ambient; % [K], initial value for simulation
T_sky = T_ini - 6; % [K]

% Time
t = 1;
Total_time = 1800; % [s], WLTC duration
timestep =  1; % [s]
time = 0:timestep:Total_time;
if (Irr==0)
    T_amb = T_ambient + 273.15; % [K]

else
    T_amb = ones(1,Total_time+1)*(T_ambient+273.15); %[K], vector
    Irr = ones(1,Total_time+1)*Irr;
end

% Radiation Properties 
sigma=0.0000000567037321; % Stefan-Boltzmann constant(sigma)[W/(m^2*K^4)]
% epsilon:emissivity[], rho=reflectivity[], tau=transmissivity[], alpha=absorptivity[]
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

% Conduction Properties 
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

% Heat transfer Convection coefficients, W((m2*K)
h_ws=5; 
h_lsw=5;
h_rsw=5;
h_rw=5;
h_ceiling=5;
h_roof=20;
h_ext=20; 
h_base = 5;

% Geometry
cabin_height =  0.5456;
cabin_width = 1.1;
cabin_roof_lenght = 1.8;
cabin_base_lenght = 1.45;
CabinVolume = (cabin_width*cabin_roof_lenght*cabin_height/2)+...
    (cabin_width*cabin_base_lenght*cabin_height/2); 

% Air properies
density_air = 1.18; % [Kg/m3] (at 25C)
cp_air = 1006; % [J/kg*K]

% Load input for validation: Golf 7 22C: LAB_20220511_06_MAC_on
t = 1;
timestep =  1; %s
lab_input='G7_22C_MACoff.mat';
struct=load(lab_input);
table=struct2table(struct);
lab_table=table2array(table);
compressorspeed_rpms=table2array(lab_table(:,1));
compressortorque_Nm=table2array(lab_table(:,2));
TC_vent=table2array(lab_table(:,3)); % C
TC_front=table2array(lab_table(:,4)); % C
TC_driver=table2array(lab_table(:,5)); % C
TC_pass=table2array(lab_table(:,6)); % C
TC_back=table2array(lab_table(:,7)); % C
TC_backpass=table2array(lab_table(:,8)); % C
TC_cell=table2array(lab_table(:,9)); % C
T_cell_OBD=table2array(lab_table(:,10)); % C
T_evap=table2array(lab_table(:,11)); % C
vehiclespeed_kmh=table2array(lab_table(:,12));

T_amb=transpose(TC_cell)+273.15; % T_amb in K
%T_amb=transpose(T_cell_OBD)+273.15; % T_amb in K

%time = 0:timestep:Total_time;%-1;
%Interpolacion: Irr = interp1(Irradiancia(:,1),Irradiancia(:,2),time);
% Ambient temperature: T_ext = 273.15 + interp1(Temp_Ext(:,1),Temp_Ext(:,2),time); Tamb=T_ext; T_ini=T_ext(1);

% Human heat % W, UNE-EN ISO 8996, Alberto Viti Corsi (IDAE)
Qhuman = ones(1,Total_time+1)*N_Humans * 120; % [W]

% Global Horizontal Irradiance (Hypothesis: same in all windows)
G_roof_inc=Irr;
G_ws_inc=Irr; 
G_rw_inc=Irr;
G_lsw_inc=Irr;
G_rsw_inc=Irr;

% Solar radiations
%ws
G_ws_r=rho_ws*G_ws_inc;
G_ws_a=alfa_ws*G_ws_inc; %bc
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

% Boundary conditions vector
Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ ...
    (A_lsw*G_lsw_t(t))+ (A_lsw*G_rsw_t(t)) );
Tbc=[T_amb(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);...
    -A_rw*G_rw_a(t);-A_lsw*G_lsw_a(t);-A_rsw*G_rsw_a(t);0;0;0;0;-Qhuman(t);0]; 

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
K11=h_ext*A_lsw;
K12=A_lsw*k_lsw/e_lsw;
K13=h_lsw*A_lsw;
% Tamb --K14-> Tw_ext_rsw --K15-> Tw_int_rsw --K16-> Tcabin_back
K14=h_ext*A_rsw;
K15=A_rsw*k_rsw/e_rsw;
K16=h_rsw*A_rsw;
% Tcabin_front --K17-> Tcabin_back
K17=h_ws*A_ws; %TBD
K=[1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    K1 -(K1+K2) K2 0 0 0 0 0 0 0 0 0 0 0;
    0 K2 -(K2+K3) 0 0 0 0 0 0 0 0 0 0 K3;
    0 0 0 -(K4) 0 0 0 0 0 0 0 0 0 K4;
    K5 0 0 0 -(K5+K6) 0 0 0 K6 0 0 0 0 0;
    K8 0 0 0 0 -(K8+K9) 0 0 0 K9 0 0 0 0;
    K11 0 0 0 0 0 -(K11+K12) 0 0 0 K12 0 0 0;
    K14 0 0 0 0 0 0 -(K14+K15) 0 0 0 K15 0 0;
    0 0 0 0 K6 0 0 0 -(K6+K7) 0 0 0 K7 0;
    0 0 0 0 0 K9 0 0 0 -(K9+K10) 0 0 0 K10;
    0 0 0 0 0 0 K12 0 0 0 -(K12+K13) 0 0 K13;
    0 0 0 0 0 0 0 K15 0 0 0 -(K15+K16) 0 K16;
    0 0 0 0 0 0 0 0 K7 0 0 0 -(K7+K17) K17;
    0 0 K3 K4 0 0 0 0 0 K10 K13 K16 K17 -(K3+K4+K10+K13+K16+K17)];

% Capacitance Matrix
C_amb=0; %b.c.
C_roof=0; %for now
C_ceiling=0; %for now
C_base=144240; %Capacitance, W
C_ws=0;
C_rw=0;
C_lsw=0;
C_rsw=0;
C_cabin_front=density_air*(CabinVolume/2)*cp_air/timestep;
C_cabin_back=density_air*(CabinVolume/2)*cp_air/timestep;
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

% %% GENERAL EQUATIONS
% % m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 
% T_air(t) = (Q_windows(t)+Q_ceiling(t)+Q_base(t)+Q_human(t))*(timestep/(CabinVolume*cp_air*density_air)) + T_air(t-1);
% Q_ws(t) = h_ws*A_ws*(T_ws(t)- T_air(t-1));
% Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
% Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));
% Q_windows(t)=Q_ws(t)+Q_lsw(t)+Q_rsw(t)+Q_rw(t);

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
temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];
delta_temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];

Q_windows(t) = 0;

for t=2:Total_time+1
    % Boundary conditions vector
    if(Irr==0)
        Qbase(t)=0;
        Qhuman(t)=N_Humans * 120;
        if(size(T_amb)==[1 1])
            Tbc=[T_amb;0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t);0];
        else % size(T_amb)==[1 1801]
            Tbc=[T_amb(t);0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t);0];
        end
    else
        Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_lsw*G_lsw_t(t)) + (A_lsw*G_rsw_t(t)) );
        Tbc=[T_amb(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);-A_lsw*G_lsw_a(t);-A_rsw*G_rsw_a(t);0;0;0;0;-Qhuman(t);0]; 
        Qhuman(t)=N_Humans * 120;
    end

    % Temperatures calculation
    temperature=(inv(K - C))*(Tbc - C*delta_temperature);
    delta_temperature=temperature;
    
    % Temperatures extraction
    Tamb(:,t)=temperature(1); %T_amb: boundary condition, Tamb: predicted
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

    Q_cv_front(:,t)=-h_ws*A_ws*(Tw_int_ws(:,t)-Tcabin_front(:,t));
    Q_cv_back(:,t)=-h_rw*A_rw*(Tw_int_rw(:,t)-Tcabin_back(:,t));
    Q_cv_left(:,t)=-h_lsw*A_lsw*(Tw_int_lsw(:,t)-Tcabin_back(:,t));
    Q_cv_right(:,t)=-h_rsw*A_rsw*(Tw_int_rsw(:,t)-Tcabin_back(:,t));
    Q_cv_top(:,t)=-h_ceiling*A_roof*(Tceiling(:,t)-Tcabin_back(:,t));

    Q_cv_cabin(:,t)=Q_cv_front(:,t)+Q_cv_back(:,t)+Q_cv_left(:,t)+Q_cv_right(:,t)+Q_cv_top(:,t);

end



%% PLOTS
figure(2)
plot(time,Q_cv_cabin,':','LineWidth',1.5)
hold on
plot(time,Q_cv_front,time,Q_cv_back,time,Q_cv_top)
hold on
plot(time,Q_cv_left)

legend('Cabin','Front','Back','Top','Left/Right','Location','best')
xt=[1000 800 1000 1000 1000];
yt = [Q_cv_front(end) Q_cv_back(end) Q_cv_left(end) Q_cv_top(end) Q_cv_cabin(end)];
str = [string(Q_cv_front(end)) string(Q_cv_back(end)) string(Q_cv_left(end)) string(Q_cv_top(end)) string(Q_cv_cabin(end))];
text(xt,yt,str)

% Deletion of first 0 value for plotting (inicialization value)
x=Tamb(2); Tamb(1)=x;

% The boundary conditon 'T_amb' needs to be a vector to be plotted
if (Irr==0) & (size(T_amb)==[1 1])
    T_amb = ones(1,Total_time+1)*(T_ambient+273.15); %[K], vector
end
% figure(1)
% plot(time,Tcabin_front-273.15,':','LineWidth',1.5)
% hold on
% plot(time,Tcabin_back-273.15,'-.','LineWidth',1.5)
% plot(time,Troof-273.15,'--','LineWidth',1.5)
% plot(time,T_amb-273.15,'LineWidth',1.5)
% hold off
% legend('T cabin front','T cabin back','T roof','T ambient','Location','best')
% grid on 
% grid minor
% xlabel('Time [s]') 
% ylabel('Temperature [ºC]')
% axis([0 1200 0 80])

figure(1)
plot(time,Tcabin_front-273.15,'-.b','LineWidth',1.5)
hold on
plot(time,Tcabin_back-273.15,':b','LineWidth',1.5)
%plot(time,Troof-273.15,'--','LineWidth',1.5)
plot(time,Tamb-273.15,'-g','LineWidth',1)
%plot(time,T_amb-273.15)
plot(time,TC_front,'-.r','LineWidth',1.5)
plot(time,TC_back,':r','LineWidth',1.5)
%plot(time,TC_cell)
% hold on
% yyaxis right
% plot(time, vehiclespeed_kmh)
hold off
legend('T front sim','T back sim','Tamb','TC front','TC back','Vehicle Speed km/h','Location','best')
grid on 
grid minor
xlabel('Time [s]') 
ylabel('Temperature [ºC]')
axis([1 1800 18 30])
