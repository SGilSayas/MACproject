% Vehicle's Cabin Model, Susana Gil-Sayas
% v2: ceiling, half base and half roof go to front node
% h_ext_over80kmh=50, lsw becomes sidewindows, rsw becomes doors
% UA_front and UA_back
clc
clear all 
close all

%% INPUTS:

% Test conditions
Irr = 0; % [W/(m2*K)] 
N_Humans = 1; % old: 1 human = 120 W, 0.66 = 80 W, 0.25 = 30 W;
heat_human = 80; %W
heat_equipment = 40 ; %40; % W
T_ambient = 22; % [C] % NOT USED WHEN IMPORTING DATA

% Initial conditions
T_ini = 273.15 + T_ambient; % [K], initial value for simulation
T_sky = T_ini - 6; % [K]
cell_length = 6; %m
cell_width = 5; %m
cell_heigh = 4; %m       
V_cell = cell_length * cell_width * cell_heigh * 1.2;

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
% Roof and ceiling
A_roof=1.8*1.1;
k_ABS=0.1; % Acrylonitrile butadiene styrene
k_steel=14.9;
k_PU = 0.022; % 0.022-0.028 W/mK polyurethane
k_ceiling=k_ABS*k_steel*k_PU;
e_ABS=0.002;
e_steel=0.0005;
e_PU = 0.025; % Polyurethane 
e_ceiling=e_ABS+e_steel+e_PU;
% Windshield
A_ws=0.63*1.3; 
k_ws=0.8; % 0.8-1.4
e_ws=0.006;
% Rear windows
A_rw=0.29*1; 
k_rw=1.4; % 1.4-1.5
e_rw=0.005;
% Side windows
A_sidewindows=2*1.45*0.29; 
k_sidewindows=1.4;
e_sidewindows=0.003;
% Doors
A_doors=2*1.45*0.29;
k_doors=k_ceiling;
e_doors=e_ceiling;
% Base
A_base=6;

A_tot = A_ws + A_rw + A_sidewindows + A_doors + A_roof;
A_front = A_ws + A_sidewindows/2 + A_roof/2 + A_doors/2 + A_base/2;
A_back = A_rw + A_sidewindows/2 + A_doors/2 + A_roof/2 + A_base/2;

% Heat transfer Convection coefficients, W((m2*K)
h_ws=5;
h_ext_over80=50;
h_lsw=5;
h_rsw=5;
h_rw=5;
h_ceiling=5;
h_roof=20;
h_ext=20; 
h_base = 5;
h_cabin=5;%10;

% Overall heat transfer coefficients, UA (W/K)
UA_front=(h_cabin*A_ws) + (k_ws*A_ws)/e_ws + (k_sidewindows*A_sidewindows/2)/e_sidewindows + (k_doors*A_doors/2)/e_doors + (k_ceiling*A_roof/2)/e_ceiling;
UA_back= (h_cabin*A_ws) + (k_rw*A_rw)/e_rw + (k_sidewindows*A_sidewindows/2)/e_sidewindows + (k_doors*A_doors/2)/e_doors + (k_ceiling*A_roof/2)/e_ceiling;

% R_front= 1/h_cabin + e_ws/(k_ws) + e_sidewindows/(k_sidewindows*2) + e_doors/(k_doors*2) + e_ceiling/(k_ceiling*2) + 1/h_ext; U_front= 1 * 1/R_front;
% R_back = 1/h_cabin + e_av_back/(k_av_back) + 1/h_ext; U_back= 1 * 1/R_back ;

% Geometry
cabin_height =  0.5456;
cabin_width = 1.1;
cabin_roof_lenght = 1.8;
cabin_base_lenght = 1.45;
CabinVolume = (cabin_width*cabin_roof_lenght*cabin_height/2)+...
    (cabin_width*cabin_base_lenght*cabin_height/2); 

% Air properies
density_air = 1.18; % [Kg/m3] (at 25C)
cp_air = 1006;% + 10*1500; % [J/ kg*K]

% Seats properties
mass_seats = 30; % 20-30 kg
cp_seats = 2000; % 500-2000 [J/kgK]
density_PE = 860; % [kg/m3] 
cp_PE = 2302.7; %polyethylene [J/ kg*K]

% Load input for validation: Golf 7 22C: LAB_20220511_06_MAC_on
lab_input='G7_35C_MACoff_TC';
struct=load(lab_input);
table=struct2table(struct);
lab_array=table2array(table);

cell_array=table2array(lab_array(:,7));
vent_array=table2array(lab_array(:,1));
front_array=table2array(lab_array(:,2));
driver_array=table2array(lab_array(:,3));
pass_array=table2array(lab_array(:,4));
back_array=table2array(lab_array(:,5)); 
backpass_array=table2array(lab_array(:,6));
n=10;
TC_cell=cell_array(1 : n : end);
TC_vent=vent_array(1 : n : end);
TC_front=front_array(1 : n : end);
TC_driver=driver_array(1 : n : end);
TC_pass=pass_array(1 : n : end);
TC_back=back_array(1 : n : end);
TC_backpass=backpass_array(1 : n : end);

t = 1; timestep =  1; % [s]
Total_time=length(TC_cell)-1; % 1800; % [s], WLTC duration
time = 0:timestep:Total_time;
if (Irr==0)
    T_amb = T_ambient + 273.15; % [K]
else
    T_amb = ones(1,Total_time+1)*(T_ambient+273.15); %[K], vector
    Irr = ones(1,Total_time+1)*Irr;
end

%lab_input='G7_XCUdata_35C_MACoff.mat';%'G7_22C_MACoff.mat';
% struct=load(lab_input);
% table=struct2table(struct);
% lab_table=table2array(table);
% compressorspeed_rpms=table2array(lab_table(:,1));
% compressortorque_Nm=table2array(lab_table(:,2));
% TC_vent=table2array(lab_table(:,3)); % C
% TC_front=table2array(lab_table(:,4)); % C
% TC_driver=table2array(lab_table(:,5)); % C
% TC_pass=table2array(lab_table(:,6)); % C
% TC_back=table2array(lab_table(:,7)); % C
% TC_backpass=table2array(lab_table(:,8)); % C
% TC_cell=table2array(lab_table(:,9)); % C
% T_cell_OBD=table2array(lab_table(:,10)); % C
% T_evap=table2array(lab_table(:,11)); % C
% vehiclespeed_kmh=table2array(lab_table(:,12));

T_amb=transpose(TC_cell)+273.15; % T_amb in K. Both TC_cell or T_cell_OBD
i=TC_front(1)+273.15-T_amb(1); 
T_ini=T_amb(1); %TC_front(1)+273.15+1; %T_amb(1)+i; 

% Correcting measurement error in thermocouples
TC_front=TC_front+(T_amb(1)-(TC_front(1)+273.15));
TC_back=TC_back+(T_amb(1)-(TC_back(1)+273.15));

% Human and equipment heat % ISO 8996, Alberto Viti Corsi (IDAE)
Qhuman = ones(1,Total_time+1)*N_Humans*heat_human; % [W]
Qequipment = ones(1,Total_time+1)*heat_equipment; % [W]

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
    (A_sidewindows*G_lsw_t(t))+ (A_sidewindows*G_rsw_t(t)) );

Tbc=[T_amb(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);...
    -A_rw*G_rw_a(t);-A_sidewindows*G_lsw_a(t);-A_doors*G_rsw_a(t);0;0;0;0;-Qhuman(t);-Qequipment(t)]; 

%% Conductances Matrix (heat from ambient to cabin)
% Tamb --K1-> Troof --K2-> Tceiling --K3-> Tcabin_back and Tcabin_front
K1=h_ext*A_roof; %h_roof?
K2=A_roof*k_ceiling/e_ceiling;
K3=h_ceiling*A_roof; %Aroof/2 to the front and Aroof/2 to the back

% Tbase --K4-> Tcabin_back
K4=h_base*A_base; %Abase/2 to the front and Abase/2 to the back

% Tamb --K5-> Tw_ext_ws --K6-> Tw_int_ws --K7-> Tcabin_front
K5=h_ext*A_ws;
K6=A_ws*k_ws/e_ws;
K7=h_ws*A_ws;

% Tamb --K8-> Tw_ext_rw --K9-> Tw_int_rw --K10-> Tcabin_back
K8=h_ext*A_rw;
K9=A_rw*k_rw/e_rw;
K10=h_rw*A_rw;

% Tamb --K11-> Twindows_ext --K12-> Twindows_int --K13-> Tcabin_back&front
K11=h_ext*A_sidewindows;
K12=A_sidewindows*k_sidewindows/e_sidewindows;
K13=h_lsw*A_sidewindows;

% Tamb --K14-> Tdoors_ext --K15-> Tdoors_int --K16-> Tcabin_front&back
K14=h_ext*A_doors;
K15=A_doors*k_doors/e_doors;
K16=h_rsw*A_doors;

% Tcabin_front --K17-> Tcabin_back
K17=h_cabin*cabin_height*cabin_width;

K= [(1) 0 0 0 0 0 0 0 0 0 0 0 0 0;
    (K1) -(K1+K2) K2 0 0 0 0 0 0 0 0 0 0 0;
    0 (K2) -(K2+K3/2+K3/2) 0 0 0 0 0 0 0 0 0 (K3/2) (K3/2);
    0 0 0 -(K4/2+K4/2) 0 0 0 0 0 0 0 0 (K4/2) (K4/2);
    (K5) 0 0 0 -(K5+K6) 0 0 0 (K6) 0 0 0 0 0;
    (K8) 0 0 0 0 -(K8+K9) 0 0 0 (K9) 0 0 0 0;
    (K11) 0 0 0 0 0 -(K11+K12) 0 0 0 (K12) 0 0 0;
    (K14) 0 0 0 0 0 0 -(K14+K15) 0 0 0 (K15) 0 0;
    0 0 0 0 (K6) 0 0 0 -(K6+K7) 0 0 0 (K7) 0;
    0 0 0 0 0 (K9) 0 0 0 -(K9+K10) 0 0 0 (K10);
    0 0 0 0 0 0 (K12) 0 0 0 -(K12+K13/2+K13/2) 0 (K13/2) (K13/2);
    0 0 0 0 0 0 0 (K15) 0 0 0 -(K15+K16/2+K16/2) (K16/2) (K16/2);
    0 0 (K3/2) (K4/2) 0 0 0 0 (K7) 0 (K13/2) (K16/2) -(K3/2+K4/2+K7+K16/2+K17) (K17);
    0 0 (K3/2) (K4/2) 0 0 0 0 0 (K10) (K13/2) (K16/2) (K17) -(K3/2+K4/2+K10+K13/2+K16/2+K17)];


%% Capacitance Matrix
C_amb=0; %b.c.
C_roof=0; %for now
C_ceiling=0; %for now
C_base= 144240 + (mass_seats*cp_seats/timestep) ; %Capacitance, W/K 
C_ws=0;
C_rw=0;
C_sidewindows=0;
C_doors=0;
C_cabin_front=(density_air*(CabinVolume/2)*cp_air/timestep);
C_cabin_back=(density_air*(CabinVolume/2)*cp_air/timestep);
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

% %% GENERAL EQUATIONS
% % m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 
% T_air(t) = (Q_windows(t)+Q_ceiling(t)+Q_base(t)+Q_human(t))*(timestep/(CabinVolume*cp_air*density_air)) + T_air(t-1);
% Q_ws(t) = h_ws*A_ws*(T_ws(t)- T_air(t-1));
% Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
% Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));
% Q_windows(t)=Q_ws(t)+Q_lsw(t)+Q_rsw(t)+Q_rw(t);

%% INI & CALC
Troof(t)=T_ini;
Tceiling(t)=T_ini;
Tbase_int(t)=T_ini;
Tws_ext(t)=T_ini;
Trw_ext(t)=T_ini;
Tsidewindows_ext(t)=T_ini;
Tdoors_ext(t)=T_ini;
Tws_int(t)=T_ini;
Trw_int(t)=T_ini;
Tsidewindows_int(t)=T_ini;
Tdoors_int(t)=T_ini;
Tcabin_front(t)=T_ini;
Tcabin_back(t)=T_ini;
temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];
delta_temperature=[T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini;T_ini];
check=0; check_f=0; check_b=0; check_t=0;
% j13(t)=TC_front(t)-T_amb(t); j14(t)=TC_back(t)-T_amb(t);

for t=2:Total_time+1

    if t>1200
        h_ext=h_ext_over80;
    end

    % Boundary conditions vector
    if(Irr==0)
        Qbase(t)=0;
        Qhuman(t) = N_Humans * heat_human;
        Qequipment(t) = heat_equipment;
        if(size(T_amb)==[1 1]) % constant value
            Tbc=[T_amb;0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t);-Qequipment(t)];
        else % size(T_amb)==[1 1801] vector from input
            Tbc=[T_amb(t);0;0;-Qbase(t);0;0;0;0;0;0;0;0;-Qhuman(t);-Qequipment(t)];
        end
    else
        Qbase(t)=alfa_base*((A_ws*G_ws_t(t)) + (A_rw*G_rw_t(t))+ (A_sidewindows*G_lsw_t(t)) + (A_sidewindows*G_rsw_t(t)) );
        Tbc=[T_amb(t);-A_roof*G_roof_a(t);0;-Qbase(t);-A_ws*G_ws_a(t);-A_rw*G_rw_a(t);-A_sidewindows*G_lsw_a(t);-A_doors*G_rsw_a(t);0;0;0;0;-Qhuman(t);-Qequipment(t)]; 
        Qhuman(t) = N_Humans * heat_human;
        Qequipment(t) = heat_equipment;
    end

    % Temperatures calculation
    temperature=(inv(K - C))*(Tbc - C*delta_temperature); % temperature=(inv(K - C))*(Tbc - C*delta_temperature);
    delta_temperature=temperature;

    if(temperature(13)>temperature(14)) 
        check_t(t)=0.5; 
    else check_t(t)=0; 
    end
  
    
  % Exchange from the front and back
    if ( temperature(13) > T_amb(t) ) % && temperature(14)<temperature(13) )
        temperature(13)= (UA_front*temperature(1) + density_air*V_cell*delta_temperature(13)/timestep) / (UA_front + density_air*V_cell/(timestep));
        delta_temperature(13)=temperature(13);
        check_f=check_f+1;
    end
    if ( temperature(14) > T_amb(t) ) % && temperature(14)>temperature(13) )
        temperature(14) = (UA_back*temperature(1) + density_air*V_cell*delta_temperature(14)/timestep) / (UA_back + density_air*V_cell/(timestep));         
        %temperature(14) = (h_cabin*cabin_height*cabin_width*temperature(13) + density_air*CabinVolume*delta_temperature(14)/timestep) / (h_cabin*cabin_height*cabin_width + density_air*CabinVolume/(timestep));
        delta_temperature(14)=temperature(14);
        check_b=check_b+1;
    end


    % Temperatures extraction
    Tamb(:,t)=temperature(1); %T_amb: boundary condition, Tamb: predicted
    Troof(:,t)=temperature(2);
    Tceiling(:,t)=temperature(3);
    Tbase_int(:,t)=temperature(4);
    Tws_ext(:,t)=temperature(5);
    Trw_ext(:,t)=temperature(6);
    Tsidewindows_ext(:,t)=temperature(7);
    Tdoors_ext(:,t)=temperature(8);
    Tws_int(:,t)=temperature(9);
    Trw_int(:,t)=temperature(10);
    Tsidewindows_int(:,t)=temperature(11);
    Tdoors_int(:,t)=temperature(12);
    Tcabin_front(:,t)=temperature(13);
    Tcabin_back(:,t)=temperature(14);
    Tcabin(:,t)= (Tcabin_front(:,t)+Tcabin_back(:,t))/2;

    %Heats calculation
    Q_cv_front(:,t)=-h_ws*A_ws*(Tws_int(:,t)-Tcabin_front(:,t));
    Q_cv_back(:,t)=-h_rw*A_rw*(Trw_int(:,t)-Tcabin_back(:,t));
    Q_cv_left(:,t)=-h_lsw*A_sidewindows*(Tsidewindows_int(:,t)-Tcabin_back(:,t));
    Q_cv_right(:,t)=-h_rsw*A_doors*(Tdoors_int(:,t)-Tcabin_back(:,t));
    Q_cv_top(:,t)=-h_ceiling*A_roof*(Tceiling(:,t)-Tcabin_back(:,t));
    Q_cv_cabin(:,t)=Q_cv_front(:,t)+Q_cv_back(:,t)+Q_cv_left(:,t)+Q_cv_right(:,t)+Q_cv_top(:,t);
end


%% PLOTS
% Deletion of first 0 value for plotting (inicialization value)
x=T_amb(2); T_amb(1)=x;
% The boundary conditon 'T_amb' needs to be a vector to be plotted
if (Irr==0) & (size(T_amb)==[1 1])
    T_amb = ones(1,Total_time+1)*(T_ambient+273.15); %[K], vector
end

figure(1)
plot(time,Tcabin_front-273.15,':b','LineWidth',1.5)
hold on
plot(time,Tcabin_back-273.15,':m','LineWidth',1.5)
plot(time,T_amb-273.15,'-k','LineWidth',1) %T amb 
plot(time,TC_front,'-b','LineWidth',1.5) % considerign ini error: +1.5
plot(time,TC_back,'-m','LineWidth',1.5) % considerign ini error: +0.7
hold off
legend('T front sim','T back sim','TC amb','TC front','TC back','Location','best')
grid minor
xlabel('Time [s]') 
ylabel('Temperature [ºC]')
axis([1 Total_time 33 38])

% figure(3)
% plot(time,Tcabin_front-273.15)
% hold on
% plot(time,Tcabin_back-273.15)
% hold on
% yyaxis right
% hold on
% plot(time,check_t)
% hold on
% ylim([0 1])

% figure(2)
% plot(time,Q_cv_cabin,':','LineWidth',1.5)
% hold on
% plot(time,Q_cv_front,time,Q_cv_back,time,Q_cv_top)
% hold on
% plot(time,Q_cv_left)
% legend('Cabin','Front','Back','Top','Left/Right','Location','best')
% xt=[1000 800 1000 1000 1000];
% yt = [Q_cv_front(end) Q_cv_back(end) Q_cv_left(end) Q_cv_top(end) Q_cv_cabin(end)];
% str = [string(Q_cv_front(end)) string(Q_cv_back(end)) string(Q_cv_left(end)) string(Q_cv_top(end)) string(Q_cv_cabin(end))];
% text(xt,yt,str)

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
