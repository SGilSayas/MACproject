%%
clear all;
close all;

%% HVAC Model
% This model calculates the HVAC system of an xEV vehicle


% Sign criteria ----->   [+] Positive heat == Heat absorbed by the refrigerant
%               ----->   [-] Negative heat == Heat rejected by the refrigerant

%% Set PATH FOR PYTHON and CoolProp Installation
% pcPythonExe = 'C:\Users\amdrben\AppData\Local\Programs\Python\Python38\python.exe';
% [ver, exec, loaded]	= pyversion(pcPythonExe);
% pyversion;
% 
% % CoolProp Instalation library instalation -python
% [v,e] = pyversion; system([e,' -m pip install --user -U CoolProp'])
 

%% Boundary conditions
% Load Inputs from cabin model

load Q_Air_Cabin_T1.mat
load Inputs06212013.mat


% GENERAL EQUATION
% m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 

N_human = 1;
Q_human = 120 * N_human;
Q_air_cabin = Q_air_cabin + Q_human;
heat = max(Q_air_cabin) ;


% Variable time
timestep = 1; % 1s
Total_time = Irradiancia(end,1);
time = 0:timestep:Total_time;

t = 1;

% Interpolacion para el vector tiempo
Irr = interp1(Irradiancia(:,1),Irradiancia(:,2),time);
% Vel = interp1(Veh_speed(:,1),Veh_speed(:,2),time);

% Sources (Temp_conf && T_ext)
T_ext = 273 + interp1(Temp_Ext(:,1),Temp_Ext(:,2),time);
Temp_conf = 23 + 273;
T_ext_const = mean(T_ext);

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


%% Initialization
Q_evap_req(t) = 0; % Heat exchanged in the evaporator [W]
Q_cond_req(t) = 0; % Heat echanged in the condenser   [W]
W_comp(t) = 0; % Work performed by the compresor [W]
COP(t) = 0; % Coefficient of performance of the HP [-]
mf(t) = 0; % Compressor refrigerant mass flow [kg/s]
eta_is(t) = 0.8; % Compressor efficiency [-]


if (mean(Q_air_cabin)>= 0)
   Mode = true;
%    Q_evap(t) = abs(Q_air_cabin(t));
   Q_evap(t) = heat;
else
   Mode = false;
%    Q_cond(t) = abs(Q_air_cabin(t));
   Q_cond(t) = heat;
end

% Summer Mode
if (Mode == true)
% Heat exchangers
T_evap(t) = Temp_conf;
% T_cond(t) = T_ext(t);
T_cond(t) = T_ext_const;
else
% Winter Mode
T_evap(t) = T_ext_const;
% T_evap(t) = T_ext(t);
T_cond(t) = Temp_conf;
end

%% Calculating Loop
t = 1 + timestep;
while (t <= Total_time)

    
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
P3(t)=P2(t); %% Liquido subenfriado (Punto 3)
T3(t)=py.CoolProp.CoolProp.PropsSI('T','P',P3(t),'Q',0,Ref);
h3(t)=py.CoolProp.CoolProp.PropsSI('H','P',P3(t),'Q',0,Ref);
s3(t)=py.CoolProp.CoolProp.PropsSI('S','P',P3(t),'Q',0,Ref);
s4(t)=s3(t);
h4(t)=py.CoolProp.CoolProp.PropsSI('H','S',s4(t),'P',P1(t),Ref);
T4(t)=py.CoolProp.CoolProp.PropsSI('T','S',s4(t),'P',P1(t),Ref);
Q4(t)=py.CoolProp.CoolProp.PropsSI('Q','S',s4(t),'P',P1(t),Ref);
%%
COP(t)=(h2(t)-h3(t))/(h2(t)-h1(t));

if (Mode == true)
    mf(t) = Q_evap_req(t)/(h1(t) - h4(t));
else
    mf(t) = Q_cond_req(t)/(h2(t) - h3(t));  
end
Q_evap(t) = mf(t)*(h1(t) - h4(t));
Q_cond(t) = mf(t)*(h3(t) - h2(t));
W_comp(t) = mf(t)*(h2(t) - h1(t));

t = t + timestep;
end

%% Ploting
plottemplate1p(time(2:end-1), COP(2:end),'COP [-]','HP COP');
tufteAxesOrig;
plottemplate1p(time(2:end-1), mf(2:end)*1000,'Mass flow [g/s]','Mass flow');
tufteAxesOrig;
plottemplate1p(time(2:end-1), W_comp(2:end),'Power [W]','Compressor Power');
tufteAxesOrig;
plottemplate1p(time(2:end-1), Q_evap(2:end),'Heat [W]','Evaporator heat');
tufteAxesOrig;
plottemplate1p(time(2:end-1), Q_cond(2:end),'Heat [W]','Condenser heat');
tufteAxesOrig;




