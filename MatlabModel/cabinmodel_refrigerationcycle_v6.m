clc;
clear;
close all;
clear py;

%%
% V6: it runs, PID working, energy calc implemented
%
% V5: check github versions comments
%
% V4_2: new version of v4 starting fresh from v3
% - 14 nodes modified to 3 nodes and it works
% - Modification of Q_engine and Q_exhaust, deleting "_av", so taking the
% cycle heat, not the average (needed for keepong 1800 s too)
% - Substitution ot T_amb for T_cell
% - Reading inputs from other folders done
% - Reset of loop end point (Total_time, instead of Total_time+1)
% - Reset initialization vectors to (Total_time, instead of Total_time+1) 
% - Reset limits for plots to keep same vector length
% - NEXT GOAL: the script out of the folder (not done)
%
% V3: it works
% - Initial heat exchanged cabin Air with the HVAC defined with T target
% - Addition of PID CONTROLLER for Compressor Speed. Input = error
% between current cabin temp and target cabin temperature. Output =
% compressor rpm, which later is used to calculate the refrigerant mass
% flow, which is used to obtain compressor's power and evaporator's heat.
%
% V2:
% - better check of compression ratio with an average minimum,
% - POINT (2) is not saturated liquid, but a saturated mix, so input
%   parameter is s2is = s1, not T2is=Tcond (COP=1.7)
% - DECIDE: mass flow fixed or dependant on Q_cabi_req?
%   fixed doesnt work but makes more physical sense
%
% V1:
% Combination cabin model and the refrigertion cycle.
% Cabin model: cabinmodelfunction4 (based on cabinmodel_v22)
% Refrigeration cycle: standard cycle (evap-comp-cond-valve),
% using Python prompt CoolProp to get refrigerant properties (p-h diagram)
% Refrigerant used: R1234yf
%
% Gil-Sayas, S., 2025
%
%% If it appears a message not finding CoolProp: Redo link python-matlab
% in matlab terminal:
% pyenv('Clear');
% pyenv('Version', 'C:\Users\susan\AppData\Python\Python310\python.exe')
% in matlab terminal: pyversion
% if all ok: check this on cmd window: pip install coolprop
% if all ok: check this on matlab terminal py.CoolProp.CoolProp.PropsSI('T','P',101325,'Q',0,'Water')

%% Input data files upload
timestep = 1;
folderPath = 'C:\Users\susan\OneDrive - UPV\Desktop_UPV\Doctorado\PhD Simulation Model\Lumped model\2025_Matlab Cabin HVAC model\thermalmodel_inputs';
% Read CSV files
csvFiles = dir(fullfile(folderPath, '*.csv'));
for i = 1:length(csvFiles)
    csvFileName = fullfile(folderPath, csvFiles(i).name);
    disp(['Reading CSV file: ', csvFileName]);
    
    % Read the entire CSV file into a table
    R1234yf_op_pres = readtable(csvFileName);
    
    % Generate a valid variable name from the filename
    csvFileBaseName = matlab.lang.makeValidName(erase(csvFiles(i).name, '.csv'));
    
    % Assign the whole table as a single variable in the workspace
    assignin('base', csvFileBaseName, R1234yf_op_pres);
   
    disp(['Extracted table: ', csvFileBaseName]);
end
% Read .mat files
matFiles = dir(fullfile(folderPath, '*.mat'));
for i = 1:length(matFiles)
    fileName = fullfile(folderPath, matFiles(i).name);
    disp(['Loading .mat file: ', fileName]);

    % Load .mat file into a struct
    matData = load(fileName);
    
    % Get all variable names in the .mat file
    fields = fieldnames(matData);
    
    % Assign each variable dynamically to the workspace
    for j = 1:length(fields)
        varName = fields{j}; % Extract variable name
        assignin('base', varName, matData.(varName)); % Assign in base workspace
        disp(['Extracted .mat variable: ', varName]);
    end
end
disp('All files have been read and extracted successfully.');

%% 
TCfiltering = 0; % 0; for no

%% TC filtering
if TCfiltering
    l=[length(TC_cell)];
    T=[fft(TC_cell(:,:)/l(1))];
    k=0.1;
    Tfil=T(:,:);
    Tfil(2:15,:)=T(2:15,:)*k; 
    Tfil(l-13:l,:)=T(l-13:l,:)*k;
    TC_fil=[ifft(Tfil(:,1).*l(1),l(1))];
    TC_cell=TC_fil(:,1);
    Total_time=length(TC_cell);
    time = 1:timestep:Total_time;
    time = time'; % Ensure it's a column vector

    % Extract the trend (data from 61 to 1800)
    time_trend = time(61:1770);
    TC_trend = TC_cell(61:1770);
    
    % Estimate frequency using FFT (optional but helps with accuracy)
    Y = fft(TC_trend - mean(TC_trend)); % Remove DC component
    p2 = abs(Y / length(TC_trend));
    p1 = p2(1:floor(length(TC_trend) / 2)); % Single-sided spectrum
    [~, idx] = max(p1(2:end)); % Find dominant frequency, ignoring DC component
    dominant_freq = (idx + 1) / length(TC_trend) * (0.4 * pi); % (2 * pi) %Convert to angular frequency
    
    % Define a sinusoidal fit model explicitly
    sinModel = fittype('a*sin(b*x + c) + d', 'independent', 'x', 'coefficients', {'a', 'b', 'c', 'd'});
    
    % Fit the data (adjust StartPoint based on estimated frequency)
    sin_fit = fit(time_trend, TC_trend, sinModel, ...
        'StartPoint', [0.5, dominant_freq, 0, mean(TC_trend)]);
    
    % Generate the first 60  and last 30 seconds using the fitted function
    time_first_60 = time(1:60);
    TC_substitute_ini = sin_fit.a * sin(sin_fit.b * time_first_60 + sin_fit.c) + sin_fit.d;
    TC_cell(1:60) = TC_substitute_ini;
    
    time_last_30 = time(1770:1800);
    TC_substitute_end = sin_fit.a * sin(sin_fit.b * time_last_30 + sin_fit.c) + sin_fit.d;
    TC_cell(1770:1800) = TC_substitute_end;
end
%% Test conditions
t = 1;
timestep = 1;
duration_WLTC = 1800; % s
Total_time = duration_WLTC;
time = 1:timestep:Total_time;
time = transpose(time);

% Cabin target temperature in C?
T_target = -10;
T_target = T_target + 273.15; %K

% Number of humans?
N_Humans = 1;

% Refrigeran fluid of MAC system
fluid = {'R1234yf'}; %{'Water'};{'R134a'};

% Vehicle fuel type?
fuel = 'PETROL'; %'DIESEL';

% Vehicle powertrain? 1 for yes, 0 for no
bev = 1;
phev_cd = 0;
phev_cs = 0;
icev = 0;

% Vehicle category?
hatchback = 0;
suv = 1;

% Setting a vehicle category factor 
if hatchback == 1    
    cat_factor = 1;
    min_mf = 20;
    max_mf = 50;
elseif suv == 1    
    cat_factor = 1.5;
    min_mf = 50;
    max_mf = 150;
else disp("Define vehicle category.")
end
min_mf = min_mf/1000; %kg/s
max_mf = max_mf/1000; 
% /small:20-50 g/s, SUVS:100-200, HDV:200-500/ mf(t)=100/1000; %kg/s


% Standard copressor's isoentropic efficiency, eta_is = (h2s-h1)/(h2-h1)
eta_is = 0.8; % av: 0.8;
compressor_on_lower_limit = 600; % W

% T between HX-sources:
% heat exchangers (evap&cond) and cold&hot sources (cabin&external air)
delta_T_evap = 7;   % Automotive: 7-12  ; HVAC: 5-10
delta_T_cond = 10;  % Automotive: 10-15 ; HVAC: 5-15

%% Ambient conditions
T_ini = TC_cell(1) + 273.15; % alternative: T_ini = TC_cabin_ini + 273.15;
T_cell = TC_cell + 273.15; %K

Irr_scalar = 0; % W/m2
Irr = ones(1,Total_time)*Irr_scalar;
% If not loading data:
    % T_cell_constant = 35; %C
    % T_cell = ones(1,Total_time)*T_cell_constant + 273.15; %K % Row vecto

%% Constant input data
% Cabin dimensions for category C
vehicle_height = 1.46;
vehicle_width = 1.8;
vehicle_lenght = 4.36;

% cabin_height = vehicle_height - 0.06 - 0.6215/2; % minus roof&base thickness and half of wheel
% cabin_width = vehicle_width - 0.06 - 0.241/2; 
% - doors thickness - mirrors=0.241, we use half bc some database width included mirrors and others did not.
% cabin_roof_lenght = vehicle_lenght*0.43;
% cabin_base_lenght = vehicle_lenght*0.35;
cabin_height =  0.5456*cat_factor;
cabin_width = 1.3*cat_factor;
cabin_roof_lenght = 1.8*cat_factor;
cabin_base_lenght = 1.45*cat_factor;

seatcm3 = 150*cat_factor; %115; % cm3
seats_volume = 5*seatcm3*10^(-6); % m3

V_cabin = (cabin_width*cabin_roof_lenght*cabin_height/2)+(cabin_width*cabin_base_lenght*cabin_height/2) - seats_volume;

% Area (m2)
A_ws=0.63*1.3*cat_factor;
A_rw=0.29*1*cat_factor;
A_roof=1.8*1.1*cat_factor;
A_sidewindows=2*1.45*0.29*cat_factor;
A_doors=2*1.45*0.29*cat_factor;
A_base = cabin_base_lenght*cabin_width*3;
A_dashboard = cabin_width*0.5; %m2
A_front = A_ws + A_dashboard; %3
A_back = A_rw; %1.5
A_side = A_sidewindows + A_doors + A_roof;
% A_front = A_ws + A_sidewindows/2 + A_doors/2 + A_roof/2;
% A_back = A_rw + A_sidewindows/2 + A_doors/2 + A_roof/2;

% Material's thickness, e (m)
e_ws=0.006;
e_rw=0.005;
e_ABS=0.002;
e_steel=0.0005;
e_PU = 0.025; % polyurethane
e_ceiling=e_ABS+e_steel+e_PU;
e_sidewindows=0.003;
e_doors=e_ceiling;

% Fuel
Vpe_petrol = 0.264; % l/kWh
CF_petrol = 2330; % gCO2/l
Vpe_diesel = 0.22; % l/kWh
CF_diesel = 2640; % gCO2/l
if strcmp(fuel, 'PETROL') || strcmp(fuel, 'PETROL/ELECTRIC')
    Vpe=Vpe_petrol;
    CF=CF_petrol;
elseif strcmp(fuel, 'DIESEL') || strcmp(fuel, 'DIESEL/ELECTRIC')
    Vpe=Vpe_diesel;
    CF=CF_diesel;
else
     disp('Specify PETROL or DIESEL as the last argument');
end

% Conduction Properties: k = thermoconductivity[W/(m*K)]
k_ws=0.8; % 0.8-1.4
k_rw=1.4; % 1.4-1.5
k_ABS=0.1; % Acrylonitrile butadiene styrene
k_steel=14.9;
k_PU = 0.022; % 0.022-0.028 W/mK polyurethane
k_ceiling=k_ABS*k_steel*k_PU;
k_sidewindows=1.4;
k_doors=k_ceiling;

% Heat transfer coefficient of interior cabin air, W/(K*m2)
h_cabin = 5; % 5-10

% Air properies
density_air = 1.18; % [Kg/m3] (at 25C)
cp_air = 1006;      % + 10*1500; % [J/ kg*K]

% Seats properties
mass_seats = 20;    % 20-30 kg
cp_seats = 2000;    % 500-2000 [J/kgK]

% Human and equipment heat % ISO 8996, Alberto Viti Corsi (IDAE)
heat_equipment = 40; % W
Qequipment = ones(1,Total_time)*heat_equipment; % W
A_skin = 1.5; % m2
T_skin = 24.648 + 273.15; % K

% Radiation Properties
sigma=0.0000000567037321; % Stefan-Boltzmann constant(sigma)[W/(m^2*K^4)]
    % epsilon:emissivity[]
epsilon_ws=0.9;
epsilon_rw=0.9;
epsilon_sidewindows=0.9;
epsilon_doors=0.9;
epsilon_roof=0.9;
    % rho=reflectivity[]
rho_ws=0.246;
rho_rw=0.1;
rho_sidewindows=0.2;
rho_doors=0.74;
rho_roof=0.04; %G7:0.74;
rho_base=0.3;
    % tau=transmissivity[]
tao_ws=0.452;
tao_rw=0.311;
tao_sidewindows=0.475;
tao_doors=0;
tao_roof=0.9;%0;
tao_base=0;
    % alpha=absorptivity[]
alpha_coef=1;
alpha_ws = alpha_coef*(1-rho_ws-tao_ws);
alpha_rw = alpha_coef*(1-rho_rw-tao_rw);
alpha_sidewindows = alpha_coef*(1-rho_sidewindows-tao_sidewindows);
alpha_doors = alpha_coef*(1-rho_doors-tao_doors);
alpha_roof = alpha_coef*(1-rho_roof-tao_roof);
alpha_base = alpha_coef*(1-rho_base-tao_base);

% Global Horizontal Irradiance (Hypothesis: same in all windows)
G_roof_inc=Irr;
G_ws_inc=Irr;
G_rw_inc=Irr;
G_sidewindows_inc=Irr;
G_doors_inc=Irr;
G_ws_r=rho_ws*G_ws_inc;
G_ws_a=alpha_ws*G_ws_inc; %bc
G_ws_t=tao_ws*G_ws_inc;
G_rw_r=rho_rw*G_rw_inc;
G_rw_a=alpha_rw*G_rw_inc;
G_rw_t=tao_rw*G_rw_inc;
G_sidewindows_r=rho_sidewindows*G_sidewindows_inc;
G_sidewindows_a=alpha_sidewindows*G_sidewindows_inc;
G_sidewindows_t=tao_sidewindows*G_sidewindows_inc;
G_doors_r=rho_doors*G_doors_inc;
G_doors_a=alpha_doors*G_doors_inc;
G_doors_t=tao_doors*G_doors_inc;
G_roof_a=alpha_roof*G_roof_inc;
% Not needed bc rad is added as G*A --> Qirr = Irr.*(A_front+A_side+A_back); % W

% Air vent volume and flow rate
vent_volumerate=0.0003; % m3/s, Faya 2013: 0.02; 0.001
leakage_volumerate = 0.0001; %% m3/s, Faya2014

% Average humidity ratio in gram of water per gram of dry air, X
Ps_22_kPa = 2.626; % kPa, Water saturation pressure at 22C
Ps_35_kPa = 5.63; % kPa, Water saturation pressure at 35C
Ps_kPa = Ps_22_kPa+(T_ini-22)/(35-22)*(Ps_35_kPa-Ps_22_kPa);
meanP_kPa_amb = mean(P_amb_kPa); % in kPa
X = 0.0732; % X = 0.62198*Ps_kPa.*RH_amb./(100.*P_amb_kPa - Ps_kPa.*RH_amb); % Faya, 2013

%% Calculation of Properties

% Enthalpies calculation. Humidity considered same in cabin and exterior
enthalpy_amb = 1006.*T_cell + (2501000 + 1770.*T_cell)*X; % J/kg, Faya 2013, % before ".*X"

% Capacitance Matrix
C_amb = 0; % boundary condition
C_base = mass_seats*cp_seats + 144240; %Capacitance, W/K
C_cabin_front = density_air*cp_air/timestep*V_cabin*0.5 + C_base/2;
C_cabin_back = density_air*cp_air/timestep*V_cabin*0.5 + C_base/2; 

C=zeros(3);
C(1,1)=C_amb;
C(2,2)=C_cabin_front;
C(3,3)=C_cabin_back;

%% Initialization of values
% Initial heat transfer coefficients, UA (W/K):
h_front(t) = 160;
h_side(t) = 0.9*h_front(t);
h_rear(t) = 0.1*h_front(t);
h_front = ones(1,Total_time)*h_front(t);
h_side = ones(1,Total_time)*h_side(t);
h_rear = ones(1,Total_time)*h_rear(t);

UA_front(t) = A_front*h_front(t) + 1*(A_side*h_side(t));
UA_back(t) = A_back*h_rear(t) + 1*(A_side*h_side(t));
UA_front = ones(1,Total_time).*UA_front(t);
UA_back = ones(1,Total_time).*UA_back(t);

%  Conductances Matrix
K1=h_side(t)*A_roof;
K2=A_roof*k_ceiling/e_ceiling;
K3=h_cabin*A_roof;
K4=h_cabin*A_base;
K5=h_front(t)*A_ws;
K6=A_ws*k_ws/e_ws;
K7=h_cabin*A_ws;
K8=h_rear(t)*A_rw;
K9=A_rw*k_rw/e_rw;
K10=h_cabin*A_rw;
K11=h_side(t)*A_sidewindows;
K12=A_sidewindows*k_sidewindows/e_sidewindows;
K13=h_cabin*A_sidewindows;
K14=h_side(t)*A_doors;
K15=A_doors*k_doors/e_doors;
K16=h_cabin*A_doors;
K17=h_cabin*cabin_width*cabin_height;

% initial cabin surfaces tmperatures
Troof(t)=T_ini;
Tceiling(t)=T_ini;
Tbase_int(t)=T_ini;
Tw_ext_ws(t)=T_ini;
Tw_ext_rw(t)=T_ini;
Tw_ext_sidewindows(t)=T_ini;
Tw_ext_doors(t)=T_ini;
Tw_int_ws(t)=T_ini;
Tw_int_rw(t)=T_ini;
Tw_int_sidewindows(t)=T_ini;
Tw_int_doors(t)=T_ini;
% Tcabin_front(t)=T_ini;
Tcabin_back(t)=T_ini;
Tcabin = ones(1,Total_time)*T_ini;
temperature=[T_ini;T_ini;T_ini];
prev_temp=[T_ini;T_ini;T_ini];

% Refrigerant cycle:
    % point 1: evaporator oulet = compressor inlet
    % point 2: compressor outlet = condenser inlet
    % point 3: condenser outlet = valve inlet
    % point 4: valve outlet = evaporator inlet
% initial parameters of refrigerant cycle followng p-h diagram
T1(t) = 0; % Temperature 1 [K]
P1(t) = 0; % Pressure 1 [Pa]
T2(t) = 0; % Temperature 2
P2(t) = 0; % Pressure 2 [Pa]
P2is(t) = 0; % Pressure 2' (isoentropic point) [Pa]
T3(t) = 0; % Temperature 3
P3(t) = 0; % Pressure 3 [Pa]
T4(t) = 0; % Temperature 4
P4(t) = 0; % Pressure 4 [Pa]
COP(t) = 0; % Coefficient of performance [-]

% initialization of heat flows
Qleakage = ones(1,Total_time);
Qvent = zeros(1,Total_time);
Qbase = zeros(1,Total_time); 
Qirr = zeros(1,Total_time);
Qhuman = ones(1,Total_time);
Qcabin_req = ones(1,Total_time);
Qcabin_received= ones(1,Total_time);
Qcv_emitted(t) = 0;
Qcv_emitted = ones(1,Total_time);
Qcabin_tot = zeros(1,Total_time);
Q_out_front = ones(1,Total_time);
Q_out_back = ones(1,Total_time);
Qtarget = ones(1,Total_time);
comp_cum_Wh = ones(1,Total_time);
Qcompressor = ones(1,Total_time).*0;
heat_human = ones(1,Total_time);
Q_evap_req = zeros(1,Total_time); % Heat exchanged in the evaporator [W]
Q_cond_req = zeros(1,Total_time); % Heat echanged in the condenser   [W]
W_comp(t) = 0; % Work performed by the compresor [W]
Q_evap(t) = 0; % Heat exchanged in the evaporator [W]
Q_cond(t) = 0; % Heat echanged in the condenser   [W]

%% v3: Initial heat exchanged of cabin Air with the HVAC, W
Q_MAC = zeros(1,Total_time);
Q_MAC(t) = cp_air*V_cabin*density_air*(T_target - T_cell(1))/timestep;

% previous versions initializations
cabin_received_cum_Wh = ones(1,Total_time);
cabin_received_Wh_WLTP = ones(1,Total_time);
mac_needed_Wh = ones(1,Total_time);
mac_needed_cum_Wh = ones(1,Total_time);
mac_needed_Wh_WLTP = ones(1,Total_time);
comp_Wh = ones(1,Total_time);
comp_cum_Wh = ones(1,Total_time);
comp_Wh_WLTP = ones(1,Total_time);
active_time(t) = 0;
active_time_simulated(t) = 0;
cooling_time_simulated(t) = 0;
heating_time_simulated(t) = 0;
compressor_time(t)=0;

% initializations of flags, PID variables, and others
check_f = zeros(1,Total_time);
check_b = zeros(1,Total_time);
temp_error_prev = 0;
error_sum = 0;
temp_error(t) = 0; % Temperature error [-]
temp_error = zeros(1,Total_time);
W_comp(t) = 0;

mf = zeros(1,Total_time); %20/1000; % Compressor refrigerant mass flow [kg/s]
comp_speed = ones(1, Total_time)*1000; % Compressor speed [rpm]
CO2 = ones(1,Total_time).*0;
cooling = zeros(1,Total_time);
heating = zeros(1,Total_time);
mac_off = zeros(1,Total_time);

% Initialize PID variables
error = zeros(1,Total_time);
error_integral = zeros(1,Total_time);
error_derivative = zeros(1,Total_time);
prev_error = 0;
prev_error_integral = 0;
PID_output = zeros(1,Total_time);

%% CALCULATION LOOP
for t = 2:Total_time
    %% 1. Determine Mode (Cooling/Heating)
    cooling(t) = (T_cell(t) > T_target);  % True if cooling
    heating(t) = (T_cell(t) < T_target);
    mac_off(t) = (T_cell(t) == T_target);

    %% 2. Heat Transfer Coefficients
    if t > 1800-323
        h_front(t) = 200;
    elseif t > 1800-323-455
        h_front(t) = 190;
    elseif t > 1800-323-455-433
        h_front(t) = 170;
    else
        h_front(t) = 160;
    end
    h_side(t) = 0.9*h_front(t);
    h_rear(t) = 0.1*h_front(t);
    UA_front(t) = A_front*h_front(t) + 1*(A_side*h_side(t));
    UA_back(t) = A_back*h_rear(t) + 1*(A_side*h_side(t));

    % Conductances Matrix
    % Tamb --K1-> Troof --K2-> Tceiling --K3-> Tcabin_back
    K1=h_side(t)*A_roof;
    K2=A_roof*k_ceiling/e_ceiling;
    K3=h_cabin*A_roof;
    % Tbase --K4-> Tcabin_back
    K4=h_cabin*A_base;
    % Tamb --K5-> Tw_ext_ws --K6-> Tw_int_ws --K7-> Tcabin_front
    K5=h_front(t)*A_ws;
    K6=A_ws*k_ws/e_ws;
    K7=h_cabin*A_ws;
    % Tamb --K8-> Tw_ext_rw --K9-> Tw_int_rw --K10-> Tcabin_back
    K8=h_rear(t)*A_rw;
    K9=A_rw*k_rw/e_rw;
    K10=h_cabin*A_rw;
    % Tamb --K11-> Tw_ext_sidewindows --K12-> Tw_int_sidewindows --K13-> Tcabin_back
    K11=h_side(t)*A_sidewindows;
    K12=A_sidewindows*k_sidewindows/e_sidewindows;
    K13=h_cabin*A_sidewindows;
    % Tamb --K14-> Tw_ext_doors --K15-> Tw_int_doors --K16-> Tcabin_back
    K14=h_side(t)*A_doors;
    K15=A_doors*k_doors/e_doors;
    K16=h_cabin*A_doors;
    % Tcabin_front --K17-> Tcabin_back
    K17=h_cabin*cabin_width*cabin_height;%A_ws; %TBD
    
    % Three-nodes approach:
    UAroof = 1/( 1/K1 + 1/K2 + 1/K3 );
    UAbase = 1/( 1/K4 );
    UAws = 1/( 1/K5 + 1/K6 + 1/K7 );
    UArs = 1/( 1/K8 + 1/K9 + 1/K10 );
    UAsw = 1/( 1/K11 + 1/K12 + 1/K13 );
    UAdoors = 1/( 1/K14 + 1/K15 + 1/K16 );
    UAcabin = 1/( 1/K17 );
   
    Kfront = UAws + UAroof/2 + UAbase/2 + UAsw/2 + UAdoors/2;
    Kback = UArs + UAroof/2 + UAbase/2 + UAsw/2 + UAdoors/2;
    Kcabin = UAcabin;
    
    K = [1 0 0;
        Kfront -(Kfront+Kcabin) Kcabin;
        Kback Kcabin -(Kback+Kcabin)];
    
    UA_cabin = 1/K1 + 1/K2 + 1/K3 + 1/K4 + 1/K5 + 1/K6 + 1/K7 + 1/K8 + ...
        1/K9 + 1/K10 + 1/K11 + 1/K12 + 1/K13 + 1/K14 + 1/K15 + 1/K16 + 1/K17;

     %% 4. Heat Flows
    enthalpy_cabin(t) = 1006*Tcabin(t-1) + (2501000 + 1770*Tcabin(t-1))*X;
    Qleakage(t) = -leakage_volumerate*density_air*(enthalpy_amb(t)-enthalpy_cabin(t));
    Qvent(t) = vent_volumerate*density_air*(enthalpy_amb(t)-enthalpy_cabin(t));
    Qhuman(t) = N_Humans*h_cabin*A_skin*abs(Tcabin(t-1)-T_skin);
    
    %% 5. Refrigeration Cycle
    % Set evaporator/condenser temps based on mode
    if cooling(t)
        T_evap(t) = T_target - delta_T_evap;
        T_cond(t) = T_cell(t) + delta_T_cond;
    elseif heating(t)
        T_evap(t) = T_cell(t) - delta_T_evap;
        T_cond(t) = T_target + delta_T_cond;
    end
    
    % Call diagram p-h on CoolProp
    Ref = char(fluid);

    % Check if T_cell is in the table
    idx = find(R1234yf_op_pres.T_ambient_C == T_cell(t));

    if ~isempty(idx) % If T_cell is found in the table
        P_LowSide_max_Pa = 1000*R1234yf_op_pres.P_LowSide_max_kPa(idx);
        P_HighSide_max_Pa = 1000*R1234yf_op_pres.P_HighSide_max_kPa(idx);
    else % If T_cell is not found, interpolate
        T_ambient = R1234yf_op_pres.T_ambient_C;
        P_LowSide_max_Pa = 1000*R1234yf_op_pres.P_LowSide_max_kPa;
        P_HighSide_max_Pa = 1000*R1234yf_op_pres.P_HighSide_max_kPa;

        P_LowSide_max_Pa = interp1(T_ambient, P_LowSide_max_Pa, T_cell, 'linear', 'extrap');
        P_HighSide_max = interp1(T_ambient, P_HighSide_max_Pa, T_cell, 'linear', 'extrap');
    end

    compression_ratio_min = 5; % Average 5.7357

        % Refrigeration Cycle Points:
        % (1): evaporator outlet, salutared vapor;
        % (2): compressor outlet, superheated vapor;
        % (3): condenser outlet, saturated fluid/liquid;
        % (4): expansion device outlet, vapor-fluid mix;

    % (1) Evaporator's outlet = Compressor's inlet, Saturated vapor Q=1
    T1(t)= T_evap(t); % T_evap_out is known
    P1(t)=py.CoolProp.CoolProp.PropsSI('P','T',T1(t),'Q',1,Ref);
    P1(t) = double(P1(t));

    % Ensure evaporator pressure does not exceed the maximum low-side pressure
    if (P1(t) >= P_LowSide_max_Pa(t))
        P1(t) = P_LowSide_max_Pa(t);
        T1(t) = py.CoolProp.CoolProp.PropsSI('T','P',P1(t),'Q',1,Ref);
    end

    h1(t)=py.CoolProp.CoolProp.PropsSI('H','T',T1(t),'Q',1,Ref);
    s1(t)=py.CoolProp.CoolProp.PropsSI('S','T',T1(t),'Q',1,Ref);

    % (2') Isoentropic compression (theoretical)
    s2is(t)=s1(t);
    P2is(t)=py.CoolProp.CoolProp.PropsSI('P','T',T_cond(t),'S',s2is(t),Ref);

    % (2) Compressor's outlet = Condenser's inlet
    P2(t)=P2is(t);

    % Actual compression ratio
    compression_ratio_actual(t)=P2(t)/P1(t);

    % Compression ratio needs to meet/exceed the minimum one
    if compression_ratio_actual(t) < compression_ratio_min
        % Adjust P2 so compression_ratio_actual = compression_ratio_min
        P2(t) = compression_ratio_min * P1(t);
        P2is(t)=P2(t);

        % Recompute T_cond based on adjusted P2
        T_cond(t) = py.CoolProp.CoolProp.PropsSI('T','P',P2(t),'Q',1,Ref);

        % Update delta_T_cond to reflect the new condenser saturation temperature
        delta_T_cond = T_cond(t) - T_cell(t);
    end

    h2is(t)=py.CoolProp.CoolProp.PropsSI('H','P',P2(t),'S',s2is(t),Ref); % falls inside of saturated mix area, so we can use T_cond
    % or: h2is(t)=py.CoolProp.CoolProp.PropsSI('H','T',T_cond(t),'P',P2(t),Ref);

    % Get enthalpy from isoenthalic efficiency
    h2(t) = (h2is(t)-h1(t))/eta_is + h1(t);
    T2(t)=py.CoolProp.CoolProp.PropsSI('T','H',h2(t),'P',P2(t),Ref);
    s2(t)=py.CoolProp.CoolProp.PropsSI('S','H',h2(t),'P',P2(t),Ref);

    % (3) Condenser's outlet = Expansion vale's inlet, Saturated liquid
    % T3(t)=T_cond(t);
    P3(t)=P2(t);
    T3(t)=py.CoolProp.CoolProp.PropsSI('T','P',P3(t),'Q',0,Ref);
    h3(t)=py.CoolProp.CoolProp.PropsSI('H','P',P3(t),'Q',0,Ref);
    s3(t)=py.CoolProp.CoolProp.PropsSI('S','P',P3(t),'Q',0,Ref);

    % (4) Valve's outlet = Evaporator's inlet
    s4(t)=s3(t);
    P4(t)= P1(t);
    h4(t)=py.CoolProp.CoolProp.PropsSI('H','S',s4(t),'P',P4(t),Ref);
    T4(t)=py.CoolProp.CoolProp.PropsSI('T','S',s4(t),'P',P4(t),Ref);

    % State of each operating point
    Q1(t)=py.CoolProp.CoolProp.PropsSI('Q','S',s1(t),'P',P1(t),Ref);
    Q4(t)=py.CoolProp.CoolProp.PropsSI('Q','S',s4(t),'P',P1(t),Ref);

    % Coefficient of performance
    COP(t)=(h1(t)-h4(t))/(h2(t)-h1(t));

    %% 6. PID Controller
    
    % Set theoretical demad
    Qcabin_req(t) = cp_air*V_cabin*density_air*(T_target - Tcabin(t-1))/t;

    % Base mass flow value (theoretical)
    mf_req(t) = (Qcabin_req(t) / (h1(t)-h4(t))); % IF FIXED:mf_req(t)=50/1000; %kg/s

    % Proportional, Derivative and Integral Error
    error(t) = T_target - Tcabin(t-1);
    error_derivative(t) = (error(t) - prev_error) / timestep;
    if (mf(t) > min_mf) && (mf(t) < max_mf) % Only integrate error when the mass flow is not saturated
        error_integral(t) = prev_error_integral + error(t) * timestep;
    end
    
    % Anti-windup
    max_integral = 1000;
    error_integral = min(max(error_integral, -max_integral), max_integral);
    
    % PID Output (mass flow adjustment in kg/s)
    Kp = 0.02;    Ki = 0.002;    Kd = 0.0001;
    PID_output(t) = Kp*error(t) + Ki*error_integral(t) + Kd*error_derivative(t);

    prev_error = error(t);
    prev_error_integral = error_integral(t);

    % Mass flow calculation with PID and main heat loads
    if cooling(t)
        % Cooling: PID should REDUCE mass flow when cabin is too cold
        mf(t) = mf_req(t) - PID_output(t);  % Base + PID adjustment   
    elseif heating(t)
        % Heating: PID should INCREASE mass flow when cabin is too cold
        mf(t) = mf_req(t) + PID_output(t);
    end

    if mf(t) > max_mf
        mf(t) = max_mf;
    elseif mf(t) < min_mf
        mf(t) = min_mf;
    end
    % mf(t) = min(max(mf(t), min_mf), max_mf);  % Clamp to [20,50] g/s for hatchback

    %% 7. Connect to Cabin Thermal Model
    Q_evap(t) = mf(t)*(h1(t)-h4(t));
    Q_cond(t) = mf(t)*(h3(t)-h2(t));
    W_comp(t) = mf(t)*(h2(t)-h1(t)); % Compressor power, J/s

    % MAC energy
    E_comp_Wh(t) = (W_comp(t)+W_comp(t-1))/2*timestep/3600; % Instantaneous energy
    cumE_comp_Wh = cumtrapz(timestep.*W_comp./3600); % Cumulative energy

    % Convert refrigeration outputs (c.v.1) to cabin inputs (control vol 2)
    if cooling(t)
        Q_MAC = -Q_evap(t);  % Heat removed from cabin
    elseif heating(t)
        Q_MAC = Q_cond(t);  % Heat  added to cabin
    end

    % Boundary condition (Tbc)
    Tbc = [T_cell(t);
           -(Qhuman(t)/2 + Qengine(t) + Qvent(t) + Q_MAC/2);
           -(Qhuman(t)/2 + Qexhaust(t) + Qleakage(t) + Q_MAC/2)];
    
    % Solve for temperature
    temperature = (K - C) \ (Tbc - C*prev_temp);
    
    % Update temperatures
    Tcabin_front(t) = temperature(2);
    Tcabin_back(t) = temperature(3);
    Tcabin(t) = (Tcabin_front(t) + Tcabin_back(t))/2;
    prev_temp = temperature;
    
    %% 8. Debug Output
    if mod(t, 100) == 0  % Print every 100 seconds
        fprintf('t=%d s: Tcabin=%.2f°C, Q_MAC= %.1fW, W_compr= %.1fW, mf= %.3f g/s\n', ...
            t,Tcabin(t) - 273.15,Q_MAC,W_comp(t),mf(t) * 1000);
    end

    %% 9. Extract values
    Names = {'P1';'T1';'h1';'T2';'P2';'h2';'T3';'P3';'h3';'T4';'P4';'h4';'COP';'mf';'W_comp';'Q_evap';'Q_cond'};
    MAC_calcs = table(P1',T1',h1',T2',P2',h2',T3',P3',h3',T4',P4',h4',COP',mf(1,1:t)',W_comp',Q_evap',Q_cond','VariableNames',Names);

end
MAC_calcs = MAC_calcs(2:end,:);

% Extract data from MAC_calcs table
h = [MAC_calcs.h1, MAC_calcs.h2, MAC_calcs.h3, MAC_calcs.h4, MAC_calcs.h1]; % kJ/kg (add first point to close cycle)
P = [MAC_calcs.P1, MAC_calcs.P2, MAC_calcs.P3, MAC_calcs.P4, MAC_calcs.P1]; % Convert MPa to Pa

%% Plots MAC components
figure(1);
hold on;
grid on;
grid minor;
plot(time(2:end), MAC_calcs.COP, 'LineWidth', 1);
xlabel('Time [s]');
ylabel('COP');

figure(2);
hold on;
grid on;
grid minor;
plot(time(2:end), MAC_calcs.mf*1000, 'LineWidth', 1); % kg/s to g/s
xlabel('Time [s]');
ylabel('Refrigerant mass flow [g/s]');

figure(3);
hold on;
grid on;
grid minor;
plot(time(2:end), MAC_calcs.W_comp, 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Compressor Work [W]');

figure(4);
hold on;
grid on;
grid minor;
plot(time(2:end), MAC_calcs.Q_evap, 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Evaporator Heat Load [W]');

figure(5);
hold on;
grid on;
grid minor;
plot(time(2:end), MAC_calcs.Q_cond, 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Condenser Heat Load [W]');

figure(6);
hold on;
grid on;
grid minor;
plot(time, Qcabin_req, 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Requested Cabin Heat Load [W]');

% figure(8)
% plot(time,Qcabin_req,':k','LineWidth',1.5)
% hold on
% plot(time,Q_MAC,'-k','LineWidth',1.5)
% grid minor
% ylabel('MAC Heat [W]')
% xlabel('Time [s]')
% yyaxis right
% ax = gca;
% ax.YColor = [1 0 0];
% hold on
% plot(time,Qcabin_tot,'-r','LineWidth',1.5)
% ylabel('Total Cabin Heat [W]','Color',[1 0 0])
% legend('Target Heat: Cabin & evap requested heat)','MAC Heat: Evap given heat','Cabin Total: present heat','Location','best')
% % saveas(figure(12),strcat(test,'_4'),'png')

figure(10);
hold on;
grid on;
grid minor;
plot(time(2:end),Tcabin_front(2:end)-273.15,':b','LineWidth',1.5);
plot(time(2:end),Tcabin_back(2:end)-273.15,':m','LineWidth',1.5);
plot(time(2:end),Tcabin(2:end)-273.15,':k','LineWidth',1.5);
plot(time(2:end),T_cell(2:end)-273.15,'-k','LineWidth',1);
legend('T front simulated','T back simulated','T cabin simulated','TC exterior');
xlabel('Time [s]');
ylabel('Temperature [ºC]');

figure(11);
hold on;
grid on;
grid minor;
plot(time(2:end),comp_speed(2:end),':k','LineWidth',1.5);
plot(time(2:end),100*MAC_calcs.mf,'-k','LineWidth',1.5);
yyaxis right;
plot(time(2:end), MAC_calcs.W_comp,'LineWidth',1);
xlabel('Time [s]');
legend('Compressor speed [rpm]','Mass flow[kg/s]*100','Compressor work (Y)[W]');

figure(13);
hold on;
grid on;
grid minor;
plot(PID_output,'-.','LineWidth',1.5);
yyaxis right;
grid on;
plot(cooling,'-k','LineWidth',1);
plot(heating,'--k','LineWidth',1);
plot(mac_off,'-.r','LineWidth',1);
ylim([-1 2]);
ylabel('Cooling/Heating Flag');
legend('PID output (kg/s)','Cooling', 'Heating','MAC off');

figure(14);
hold on;
grid on;
grid minor;
plot(error,':','LineWidth',1.5);
plot(error_integral,'--','LineWidth',1.5);
plot(error_derivative,'-.','LineWidth',1.5);
legend('proportional error','integral error', 'derivative error')

figure(15);
hold on;
grid on;
plot(time,cumE_comp_Wh,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('MAC Comulative Energy (Wh)');
last_y = cumE_comp_Wh(end);
annotation(figure(15),'textarrow',[0.717857142857142 0.9],...
    [0.890476190476193 0.895238095238096],'String',sprintf('Last value: %.1f', last_y));
grid minor;

 %% p-h diagram
% figure(12);
% hold on;
% grid on;
% set(gca, 'YScale', 'log'); % Logarithmic pressure scale for better visualization
% 
% % Plot the cycle with different colors for each process
% plot([h(1,4) h(1,1)], [P(1,4) P(1,1)], 'b-', 'LineWidth', 1, 'DisplayName', '4→1: Evaporator');
% plot([h(1,1) h(1,2)], [P(1,1) P(1,2)], 'r-', 'LineWidth', 1, 'DisplayName', '1→2: Compressor');
% plot([h(1,2) h(1,3)], [P(1,2) P(1,3)], 'g-', 'LineWidth', 1, 'DisplayName', '2→3: Condenser');
% plot([h(1,3) h(1,4)], [P(1,3) P(1,4)], 'm-', 'LineWidth', 1, 'DisplayName', '3→4: Expansion');
% 
% % Mark the points
% plot(h(1,1:4), P(1,1:4), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
% 
% % Add labels for points
% text(h(1,1), P(1,1), ' 1', 'VerticalAlignment', 'bottom', 'FontSize', 8);
% text(h(1,2), P(1,2), ' 2', 'FontSize', 8);
% text(h(1,3), P(1,3), ' 3', 'FontSize', 8);
% text(h(1,4), P(1,4), ' 4', 'VerticalAlignment', 'bottom', 'FontSize', 8);
% 
% % Add saturation lines (optional - requires CoolProp)
% try
%     % Get saturation properties
%     h_f = py.CoolProp.CoolProp.PropsSI('H', 'Q', 0, 'P', P(1,1), Ref);
%     h_g = py.CoolProp.CoolProp.PropsSI('H', 'Q', 1, 'P', P(1,1), Ref);
%     P_sat = linspace(min(P(:)), max(P(:)), 100);
%     hf_sat = arrayfun(@(p) py.CoolProp.CoolProp.PropsSI('H', 'Q', 0, 'P', p, Ref), P_sat);
%     hg_sat = arrayfun(@(p) py.CoolProp.CoolProp.PropsSI('H', 'Q', 1, 'P', p, Ref), P_sat);
% 
%     plot(hf_sat, P_sat, 'k--', 'DisplayName', 'Saturated Liquid');
%     plot(hg_sat, P_sat, 'k--', 'DisplayName', 'Saturated Vapor');
% catch
%     warning('Could not plot saturation lines - verify CoolProp installation');
% end
% 
% % Formatting
% xlabel('Specific Enthalpy (kJ/kg)');
% ylabel('Pressure (MPa) - Logarithmic Scale');
% title('Refrigeration Cycle on p-h Diagram');
% legend('Location', 'best');
% set(gca);
% hold off;

