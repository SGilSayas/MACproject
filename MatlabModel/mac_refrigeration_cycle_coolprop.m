clc;
clear;
close all;
clear py;

%% Redo link python-matlab
% in matlab terminal: 
% pyenv('Clear');
% pyenv('Version', 'C:\Users\susan\AppData\Python\Python310\python.exe')
% in matlab terminal: pyversion
% if all ok: check this on cmd window: pip install coolprop
% if all ok: check this on matlab terminal py.CoolProp.CoolProp.PropsSI('T','P',101325,'Q',0,'Water')

%% Load Inputs from cabin model
% Timestep series
Total_time = 1800;
timestep = 1;
time = 1:timestep:Total_time;

Tmax = 35 + 273.15; 
Tmin = 30 + 273.15;
temp_step = (Tmax-Tmin)/(Total_time/timestep);
T_amb = Tmin:temp_step:Tmax; % example

Delta_T = 5; % Delta T between the heat exchangers (Evap. + Cond.) and the sources (Cold and Hot)
Q_air_cabin = ones(1,Total_time)*15000; % W, heat requested from the MAC system (evaporator's or condenser's work) --> In cabin model: Q_mac
T_target = 22 + 273.15;

% Load working fluid operation pressures
fluid = {'Water'}; %{'R1234yf'}; %{'R134a'};
R1234yf_op_pres = readtable('R1234yf_operating_pressures.xlsx','VariableNamesRow',1);

% Display the table to verify its contents: disp(R1234yf_op_pres);

%% LOOP for each time step need to start here
t = 1;

% MAC properties
eta_is=0.8; % Standard copressor's isoentropic efficiency, eta_is = (h2s-h1)/(h2-h1)

% Initialization
Q_evap_req = zeros(1,Total_time); % Heat exchanged in the evaporator [W]
Q_cond_req = zeros(1,Total_time); % Heat echanged in the condenser   [W]
W_comp(t) = 0; % Work performed by the compresor [W]
COP(t) = 0; % Coefficient of performance of the HP [-]
mf(t) = 0; % Compressor refrigerant mass flow [kg/s]
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
Q_MAC(t) = 0; % Heat exchanged with cabin Air with the HVAC [W]

for t=2:Total_time+1 % better for for the mf
%while (t <= Total_time)

    %
    % Here the cabin model already done
    %
    %

    % Check if T_amb is in the table
    idx = find(R1234yf_op_pres.T_ambient_C == T_amb(t));
    
    if ~isempty(idx) % If T_amb is found in the table
        P_LowSide_max = 1000*R1234yf_op_pres.P_LowSide_max_kPa(idx);
        P_HighSide_max_Pa = 1000*R1234yf_op_pres.P_HighSide_max_kPa(idx);
    else % If T_amb is not found, interpolate
        T_ambient = R1234yf_op_pres.T_ambient_C;
        P_LowSide_max_Pa = 1000*R1234yf_op_pres.P_LowSide_max_kPa;
        P_HighSide_max_Pa = 1000*R1234yf_op_pres.P_HighSide_max_kPa;
        
        P_LowSide_max = interp1(T_ambient, P_LowSide_max_Pa, T_amb, 'linear', 'extrap');
        P_HighSide_max = interp1(T_ambient, P_HighSide_max_Pa, T_amb, 'linear', 'extrap');
    end
    
    Comp_ratio_min = 5.7357; % " Average"

    % Define cooling or heating mode
    if (mean(Q_MAC)>= 0)
       cooling = true; % Cooling
       Q_evap_req(t) = abs(Q_air_cabin(t-1));
    else
       cooling = false; % Heating
       Q_cond_req(t) = abs(Q_air_cabin(t-1));
    end
    
    % Define evaporator's and condenser's temperatures depending on MAC mode
    if (cooling == true) % Cooling
        % Heat exchangers
        T_evap(t) = T_target - Delta_T ;
        T_cond(t) = T_amb(t) + Delta_T;
    else % Heating
        T_evap(t) = T_ext(t) - Delta_T ;
        T_cond(t) = Temp_conf + Delta_T;
    end
    
    % Code Diagram p-h 
    Ref = char(fluid);
    
    %% Refrigeration Cycle Points:
    % (1): evaporator outlet, salutared vapor; 
    % (2): compressor outlet, superheated vapor;
    % (3): condenser outlet, saturated fluid/liquid; 
    % (4): expansion device outlet, vapor-fluid mix;
    
    % (1) Evaporator's outlet = Compressor's inlet, Saturated vapor Q=1
    T1(t)= T_evap(t); % T_evap_out is known
    P1(t)=py.CoolProp.CoolProp.PropsSI('P','T',T1(t),'Q',1,Ref);
    P1(t) = double(P1(t));
    if (P1(t) >= P_LowSide_max)
        P1(t) = P_LowSide_max;
        T1(t) = py.CoolProp.CoolProp.PropsSI('T','P',P1(t),'Q',1,Ref);
    end
    h1(t)=py.CoolProp.CoolProp.PropsSI('H','T',T1(t),'Q',1,Ref);     
    s1(t)=py.CoolProp.CoolProp.PropsSI('S','T',T1(t),'Q',1,Ref);
    
    % (2) Compressor's outlet = Condenser's inlet
    T2sat(t) = T_cond(t);  % T_cond_in is known
    P2(t)=py.CoolProp.CoolProp.PropsSI('P','T',T2sat(t),'Q',1,Ref);    
    if (P2(t) <= Comp_ratio_min*P1(t))
        P2(t) = Comp_ratio_min*P1(t);
    end
    beta(t)=P2(t)/P1(t);   
    s2is(t)=s1(t); 
    h2is(t)=py.CoolProp.CoolProp.PropsSI('H','P',P2(t),'S',s2is(t),Ref);    
    h2(t)=(h2is(t)-h1(t))/eta_is + h1(t);
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
    h4(t)=py.CoolProp.CoolProp.PropsSI('H','S',s4(t),'P',P1(t),Ref);
    T4(t)=py.CoolProp.CoolProp.PropsSI('T','S',s4(t),'P',P1(t),Ref);
    Q4(t)=py.CoolProp.CoolProp.PropsSI('Q','S',s4(t),'P',P1(t),Ref);
    P4(t)= P1(t);

    % Coefficient of performance
    COP(t)=(h2(t)-h3(t))/(h2(t)-h1(t));

    % Mass flow calculation
    if (cooling == true)
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

    % Components' heat loads calculation
    Q_evap(t) = mf(t)*(h1(t) - h4(t));
    Q_cond(t) = mf(t)*(h3(t) - h2(t));
    W_comp(t) = mf(t)*(h2(t) - h1(t));

    if (cooling == true)
        Q_MAC(t) = Q_evap(t);
    else
        Q_MAC(t) = Q_cond(t);  
    end

    % Extract values
    Names = {'P1';'T1';'T2';'P2';'T3';'P3';'T4';'P4';'COP';'mf';'W_comp';'Q_evap';'Q_cond'};
    MAC_calcs = table(P1',T1',T2',P2',T3',P3',T4',P4',COP',mf',W_comp',Q_evap',Q_cond','VariableNames',Names);

    % GENERAL EQUATION
    % m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 
    % Tcabin(t) = (Q_windows(t)+Q_ceiling(t)+Qbase(t)+Qhuman(t) - Q_HVAC(t) )*(timestep/(V_cabin*cp_air*rho_air)) + Tcabin(t-1);

    % Add a timestep
    t = t + timestep;

end
MAC_calcs = MAC_calcs(2:end,:);

% Plots
figure(1)
plot(time, MAC_calcs.COP, 'LineWidth', 1);
xlabel('Time (s)')
ylabel('COP')
grid on

figure(2)
plot(time, MAC_calcs.mf*1000, 'LineWidth', 1);
xlabel('Time (s)')
ylabel('Refrigerant mass flow, g/s')
grid on

figure(3)
plot(time, MAC_calcs.W_comp*1000, 'LineWidth', 1);
xlabel('Time (s)')
ylabel('Compressor Work (kW)')
grid on

figure(4)
plot(time, MAC_calcs.Q_evap*1000, 'LineWidth', 1);
xlabel('Time (s)')
ylabel('Evaporator Heat Load (kW)')
grid on

figure(5)
plot(time, MAC_calcs.Q_cond*1000, 'LineWidth', 1);
xlabel('Time (s)')
ylabel('Condenser Heat Load (kW)')
grid on

figure(6)
plot(time, Q_air_cabin*1000, 'LineWidth', 1);
xlabel('Time (s)')
ylabel('Cabin Heat Load (kW)')
grid on