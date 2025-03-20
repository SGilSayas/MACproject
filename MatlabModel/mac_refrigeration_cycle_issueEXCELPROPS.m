clc;
clear all;
close all;

%% Define variable names
var_names = {'T_C', 'p_liq_kPa', 'p_vap_kPa', 'v_liq_m3_kg', 'v_vap_m3_kg', 'rho_liq_kg_m3', 'rho_vap_kg_m3', 'H_liq_kJ_kg', 'H_latent_kJ_kg', 'H_vap_kJ_kg', 's_liq_kJ_kgK', 's_vap_kJ_kgK'};

%% Load data from Excel file
data = readtable('R1234yf_properties.xlsx', 'HeaderLines', 2); %'VariableNamesRow',1);
data.Properties.VariableNames = var_names;

%% Convert data to struct array
data_struct = struct();
for i = 1:length(var_names)
    data_struct.(var_names{i}) = data.(var_names{i});
end

%% Define inputs
evaporator_work = 1000;  % W
evaporator_temperature = 10;  % Celsius

%% Convert temperature to Kelvin
evaporator_temperature_K = evaporator_temperature + 273.15;

%% Find saturated liquid and vapor properties at evaporator temperature
temp_idx = find(data_struct.T_C == evaporator_temperature);
if isempty(temp_idx)
    % Interpolate between closer values
    idx = find(data_struct.T_C < evaporator_temperature, 1, 'last');
    idx2 = find(data_struct.T_C > evaporator_temperature, 1, 'first');
    if isempty(idx)
        idx = 1;
    end
    if isempty(idx2)
        idx2 = length(data_struct.T_C);
    end
    p_liq = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.p_liq_kPa(idx) data_struct.p_liq_kPa(idx2)], evaporator_temperature);
    p_vap = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.p_vap_kPa(idx) data_struct.p_vap_kPa(idx2)], evaporator_temperature);
    v_liq = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.v_liq_m3_kg(idx) data_struct.v_liq_m3_kg(idx2)], evaporator_temperature);
    v_vap = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.v_vap_m3_kg(idx) data_struct.v_vap_m3_kg(idx2)], evaporator_temperature);
    rho_liq = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.rho_liq_kg_m3(idx) data_struct.rho_liq_kg_m3(idx2)], evaporator_temperature);
    rho_vap = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.rho_vap_kg_m3(idx) data_struct.rho_vap_kg_m3(idx2)], evaporator_temperature);
    H_liq = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.H_liq_kJ_kg(idx) data_struct.H_liq_kJ_kg(idx2)], evaporator_temperature);
    H_latent = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.H_latent_kJ_kg(idx) data_struct.H_latent_kJ_kg(idx2)], evaporator_temperature);
    H_vap = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.H_vap_kJ_kg(idx) data_struct.H_vap_kJ_kg(idx2)], evaporator_temperature);
    s_liq = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.s_liq_kJ_kgK(idx) data_struct.s_liq_kJ_kgK(idx2)], evaporator_temperature);
    s_vap = interp1([data_struct.T_C(idx) data_struct.T_C(idx2)], [data_struct.s_vap_kJ_kgK(idx) data_struct.s_vap_kJ_kgK(idx2)], evaporator_temperature);
else
    p_liq = data_struct.p_liq_kPa(temp_idx);
    p_vap = data_struct.p_vap_kPa(temp_idx);
    v_liq = data_struct.v_liq_m3_kg(temp_idx);
    v_vap = data_struct.v_vap_m3_kg(temp_idx);
    rho_liq = data_struct.rho_liq_kg_m3(temp_idx);
    rho_vap = data_struct.rho_vap_kg_m3(temp_idx);
    H_liq = data_struct.H_liq_kJ_kg(temp_idx);
    H_latent = data_struct.H_latent_kJ_kg(temp_idx);
    H_vap = data_struct.H_vap_kJ_kg(temp_idx);
    s_liq = data_struct.s_liq_kJ_kgK(temp_idx);
    s_vap = data_struct.s_vap_kJ_kgK(temp_idx);
end

%% Calculate compressor work
s_comp = s_vap;  % Constant entropy
H_comp = H_vap;  % Initial enthalpy
p_comp = p_vap;  % Initial pressure

% Find the temperature at which the entropy is equal to s_comp
% and the pressure is higher than p_comp (condenser pressure)
condenser_pressure_idx = find(data_struct.p_vap_kPa > p_comp & data_struct.s_vap_kJ_kgK == s_comp);
if isempty(condenser_pressure_idx)
    % Interpolate between closer values
    idx = find(data_struct.s_vap_kJ_kgK < s_comp, 1, 'last');
    idx2 = find(data_struct.s_vap_kJ_kgK > s_comp, 1, 'first');
    if isempty(idx)
        idx = 1;
    end
    if isempty(idx2)
        idx2 = length(data_struct.s_vap_kJ_kgK);
    end
    H_cond = interp1([data_struct.s_vap_kJ_kgK(idx) data_struct.s_vap_kJ_kgK(idx2)], [data_struct.H_vap_kJ_kg(idx) data_struct.H_vap_kJ_kg(idx2)], s_comp);
    p_cond = interp1([data_struct.s_vap_kJ_kgK(idx) data_struct.s_vap_kJ_kgK(idx2)], [data_struct.p_vap_kPa(idx) data_struct.p_vap_kPa(idx2)], s_comp);
else
    H_cond = data_struct.H_vap_kJ_kg(condenser_pressure_idx);
    p_cond = data_struct.p_vap_kPa(condenser_pressure_idx);
end

%% Calculate condenser heat transfer
Q_cond = H_vap - H_cond;

%% Calculate expansion valve enthalpy
H_exp = H_cond;  % Constant enthalpy

% Find the pressure at which the enthalpy is equal to H_exp
expansion_pressure_idx = find(data_struct.H_vap_kJ_kg == H_exp);
if isempty(expansion_pressure_idx)
    % Interpolate between closer values
    idx = find(data_struct.H_vap_kJ_kg < H_exp, 1, 'last');
    idx2 = find(data_struct.H_vap_kJ_kg > H_exp, 1, 'first');
    if isempty(idx)
        idx = 1;
    end
    if isempty(idx2)
        idx2 = length(data_struct.H_vap_kJ_kg);
    end
    p_exp = interp1([data_struct.H_vap_kJ_kg(idx) data_struct.H_vap_kJ_kg(idx2)], [data_struct.p_vap_kPa(idx) data_struct.p_vap_kPa(idx2)], H_exp);
else
    p_exp = data_struct.p_vap_kPa(expansion_pressure_idx);
end

%% Plot refrigeration cycle
figure;
plot(data_struct.H_liq_kJ_kg, data_struct.p_liq_kPa, 'b-');  % Saturated liquid line
hold on;
plot(data_struct.H_vap_kJ_kg, data_struct.p_vap_kPa, 'r-');  % Saturated vapor line
plot([H_liq H_vap], [p_vap p_vap], 'g-');  % Evaporator line
plot([H_vap H_cond], [p_vap p_cond], 'b-');  % Compressor line
plot([H_cond H_cond], [0 p_cond], 'k-');  % Condenser line
plot([H_cond H_exp], [p_cond p_exp], '-');  % Expansion valve line
xlabel('Enthalpy (kJ/kg)');
ylabel('Pressure (kPa)');
legend('Saturated liquid', 'Saturated vapor', 'Evaporator', 'Compressor', 'Condenser', 'Expansion valve');

%% Display compressor work
W_comp = H_vap - H_cond;
fprintf('Compressor work: %.2f kJ/kg\n', W_comp);