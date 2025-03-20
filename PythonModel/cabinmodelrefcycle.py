import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import scipy.io as sio

# Test conditions
#t = 1
timestep = 1
duration_WLTC = 1800  # s
Total_time = duration_WLTC
time = np.arange(0, Total_time, timestep)

N_Humans = 1

fuel = 'PETROL'  # or 'DIESEL'

# Cabin target temperature
T_target = 22 + 273.15  # K

# Standard compressor's isoentropic efficiency, eta_is = (h2s-h1)/(h2-h1)
eta_is = 0.8  # av: 0.8;
compressor_on_lower_limit = 600  # W

# T between HX-sources:
# heat exchangers (evap&cond) and cold&hot sources (cabin&external air)
delta_T_evap = 7  # Automotive: 7-12  ; HVAC: 5-10
delta_T_cond = 10  # Automotive: 10-15 ; HVAC: 5-15

# Ambient conditions
Irr_scalar = 0  # W/m2
Irr = np.ones((Total_time + 1,1)) * Irr_scalar

# Load ambient temperature
# T_amb = np.load('G8_35_cell_front_back.npy') + 273.15  # K
#T_amb = np.ones((Total_time-1,1))
T_amb = sio.loadmat('G8_35_cell_front_back.mat')['TC_cell'] + 273.15  # K
#T_amb = np.squeeze(T_amb)
#T_amb = np.append(T_amb, T_amb[-1])
#T_amb[Total_time,0] = T_amb[-1,0]
#T_amb[Total_time+1,0] = T_amb[Total_time-1,0]

T_amb = np.pad(T_amb, (0, 2), mode='edge')

T_ini = T_amb[0]  # alternative: T_ini = TC_cabin_ini + 273.15;

# Select and load working fluid operation pressures
fluid = 'R1234yf'
R1234yf_op_pres = pd.read_excel('R1234yf_operating_pressures.xlsx')

# Constant input data

# Load input cell humidity (RH) and pressure (P)
RH_amb_av = np.mean(sio.loadmat('G7_35_off_RH_P_amb.mat')['RH_amb'])
P_amb_kPa_av = np.mean(sio.loadmat('G7_35_off_RH_P_amb.mat')['P_amb_kPa'])

# Load input Engine and exhaust heat flows
Qengine_av = np.mean(sio.loadmat('heatflow_engine_exhaust_W.mat')['Qengine'])
Qexhaust_av = np.mean(sio.loadmat('heatflow_engine_exhaust_W.mat')['Qexhaust'])

# TO DO: make vector chape with same av value:
RH_amb_av = np.ones((Total_time + 1,1)) * RH_amb_av
P_amb_kPa_av = np.ones((Total_time + 1,1)) * P_amb_kPa_av
Qengine_av = np.ones((Total_time + 1,1)) * Qengine_av
Qexhaust_av = np.ones((Total_time + 1,1)) * Qexhaust_av

# Cabin dimensions for category C
vehicle_height = 1.46
vehicle_width = 1.8
vehicle_lenght = 4.36
cat_factor = 1.1
cabin_height = vehicle_height - 0.06 - 0.6215 / 2  # minus roof&base thickness and half of wheel
cabin_width = vehicle_width - 0.06 - 0.241 / 2  # minus doors thickness % mirrors=0.241, we use half bc some database width included mirrors and others did not.
cabin_roof_lenght = vehicle_lenght * 0.43
cabin_base_lenght = vehicle_lenght * 0.35
seatcm3 = 150 * cat_factor  # 115; % cm3
seats_volume = 5 * seatcm3 * 10 ** (-6)  # m3
V_cabin = (cabin_width * cabin_roof_lenght * cabin_height / 2) + (
            cabin_width * cabin_base_lenght * cabin_height / 2) - seats_volume
V_cabin = V_cabin * 1

# A=area[m2]
A_ws = 0.63 * 1.3 * cat_factor
A_rw = 0.29 * 1 * cat_factor
A_roof = 1.8 * 1.1 * cat_factor
A_sidewindows = 2 * 1.45 * 0.29 * cat_factor
A_doors = 2 * 1.45 * 0.29 * cat_factor
A_base = cabin_base_lenght * cabin_width * 3
A_dashboard = cabin_width * 0.5  # m2
A_front = A_ws + A_dashboard  # 3
A_back = A_rw  # 1.5
A_side = A_sidewindows + A_doors + A_roof
# A_front = A_ws + A_sidewindows/2 + A_doors/2 + A_roof/2
# A_back = A_rw + A_sidewindows/2 + A_doors/2 + A_roof/2

# e=thickness[m]
e_ws = 0.006
e_rw = 0.005
e_ABS = 0.002
e_steel = 0.0005
e_PU = 0.025  # polyurethane
e_ceiling = e_ABS + e_steel + e_PU
e_sidewindows = 0.003
e_doors = e_ceiling

# Fuel
Vpe_petrol = 0.264  # l/kWh
CF_petrol = 2330  # gCO2/l
Vpe_diesel = 0.22  # l/kWh
CF_diesel = 2640  # gCO2/l
if fuel == 'PETROL' or fuel == 'PETROL/ELECTRIC':
    Vpe = Vpe_petrol
    CF = CF_petrol
elif fuel == 'DIESEL' or fuel == 'DIESEL/ELECTRIC':
    Vpe = Vpe_diesel
    CF = CF_diesel
else:
    print('Specify PETROL or DIESEL as the last argument')

# Conduction Properties: k = thermoconductivity[W/(m*K)]
k_ws = 0.8  # 0.8-1.4
k_rw = 1.4  # 1.4-1.5
k_ABS = 0.1  # Acrylonitrile butadiene styrene
k_steel = 14.9
k_PU = 0.022  # 0.022-0.028 W/mK polyurethane
k_ceiling = k_ABS * k_steel * k_PU
k_sidewindows = 1.4
k_doors = k_ceiling

# Heat transfer coefficient of interior cabin air, W/(K*m2)
h_cabin = 5  # 5-10

# Air properies
density_air = 1.18  # [Kg/m3] (at 25C)
cp_air = 1006  # + 10*1500; % [J/ kg*K]

# Seats properties
mass_seats = 20  # 20-30 kg
cp_seats = 2000  # 500-2000 [J/kgK]

# Human and equipment heat % ISO 8996, Alberto Viti Corsi (IDAE)
heat_equipment = 40  # W
Qequipment = np.ones((Total_time + 1,1)) * heat_equipment  # W
A_skin = 1.5  # m2
T_skin = 24.648 + 273.15  # K

# Radiation Properties
sigma = 0.0000000567037321  # Stefan-Boltzmann constant(sigma)[W/(m^2*K^4)]
# epsilon:emissivity[]
epsilon_ws = 0.9
epsilon_rw = 0.9
epsilon_sidewindows = 0.9
epsilon_doors = 0.9
epsilon_roof = 0.9
# rho=reflectivity[]
rho_ws = 0.246
rho_rw = 0.1
rho_sidewindows = 0.2
rho_doors = 0.74
rho_roof = 0.04  # G7:0.74;
rho_base = 0.3
# tao=transmissivity[]
tao_ws = 0.452
tao_rw = 0.311
tao_sidewindows = 0.475
tao_doors = 0
tao_roof = 0.9  # 0;
tao_base = 0
# alpha=absorptivity[]
alpha_coef = 1
alpha_ws = alpha_coef * (1 - rho_ws - tao_ws)
alpha_rw = alpha_coef * (1 - rho_rw - tao_rw)
alpha_sidewindows = alpha_coef * (1 - rho_sidewindows - tao_sidewindows)
alpha_doors = alpha_coef * (1 - rho_doors - tao_doors)
alpha_roof = alpha_coef * (1 - rho_roof - tao_roof)
alpha_base = alpha_coef * (1 - rho_base - tao_base)

# Global Horizontal Irradiance (Hypothesis: same in all windows)
G_roof_inc = Irr
G_ws_inc = Irr
G_rw_inc = Irr
G_sidewindows_inc = Irr
G_doors_inc = Irr
G_ws_r = rho_ws * G_ws_inc
G_ws_a = alpha_ws * G_ws_inc  # bc
G_ws_t = tao_ws * G_ws_inc
G_rw_r = rho_rw * G_rw_inc
G_rw_a = alpha_rw * G_rw_inc
G_rw_t = tao_rw * G_rw_inc
G_sidewindows_r = rho_sidewindows * G_sidewindows_inc
G_sidewindows_a = alpha_sidewindows * G_sidewindows_inc
G_sidewindows_t = tao_sidewindows * G_sidewindows_inc
G_doors_r = rho_doors * G_doors_inc
G_doors_a = alpha_doors * G_doors_inc
G_doors_t = tao_doors * G_doors_inc
G_roof_a = alpha_roof * G_roof_inc

# Air vent volume and flow rate
vent_volumerate = 0.0003  # m3/s, Faya 2013: 0.02; 0.001
leakage_volumerate = 0.0001  # m3/s, Faya2014

# Average humidity ratio in gram of water per gram of dry air, X
Ps_22_kPa = 2.626  # kPa, Water saturation pressure at 22C
Ps_35_kPa = 5.63  # kPa, Water saturation pressure at 35C
Ps_kPa = Ps_22_kPa + (T_ini - 22) / (35 - 22) * (Ps_35_kPa - Ps_22_kPa)
X = 0.0732  # X = 0.62198*Ps_kPa.*RH_amb./(100.*P_amb_kPa - Ps_kPa.*RH_amb); % Faya, 2013

# Calculation of Properties

# Enthalpies calculation. Humidity considered same in cabin and exterior
e_amb = 1006 * T_amb + (2501000 + 1770 * T_amb) * X  # J/kg, Faya 2013, % before ".*X"

# Capacitance Matrix
C_amb = 0  # boundary condition
C_roof = 0  # for now
C_ceiling = 0  # for now
C_base = mass_seats * cp_seats + 144240  # Capacitance, W/K
C_ws = 0
C_rw = 0
C_sidewindows = 0
C_doors = 0
C_cabin_front = density_air * cp_air / timestep * V_cabin * 0.8  # front air volume = 80% of total volume % other hypothesis: 0.9;
C_cabin_back = density_air * cp_air / timestep * V_cabin * 0.2  # back air volume = 20% of total cabin volume % other hypothesis: 0.001;
C = np.zeros((14, 14))
C[0, 0] = C_amb
C[1, 1] = C_roof
C[2, 2] = C_ceiling
C[3, 3] = C_base
C[4, 4] = C_ws
C[5, 5] = C_rw
C[6, 6] = C_sidewindows
C[7, 7] = C_doors
C[8, 8] = C_ws
C[9, 9] = C_rw
C[10, 10] = C_sidewindows
C[11, 11] = C_doors
C[12, 12] = C_cabin_front
C[13, 13] = C_cabin_back

# Initialization of values
# Initial heat transfer coefficients, UA (W/K):
h_front = np.ones((Total_time + 1,1)) * 160
h_side = np.ones((Total_time + 1,1)) * 0.9 * 160
h_rear = np.ones((Total_time + 1,1)) * 0.1 * 160
UA_front = np.ones((Total_time + 1,1)) * (A_front * h_front[0,0] + 1 * (A_side * h_side[0,0]))
UA_back = np.ones((Total_time + 1,1)) * (A_back * h_rear[0,0] + 1 * (A_side * h_side[0,0]))

#  Conductances Matrix
K1 = h_side[0,0] * A_roof
K2 = A_roof * k_ceiling / e_ceiling
K3 = h_cabin * A_roof
# Tbase --K4-> Tcabin_back
K4 = h_cabin * A_base
# Tamb --K5-> Tw_ext_ws --K6-> Tw_int_ws --K7-> Tcabin_front
K5 = h_front[0,0] * A_ws
K6 = A_ws * k_ws / e_ws
K7 = h_cabin * A_ws
# Tamb --K8-> Tw_ext_rw --K9-> Tw_int_rw --K10-> Tcabin_back
K8 = h_rear[0,0] * A_rw
K9 = A_rw * k_rw / e_rw
K10 = h_cabin * A_rw
# Tamb --K11-> Tw_ext_sidewindows --K12-> Tw_int_sidewindows --K13-> Tcabin_back
K11 = h_side[0,0] * A_sidewindows
K12 = A_sidewindows * k_sidewindows / e_sidewindows
K13 = h_cabin * A_sidewindows
# Tamb --K14-> Tw_ext_doors --K15-> Tw_int_doors --K16-> Tcabin_back
K14 = h_side[0,0] * A_doors
K15 = A_doors * k_doors / e_doors
K16 = h_cabin * A_doors
# Tcabin_front --K17-> Tcabin_back
K17 = h_cabin * cabin_width * cabin_height  # A_ws; % TBD
K = np.zeros((14, 14))
K[0, 0] = 1
K[1, 0] = K1
K[1, 1] = - (K1 + K2)
K[1, 2] = K2
K[2, 1] = K2
K[2, 2] = - (K2 + K3 / 2 + K3 / 2)
K[2, 12] = K3 / 2
K[2, 13] = K3 / 2
K[3, 2] = K4 / 2
K[3, 3] = - (K4 / 2 + K4 / 2)
K[3, 12] = K4 / 2
K[3, 13] = K4 / 2
K[4, 0] = K5
K[4, 4] = - (K5 + K6)
K[4, 8] = K6
K[5, 0] = K8
K[5, 5] = - (K8 + K9)
K[5, 9] = K9
K[6, 0] = K11
K[6, 6] = - (K11 + K12)
K[6, 10] = K12
K[7, 0] = K14
K[7, 7] = - (K14 + K15)
K[7, 11] = K15
K[8, 4] = K6
K[8, 8] = - (K6 + K7)
K[8, 12] = K7
K[9, 5] = K9
K[9, 9] = - (K9 + K10)
K[9, 13] = K10
K[10, 6] = K12
K[10, 10] = - (K12 + K13 / 2 + K13 / 2)
K[10, 12] = K13 / 2
K[10, 13] = K13 / 2
K[11, 7] = K15
K[11, 11] = - (K15 + K16 / 2 + K16 / 2)
K[11, 12] = K16 / 2
K[11, 13] = K16 / 2
K[12, 2] = K3 / 2
K[12, 3] = K4 / 2
K[12, 8] = K7
K[12, 10] = K13 / 2
K[12, 11] = K16 / 2
K[12, 12] = - (K3 / 2 + K4 / 2 + K7 + K16 / 2 + K17)
K[12, 13] = K17
K[13, 2] = K3 / 2
K[13, 3] = K4 / 2
K[13, 9] = K10
K[13, 10] = K13 / 2
K[13, 11] = K16 / 2
K[13, 12] = K17
K[13, 13] = - (K3 / 2 + K4 / 2 + K10 + K13 / 2 + K16 / 2 + K17)
t=0
# initial cabin surfaces tmperatures
Troof = np.ones((Total_time + 1,1)) * T_ini
Tceiling = np.ones((Total_time + 1,1)) * T_ini
Tbase_int = np.ones((Total_time + 1,1)) * T_ini
Tw_ext_ws = np.ones((Total_time + 1,1)) * T_ini
Tw_ext_rw = np.ones((Total_time + 1,1)) * T_ini
Tw_ext_sidewindows = np.ones((Total_time + 1,1)) * T_ini
Tw_ext_doors = np.ones((Total_time + 1,1)) * T_ini
Tw_int_ws = np.ones((Total_time + 1,1)) * T_ini
Tw_int_rw = np.ones((Total_time + 1,1)) * T_ini
Tw_int_sidewindows = np.ones((Total_time + 1,1)) * T_ini
Tw_int_doors = np.ones((Total_time + 1,1)) * T_ini
Tcabin_front = np.ones((Total_time + 1,1)) * T_ini
Tcabin_back = np.ones((Total_time + 1,1)) * T_ini
Tcabin = (Tcabin_front + Tcabin_back) / 2
temperature = np.ones((14,1)) * T_ini

prev_temp = temperature.copy()

# Refrigerant cycle:
# point 1: evaporator oulet = compressor inlet
# point 2: compressor outlet = condenser inlet
# point 3: condenser outlet = valve inlet
# point 4: valve outlet = evaporator inlet
# initial parameters of refrigerant cycle followng p-h diagram
T1 = np.zeros((Total_time + 1,1))
P1 = np.zeros((Total_time + 1,1))
h1 = np.zeros((Total_time + 1,1))
s1 = np.zeros((Total_time + 1,1))
T2 = np.zeros((Total_time + 1,1))
P2 = np.zeros((Total_time + 1,1))
h2 = np.zeros((Total_time + 1,1))
s2 = np.zeros((Total_time + 1,1))
P2is = np.zeros((Total_time + 1,1))
h2is = np.zeros((Total_time + 1,1))
s2is = np.zeros((Total_time + 1,1))
T3 = np.zeros((Total_time + 1,1))
P3 = np.zeros((Total_time + 1,1))
h3 = np.zeros((Total_time + 1,1))
s3 = np.zeros((Total_time + 1,1))
T4 = np.zeros((Total_time + 1,1))
P4 = np.zeros((Total_time + 1,1))
h4 = np.zeros((Total_time + 1,1))
s4 = np.zeros((Total_time + 1,1))
COP = np.zeros((Total_time + 1,1))

# initialization of heat flows
Qcv_received = np.ones((Total_time + 1,1))
Qleakage = np.ones((Total_time + 1,1))
Qvent = np.ones((Total_time + 1,1))
Qhuman = np.ones((Total_time + 1,1))
Qbase = np.ones((Total_time + 1,1))
Qcabin_req = np.ones((Total_time + 1,1))
Qcabin_received = np.ones((Total_time + 1,1))
Qcv_emitted = np.zeros((Total_time + 1,1))
Qcabin_tot = np.zeros((Total_time + 1,1))
Q_out_front = np.ones((Total_time + 1,1))
Q_out_back = np.ones((Total_time + 1,1))
Qtarget = np.ones((Total_time + 1,1))
comp_cum_Wh = np.ones((Total_time + 1,1))
Qcompressor = np.zeros((Total_time + 1,1))
heat_human = np.ones((Total_time + 1,1))
Q_evap_req = np.zeros((Total_time + 1,1))
Q_cond_req = np.zeros((Total_time + 1,1))
W_comp = np.zeros((Total_time + 1,1))  # Work performed by the compressor [W]
Q_evap = np.zeros((Total_time + 1,1))  # Heat exchanged in the evaporator [W]
Q_cond = np.zeros((Total_time + 1,1))  # Heat echanged in the condenser   [W]

e_cabin = np.zeros((Total_time + 1,1))
error_sum = 0
mf = np.zeros((Total_time + 1,1))

# CALCULATION LOOP
for t in range(1, Total_time + 1):
    print(f"loop index: {t}")
    # Overall heat transfer coefficients, UA W/(Km2):
    if t > 1800 - 323:  # extra high
        h_front[t,0] = 200
    elif t > 1800 - 323 - 455:  # high
        h_front[t,0] = 190
    elif t > 1800 - 323 - 455 - 433:  # medium
        h_front[t,0] = 170
    else:  # low
        h_front[t,0] = 160
    h_side[t,0] = 0.9 * h_front[t,0]
    h_rear[t,0] = 0.1 * h_front[t,0]
    UA_front[t,0] = A_front * h_front[t,0] + 1 * (A_side * h_side[t,0])
    UA_back[t,0] = A_back * h_rear[t,0] + 1 * (A_side * h_side[t,0])

    # Conductances Matrix
    K1 = h_side[t,0] * A_roof
    K2 = A_roof * k_ceiling / e_ceiling
    K3 = h_cabin * A_roof
    # Tbase --K4-> Tcabin_back
    K4 = h_cabin * A_base
    # Tamb --K5-> Tw_ext_ws --K6-> Tw_int_ws --K7-> Tcabin_front
    K5 = h_front[t,0] * A_ws
    K6 = A_ws * k_ws / e_ws
    K7 = h_cabin * A_ws
    # Tamb --K8-> Tw_ext_rw --K9-> Tw_int_rw --K10-> Tcabin_back
    K8 = h_rear[t,0] * A_rw
    K9 = A_rw * k_rw / e_rw
    K10 = h_cabin * A_rw
    # Tamb --K11-> Tw_ext_sidewindows --K12-> Tw_int_sidewindows --K13-> Tcabin_back
    K11 = h_side[t,0] * A_sidewindows
    K12 = A_sidewindows * k_sidewindows / e_sidewindows
    K13 = h_cabin * A_sidewindows
    # Tamb --K14-> Tw_ext_doors --K15-> Tw_int_doors --K16-> Tcabin_back
    K14 = h_side[t,0] * A_doors
    K15 = A_doors * k_doors / e_doors
    K16 = h_cabin * A_doors
    # Tcabin_front --K17-> Tcabin_back
    K17 = h_cabin * cabin_width * cabin_height  # A_ws; % TBD
    K = np.zeros((14, 14))
    K[0, 0] = 1
    K[1, 0] = K1
    K[1, 1] = - (K1 + K2)
    K[1, 2] = K2
    K[2, 1] = K2
    K[2, 2] = - (K2 + K3/2 + K3/2)
    K[2, 12] = K3/2
    K[2, 13] = K3/2
    K[3, 2] = K4/2
    K[3, 3] = - (K4/2 + K4/2)
    K[3, 12] = K4/2
    K[3, 13] = K4/2
    K[4, 0] = K5
    K[4, 4] = - (K5 + K6)
    K[4, 8] = K6
    K[5, 0] = K8
    K[5, 5] = - (K8 + K9)
    K[5, 9] = K9
    K[6, 0] = K11
    K[6, 6] = - (K11 + K12)
    K[6, 10] = K12
    K[7, 0] = K14
    K[7, 7] = - (K14 + K15)
    K[7, 11] = K15
    K[8, 4] = K6
    K[8, 8] = - (K6 + K7)
    K[8, 12] = K7
    K[9, 5] = K9
    K[9, 9] = - (K9 + K10)
    K[9, 13] = K10
    K[10, 6] = K12
    K[10, 10] = - (K12 + K13/2 + K13/2)
    K[10, 12] = K13/2
    K[10, 13] = K13/2
    K[11, 7] = K15
    K[11, 11] = - (K15 + K16/2 + K16/2)
    K[11, 12] = K16/2
    K[11, 13] = K16/2
    K[12, 2] = K3/2
    K[12, 3] = K4/2
    K[12, 8] = K7
    K[12, 10] = K13/2
    K[12, 11] = K16/2
    K[12, 12] = - (K3/2 + K4/2 + K7 + K16/2 + K17)
    K[12, 13] = K17
    K[13, 2] = K3/2
    K[13, 3] = K4/2
    K[13, 9] = K10
    K[13, 10] = K13/2
    K[13, 11] = K16/2
    K[13, 12] = K17
    K[13, 13] = - (K3/2 + K4/2 + K10 + K13/2 + K16/2 + K17)

    # Heat Flows
    e_cabin[t,0] = 1006 * Tcabin[t,0] + (2501000 + 1770 * Tcabin[t,0]) * X  # J/kg, Faya 2013 % Cabin air enthalpy
    Qleakage[t,0] = leakage_volumerate * density_air * (e_amb[t,0] - e_cabin[t,0])  # W, J/s
    Qvent[t,0] = vent_volumerate * density_air * (e_amb[t,0] - e_cabin[t,0])  # W, J/s
    Qbase[t,0] = alpha_base * (A_ws * G_ws_t[t,0] + A_rw * G_rw_t[t,0] + A_sidewindows * G_sidewindows_t[t,0] + A_sidewindows * G_doors_t[t,0])
    Qhuman[t,0] = N_Humans * h_cabin * A_skin * abs(temperature[12, 0] - T_skin)
    Qequipment[t,0] = heat_equipment

    # Boundary conditions vector for scalar or vector T amb inputs
    Tbc = np.zeros((14,1))
    Tbc[0,0] = T_amb[t,0]
    Tbc[1,0] = - A_roof * G_roof_a[t,0]
    Tbc[2,0] = 0
    Tbc[3,0] = - Qbase[t,0]
    Tbc[4,0] = - A_ws * G_ws_a[t,0]
    Tbc[5,0] = - A_rw * G_rw_a[t,0]
    Tbc[6,0] = - A_sidewindows * G_sidewindows_a[t,0]
    Tbc[7,0] = - A_doors * G_doors_a[t,0]
    Tbc[8,0] = 0
    Tbc[9,0] = 0
    Tbc[10,0] = 0
    Tbc[11,0] = 0
    Tbc[12,0] = - Qhuman[t,0] - Qengine_av[t,0] - Qvent[0,0] + Q_evap_req[t-0,0]/2 + Qcv_emitted[t-0,0]/2
    Tbc[13,0] = - Qequipment[t,0] - Qexhaust_av[t,0] + Qleakage[t,0] + Q_evap_req[t-0,0]/2 + Qcv_emitted[t-0,0]/2

    # Temperatures calculation
    #temperature[:,t] = np.linalg.solve( K - C, Tbc - np.dot(C, prev_temp[:,t-1]) )
    rhs = Tbc - np.dot(C, prev_temp[:,0])
    temperature = np.linalg.solve(K - C, rhs)

    # Exchange from the front and back:
    if temperature[12,0] > T_amb[t,0]:
        temperature[12,0] = (density_air * V_cabin/2 * cp_air * prev_temp[12,0]/timestep + UA_front[t,0] *
                              T_amb[t,0]) / (density_air * V_cabin/2 * cp_air/timestep + UA_front[t,0])
        prev_temp[12,0] = temperature[12,0]
        Qcv_emitted[t,0] = (temperature[12,0] - T_amb[t,0]) * UA_front[t,0]
        Qcv_received[t,0] = 0
    else:
        Qcv_emitted[t,0] = 0
        Qcv_received[t,0] = (T_amb[t,0] - temperature[12,0]) * UA_front[t,0]

    if temperature[13,0] > T_amb[t,0]:
        temperature[13,0] = (density_air * V_cabin/2 * cp_air * prev_temp[13,0]/timestep + UA_back[t,0] * T_amb[
            t]) / (density_air * V_cabin/2 * cp_air/timestep + UA_back[t,0])
        prev_temp[13,0] = temperature[13,0]
        Qcv_emitted[t,0] = (temperature[13,0] - T_amb[t,0]) * UA_back[t,0]
    else:
        Qcv_emitted[t,0] = 0
        Qcv_received[t,0] = (T_amb[t,0] - temperature[13, 0]) * UA_back[t,0]

    # Heat requested from the MAC system (evaporator's or condenser's work), Q cabin
    Qcabin_received[t,0] = Qhuman[t,0] + Qequipment[t,0] + Qengine_av[t,0] + Qexhaust_av[t,0] + Qvent[t,0] + Qleakage[t,0] + Qcv_received[t,0]
    Qcabin_tot[t,0] = Qcabin_received[t,0] + Qcv_emitted[t,0]

    # Define cooling or heating mode & evap/condenser temps
    if T_amb[0,0] >= T_target:
        cooling = True
        Q_evap_req[t,0] = abs(Qcabin_req[t-1,0])
        T_evap = T_target - delta_T_evap
        T_cond = T_amb[t,0] + delta_T_cond
    else:  # Heating
        cooling = False
        Q_cond_req[t,0] = abs(Qcabin_req[t-1,0])
        T_evap = T_amb[t,0] - delta_T_evap
        T_cond = T_target + delta_T_cond

    # Call diagram p-h on CoolProp
    Ref = fluid

    # Check if T_amb is in the table
    idx = np.where(R1234yf_op_pres['T_ambient_C'] == T_amb[t, 0] - 273.15)[0]
    #print(f"idx = {idx} with size: {idx.size}")

    if idx.size > 0:  # If T_amb exists in the table, use exact value #before: len(idx) > 0:
        P_LowSide_max_Pa = 1000 * R1234yf_op_pres['P_LowSide_max_kPa'].iloc[idx[0]]
        P_HighSide_max_Pa = 1000 * R1234yf_op_pres['P_HighSide_max_kPa'].iloc[idx[0]]
        print(f"No interpolation. P_LowSide_max_Pa at T_ambient = {T_amb[t, 0]-273.15}°C: {P_LowSide_max_Pa} Pa")
    else:  # (meaning idx.size = 0) If T_amb is not found, interpolate
        # Extract column values
        T_ambient = R1234yf_op_pres['T_ambient_C'].values
        P_LowSide_max_kPa = R1234yf_op_pres['P_LowSide_max_kPa'].values
        P_HighSide_max_kPa = R1234yf_op_pres['P_HighSide_max_kPa'].values

        # Ensure T_ambient is sorted for np.interp
        sorted_indices = np.argsort(T_ambient)
        T_ambient = T_ambient[sorted_indices]
        P_LowSide_max_kPa = P_LowSide_max_kPa[sorted_indices]
        P_HighSide_max_kPa = P_HighSide_max_kPa[sorted_indices]

        # Perform interpolation
        P_LowSide_max_Pa = np.interp(T_amb[t, 0], T_ambient, P_LowSide_max_kPa) * 1000
        P_HighSide_max_Pa = np.interp(T_amb[t, 0], T_ambient, P_HighSide_max_kPa) * 1000
        #print(f"Interpolated P_LowSide_max_Pa at T_ambient = {T_amb[t, 0]-273.15}°C: {P_LowSide_max_Pa} Pa")
    #print(R1234yf_op_pres.head())  # Display first few rows

    compression_ratio_min = 5.7357  # Average

    # (1): evaporator outlet, saturated vapor;
    T1[t,0] = T_evap  # T_evap_out is known
    P1[t,0] = PropsSI('P', 'T', T1[t,0], 'Q', 1, Ref)
    # Ensure evaporator pressure does not exceed the maximum low-side pressure
    if P1[t,0] >= P_LowSide_max_Pa:
        P1[t,0] = P_LowSide_max_Pa
        T1[t,0] = PropsSI('T', 'P', P1[t,0], 'Q', 1, Ref)
    h1[t,0] = PropsSI('H', 'T', T1[t,0], 'Q', 1, Ref)
    s1[t,0] = PropsSI('S', 'T', T1[t,0], 'Q', 1, Ref)

    # (2') Isoentropic compression (theoretical)
    s2is[t,0] = s1[t,0]
    # (2) Compressor's outlet = Condenser's inlet
    P2is[t,0] = PropsSI('P', 'T', T_cond, 'S', s2is[t,0], Ref)

    # Actual compression ratio
    P2[t,0] = P2is[t,0]

    # Compression ratio needs to meet/exceed the minimum one
    if P2[t,0] / P1[t,0] < compression_ratio_min:
        # Adjust P2 so compression_ratio_actual = compression_ratio_min
        P2[t,0] = compression_ratio_min * P1[t,0]
        P2is[t,0] = P2[t,0]

        # Recompute T_cond based on adjusted P2
        T_cond = PropsSI('T', 'P', P2[t,0], 'Q', 1, Ref)

        # Update delta_T_cond to reflect the new condenser saturation temperature
        delta_T_cond = T_cond - T_amb[t,0]

    h2is[t,0] = PropsSI('H', 'P', P2[t,0], 'S', s2is[t,0], Ref)
    h2[t,0] = (h2is[t,0] - h1[t,0]) / eta_is + h1[t,0]
    T2[t,0] = PropsSI('T', 'H', h2[t,0], 'P', P2[t,0], Ref)
    s2[t,0] = PropsSI('S', 'H', h2[t,0], 'P', P2[t,0], Ref)

    # (3) Condenser's outlet = Expansion vale's inlet, Saturated liquid
    T3[t,0] = PropsSI('T', 'P', P2[t,0], 'Q', 0, Ref)
    h3[t,0] = PropsSI('H', 'P', P2[t,0], 'Q', 0, Ref)
    s3[t,0] = PropsSI('S', 'P', P2[t,0], 'Q', 0, Ref)

    # (4) Valve's outlet = Evaporator's inlet
    s4[t,0] = s3[t,0]
    P4[t,0] = P1[t,0]
    h4[t,0] = PropsSI('H', 'S', s4[t,0], 'P', P4[t,0], Ref)
    T4[t,0] = PropsSI('T', 'S', s4[t,0], 'P', P4[t,0], Ref)

    # Coefficient of performance
    COP[t,0] = (h1[t,0] - h4[t,0]) / (h2[t,0] - h1[t,0])

    # PID CONTROLLER for Compressor Speed
    temp_error = Tcabin[t,0] - T_target

    # PID gains
    Kp = 0.05
    Ki = 0.005
    Kd = 0
    # Proportional, integral and derivative terms
    P = Kp * temp_error
    error_sum = error_sum + temp_error * timestep
    I = Ki * error_sum
    D = Kd * (temp_error - temp_error) / timestep

    # Total control output
    comp_speed = P + I + D

    # Relating Speed to Mass Flow
    k = 0.005
    mf[t,0] = k * comp_speed

    # Components' heat loads calculation
    Q_evap[t,0] = mf[t,0] * (h1[t,0] - h4[t,0])
    Q_cond[t,0] = mf[t,0] * (h3[t,0] - h2[t,0])
    W_comp[t,0] = mf[t,0] * (h2[t,0] - h1[t,0])

    # Temperature recalculation after the MAC load:
    Q_MAC_zone = Q_evap[t,0] / 2
    temperature[12, 0] = (Qcabin_tot[t,0] / 2 - Q_MAC_zone) * (timestep / (V_cabin / 2 * cp_air * density_air)) + \
                         prev_temp[12, 0]
    temperature[13, 0] = (Qcabin_tot[t,0] / 2 - Q_MAC_zone) * (timestep / (V_cabin / 2 * cp_air * density_air)) + \
                         prev_temp[13, 0]

    # Save this timestep temperatures to be used in the following timestep
    prev_temp[:, 0] = temperature[:, 0]

# Plots MAC components
# Plots MAC components
plt.figure(1)
plt.plot(time, COP[1:, 0], linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('COP')
plt.grid(True)

plt.figure(2)
plt.plot(time, mf[1:] * 1000, linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Refrigerant mass flow, g/s')
plt.grid(True)

plt.figure(3)
plt.plot(time, W_comp[1:], linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Compressor Work (W)')
plt.grid(True)

plt.figure(4)
plt.plot(time, Q_evap[1:], linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Evaporator Heat Load (W)')
plt.grid(True)

plt.figure(5)
plt.plot(time, Q_cond[1:], linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Condenser Heat Load (W)')
plt.grid(True)

plt.figure(6)
plt.plot(time, Qcabin_received[:len(time)].flatten(), linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Cabin Heat Load (W)')
plt.grid(True)

plt.show()  # Display all plots