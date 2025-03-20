# Vehicle HVAC and Cabin Model
# sgs, 2021

from scipy.io import loadmat
from scipy.interpolate import interp1d
import numpy as np
import CoolProp.CoolProp as CP
#import matplotlib.pyplot as plt
import pandas as pd

### INPUTS ###

# Radiation Properties
# epsilon=emissivity[-], rho=reflectivity[-], tau=transmissivity[-], alpha=absorptivity[-]
# ws: Windshield, rw: Rear Window, lsw: Left Side Windows, rsw: Right Side Windows
sigma = 5.67037321 * 10 ** (-8)  # W/(m^2*K^4) #Stefan-Boltzmann constant
epsilon_ws = 0.9
rho_ws = 0.246
tao_ws = 0.452
alfa_ws = 1 - (rho_ws + tao_ws)
epsilon_rw = 0.9
rho_rw = 0.1
tao_rw = 0.311
alfa_rw = 1 - rho_rw - tao_rw
epsilon_lsw = 0.9
epsilon_rsw = 0.9
rho_lsw = 0.2
tao_lsw = 0.475
alfa_lsw = 1 - rho_lsw - tao_lsw
rho_rsw = rho_lsw
tao_rsw = tao_lsw
alfa_rsw = alfa_lsw
epsilon_roof = 0.9
rho_roof = 0.74
tao_roof = 0
alfa_roof = 1 - rho_roof - tao_roof
rho_base = 0.3
alfa_base = 0.7

# Conduction Properties
# A=area[m2], k=thermoconductivity[W/(m*K)], e=thickness[m]
# ws: Windshield, rw: Rear Window, lsw: Left Side Windows, rsw: Right Side Windows
A_ws = 0.63 * 1.3
k_ws = 1.4
e_ws = 0.006
A_rw = 0.29 * 1
k_rw = 1.4
e_rw = 0.005
A_lsw = 1.45 * 0.29
A_rsw = 1.45 * 0.29
k_lsw = 1.4
k_rsw = 1.4
e_lsw = 0.003
e_rsw = 0.003
A_roof = 1.8 * 1.1
k_ceiling = 14.9  # ONLY STEEL LAYER
e_ceiling = 0.0005  # ONLY STEEL LAYER
A_base = 6
C_base = 144240  # Capacitance, W

# Thermal Coefficients # TODO: add correlations
h_ws = 5
h_lsw = 5
h_rsw = 5
h_rw = 5
h_ceiling = 5
h_roof = 20
h_ext = 20
h_base = 5

# Cabin Geometry aprox
cabin_h = 0.5456  # m, heigh
cabin_w = 1.1  # m, width
cabin_r_l = 1.8  # m, rear lenght
cabin_b_l = 1.45  # m, ?? lenght
V_cabin = cabin_w * (cabin_r_l + cabin_b_l) * cabin_h  # = 1.45*1.3*1
S_cabin = cabin_w + cabin_r_l + cabin_b_l + cabin_h
Temp_comfort = 23 + 273.15  # K

# Air properies
rho_air = 1.18  # kg/m3 (25ºC)
cp_air = 1006  # J/kg*K  "
volume_flow_HVAC = 111.1  # l/s, Air flow volume rate (Lee et al., 2015)
mass_flow_HVAC = volume_flow_HVAC * rho_air  # =131.098 g/s  0.1311 kg/s
mass_flow_HVAC_heating = rho_air * 97.2
T_sky = 293  # K  #T_sky=T_ini-6

# Test case inputs
Total_time_h = 1  # h
t = 1  # s
timestep = 1  # s
Tmax_C = 50  # ºC
Tmin_C = 50  # ºC
N_Humans = 1
# irradiance="false"
irradiance = "true"

if irradiance == "true":  # Load input IRRADIANCE #TODO: actualizar funciones a python
    inputs_mat = loadmat(r"C:\Users\susan\OneDrive\Escritorio\Python_Projects/Inputs.mat")  # scipy.io.loadmat
    Total_time = inputs_mat.get("Irradiancia")[:,0] #inputs_mat['Irradiancia'][:,0]
    Irr_raw = inputs_mat.get("Irradiancia")[:,1]     # .T (transposada) al final per columna, imp al multiplicar
    time = np.arange(0, Total_time[-1], timestep)
    Irr = interp1d(Total_time, Irr_raw)(time)

    temp_mat = loadmat(r"C:\Users\susan\OneDrive\Escritorio\Python_Projects/T_Air_mod_Ref.mat")
    time_temp_ext_raw = temp_mat.get("T_Air_Mod_Ref")[:, 0]
    Temp_ext_raw = temp_mat.get("T_Air_Mod_Ref")[:,1]
    time2 = np.arange(0, time_temp_ext_raw[-1], timestep)
    T_ext = 273.15 + interp1d(time_temp_ext_raw, Temp_ext_raw)(time2)
else:
    Total_time = Total_time_h * 3600  # s
    T_ext_C = np.linspace(Tmin_C, Tmax_C, Total_time + 1)  # start, end, no_samples
    T_ext = 273.15 + T_ext_C
    time = np.arange(0, Total_time, timestep)

T_ini = T_ext[1]  # () for functions, [] for elements in array

if irradiance == "true":  # Global Horizontal Irradiance (Hipótesis: same all windows)
    G_roof_inc = Irr
    G_ws_inc = Irr
    G_rw_inc = Irr
    G_lsw_inc = Irr
    G_rsw_inc = Irr
    G_ws_r = rho_ws * G_ws_inc
    G_ws_a = alfa_ws * G_ws_inc
    G_ws_t = tao_ws * G_ws_inc
    G_rw_r = rho_rw * G_rw_inc
    G_rw_a = alfa_rw * G_rw_inc
    G_rw_t = tao_rw * G_rw_inc
    G_lsw_r = rho_lsw * G_lsw_inc
    G_lsw_a = alfa_lsw * G_lsw_inc
    G_lsw_t = tao_lsw * G_lsw_inc
    G_rsw_r = rho_rsw * G_rsw_inc
    G_rsw_a = alfa_rsw * G_rsw_inc
    G_rsw_t = tao_rsw * G_rsw_inc
    G_roof_a = alfa_roof * G_roof_inc

    # Conductance Matrix
    # T_ext --K1-> Troof --K2-> Tceiling --K3-> Tcabin
    K1 = h_ext * A_roof  # h_roof por qué existe?
    K2 = A_roof * k_ceiling / e_ceiling
    K3 = h_ceiling * A_roof
    # Tbase --K4-> Tcabin_back
    K4 = h_base * A_base
    # T_ext --K5-> Tw_ext_ws --K6-> Tw_int_ws --K7-> Tcabin
    K5 = h_ext * A_ws
    K6 = A_ws * k_ws / e_ws
    K7 = h_ws * A_ws
    # T_ext --K8-> Tw_ext_rw --K9-> Tw_int_rw --K10-> Tcabin
    K8 = h_ext * A_rw
    K9 = A_rw * k_rw / e_rw
    K10 = h_rw * A_rw
    # T_ext --K11-> Tw_ext_lsw --K12-> Tw_int_lsw --K13-> Tcabin
    K11 = h_ext * A_lsw
    K12 = A_lsw * k_lsw / e_lsw
    K13 = h_lsw * A_lsw
    # T_ext --K14-> Tw_ext_rsw --K15-> Tw_int_rsw --K16-> Tcabin
    K14 = h_ext * A_rsw
    K15 = A_rsw * k_rsw / e_rsw
    K16 = h_rsw * A_rsw
    K = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0...
...K1, -K1-K2, K2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
    print(K)
    """
        0, K2, -K2-K3, 0, 0, 0, 0, 0, 0, 0, 0, 0, K3
         0, 0, 0, -K4, 0, 0, 0, 0, 0, 0, 0, 0, K4
         K5, 0, 0, 0, -K5-K6, 0, 0, 0, K6, 0, 0, 0, 0
         K8, 0, 0, 0, 0, -K8-K9, 0, 0, 0, K9, 0, 0, 0
         K11, 0, 0, 0, 0, 0, -K11-K12, 0, 0, 0, K12, 0, 0
         K14, 0, 0, 0, 0, 0, 0, -K14-K15, 0, 0, 0, K15, 0
         0, 0, 0, 0, K6, 0, 0, 0, -K6-K7, 0, 0, 0, K7
         0, 0, 0, 0, 0, K9, 0, 0, 0, -K9-K10, 0, 0, K10
         0, 0, 0, 0, 0, 0, K12, 0, 0, 0, -K12-K13, 0, K13
         0, 0, 0, 0, 0, 0, 0, K15, 0, 0, 0, -K15-K16, K16
         0, 0, K3, K4, 0, 0, 0, 0, K7, K10, K13, K16, -K3-K4-K7-K10-K13-K16]

    # Capacitance Matrix
    C_amb = 0  # b.c. #todo: unexpected indent?
    C_roof = 0  # de momento
    C_ceiling = 0  # de momento
    C_base = 1442.40  # Capacitance, W
    C_ws = 0
    C_rw = 0
    C_lsw = 0
    C_rsw = 0
    C_air = rho_air * V_cabin * cp_air / timestep
    C_cabin = C_air
    C = zeros[13] #in pyhton like this?
    C[0,0]= Camb
    C[1, 1]= C_roof
    C[2, 2]= C_ceiling
    C[3, 3]= C_base
    C[4, 4]= C_ws
    C[5, 5]= C_rw
    C[6, 6]= C_lsw
    C[7, 7] = C_rsw
    C[8, 8] = C_ws
    C[9, 9]= C_rw
    C[10, 10] = C_lsw
    C[11, 11] = C_rsw
    C[12, 12] = C_cabin

    ## INI
    Troof[0] = T_ini
    Tceiling[0] = T_ini
    Tbase_int[0] = T_ini
    Tw_ext_ws[0] = T_ini
    Tw_ext_rw[0] = T_ini
    Tw_ext_lsw[0] = T_ini
    Tw_ext_rsw[0] = T_ini
    Tw_int_ws[0] = T_ini
    Tw_int_rw[0] = T_ini
    Tw_int_lsw[0] = T_ini
    Tw_int_rsw[0] = T_ini
    Tcabin[0] = T_ini
    temperature = [T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini]
    delta_temperature = [T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini T_ini + 1]
    Q_target_cum[0] = rho_air * V_cabin * cp_air * (temperature[12] - Temp_conf) / timestep  # ini
    Q_HVAC_cum[0] = Q_target_cum[0]
    Q_target[0] = rho_air * V_cabin * cp_air * (temperature[12] - Temp_conf) / timestep  # Tcabin(:,t-1)=324.7876
"""