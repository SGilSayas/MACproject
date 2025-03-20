% FIRS VERSION
while (t <= Total_time)
%% CALCULO DE COEFICIENTES DE CONVECCION

% T_Sky set to Ambient Temperature

T_sky = T_ext(t);
%% WINDSHIELD

% Q_cd(t)=Q_ws(t); %2nd Eq

G_ws_r(t)=ro_ws*G_ws_inc(t);
G_ws_a(t)=alfa_ws*G_ws_inc(t); %añadir como bc
G_ws_t(t)=tao_ws*G_ws_inc(t);

G_ws_inc(t)=G_ws_r(t)+G_ws_a(t)+G_ws_t(t);

% Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
% Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
% Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));

% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_ws*G_ws_a(t); %1er Eq

% 1er Eq (Sustitucion)
% h_ext*A_ws*(T_ws_ext(t)-T_ext(t))+epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4) = (A_ws*k_ws/e_ws)*(T_ws_ext(t) - T_ws(t)) +  A_ws*G_ws_a(t);

% 2nd Eq 
% T_ws_ext(t) = (h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t);

% Iteracion
 It_max = 1000;
 err_prev = 0;
 err = 0;
 T_ws(t) = T_ws(t-timestep);
 rest  = false;
 DT = 0.01;
 for i = 0:0.02:It_max
     
     Ter1 = h_ext*A_ws*(((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))-T_ext(t))+epsilon_ws*sigma*A_ws*(((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))^4-(T_sky)^4);
     Ter2 = (A_ws*k_ws/e_ws)*(T_ws(t) - ((h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t))) +  A_ws*G_ws_a(t);
     err = abs(Ter2 - Ter1);
     
     if (err < err_prev)&&(rest == true)
         T_ws(t) = T_ws(t) - DT;
         rest = true; 
     elseif (err < err_prev)&&(rest == false)
         T_ws(t) = T_ws(t) + DT;
         rest = false;
     elseif (err >= err_prev)&&(rest == true)
         T_ws(t) = T_ws(t) + DT;
         rest = false; 
     elseif (err >= err_prev)&&(rest == false)
         T_ws(t) = T_ws(t) - DT;
         rest = true;
     end
     err_prev = err;
     if err < 0.01
         i = It_max;
     end
 end
 
%Calcula T_ws_ext(t)
T_ws_ext(t) = (h_ws*A_ws*(T_ws(t)- T_air(t-1)))/((A_ws*k_ws/e_ws)) + T_ws(t);
% Comprobacion del balance
% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_ws*G_ws_a(t); %1er Eq
% Q_cd(t)=Q_ws(t); %2nd Eq
Q_ws(t) = h_ws*A_ws*(T_ws(t)- T_air(t-1));
Q_cv_ext(t)=h_ext*A_ws*(T_ws_ext(t)-T_ext(t)); %Exterior convection
Q_rd_ext(t)=epsilon_ws*sigma*A_ws*((T_ws_ext(t))^4-(T_sky)^4); %Exterior radiation
Q_cd(t)=((k_ws*A_ws)/e_ws)*(T_ws_ext(t)-T_ws(t));

Conv_Rad_ws(t) = Q_cv_ext(t) + Q_rd_ext(t);
Cond_RadAbs_ws(t) = - Q_cd(t) + A_ws*G_ws_a(t);
Cond_ws(t) = Q_cd(t);
Convin_ws(t) = Q_cd(t);
Error_WS(t) = err;
%% LEFT SIDE WINDOW

G_lsw_r(t)=ro_lsw*G_lsw_inc(t);
G_lsw_a(t)=alfa_lsw*G_lsw_inc(t);
G_lsw_t(t)=tao_lsw*G_lsw_inc(t);

G_lsw_inc(t)=G_lsw_r(t)+G_lsw_a(t)+G_lsw_t(t);

% Iteracion
 It_max = 1000;
 err_prev = 0;
 err = 0;
 T_lsw(t) = T_lsw(t-timestep);
 rest  = false;
 DT = 0.01;
 for i = 0:0.02:It_max
     
     Ter1 = h_ext*A_lsw*(((h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t))-T_ext(t))+epsilon_lsw*sigma*A_lsw*(((h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t))^4-(T_sky)^4);
     Ter2 = (A_lsw*k_lsw/e_lsw)*(T_lsw(t) - ((h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t))) +  A_lsw*G_lsw_a(t);
     err = abs(Ter2 - Ter1);
     
     if (err < err_prev)&&(rest == true)
         T_lsw(t) = T_lsw(t) - DT;
         rest = true; 
     elseif (err < err_prev)&&(rest == false)
         T_lsw(t) = T_lsw(t) + DT;
         rest = false;
     elseif (err >= err_prev)&&(rest == true)
         T_lsw(t) = T_lsw(t) + DT;
         rest = false; 
     elseif (err >= err_prev)&&(rest == false)
         T_lsw(t) = T_lsw(t) - DT;
         rest = true;
     end
     err_prev = err;
     if err < 0.01
         i = It_max;
     end
 end
 
%Calcula T_lsw_ext(t)
 T_lsw_ext(t) = (h_lsw*A_lsw*(T_lsw(t)- T_air(t-1)))/((A_lsw*k_lsw/e_lsw)) + T_lsw(t);
% Comprobacion del balance
% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_lsw*G_lsw_a(t); %1er Eq
% Q_cd(t)=Q_lsw(t); %2nd Eq
Q_lsw(t) = h_lsw*A_lsw*(T_lsw(t)- T_air(t-1));
Q_cv_ext_lsw(t)=h_ext*A_lsw*(T_lsw_ext(t)-T_ext(t)); %Exterior convection
Q_rd_ext_lsw(t)=epsilon_lsw*sigma*A_lsw*((T_lsw_ext(t))^4-(T_sky)^4); %Exterior radiation
Q_cd_lsw(t)=((k_lsw*A_lsw)/e_lsw)*(T_lsw_ext(t)-T_lsw(t));

Conv_Rad_lsw(t) = Q_cv_ext_lsw(t) + Q_rd_ext_lsw(t);
Cond_RadAbs_lsw(t) = - Q_cd_lsw(t) + A_lsw*G_lsw_a(t);
Cond_lsw(t) = Q_cd_lsw(t);
Convin_lsw(t) = Q_cd_lsw(t);

Error_LSW(t) = err;
%% RIGHT SIDE WINDOW

G_rsw_r(t)=ro_lsw*G_rsw_inc(t);
G_rsw_a(t)=alfa_rsw*G_rsw_inc(t);
G_rsw_t(t)=tao_lsw*G_rsw_inc(t);

G_rsw_inc(t)=G_rsw_r(t)+G_rsw_a(t)+G_rsw_t(t);

% Iteracion
 It_max = 1000;
 err_prev = 0;
 err = 0;
 T_rsw(t) = T_rsw(t-timestep);
 rest  = false;
 DT = 0.01;
 for i = 0:0.02:It_max
     
     Ter1 = h_ext*A_rsw*(((h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t))-T_ext(t))+epsilon_rsw*sigma*A_rsw*(((h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t))^4-(T_sky)^4);
     Ter2 = (A_rsw*k_rsw/e_rsw)*(T_rsw(t) - ((h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t))) +  A_rsw*G_rsw_a(t);
     err = abs(Ter2 - Ter1);
     
     if (err < err_prev)&&(rest == true)
         T_rsw(t) = T_rsw(t) - DT;
         rest = true; 
     elseif (err < err_prev)&&(rest == false)
         T_rsw(t) = T_rsw(t) + DT;
         rest = false;
     elseif (err >= err_prev)&&(rest == true)
         T_rsw(t) = T_rsw(t) + DT;
         rest = false; 
     elseif (err >= err_prev)&&(rest == false)
         T_rsw(t) = T_rsw(t) - DT;
         rest = true;
     end
     err_prev = err;
     if err < 0.01
         i = It_max;
     end
 end
 
%Calcula T_rsw_ext(t)
 T_rsw_ext(t) = (h_rsw*A_rsw*(T_rsw(t)- T_air(t-1)))/((A_rsw*k_rsw/e_rsw)) + T_rsw(t);
% Comprobacion del balance
% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_rsw*G_rsw_a(t); %1er Eq
% Q_cd(t)=Q_rsw(t); %2nd Eq
Q_rsw(t) = h_rsw*A_rsw*(T_rsw(t)- T_air(t-1));
Q_cv_ext_rsw(t)=h_ext*A_rsw*(T_rsw_ext(t)-T_ext(t)); %Exterior convection
Q_rd_ext_rsw(t)=epsilon_rsw*sigma*A_rsw*((T_rsw_ext(t))^4-(T_sky)^4); %Exterior radiation
Q_cd_rsw(t)=((k_rsw*A_rsw)/e_rsw)*(T_rsw_ext(t)-T_rsw(t));

Conv_Rad_rsw(t) = Q_cv_ext_rsw(t) + Q_rd_ext_rsw(t);
Cond_RadAbs_rsw(t) = - Q_cd_rsw(t) + A_rsw*G_rsw_a(t);
Cond_rsw(t) = Q_cd_rsw(t);
Convin_rsw(t) = Q_cd_rsw(t);

Error_RSW(t) = err;

%% REAR WINDOW

G_rw_r(t)=ro_rw*G_rw_inc(t);
G_rw_a(t)=alfa_rw*G_rw_inc(t);
G_rw_t(t)=tao_rw*G_rw_inc(t);

G_rw_inc(t)=G_rw_r(t)+G_rw_a(t)+G_rw_t(t);
% Iteracion
 It_max = 1000;
 err_prev = 0;
 err = 0;
 T_rw(t) = T_rw(t-timestep);
 rest  = false;
 DT = 0.01;
 for i = 0:0.02:It_max
     
     Ter1 = h_ext*A_rw*(((h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t))-T_ext(t))+epsilon_rw*sigma*A_rw*(((h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t))^4-(T_sky)^4);
     Ter2 = (A_rw*k_rw/e_rw)*(T_rw(t) - ((h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t))) +  A_rw*G_rw_a(t);
     err = abs(Ter2 - Ter1);
     
     if (err < err_prev)&&(rest == true)
         T_rw(t) = T_rw(t) - DT;
         rest = true; 
     elseif (err < err_prev)&&(rest == false)
         T_rw(t) = T_rw(t) + DT;
         rest = false;
     elseif (err >= err_prev)&&(rest == true)
         T_rw(t) = T_rw(t) + DT;
         rest = false; 
     elseif (err >= err_prev)&&(rest == false)
         T_rw(t) = T_rw(t) - DT;
         rest = true;
     end
     err_prev = err;
     if err < 0.01
         i = It_max;
     end
 end
 
%Calcula T_rw_ext(t)
 T_rw_ext(t) = (h_rw*A_rw*(T_rw(t)- T_air(t-1)))/((A_rw*k_rw/e_rw)) + T_rw(t);
% Comprobacion del balance
% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_rw*G_rw_a(t); %1er Eq
% Q_cd(t)=Q_rw(t); %2nd Eq
Q_rw(t) = h_rw*A_rw*(T_rw(t)- T_air(t-1));
Q_cv_ext_rw(t)=h_ext*A_rw*(T_rw_ext(t)-T_ext(t)); %Exterior convection
Q_rd_ext_rw(t)=epsilon_rw*sigma*A_rw*((T_rw_ext(t))^4-(T_sky)^4); %Exterior radiation
Q_cd_rw(t)=((k_rw*A_rw)/e_rw)*(T_rw_ext(t)-T_rw(t));

Conv_Rad_rw(t) = Q_cv_ext_rw(t) + Q_rd_ext_rw(t);
Cond_RadAbs_rw(t) = - Q_cd_rw(t) + A_rw*G_rw_a(t);
Cond_rw(t) = Q_cd_rw(t);
Convin_rw(t) = Q_cd_rw(t);

Error_RW(t) = err;


%% GENERAL WINDOW EQUATION

Q_windows(t)=Q_ws(t)+Q_lsw(t)+Q_rsw(t)+Q_rw(t); %% Listo

%% ROOF AND CEILING

G_roof_a(t)=alfa_roof*G_roof_inc(t);

% Iteracion
 It_max = 1000;
 err_prev = 0;
 err = 0;
 T_ceiling(t) = T_ceiling(t-timestep);
 rest  = false;
 DT = 0.01;
 for i = 0:0.02:It_max
     
     Ter1 = h_roof*A_roof*(((h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t))-T_ext(t))+epsilon_roof*sigma*A_roof*(((h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t))^4-(T_sky)^4);
     Ter2 = (A_roof*k_ceiling/e_ceiling)*(T_ceiling(t) - ((h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t))) +  A_roof*G_roof_a(t);
     err = abs(Ter2 - Ter1);
     
     if (err < err_prev)&&(rest == true)
         T_ceiling(t) = T_ceiling(t) - DT;
         rest = true; 
     elseif (err < err_prev)&&(rest == false)
         T_ceiling(t) = T_ceiling(t) + DT;
         rest = false;
     elseif (err >= err_prev)&&(rest == true)
         T_ceiling(t) = T_ceiling(t) + DT;
         rest = false; 
     elseif (err >= err_prev)&&(rest == false)
         T_ceiling(t) = T_ceiling(t) - DT;
         rest = true;
     end
     err_prev = err;
     if err < 0.01
         i = It_max;
     end
 end
 
%Calcula T_rw_ext(t)
 T_roof(t) = (h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1)))/((A_roof*k_ceiling/e_ceiling)) + T_ceiling(t);
% Comprobacion del balance
% Q_cv_ext(t)+Q_rd_ext(t) = Q_cd(t)+A_roof*G_roof_a(t); %1er Eq
% Q_cd(t)=Q_ceiling(t); %2nd Eq
Q_ceiling(t) = h_ceiling*A_roof*(T_ceiling(t)- T_air(t-1));
Q_cv_ext_roof(t)=h_ext*A_roof*(T_roof(t)-T_ext(t)); %Exterior convection
Q_rd_ext_roof(t)=epsilon_roof*sigma*A_roof*((T_roof(t))^4-(T_sky)^4); %Exterior radiation
Q_cd_ceiling(t)=((k_ceiling*A_roof)/e_ceiling)*(T_roof(t)-T_ceiling(t));

Conv_Rad_ceiling(t) = Q_cv_ext_roof(t) + Q_rd_ext_roof(t);
Cond_RadAbs_ceiling(t)= - Q_cd_ceiling(t) + A_roof*G_roof_a(t);
Cond_ceiling(t) = Q_cd_ceiling(t);
Convin_ceiling(t) = Q_ceiling(t);

Error_CR(t) = err;



%% BASE
% 1a Eq: (C_base*(T_base(t)-T_base(t-1))/(timestep)) = alfa_base*((A_ws*G_ws_t)+(A_rw*G_rw_t)+(A_lsw*G_lsw_t)+(A_lsw*G_rsw_t))-Q_base;

% 2a Eq:   Q_base=h_base*A_base*(T_base(t)-T_air(t-1));

% Sustituyendo 2a en 1a
T_base(t) = (timestep/(C_base+h_base*A_base)) * (alfa_base*((A_ws*G_ws_t(t))+(A_rw*G_rw_t(t))+(A_lsw*G_lsw_t(t))+(A_lsw*G_rsw_t(t))) + C_base*T_base(t-1) + h_base*A_base*T_air(t-1));
Q_base(t)=h_base*A_base*(T_base(t)-T_air(t-1));

%% BEINGS
Q_human(t)=N_Humans * 120;

%% GENERAL EQUATION
% m_air*Cp_air*(T_air(t)-T_air(t-1))/timestep=Q_windows+Q_ceiling+Q_base+Q_human; 

T_air(t) = (Q_windows(t)+Q_ceiling(t)+Q_base(t)+Q_human(t))*(timestep/(CabinVolume*cp_air*density_air))  + T_air(t-1);

t = t +timestep;
end

PLOT

% Temp Air
plot(time(2:end),T_air - 273,'b',T_Air_Exp(:,1),T_Air_Exp(:,2),'r',T_Air_Mod_Ref(:,1),T_Air_Mod_Ref(:,2),'g');

