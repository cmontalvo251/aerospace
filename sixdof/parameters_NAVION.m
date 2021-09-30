% Aircraft data for the NAVION General Aviation Aircraft

%% Atmosphere data:

rho = 1.22566578494891*(1+0.0000225696709*z)^4.258;
rho = rho*(0.002378/1.22566578494891); %%convert to SI units

%%Rho0
rho0 = 1.22566578494891*(1-0.0000225696709*200)^4.258;
rho0 = rho0*(0.002378/1.22566578494891); %%convert to SI units

V0 = 175; %%trim velocity

% Mass Moment of Inertia Matrix in Body Axis (slug.ft^2)
Ix = 1048;
Iy = 3000;
Iz = 3530;
I_b = diag([Ix Iy Iz]);
I_b_inv = diag(1./[Ix Iy Iz]);

%% Longitudinal Coefficients

C_L_0 = 0.41;
C_D_0 = 0.05;
C_L_alpha = 4.44;
C_D_alpha = 0.33;
C_m_alpha = -0.683;
C_m_alpha_dot = -4.36;
C_L_q = 3.8;
C_m_q = -9.96;
C_L_de = -0.355;
C_m_de = 0.923;
C_x_dT = 10;

%% Lateral Coefficients

C_y_beta = -0.564;
C_l_beta = -0.074;
C_n_beta = -0.071;
C_l_p = -0.410;
C_n_p = -0.0575;
C_l_r = 0.107;
C_n_r = -0.125;
C_l_da = -0.134;
C_n_da = -0.0035;
C_y_dr = 0.157;
C_l_dr = 0.107;
C_n_dr = -0.072;

%% Geometry and Inertias

% Reference Area (ft^2)
S = 184;

% Wingspan (ft)
b = 33.4;

% Length
len = 5.8/7.2*b;

% Height
hei = 1.5/7.2*b;

% Mean Chord (ft)
c_bar = 5.7;

% Gross Weight (lbs)
GW = abs(C_L_0*1/2*rho0*V0^2*S); %%this makes sure that W=L
mass = GW/32.2;

%% Thrust information

T0 = abs(C_D_0*1/2*rho*V0^2*S); %%this just makes sure that T=D


% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
