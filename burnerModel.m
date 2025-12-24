function burner = burnerModel(x, N_spools, data, inlet, comp)
% ---------------------------------------------------------
% burnerModel
% ---------------------------------------------------------
% x(1) = pi_c
% x(2) = Tt4_target
%
% Includes:
%  - Fuel equilibrium via fuelEquilibrium.m
%  - Burner pressure loss
%  - Annular diffuser at burner exit (Farokhi Ch. 6)
% ---------------------------------------------------------

%% Unpack design variables
pi_c = x(1);
Tt4_target = x(2);

%% Unpack constants
gamma_03 = data.gamma_03;
Cp_03    = data.Cp_03;
gamma_49 = data.gamma_49;
Cp_49    = data.Cp_49;

T0 = data.T0;
P0 = data.P0;
V0 = data.V0;

pi_b = data.pi_b;
f_st = data.f_st;

%% ---------------- Inlet → compressor exit (station 3) ----------------
R_03 = Cp_03*(gamma_03-1)/gamma_03;
a0   = sqrt(gamma_03*R_03*max(T0,1));
M0   = V0/a0;

Tt0  = T0*(1 + (gamma_03-1)/2*M0^2);
Pt0  = P0*(1 + (gamma_03-1)/2*M0^2)^(gamma_03/(gamma_03-1));

if isfield(data,'eta_c'), eta_c = data.eta_c; else, eta_c = 0.9; end

Pt3 = comp.Pt3;
Tt3 = comp.Tt3;

% Station 3 static (assumed)
M3 = 0.15;

tau3 = 1 + (gamma_03-1)/2*M3^2;
T3   = Tt3/tau3;
P3   = Pt3/tau3^(gamma_03/(gamma_03-1));
rho3 = P3/(R_03*T3);
V3   = M3*sqrt(gamma_03*R_03*max(T3,1));

mdot3 = comp.mdot3;
A3 = mdot3/(rho3*V3);
burner.A3 = A3;

%% ---------------- Fuel equilibrium (Farokhi Ch. 7) ----------------
% Use equilibrium table
phi_guess = (Cp_49*(Tt4_target-Tt3))/(data.LHV*f_st);
[Tt4_eq, X_CO, Q_R] = fuelEquilibrium(Tt3, Pt3, phi_guess);

% Enforce Tt4 from design variable
Tt4 = Tt4_target;
f   = phi_guess*f_st;

%% ---------------- Burner pressure loss ----------------
Pt4 = pi_b*Pt3;

%% ---------------- Burner exit static (before diffuser) ----------------
R_49 = Cp_49*(gamma_49-1)/gamma_49;
M4 = 0.25;

tau4 = 1 + (gamma_49-1)/2*M4^2;
T4   = Tt4/tau4;
P4   = Pt4/tau4^(gamma_49/(gamma_49-1));
rho4 = P4/(R_49*T4);
V4   = M4*sqrt(gamma_49*R_49*max(T4,1));

mdot4 = mdot3*(1+f);
A4    = mdot4/(rho4*V4);
AR_diff = inlet.A2 / inlet.A1;
A4d = AR_diff*A4;

phi_bd_deg = data.burner_phi_deg;
phi_bd = phi_bd_deg*pi/180;

rbar = 0.5*(inlet.Ro2 + inlet.Ri2);
w4  = A4 /(2*pi*rbar);
w4d = A4d/(2*pi*rbar);

% outer radius change (inner is symmetric in this simple model)
dRo = 0.5*(w4d - w4);

burner.l_bd = abs(dRo)/max(tan(phi_bd),1e-9);
Lb = 1.1;
burner.l_b  = Lb;


%% ---------------- Annular diffuser (Farokhi Ch. 6) ----------------
% Geometry assumptions
pi_d    = 0.985;       % total pressure recovery

Pt4d = pi_d*Pt4;
Tt4d = Tt4;


% Solve diffuser exit Mach iteratively
M = 0.15;
for k=1:30
    tau = 1+(gamma_49-1)/2*M^2;
    T   = Tt4d/tau;
    P   = Pt4d/tau^(gamma_49/(gamma_49-1));
    rho = P/(R_49*T);
    V   = M*sqrt(gamma_49*R_49*max(T,1));
    mdot_guess = rho*V*A4d;
    M = max(0.02, M*(mdot4/mdot_guess)^0.5);
end

%% ---------------- Residence time ----------------

tau_res = Lb/V4;

%% ---------------- Pack outputs ----------------
burner.Tt3 = Tt3;
burner.Pt3 = Pt3;
burner.mdot3 = mdot3;

burner.Tt4 = Tt4;
burner.Pt4 = Pt4;
burner.Tt4d = Tt4d;
burner.Pt4d = Pt4d;

burner.f = f;
burner.phi = f/f_st;
burner.X_CO = X_CO;
burner.Q_R = Q_R;

burner.mdot4 = mdot4;
burner.A4 = A4;
burner.A4d = A4d;
burner.tau_res = tau_res;

end
