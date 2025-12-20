function turbine = turbineModel(x, N_spools, data)
% ---------------------------------------------------------
% turbineModel
% ---------------------------------------------------------
% Uses:
%  - Burner output
%  - Compressor shaft power
%  - Ainley profile loss (reaction blades assumed)
%  - Cooling per project spec
% ---------------------------------------------------------

%% Get burner state
burner = burnerModel(x, N_spools, data);

Tt4 = burner.Tt4d;
Pt4 = burner.Pt4d;
mdot_g = burner.mdot4;

gamma = data.gamma_49;
cp    = data.Cp_49;
R     = cp*(gamma-1)/gamma;

%% Compressor power (placeholder until compressorModel done)
Cp_03 = data.Cp_03;
W_comp = mdot_g*Cp_03*(burner.Tt3 - data.T0);

%% Turbine mass flow including coolant
beta_init = 0.02;
mdot_t = mdot_g*(1+beta_init);

%% Turbine power balance
eta_m = data.eta_m;
DeltaTt = W_comp/(eta_m*mdot_t*cp);
Tt5 = Tt4 - DeltaTt;

%% Ainley profile loss → uncooled efficiency
% Reaction blade, s/c ≈ 0.7, beta3 ≈ 60°
Yp = 0.035;   % from Fig. 10.16 (typical)
eta_uncooled = 1 - Yp;

%% Cooling calculation (Farokhi Ch. 7.5)
c = data.c;
r_h = 0.25; r_t = 0.45;
Ag = 2*c*(r_t-r_h);
hg = 7000;

gamma_c = data.gamma_cool;
cp_c    = data.Cp_cool;
mu_c    = data.mu_cool;

Tt3 = burner.Tt3;
Pt3 = burner.Pt3;

M_c = 0.1;
tau_c = 1 + (gamma_c-1)/2*M_c^2;
Tc = Tt3/tau_c;
Pc = Pt3/tau_c^(gamma_c/(gamma_c-1));
rho_c = Pc/(cp_c*(gamma_c-1)/gamma_c*Tc);
Vc = M_c*sqrt(gamma_c*(cp_c*(gamma_c-1)/gamma_c)*Tc);

xchar = (r_t-r_h)/2;
Re_x = rho_c*Vc*xchar/mu_c;
St = 0.0296/(Re_x^0.2);
hc = St*rho_c*cp_c*Vc;

h_eff = 1/(1/hg + data.blade_tw/data.blade_kw + 1/hc);

% Gas-side static temp approx
Tg = Tt4;

qpp = h_eff*(Tg-Tc);
Twg = Tg - qpp/hg;

Qdot = qpp*Ag;
mdot_c = Qdot/(cp_c*(1000-Tc));
beta_c = mdot_c/mdot_t;

eta_pen = 0.032*(100*beta_c);
eta_t = eta_uncooled*(1-eta_pen);

%% Turbine exit pressure
Tt5s = Tt4 - (Tt4-Tt5)/eta_t;
Pt5 = Pt4*(Tt5s/Tt4)^(gamma/(gamma-1));

%% Tip Mach and centrifugal stress
omega = data.omega;
a5 = sqrt(gamma*R*Tt5);
U_tip = omega*r_t;

M_tip = U_tip/a5;
sigma_c = data.blade_density*omega^2*((r_h+r_t)/2)^2;

%% Pack outputs
turbine.Tt5 = Tt5;
turbine.Pt5 = Pt5;
turbine.Twg = Twg;
turbine.T_cool = Tc + Qdot/(mdot_c*cp_c);
turbine.beta_cool = beta_c;

turbine.eta_t = eta_t;
turbine.M_tip_turb = M_tip;
turbine.sigma_c = sigma_c;

end
