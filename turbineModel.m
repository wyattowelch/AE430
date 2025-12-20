function turbine = turbineModel(x, N_spools, data)
% turbineModel
% ------------------------------------------------------------
% Group interface: turbine = turbineModel(x, N_spools, data)
%
% This turbine model is designed to work with the burnerModel above.
% Since engineModel currently doesn't pass burner/comp structs into turbineModel,
% turbineModel recomputes burner internally (same x, data) to get station 4.
%
% Required outputs for constraints (minimum):
%   turbine.Twg, turbine.T_cool, turbine.sigma_c, turbine.M_tip_turb
%   turbine.Tt5, turbine.Pt5, turbine.eta_t
%
% Uses Farokhi Ch.7 turbine relations:
%   Power balance:  eta_m * mdot_t * cp * (Tt4 - Tt5) = W_req
%   Turbine eff:    eta_t = (Tt4 - Tt5)/(Tt4 - Tt5s)
%   Isentropic pt:  Pt5 = Pt4 * (Tt5s/Tt4)^(gamma/(gamma-1))
%
% Cooling model: simplified Stanton-based correlation per class spec.

%% -------------------- Upstream: get station 4 from burner --------------------
burner = burnerModel(x, N_spools, data);
Tt4 = burner.Tt4;
Pt4 = burner.Pt4;

% Hot gas properties (49)
gamma = data.gamma_49;
cp    = data.Cp_49;

% Turbine mass flow includes coolant
mdot_core = burner.mdot4;

% If you have a cooling mass flow in data, use it; otherwise default 5% core.
if isfield(data,'mdot_cool_total')
    mdot_cool_total = data.mdot_cool_total;
else
    mdot_cool_total = 0.05 * mdot_core;
end
mdot_t = mdot_core + mdot_cool_total;

%% -------------------- Compressor power requirement (simple until compModel exists) --------------------
% Approx compressor shaft power from inlet -> compressor exit:
%   W_comp = mdot3 * Cp_03 * (Tt3 - Tt0)
% turbine must supply this (plus any fan, not yet modeled).
Cp_03 = data.Cp_03;

% Recompute Tt0 and Tt3 similarly to burnerModel
gamma_03 = data.gamma_03;
T0 = data.T0; P0 = data.P0; V0 = data.V0;

R_03 = Cp_03 * ((gamma_03 - 1)/gamma_03);
a0   = sqrt(gamma_03 * R_03 * T0);
M0   = V0 / a0;

Tt0  = T0 * (1 + ((gamma_03 - 1)/2) * M0^2);
Pt0  = P0 * (1 + ((gamma_03 - 1)/2) * M0^2)^(gamma_03/(gamma_03 - 1));

pi_c = x(1);
if isfield(data,'eta_c'), eta_c = data.eta_c; else, eta_c = 0.90; end

Tt3 = Tt0 * pi_c^((gamma_03 - 1)/(gamma_03 * eta_c));

W_comp = burner.mdot3 * Cp_03 * (Tt3 - Tt0);

% Fan not in current framework
W_req = W_comp;

%% -------------------- Turbine power balance (Farokhi Ch.7) --------------------
eta_m = data.eta_m;

DeltaTt_total = W_req / (eta_m * mdot_t * cp);
Tt5 = Tt4 - DeltaTt_total;

%% -------------------- Turbine efficiency -> pressure ratio --------------------
if isfield(data,'eta_t'), eta_t = data.eta_t; else, eta_t = 0.90; end

% For turbine: eta_t = (Tt4 - Tt5)/(Tt4 - Tt5s)  =>  Tt5s = Tt4 - (Tt4 - Tt5)/eta_t
Tt5s = Tt4 - (Tt4 - Tt5) / eta_t;

Pt5 = Pt4 * (Tt5s/Tt4)^(gamma/(gamma - 1));

%% -------------------- Exit static (assume M5 target) --------------------
if isfield(data,'M5'), M5 = data.M5; else, M5 = 0.50; end
R_49 = cp * ((gamma - 1)/gamma);

tau5 = 1 + (gamma - 1)/2 * M5^2;
T5 = Tt5 / tau5;
P5 = Pt5 / tau5^(gamma/(gamma - 1));
a5 = sqrt(gamma * R_49 * T5);

%% -------------------- Cooling model for Twg and T_cool --------------------
% Use given blade / coolant parameters from AE430_Final.m
% Gas-side HTC is specified in project; if not in data, use 7000 W/m^2-K
if isfield(data,'h_g'), h_g = data.h_g; else, h_g = 7000; end

% Cooling inputs
c_blade = data.c;
k_w     = data.blade_kw;
t_w     = data.blade_tw;
rho_blade = data.blade_density;

Cp_cool = data.Cp_cool;
Pr      = data.Pr;
mu_cool = data.mu_cool;

if isfield(data,'T_wg_max'), T_wg_max = data.T_wg_max; else, T_wg_max = 1300; end
if isfield(data,'T_cool_max'), T_cool_max = data.T_cool_max; else, T_cool_max = 1000; end

% Geometry for cooling calc
if isfield(data,'r_h'), r_h = data.r_h; else, r_h = 0.25; end
if isfield(data,'r_t'), r_t = data.r_t; else, r_t = 0.45; end
r_m = 0.5*(r_h + r_t);

A_ht = 2 * c_blade * (r_t - r_h);     % class spec
xchar = 0.5 * (r_t - r_h);

% Coolant inlet stagnation assumed from station 3 totals
Tt_c_in = Tt3;
Pt_c_in = burner.Pt3;

% Assume a small coolant Mach (can be tuned)
if isfield(data,'M_cool'), M_cool = data.M_cool; else, M_cool = 0.10; end
gamma_c = data.gamma_cool;
R_c = Cp_cool * ((gamma_c - 1)/gamma_c);

tau_c = 1 + (gamma_c - 1)/2 * M_cool^2;
T_c_in = Tt_c_in / tau_c;
P_c_in = Pt_c_in / tau_c^(gamma_c/(gamma_c - 1));
rho_c  = P_c_in / (R_c * T_c_in);
a_c    = sqrt(gamma_c * R_c * T_c_in);
V_c    = M_cool * a_c;

% Iterate coolant fraction beta until wall temp <= limit
beta = 0.02;  % initial guess (2% of turbine flow)
T_wg = NaN;
T_c_out = NaN;

for it = 1:25
    mdot_c = beta * mdot_t;

    Re_x = rho_c * V_c * xchar / mu_cool;
    St  = 0.0296 / (Re_x^0.2);
    h_c = St * rho_c * Cp_cool * V_c;

    h_eff = 1 / ( (1/h_g) + (t_w/k_w) + (1/h_c) );

    qpp_un = h_eff * (T5 - T_c_in);     % W/m^2
    Q_un   = qpp_un * A_ht;

    Q_max  = mdot_c * Cp_cool * (T_cool_max - T_c_in);

    if Q_un > Q_max
        qpp = Q_max / A_ht;
        T_c_out = T_cool_max;
    else
        qpp = qpp_un;
        T_c_out = T_c_in + Q_un/(mdot_c * Cp_cool);
    end

    T_wg = T5 - qpp / h_g;

    if T_wg > T_wg_max
        beta = min(0.10, beta * 1.3);
    else
        beta = max(0.0,  beta * 0.9);
    end
end

%% -------------------- Tip Mach and centrifugal stress --------------------
if isfield(data,'omega'), omega = data.omega; else, omega = 800; end
U_tip = omega * r_t;
M_tip_turb = U_tip / a5;

% Very simplified centrifugal stress measure
sigma_c = rho_blade * omega^2 * r_m^2;

%% -------------------- Pack outputs --------------------
turbine = struct();

% Stations
turbine.Tt4 = Tt4; turbine.Pt4 = Pt4;
turbine.Tt5 = Tt5; turbine.Pt5 = Pt5;
turbine.T5  = T5;  turbine.P5  = P5;
turbine.M5  = M5;

% Performance
turbine.W_req = W_req;
turbine.DeltaTt_total = DeltaTt_total;
turbine.eta_t = eta_t;

% Cooling / constraints outputs
turbine.Twg = T_wg;
turbine.T_cool = T_c_out;
turbine.beta_cool = beta;

turbine.M_tip_turb = M_tip_turb;
turbine.sigma_c = sigma_c;

end
