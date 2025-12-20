function burner = burnerModel(x, N_spools, data)
% burnerModel
% ------------------------------------------------------------
% Group interface: burner = burnerModel(x, N_spools, data)
%
% Design vars assumed from your current AE430_Final.m:
%   x(1) = pi_c   (compressor total pressure ratio)
%   x(2) = Tt4    (burner exit stagnation temperature target) [K]
%
% Outputs (minimum):
%   burner.Tt4, burner.Pt4, burner.f, burner.phi, burner.X_CO, burner.Q_R
%   burner.tau_res, burner.Twg (placeholder), burner.mdot4
%   plus station-3 values used downstream: burner.Tt3, burner.Pt3, burner.mdot3
%
% Uses Farokhi Ch.7 combustor relations:
%   Pt4 = pi_b * Pt3
%   Fuel-air ratio from energy balance (approx):
%     f = cp*(Tt4 - Tt3) / (eta_b*LHV - cp*Tt4)
%
% If fuelEquilibrium + dieselEquilTable.mat exist, it uses them to estimate X_CO.

%% -------------------- Unpack design variables --------------------
pi_c = x(1);
Tt4_target = x(2);

%% -------------------- Unpack needed constants --------------------
% Cold side (03)
gamma_03 = data.gamma_03;
Cp_03    = data.Cp_03;   % J/kg-K
T0       = data.T0;
P0       = data.P0;
V0       = data.V0;

% Hot side (49) used for burner exit properties
gamma_49 = data.gamma_49;
Cp_49    = data.Cp_49;

% Burner loss + fuel
pi_b   = data.pi_b;     % combustor total pressure ratio (<=1)
LHV    = data.LHV;      % J/kg fuel
f_st   = data.f_st;     % stoichiometric FAR
if isfield(data,'eta_b'), eta_b = data.eta_b; else, eta_b = 0.99; end

% Burner geometric / constraint knobs (defaults if not yet in data)
if isfield(data,'Lb'), Lb = data.Lb; else, Lb = 1.10; end         % m
if isfield(data,'M4'), M4 = data.M4; else, M4 = 0.25; end         % -
if isfield(data,'mdot_bleed_burner'), mdot_bleed = data.mdot_bleed_burner;
else, mdot_bleed = 0.0; end

% If no inlet area is specified, pick a reasonable placeholder.
if isfield(data,'A3'), A3 = data.A3; else, A3 = 0.30; end         % m^2

%% -------------------- Station 0 -> 00 (inlet totals) --------------------
% Use the same approach as inletModel.m, but self-contained:
R_03 = Cp_03 * ((gamma_03 - 1)/gamma_03);
a0   = sqrt(gamma_03 * R_03 * T0);
M0   = V0 / a0;

Tt0  = T0 * (1 + ((gamma_03 - 1)/2) * M0^2);
Pt0  = P0 * (1 + ((gamma_03 - 1)/2) * M0^2)^(gamma_03/(gamma_03 - 1));

%% -------------------- Station 2/3 (compressor exit totals, simple) --------------------
% Until compressorModel is complete, approximate compressor exit totals.
% Farokhi-style: Tt3 = Tt0 * pi_c^((gamma-1)/(gamma*eta_c))
if isfield(data,'eta_c'), eta_c = data.eta_c; else, eta_c = 0.90; end

Pt3 = Pt0 * pi_c;
Tt3 = Tt0 * pi_c^((gamma_03 - 1)/(gamma_03 * eta_c));

% Estimate mdot at station 3 from (static) using assumed M3, A3
if isfield(data,'M3'), M3 = data.M3; else, M3 = 0.15; end

tau3 = 1 + (gamma_03 - 1)/2 * M3^2;
T3   = Tt3 / tau3;
P3   = Pt3 / tau3^(gamma_03/(gamma_03 - 1));
rho3 = P3 / (R_03 * T3);
a3   = sqrt(gamma_03 * R_03 * T3);
V3   = M3 * a3;

mdot3 = rho3 * V3 * A3;

% Apply any bleed before burner
mdot_burn_in = mdot3 - mdot_bleed;

%% -------------------- Burner pressure loss --------------------
Pt4 = pi_b * Pt3;

%% -------------------- Fuel-air ratio from energy balance --------------------
% Farokhi Ch.7 combustor energy balance (simplified constant cp form):
%   (1+f) * cp_hot * Tt4 = cp_cold * Tt3 + f * eta_b * LHV
% Solve for f:
%   f = cp_hot*(Tt4 - Tt3) / (eta_b*LHV - cp_hot*Tt4)
%
% Use Cp_49 on hot side.
numer = Cp_49 * (Tt4_target - Tt3);
denom = (eta_b * LHV - Cp_49 * Tt4_target);

if denom <= 0
    f = NaN;
else
    f = numer / denom;
end

phi = f / f_st;

% Fuel heat release per kg fuel: report LHV (or Q_R if using table)
Q_R = LHV;

%% -------------------- Optional equilibrium-based CO estimate --------------------
% If you have fuelEquilibrium.m and dieselEquilTable.mat on path, use it.
X_CO = NaN;
try
    if exist('fuelEquilibrium','file') == 2
        % fuelEquilibrium expects (Tt3, Pt3, phi)
        [Tt4_eq, X_CO_eq, Q_R_eq] = fuelEquilibrium(Tt3, Pt3, phi); %#ok<ASGLU>
        X_CO = X_CO_eq;
        Q_R  = Q_R_eq;
    end
catch
    % leave X_CO as NaN if not available
end

%% -------------------- Station 4 static state from M4 --------------------
R_49 = Cp_49 * ((gamma_49 - 1)/gamma_49);

tau4 = 1 + (gamma_49 - 1)/2 * M4^2;
T4   = Tt4_target / tau4;
P4   = Pt4 / tau4^(gamma_49/(gamma_49 - 1));
rho4 = P4 / (R_49 * T4);
a4   = sqrt(gamma_49 * R_49 * T4);
V4   = M4 * a4;

mdot4 = (1 + f) * mdot_burn_in;
A4    = mdot4 / (rho4 * V4);

%% -------------------- Residence time constraint --------------------
tau_res = Lb / V4;

%% -------------------- Pack outputs --------------------
burner = struct();

% Station 3 totals for downstream use
burner.Tt3   = Tt3;
burner.Pt3   = Pt3;
burner.mdot3 = mdot3;
burner.M3    = M3;
burner.A3    = A3;

% Station 4 totals and statics
burner.Tt4   = Tt4_target;
burner.Pt4   = Pt4;
burner.T4    = T4;
burner.P4    = P4;
burner.M4    = M4;
burner.rho4  = rho4;
burner.V4    = V4;
burner.A4    = A4;
burner.mdot4 = mdot4;

% Combustor performance
burner.f     = f;
burner.phi   = phi;
burner.X_CO  = X_CO;
burner.Q_R   = Q_R;

% Constraint support
burner.tau_res = tau_res;

% Placeholder for Twg (wall temperature) – computed in turbine cooling
burner.Twg = NaN;

end
