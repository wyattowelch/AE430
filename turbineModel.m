function turbine = turbineModel(x, N_spools, data)
% turbineModel(x, N_spools, data)
% ============================================================
% AE 430 Final Project – Consolidated Turbine Model (single file)
%
% What this file includes (merged + cleaned):
%   1) Stage-by-stage total temperature / pressure march
%   2) Mechanical efficiency eta_m = 0.98 (power delivered = eta_m * turbine_shaft_power)
%   3) Uncooled turbine efficiency from (surrogate) Ainley profile loss only
%   4) Cooling model per project sheet:
%        - Ag = Aw = Ac = 2*c*(rt-rh)
%        - hg = 7000 W/m^2-K
%        - hc from flat-plate Stanton: St = 0.0296 Re_x^-0.2 Pr^-2/3
%        - Qdot_cool = cp_c * mdot_c * (Tc_max - Tc_in)
%        - Wall constraint: Twg <= 1300 K
%        - Coolant constraint: Tc_out <= 1000 K
%        - Turbine efficiency penalty: eta_t drops by 3.2% per 1% coolant fraction
%   5) Constraints checks:
%        - Relative tip Mach <= 1.2
%        - Centrifugal stress <= 70 MPa
%        - Optional rotor exit beta guideline (>= -70 deg)
%
% Notes:
% - This model is written to be robust to different project codebases.
% - It assumes you have compressorModel(...) and burnerModel(...) in your path.
% - If your burnerModel supports a "bleed" input, this model will attempt to pass it.
%
% Output "turbine" is a struct with commonly-used fields:
%   Tt5, Pt5, T5, P5, mdot5, beta_cool, mdot_cool, eta_uncooled, eta_eff,
%   stage(i).{Tt_in,Tt_out,Pt_in,Pt_out,DeltaTt,pi,alpha1,beta2,Mtip_r,sigma_c,R_est}
%
% ============================================================

%% -----------------------------
% 0) Basic required data fields
% ------------------------------
must = {'N_turb_stages','c','r_h','r_t','gamma_g','cp_g','R_g'};
for k = 1:numel(must)
    if ~isfield(data, must{k})
        error('turbineModel:MissingData','Missing data.%s', must{k});
    end
end

% mechanical efficiency (project)
eta_m = 0.98;

% geometry
N_stages = data.N_turb_stages;
c = data.c;              % chord (core blades: 0.085 m typical)
r_h = data.r_h;
r_t = data.r_t;
r_m = 0.5*(r_h + r_t);

% gas properties (after burner)
gamma_g = data.gamma_g;
cp_g    = data.cp_g;
Rg      = data.R_g;

% blade / stress properties
rho_blade = getOr(data,'blade_density',8000); % kg/m^3
sigma_max = getOr(data,'sigma_max',70e6);     % Pa

% optional guideline
beta2_min_deg = getOr(data,'beta2_min_deg',-70); % deg guideline (not hard constraint)

%% -----------------------------
% 1) Upstream states
% ------------------------------
% Compressor (for spool speed and required power)
if exist('compressorModel','file') ~= 2
    error('turbineModel:MissingFunction','compressorModel.m not found on MATLAB path.');
end
compressor = compressorModel(x, N_spools, data);

% Burner / combustor (turbine inlet)
if exist('burnerModel','file') ~= 2
    error('turbineModel:MissingFunction','burnerModel.m not found on MATLAB path.');
end

burner = burnerModel(x, N_spools, data);

% turbine inlet totals (station 4)
Tt4 = burner.Tt4;
Pt4 = burner.Pt4;
mdot4_nom = burner.mdot4;   % hot-gas flow leaving combustor (no turbine coolant added yet)

% spool speed: try compressor output first, else data
omega = [];
if isfield(compressor,'omega') && ~isempty(compressor.omega)
    omega = compressor.omega;
elseif isfield(data,'omega') && ~isempty(data.omega)
    omega = data.omega;
elseif isfield(data,'omega_turb') && ~isempty(data.omega_turb)
    omega = data.omega_turb;
else
    error('turbineModel:MissingOmega','Need omega from compressor.omega or data.omega/data.omega_turb.');
end

U_m   = omega * r_m;
U_tip = omega * r_t;

% axial velocity (constant Cz)
Cz = getOr(data,'Cz_turb',150); % if not specified, pick a reasonable design value

%% -----------------------------
% 2) Required turbine shaft work
% ------------------------------
% Your codebase may store required compressor power in different fields.
% Priority:
%   data.P_required_shaft (W) if provided
%   compressor.P_required_shaft
%   compressor.W_comp * mdot3 (if available)
P_required_shaft = [];
if isfield(data,'P_required_shaft') && ~isempty(data.P_required_shaft)
    P_required_shaft = data.P_required_shaft;
elseif isfield(compressor,'P_required_shaft') && ~isempty(compressor.P_required_shaft)
    P_required_shaft = compressor.P_required_shaft;
elseif isfield(compressor,'W_comp') && isfield(compressor,'mdot3')
    P_required_shaft = compressor.W_comp * compressor.mdot3;
elseif isfield(compressor,'P_comp')
    P_required_shaft = compressor.P_comp;
else
    error(['turbineModel:MissingPowerRequirement ' ...
           'Provide data.P_required_shaft or compressor.{P_required_shaft, W_comp & mdot3, or P_comp}.']);
end

% turbine shaft power needed BEFORE mech losses
P_turb_shaft_needed = P_required_shaft / eta_m;

%% -----------------------------
% 3) Cooling model (solves beta_c)
% ------------------------------
cool = struct();
beta_c = 0.0;
mdot_c = 0.0;

cool.enable = getOr(data,'cooling_enable',true);

if cool.enable
    % Required cooling inputs (project defaults if not supplied)
    cool.hg      = getOr(data,'hg',7000);      % W/m^2-K
    cool.kw      = getOr(data,'kw',50);        % W/m-K
    cool.tw      = getOr(data,'tw',0.002);     % m
    cool.Twg_max = getOr(data,'Twg_max',1300); % K
    cool.Tc_max  = getOr(data,'Tc_max',1000);  % K

    cool.gamma_c = getOr(data,'gamma_c',1.4);
    cool.cp_c    = getOr(data,'cp_c',1000);    % J/kg-K
    cool.Pr      = getOr(data,'Pr_c',0.7);
    cool.mu      = getOr(data,'mu_c',4e-5);    % kg/(m*s)

    % Coolant comes from compressor exit (station 3)
    % Expect compressor.Tt3, compressor.Pt3, compressor.mdot3
    if ~isfield(compressor,'Tt3') || ~isfield(compressor,'Pt3')
        error('turbineModel:MissingCompressorCoolingState','Need compressor.Tt3 and compressor.Pt3 for cooling.');
    end
    cool.Tt3 = compressor.Tt3;
    cool.Pt3 = compressor.Pt3;

    % Need an estimate of coolant velocity in channels. If not provided, assume ~50 m/s.
    cool.Vc = getOr(data,'Vc_cool',50);

    % Determine beta_c such that Twg <= Twg_max AND Tc_out <= Tc_max.
    % Use the turbine inlet temperature as the hot-gas driving temp.
    Tg_for_cooling = Tt4;

    [beta_c, mdot_c, cool] = solveCoolingBeta( ...
        mdot4_nom, Tg_for_cooling, r_h, r_t, c, cool);

    % If your project intends coolant to be subtracted from combustor inlet mdot,
    % optionally rerun burner with a bleed if burnerModel supports it.
    % We'll attempt it safely: if it errors, we fall back to mdot4_nom.
    if isfield(data,'rerun_burner_with_bleed') && data.rerun_burner_with_bleed
        try
            data2 = data;
            data2.mdot_bleed_cool = mdot_c; %#ok<STRNU>
            burner2 = burnerModel(x, N_spools, data2);
            Tt4 = burner2.Tt4;
            Pt4 = burner2.Pt4;
            mdot4_nom = burner2.mdot4;
            Tg_for_cooling = Tt4;
            % re-solve cooling once with updated gas temp / flow
            [beta_c, mdot_c, cool] = solveCoolingBeta( ...
                mdot4_nom, Tg_for_cooling, r_h, r_t, c, cool);
        catch
            % do nothing (keep the first-pass values)
        end
    end
end

% turbine gas massflow AFTER accounting for cooling bleed from combustor inlet
% (per project note: subtract coolant from burner inlet mdot)
mdot_g = max(mdot4_nom - mdot_c, 1e-6);

%% -----------------------------
% 4) Turbine efficiency model (Ainley only + cooling penalty)
% ------------------------------
% Uncooled loss from Ainley profile loss surrogate (your own curve-fit / lookup)
s_over_c  = getOr(data,'turb_s_over_c',0.8);
beta3_deg = getOr(data,'turb_beta3_deg',60);
bladeType = getOr(data,'turb_bladeType','reaction'); % 'reaction' or 'impulse'

Yp = ainleyProfileLossSurrogate(bladeType, s_over_c, beta3_deg);

% "small loss" approximation: eta_uncooled ≈ 1 - Yp, clamped
eta_uncooled = max(0.70, min(0.95, 1 - Yp));

% Cooling penalty: eta drops by 3.2% per 1% coolant fraction (fraction of turbine flow)
% beta_c is mdot_c / mdot4_nom (fraction of turbine inlet flow BEFORE subtracting from combustor)
eta_eff = eta_uncooled * (1 - 3.2*beta_c);
eta_eff = max(0.05, eta_eff);

%% -----------------------------
% 5) Determine required total ΔTt across turbine
% ------------------------------
w_needed_actual = P_turb_shaft_needed / mdot_g;       % J/kg hot gas through turbine
DeltaTt_actual_total = w_needed_actual / cp_g;        % K (actual)
DeltaTt_is_total = DeltaTt_actual_total / eta_eff;    % isentropic equivalent total drop

DeltaTt_is_stage = DeltaTt_is_total / N_stages;

%% -----------------------------
% 6) Stage-by-stage march
% ------------------------------
Tt_stage = zeros(N_stages+1,1);
Pt_stage = zeros(N_stages+1,1);

Tt_stage(1) = Tt4;
Pt_stage(1) = Pt4;

stage = repmat(struct(), N_stages, 1);

MtipOK = true; sigmaOK = true; reactionOK = true; beta2OK = true;

for i = 1:N_stages
    Tt1 = Tt_stage(i);
    Pt1 = Pt_stage(i);

    % Isentropic exit total temperature for this stage
    Tt2s = Tt1 - DeltaTt_is_stage;

    % Actual exit total temperature using turbine efficiency definition:
    % eta_t = (Tt1 - Tt2)/(Tt1 - Tt2s)
    Tt2 = Tt1 - eta_eff*(Tt1 - Tt2s);

    % Stage total pressure ratio (turbine: <1)
    pi_stage = (Tt2s/Tt1)^(gamma_g/(gamma_g-1));
    Pt2 = Pt1 * pi_stage;

    % Store
    Tt_stage(i+1) = Tt2;
    Pt_stage(i+1) = Pt2;

    % --- Meanline velocity triangle (design-level) ---
    % Euler turbine: Δh = U*(Cθ1 - Cθ2). Assume Cθ1 ~ 0 to solve Cθ2.
    Deltaht = cp_g*(Tt1 - Tt2);   % actual enthalpy drop J/kg
    Ctheta1 = 0;
    Ctheta2 = Ctheta1 - Deltaht/max(U_m,1e-9);

    alpha1 = atan2(Ctheta1, Cz);
    alpha2 = atan2(Ctheta2, Cz);

    Wtheta2 = Ctheta2 - U_m;
    beta2 = atan2(Wtheta2, Cz);

    % Degree of reaction estimate (common axial-stage approx)
    R_est = 0.5 + (Ctheta2 + Ctheta1)/(2*max(U_m,1e-9));

    if (R_est < -1e-6) || (R_est > 1+1e-6)
        reactionOK = false;
    end

    if rad2deg(beta2) < beta2_min_deg
        beta2OK = false;
    end

    % Tip relative Mach check per equation sheet form
    a_guess = sqrt(gamma_g*Rg*Tt1);  % conservative using total temp
    Mz = Cz/max(a_guess,1e-9);
    MT_tip = U_tip/max(a_guess,1e-9);
    Mtip_r = sqrt(Mz^2 + (MT_tip - Mz*tan(alpha1))^2);
    if Mtip_r > 1.2 + 1e-9
        MtipOK = false;
    end

    % Centrifugal stress (order-of-mag ring/blade-root estimate)
    sigma_c = rho_blade * omega^2 * (r_t^2 - r_h^2)/2;
    if sigma_c > sigma_max + 1e-6
        sigmaOK = false;
    end

    stage(i).Tt_in  = Tt1;
    stage(i).Tt_out = Tt2;
    stage(i).Pt_in  = Pt1;
    stage(i).Pt_out = Pt2;
    stage(i).DeltaTt = (Tt1 - Tt2);
    stage(i).pi = pi_stage;
    stage(i).Ctheta1 = Ctheta1;
    stage(i).Ctheta2 = Ctheta2;
    stage(i).alpha1_deg = rad2deg(alpha1);
    stage(i).alpha2_deg = rad2deg(alpha2);
    stage(i).beta2_deg  = rad2deg(beta2);
    stage(i).R_est = R_est;
    stage(i).Mtip_r = Mtip_r;
    stage(i).sigma_c = sigma_c;
end

Tt5 = Tt_stage(end);
Pt5 = Pt_stage(end);

%% -----------------------------
% 7) Static exit estimate (if you specify exit Mach)
% ------------------------------
M5 = getOr(data,'M5_turb_exit',0.4);
tau5 = 1 + (gamma_g-1)/2*M5^2;
T5 = Tt5 / tau5;
P5 = Pt5 / tau5^(gamma_g/(gamma_g-1));

%% -----------------------------
% 8) Pack outputs
% ------------------------------
turbine = struct();
turbine.Tt5 = Tt5;
turbine.Pt5 = Pt5;
turbine.T5  = T5;
turbine.P5  = P5;

turbine.mdot5 = mdot_g;

turbine.eta_m = eta_m;
turbine.eta_uncooled = eta_uncooled;
turbine.eta_eff = eta_eff;

turbine.beta_cool = beta_c;
turbine.mdot_cool = mdot_c;
turbine.cool = cool;

turbine.P_required_shaft = P_required_shaft;
turbine.P_turb_shaft_needed = P_turb_shaft_needed;
turbine.P_delivered = eta_m * (mdot_g * cp_g * (Tt4 - Tt5)); % delivered to compressor

turbine.stage = stage;

% Constraints
turbine.constraints = struct();
turbine.constraints.tipMach_le_1p2 = MtipOK;
turbine.constraints.sigma_le_70MPa = sigmaOK;
turbine.constraints.reaction_0_to_1 = reactionOK;
turbine.constraints.beta2_guideline = beta2OK;
turbine.constraints.coolant_T_le_1000K = (~cool.enable) || (cool.Tc_out <= cool.Tc_max + 1e-6);
turbine.constraints.wall_T_le_1300K = (~cool.enable) || (cool.Twg <= cool.Twg_max + 1e-6);

turbine.success = all(struct2array(turbine.constraints));

end

%% ============================================================
% Helper: safe get with default
% ============================================================
function v = getOr(s, field, default)
if isfield(s, field) && ~isempty(s.(field))
    v = s.(field);
else
    v = default;
end
end

%% ============================================================
% Cooling solver: find coolant fraction beta_c to satisfy constraints
% ============================================================
function [beta_c, mdot_c, cool] = solveCoolingBeta(mdot4, Tg, r_h, r_t, c, cool)
% Returns:
%   beta_c = mdot_c / mdot4
%   mdot_c = coolant massflow [kg/s]
%   cool struct updated with Twg, Tc_out

Ag = 2*c*(r_t - r_h);
Aw = Ag;
Ac = Ag;

hg = cool.hg;
kw = cool.kw;
tw = cool.tw;

% characteristic length for Re_x per project: x = (rt-rh)/2
x = (r_t - r_h)/2;

% Use an iterative approach on beta_c
beta_c = max(0, min(0.15, getOr(cool,'beta_guess',0.01)));

for it = 1:12
    mdot_c = beta_c * mdot4;

    % coolant inlet static temperature: T = Tt - V^2/(2cp)
    Tc_in = cool.Tt3 - cool.Vc^2/(2*cool.cp_c);
    Tc_in = max(1, Tc_in);

    % coolant-side convection coefficient hc from Stanton number (flat plate)
    % Approximate rho using ideal gas with R~287 if no R provided
    R_c = getOr(cool,'R_c',287);
    rho_c = cool.Pt3 / (R_c * cool.Tt3);
    Re_x = rho_c * cool.Vc * x / max(cool.mu, 1e-12);

    St = 0.0296 * Re_x^(-0.2) * cool.Pr^(-2/3);
    hc = St * rho_c * cool.Vc * cool.cp_c;

    % Overall UA (series resistances)
    UA = 1 / ( 1/(hg*Ag) + tw/(kw*Aw) + 1/(hc*Ac) );

    % coolant heat capacity cap per project
    Qdot_cap = cool.cp_c * max(mdot_c,1e-12) * max(cool.Tc_max - Tc_in, 0);

    % maximum heat transfer available from UA and temperature difference
    Qdot_UA = UA * max(Tg - Tc_in, 0);

    % actual removable heat
    Qdot = min(Qdot_cap, Qdot_UA);

    % estimate wall temperature on gas side
    Twg = Tg - Qdot*(1/(hg*Ag) + tw/(kw*Aw));
    Tc_out = Tc_in + Qdot/(max(mdot_c,1e-12)*cool.cp_c);

    cool.Twg = Twg;
    cool.Tc_out = Tc_out;

    okWall = (Twg <= cool.Twg_max);
    okCool = (Tc_out <= cool.Tc_max);

    if okWall && okCool
        % try to reduce beta slightly to tighten
        beta_c = beta_c * 0.85;
    else
        % increase beta
        beta_c = min(0.25, beta_c*1.30 + 1e-4);
    end
end

% finalize
beta_c = max(0, beta_c);
mdot_c = beta_c * mdot4;

end

%% ============================================================
% Ainley profile loss surrogate (single-file stand-in)
% ============================================================
function Yp = ainleyProfileLossSurrogate(bladeType, s_over_c, beta3_deg)
% Very compact surrogate (not a full Ainley method):
% - Captures: higher loss at very low/high s/c, and loss increasing with beta3.
% - Keeps your project code from crashing; you can later replace this with
%   your exact Fig 10.16–10.18 digitization.

s = s_over_c;
beta3 = abs(beta3_deg);

if strcmpi(bladeType,'impulse')
    % Impulse-like blades: generally higher losses
    Ymin = interp1([40 50 60 70],[0.030 0.040 0.055 0.070],beta3,'linear','extrap');
    sopt = 0.85;
    Y30 = 0.20;
else
    % Reaction blades
    Ymin = interp1([40 50 60 70],[0.015 0.020 0.030 0.040],beta3,'linear','extrap');
    sopt = 0.80;
    Y30 = 0.065;
end

k = (Y30 - Ymin) / (0.3 - sopt)^2;
Yp = Ymin + k*(s - sopt)^2;

Yp = max(0.005, min(0.25, Yp)); % clamp
end
