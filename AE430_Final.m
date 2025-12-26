%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% AE430 Final Project
% Wyatt, Ethan, Nick, Jayden, Luis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
% Define set constants
data = struct();

data.F_req = 350e3; % N


    % Blades
data.blade_density = 8e3; % kg/m^3
data.blade_kw = 50; % W/mK
data.blade_tw = 2e-3; % m
data.c = 0.085; % m

    % Cold air
data.gamma_03 = 1.4;
data.Cp_03 = 1000; % J/kgK
data.T0 = 220; % K
data.P0 = 23000; % Pa
data.V0 = 250; % m/s

    % Hot air
data.gamma_49 = 1.3;
data.Cp_49 = 1250; % J/kgK

    % Coolant
data.gamma_cool = data.gamma_03;
data.Cp_cool = data.Cp_03; % J/kgK
data.Pr = 0.7;
data.mu_cool = 4e-5; % kg/ms
% Use Pt3 and Tt3 as stagnation pressure and temperature

    % Fuel
data.LHV = 43.3e6; %J/kg
data.f_st = 0.068;

    % Efficiencies
data.pi_b = 0.98;
data.eta_m = 0.98;
data.eta_n = 0.95;
data.epsilon = 0.55 * 9.81; % Tip clearance parameter is epsilon/g = 0.55

% Repeat data (lazy)
% --- Derived gas constants ---
data.R_03 = data.Cp_03*(data.gamma_03-1)/data.gamma_03;
data.R_49 = data.Cp_49*(data.gamma_49-1)/data.gamma_49;

% --- Alias to match turbineModel "must" fields ---
data.gamma_g = data.gamma_49;
data.cp_g    = data.Cp_49;
data.R_g     = data.R_49;

% Cooling aliases expected by turbineModel
data.gamma_c = data.gamma_cool;
data.cp_c    = data.Cp_cool;
data.Pr_c    = data.Pr;
data.mu_c    = data.mu_cool;

% blade thermal aliases expected by turbineModel
data.kw = data.blade_kw;
data.tw = data.blade_tw;

% stress limit (project)
data.sigma_max = 70e6;

data.N_turb_stages = 3;

% Constant axial velocities (recommended constant Cz)
data.Cz_comp = 150;   % m/s
data.Cz_turb = 150;   % m/s

% Pick an annulus (hub/tip). Start reasonable; later you should tie this to flow areas.
data.r_h_comp = 0.25;  % m
data.r_t_comp = 0.50;  % m
data.r_h      = data.r_h_comp;  % turbineModel uses data.r_h/data.r_t
data.r_t      = data.r_t_comp;

% Pick a shaft speed. Start conservative; you can also compute from the tip-Mach limit (below).
data.omega_comp = 600;  % rad/s
data.omega_turb = 600;  % rad/s
data.omega      = data.omega_comp;  % turbineModel will use data.omega

a_ref = sqrt(data.gamma_03 * data.R_03 * data.T0); % conservative (coldest)
Mz    = data.Cz_comp / a_ref;
omega_max = (a_ref / data.r_t_comp) * sqrt(max(1.2^2 - Mz^2, 0.01));
data.omega_comp = 0.9 * omega_max;
data.omega_turb = data.omega_comp;
data.omega      = data.omega_comp;
data.eta_c = .85;
data.burner_phi_deg = 7;

% Define design variables

%x = [pi_c, Tt5];

x = [1 2];
x0 = [20; 1700];
    % pi_c; Tt5


    % Bounds

    lb = [5; 1200];
    ub = [60; 1900];

    % Linear constraints

    A = []; 
    b = [];

    Aeq = [];
    beq = [];


% Calculate minimum f_c for a set of N_spools

best_f  = Inf;
best_x  = [];
best_N  = [];

x0 = [10; 1750];
lb = [5; 1200];
ub = [60; 1900];

[final, results] = engineModel(x, N_spools, data);
finalT = struct2table(final);
disp(finalT)

for N_spools = [1 2 3]

    % dbstop if error
    fun     = @(x) costFun(x, N_spools, data);
    nonlcon = @(x) constraints(x, N_spools, data);

    % fun(x0)          % <-- this will show exactly where it fails
    % nonlcon(x0)      % <-- then check constraints


    % [x_opt, fval] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon);

    if fval < best_f
        best_f = fval;
        best_x = x_opt;
        best_N = N_spools;
    end


end

