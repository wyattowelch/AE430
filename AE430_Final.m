%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% AE430 Final Project
% Wyatt, Ethan, Nick, Jayden, Luis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define set constants

data = struct()

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
data.gamma_cool = gamma.gamma_03;
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


% Define design variables

x = [pi_c, Tt5];

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

x0 = [20; 1700];
lb = [5; 1200];
ub = [60; 1900];

for N_spools = [1 2 3]
    fun     = @(x) costFun(x, N_spools, data);
    nonlcon = @(x) constraints(x, N_spools, data);
    [x_opt, fval] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

    if fval < best_f
        best_f = fval;
        best_x = x_opt;
        best_N = N_spools;
    end
end
