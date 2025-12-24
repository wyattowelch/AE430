function nozzle = nozzleModel(x, N_spools, data, inlet, comp, burner, turbine)
    
    % nozzle = struct containing full nozzle solution

%Ambient conditions
Pa = data.P0;
Ta = data.T0;

% Nozzle inlet (station 9, turbine exit)
Tt9 = turbine.Tt5; % K
Pt9 = turbine.Pt5; % Pa
mdot = turbine.mdot5; % kg/s

gamma = data.gamma_49;
Cp = data.Cp_49;
R = Cp*(gamma-1)/gamma;

eta_n = data.eta_n;

% Critical pressure ratio (choking check)
PR_crit = (2/(gamma+1))^(gamma/(gamma-1));

% Assume isentropic nozzle efficiency
Pt9_eff = Pt9 * eta_n;

% --- hard guard: if total pressure <= ambient (or non-finite), no real expansion ---
if ~isfinite(Pt9_eff) || ~isreal(Pt9_eff) || Pt9_eff <= Pa*(1+1e-6)
    nozzle.choked = false;

    M_e = 0;
    Te  = Tt9;
    Pe  = Pa;
    Ve  = 0;

    % big areas to penalize A_max and length (keeps objective real/finite)
    At = 1e4;
    Ae = 1e4;
    return;   % IMPORTANT: prevents later code from overwriting these safe values
end

% -------------------------
% Check for choking (now safe because Pt9_eff > Pa)
% -------------------------
if (Pa/Pt9_eff <= PR_crit)
    % -------- CHOKED FLOW --------
    nozzle.choked = true;

    M_e = 1.0;

    Te = Tt9*(2/(gamma+1));
    Pe = Pt9_eff*PR_crit;
    Ve = sqrt(gamma*R*max(Te,1));

    % Throat area
    At = mdot*sqrt(max(Tt9,1))/(Pt9_eff*sqrt(gamma/R)) * ...
        ((gamma+1)/2)^((gamma+1)/(2*(gamma-1)));

    Ae = At;

else
    % -------- UNCHOKED FLOW --------
    nozzle.choked = false;

    % Guard against tiny negative due to roundoff
    term = (Pt9_eff/Pa)^((gamma-1)/gamma) - 1;
    term = max(term, 0);

    M_e = sqrt((2/(gamma-1)) * term);

    % Avoid divide-by-zero if M_e is extremely small
    M_eff = max(M_e, 1e-6);

    Te = Tt9/(1 + (gamma-1)/2*M_e^2);
    Ve = M_eff*sqrt(gamma*R*max(Te,1));
    Pe = Pa;

    % Exit area
    Ae = mdot*sqrt(max(Te,1))/(Pe*M_eff*sqrt(gamma/R));

    % Throat (virtual)
    At = Ae / ((1/M_eff) * ...
        ((2/(gamma+1))*(1+(gamma-1)/2*M_e^2))^((gamma+1)/(2*(gamma-1))));
end


% Length
theta_deg = 7;
theta = theta_deg*pi/180;

rt = sqrt(max(At,1)/pi);
re = sqrt(max(Ae,1)/pi);

nozzle.l_n = abs(re - rt)/max(tan(theta),1e-9);


% Thrust calculation
F = mdot*Ve + (Pe - Pa)*Ae;

% Populate output structure
nozzle.gamma = gamma;
nozzle.R = R;

nozzle.M_e = M_e;
nozzle.T_e = Te;
nozzle.P_e = Pe;
nozzle.V_e = Ve;

nozzle.A_t = At;
nozzle.A_e = Ae;

nozzle.mdot = mdot;
nozzle.thrust = F;
nozzle.thrust_error = F - data.F_req;

end