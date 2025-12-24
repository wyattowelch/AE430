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

% Check for choking
if (Pa/Pt9_eff <= PR_crit)
% -------- CHOKED FLOW --------
nozzle.choked = true;

M_e = 1.0;

Te = Tt9*(2/(gamma+1));
Pe = Pt9_eff*PR_crit;
Ve = sqrt(gamma*R*Te);

% Throat area
At = mdot*sqrt(Tt9)/(Pt9_eff*sqrt(gamma/R)) * ...
    ((gamma+1)/2)^((gamma+1)/(2*(gamma-1)));

Ae = At;
else
% -------- UNCHOKED FLOW --------
nozzle.choked = false;

% Exit Mach from pressure ratio
M_e = sqrt((2/(gamma-1))* ...
    ((Pt9_eff/Pa)^((gamma-1)/gamma)-1));

Te = Tt9/(1 + (gamma-1)/2*M_e^2);
Ve = M_e*sqrt(gamma*R*Te);
Pe = Pa;

% Exit area
Ae = mdot*sqrt(Te)/(Pe*M_e*sqrt(gamma/R));

% Throat (virtual)
At = Ae / ((1/M_e)* ...
    ((2/(gamma+1))*(1+(gamma-1)/2*M_e^2))^ ...
    ((gamma+1)/(2*(gamma-1))));
end

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