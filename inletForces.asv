function [D_add, F_lip] = inletForces(mdot0, data, inlet)

gamma = data.gamma_03;
R     = data.R_03;

T0 = data.T0;  P0 = data.P0;  V0 = data.V0;
A1 = inlet.A1;

% freestream stagnation
a0  = sqrt(gamma*R*T0);
M0  = V0/a0;
Tt0 = T0*(1 + (gamma-1)/2*M0^2);
Pt0 = P0*(1 + (gamma-1)/2*M0^2)^(gamma/(gamma-1));

% solve for highlight Mach from mdot relation (subsonic)
mfp = @(M) M .* (1 + (gamma-1)/2*M.^2).^(-(gamma+1)/(2*(gamma-1)));
K   = Pt0*A1/sqrt(max(Tt0,1)) * sqrt(gamma/R);
fun = @(M) K*mfp(M) - mdot0;
M1 = fzero(fun,[1e-6 0.99]);

P1 = Pt0 / (1 + (gamma-1)/2*M1^2)^(gamma/(gamma-1));

% lip suction
F_lip = (P0 - P1) * A1;

% spillage/additive drag
rho0 = P0/(R*T0);
q0   = 0.5*rho0*V0^2;
A_cap = mdot0/(rho0*V0);

Cspill = 1.0;
D_add = max(0, Cspill*q0*(A_cap - A1));

end
