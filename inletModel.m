function inlet = inletModel(x, N_spools, data, nozzle)
    
    % Simplify variables
    gam = data.gamma_03;

    % Initial Free-air calcs
    R = data.Cp_03 * ((gam - 1) / gam); % J/kgK
    a0 = sqrt(gam * R * data.T0); % m/s
    M0 = data.V0 / a0; 
    inlet.Cv_03 = data.Cp_03 - R;

    % Stagnation Pressure / Temperature
    inlet.Tt0 = data.T0 * (1 + ((gam - 1) / 2) * M0 ^ 2);
    inlet.Pt0 = data.P0 * (1 + ((gam - 1) / 2) * M0 ^ 2) ^ (gam / (gam - 1));

    % Stage 1/2 Vars
    inlet.Tt1 = inlet.Tt0;
    inlet.T1 = inlet.T0;
    inlet.Pt1 = inlet.Pt0;
    inlet.P1 = inlet.P0;

    inlet.Tt2 = inlet.Tt1;
    inlet.T2 = inlet.T1;
    inlet.Pt2 = inlet.Pt2;
    inlet.P2 = inlet.P2;

    % m' and A1



end