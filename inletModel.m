function inlet = inletModel(x, N_spools, data)
    
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



end