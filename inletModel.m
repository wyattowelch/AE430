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

    % Stage 1/2 Vars
    inlet.Tt1 = inlet.Tt0;
    inlet.T1 = data.T0;
    inlet.Pt1 = inlet.Pt0;
    inlet.P1 = data.P0;
    M1 = M0;

    inlet.Tt2 = inlet.Tt1;
    inlet.T2 = inlet.T1;
    inlet.Pt2 = inlet.Pt1;
    inlet.P2 = inlet.P1;


    % Geometry assumptions % Fig 6.11c
    AR = 2;
    if isfield(data, 'inlet_AR')
        AR = data.inlet_AR;
    end
    if AR <= 1
        AR = 1.001;
    end

    phi_deg = 7;
    if isfield(data, 'inlet_phi_deg')
        phi_deg = data.inlet_phi_deg
    end
    phi = phi_deg * pi / 180;
    if phi <= 0
        phi = 1e-6;
    end

    w1 = .25; 
    if isfield(data, 'inlet_w1')
        w1 = data.inlet_w1;
    end
    if w1 <= 0
        w1 = 1e-6;
    end

    rbar = .75;
    if isfield(data, 'inlet_rbar')
        rbar = data.inlet_rbar;
    end

        % Height at diffuser exit
    w2 = AR * w1;
    dw = w2 - w1;

        % radii at 1 and 2
    Ri1 = rbar - w1 / 2;
    Ro1 = rbar + w1 / 2;

    Ri2 = rbar - w2 / 2;
    Ro2 = rbar + w2 / 2;

    DeltaR1 = Ro1 - Ri1;
    DeltaR2 = Ro2 - Ri2;

    dRo = Ro2 - Ro1;
    dRi = Ri2 - Ri1;
    
        % axial length from half angle
    L_ax = abs(dRo) / tan(phi);

        % Wall length
    Lo = sqrt(L_ax ^ 2 + dRo ^ 2);
    Li = sqrt(L_ax ^ 2 + dRi ^ 2);

    Lbar = .5 * (Li + Lo);

        % Areas
    A1 = pi * (Ro1 ^ 2 - Ri1 ^ 2);
    A2 = pi * (Ro2 ^ 2 - Ri2 ^ 2);


    x611 = Lbar / DeltaR1;
    y611 = AR - 1;

    % Geometry outputs
    inlet.diffuserType = 'annular';
    
    inlet.AR = AR;
    inlet.AR_minus_1 = y611;
    
    inlet.w1 = w1;
    inlet.w2 = w2;
    
    inlet.Ri1 = Ri1; inlet.Ro1 = Ro1;
    inlet.Ri2 = Ri2; inlet.Ro2 = Ro2;
    
    inlet.DeltaR1 = DeltaR1;
    
    inlet.L_ax  = L_ax;
    inlet.Li    = Li;
    inlet.Lo    = Lo;
    inlet.Lbar  = Lbar;
    
    inlet.A1 = A1;
    inlet.A2 = A2;
    
    inlet.fig611c_x = x611;
    inlet.fig611c_y = y611;
    
    inlet.phi_deg = phi_deg;

    % Diffuser Performance
        % Cp from fig 6.11c is found to be 0.6

    Cp_rec = .6;
    if isfield(data, 'inlet_Cp')
        Cp_rec = data.inlet_Cp;
    end
    inlet.Cp_rec = Cp_rec;

    if ~isnan(Cp_rec)
        P2 = inlet.P1 + Cp_rec * (inlet.Pt1 - inlet.P1);
    else
        P2 = inlet.P1;
    end

    inlet.V1 = data.V0;
    V2 = inlet.V1 / AR;

    Tt2 = inlet.Tt1;

    T2 = Tt2 - V2 ^ 2 / (2 * data.Cp_03);
    if T2 <= 1
        T2 = 1;
    end

    a2 = sqrt(gam * R  * T2);
    M2 = V2 / a2;

    % Packaging

    inlet.V2  = V2;
    inlet.a2  = a2;
    inlet.M2  = M2;

    inlet.pi_d = inlet.Pt2 / inlet.Pt0;
    
    inlet.D_add = 0;
    inlet.F_lip = 0;

end