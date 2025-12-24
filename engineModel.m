function results = engineModel(x, N_spools, data)

    % Stage Calculations
    
    inlet = inletModel(x, N_spools, data);
    comp = compressorModel(x, N_spools, data, inlet);
    burner = burnerModel(x, N_spools, data, inlet, comp);
    turbine = turbineModel(x, N_spools, data, inlet, comp, burner);
    nozzle = nozzleModel(x, N_spools, data, inlet, comp, burner, turbine);

    mdot0 = burner.mdot3;                 % first pass
    [D_add, F_lip] = inletForces(mdot0, data, inlet);
    

    
    % Compute outputs & Unpack   

    results.A_max_eng = max([inlet.A1, inlet.A2, burner.A3, nozzle.A_t, nozzle.A_e]);

    results.D_add = D_add;
    results.F_lip = F_lip;

    results.F_un = nozzle.thrust - inlet.mdot0 * data.V0;
    results.F_in = results.F_un - D_add + F_lip;

    mdot_fuel = burner.f * burner.mdot3;
    results.eta_0 = results.F_in * data.V0 / (mdot_fuel * data.LHV);

    results.Tt4 = burner.Tt4;
    results.X_CO = burner.X_CO;
    results.tau_res = burner.tau_res;

    results.Twg = turbine.cool.Twg;
    results.T_cool = turbine.cool.Tc_out;

    results.sigma_c = turbine.sigma_c;
    results.M_tip_turb = turbine.M_tip_turb;
    results.M_tip_comp = comp.M_tip_comp;
    results.reactionVals = [comp.reactionVals; turbine.reactionVals];
    results.stalls = comp.stalls;

    l_id = inlet.L_ax;
    l_c = comp.l_c;
    l_b = burner.l_b;
    l_bd = burner.l_bd;
    l_t = turbine.l_t;
    l_n = nozzle.l_n;

    results.l_tot = l_id + l_c + l_b + l_bd + l_t + l_n;

end