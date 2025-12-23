function results = engineModel(x, N_spools, data)

    % Unpack Variables
        % Constants
    
        V0 = data.V0;
        gamma = data.gamma;
        F_req = data.F_req;
    
        % Design Vars
    
        pi_c = x(1);
        % ...
        % ...
    

    % Stage Calculations
    
    inlet = inletModel(x, N_spools, data);
    comp = compressorModel(x, N_spools, data, inlet);
    burner = burnerModel(x, N_spools, data, inlet, comp);
    turbine = turbineModel(x, N_spools, data, inlet, comp, burner);
    nozzle = nozzleModel(x, N_spools, data, inlet, comp, burner, turbine);
    

    % Compute outputs    

end