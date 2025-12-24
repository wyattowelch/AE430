function [c, ceq] = constraints(x, N_spools, data)
    
    results = engineModel(x, N_spools, data);
    % Returns (at minimum):
       % Tt4
       % Twg
       % X_CO
       % sigma_c
       % tau_res
       % T_cool
       % M_tip_comp
       % M_tip_turb
       % stalls
       % reactionVals
    

    c = [];
    ceq = []; 

    
    % Inequality Constraints
    
    c(end+1) = results.Tt4 - 1900; % K
    c(end+1) = results.Twg - 1300; % K
    c(end+1) = results.X_CO - 2E-4;
    c(end+1) = results.sigma_c - 70E6; % Pa
    c(end+1) = 5e-3 - results.tau_res; % s
    c(end+1) = results.T_cool - 1000; % K
    c(end+1) = results.M_tip_comp - 1.2;
    c(end+1) = results.M_tip_turb - 1.2;

    % rh, rm, rt cannot ctall
    for k = 1:numel(results.stalls)
        c(end+1) = results.stalls(k);
    end

    % degree of reaction between 1 and 0 
    for k = 1:numel(results.reactionVals)
        R = results.reactionVals(k);
        c(end+1) = -R;     % R >= 0
        c(end+1) = R - 1;  % R <= 1
    end

    % Equality Constraints:

    ceq(end+1) = results.F_in - data.F_req;


end