function f_c = costFun(x, N_spools, data)    
    

    results = engineModel(x, N_spools, data);
    % Returns (at minimum):
        % A_max_eng
        % tot_l
        % eta_0
    
    
    % f_c function!
    
    f_c = sqrt(results.A_max_eng) * sqrt(results.tot_l) / (results.eta_0^2) * (3 + N_spools);

end