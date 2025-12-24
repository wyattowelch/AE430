function f_c = costFun(x, N_spools, data)

    try
        results = engineModel(x, N_spools, data);

        A   = results.A_max_eng;
        l   = results.l_tot;
        eta = results.eta_0;

        if ~isfinite(A) || ~isreal(A) || A <= 0 || ...
           ~isfinite(l) || ~isreal(l) || l <= 0 || ...
           ~isfinite(eta) || ~isreal(eta) || eta <= 0
            f_c = 1e30; 
            return
        end

        eta = max(eta, 1e-6);
        f_c = sqrt(A) * sqrt(l) / (eta^2) * (3 + N_spools);

    catch
        f_c = 1e30;
    end
end
