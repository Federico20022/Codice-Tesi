function PrezziCalcolati = PrezzaSuperficieSINC_yield_new(S0_input, T_vec, Strikes_vec, R_vec, Q_vec, model, params, xi, j, Xc, Xm, pow, CP_input)
% PREZZASUPERFICIESINC_YIELD_NEW  Come PrezzaSuperficieSINC_yield, con
% sottostante e flag call/put anche vettoriali (uno per scadenza).
    PrezziCalcolati = zeros(length(Strikes_vec), 1);

    ScadenzeUniche = unique(T_vec);

    for i = 1:length(ScadenzeUniche)
        t_corrente = ScadenzeUniche(i);

        idx = (T_vec == t_corrente);
        K_correnti = Strikes_vec(idx);

        % Sottostante vettoriale (es. F_surf): valore della scadenza
        if length(S0_input) > 1
            S0_corrente = S0_input(find(idx, 1));
        else
            S0_corrente = S0_input;
        end

        % Idem per il flag call/put
        if length(CP_input) > 1
            CP_corrente = CP_input(find(idx, 1));
        else
            CP_corrente = CP_input;
        end

        r_corrente = R_vec(find(idx, 1));
        q_corrente = Q_vec(find(idx, 1));

        prezzi_T = SINC_fastFT(S0_corrente, t_corrente, K_correnti, r_corrente, q_corrente, model, params, xi, j, Xc, Xm, pow, CP_corrente);

        PrezziCalcolati(idx) = prezzi_T;
    end
end
