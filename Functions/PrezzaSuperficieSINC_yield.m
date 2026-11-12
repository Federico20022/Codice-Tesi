% PREZZASUPERFICIESINC_YIELD  Prezza con SINC_fastFT tutti gli strike,
% scadenza per scadenza, con tassi R_vec e dividendi Q_vec.
function PrezziCalcolati = PrezzaSuperficieSINC_yield(S0, T_vec, Strikes_vec, R_vec, Q_vec, model, params, xi, j, Xc, Xm, pow, CP)
    PrezziCalcolati = zeros(length(Strikes_vec), 1);
    ScadenzeUniche = unique(T_vec);

    for i = 1:length(ScadenzeUniche)
        t_corrente = ScadenzeUniche(i);
        idx = (T_vec == t_corrente);
        K_correnti = Strikes_vec(idx);

        r_corrente = R_vec(find(idx, 1));
        q_corrente = Q_vec(find(idx, 1));  % Estrai il dividendo corretto!

        prezzi_T = SINC_fastFT(S0, t_corrente, K_correnti, r_corrente, q_corrente, model, params, xi, j, Xc, Xm, pow, CP);
        PrezziCalcolati(idx) = prezzi_T;
    end
end
