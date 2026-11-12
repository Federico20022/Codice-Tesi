clearvars; clc;

tic
% ROUGH HAWKES HESTON - SPX 2026-05-15 (Bondi, Pulido, Scotti 2022).
%   params = [alpha, b, c, rho, Lambda, sigma0_sq, beta]
% La curva iniziale (sigma0_sq, beta) si calibra: con i salti il legame con
% il variance swap passa per il risolvente.
% Stesso campione, filtri e pesi degli altri modelli. CF: CharFunRHH.
REDUCED_GRID = true;  % true = ~12 scadenze (primo run) | false = tesi

fprintf('\nCaricamento dataset Refinitiv SPX 2026-05-15 (merged)...\n');

CSV_PATH   = '20260515SPX_merged.csv';
QUOTE_DATE = datenum('2026-05-15');

raw = readtable(CSV_PATH, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
rng(42);

% Riga spot (Root = SPOT)
is_spot = strcmp(raw.Root, 'SPOT');
S0_feed = raw.CF_LAST(is_spot);
assert(isscalar(S0_feed), ...
    'Attesa una sola riga SPOT nel CSV, trovate %d.', numel(S0_feed));
fprintf('Livello indice da feed (non sincrono): %.2f | REDUCED_GRID = %d\n', S0_feed, REDUCED_GRID);
raw     = raw(~is_spot, :);

% Flag call/put dalla colonna pre-calcolata
is_call = strcmp(raw.PutCall, 'C');
keep    = ~isnan(raw.STRIKE_PRC) & ...
          raw.CF_BID > 0 & raw.CF_ASK > raw.CF_BID & ...
          raw.CF_BID >= 0.375;
raw     = raw(keep, :);
is_call = is_call(keep);

K_all   = raw.STRIKE_PRC;
mid_all = 0.5*(raw.CF_BID + raw.CF_ASK);
T_all   = (datenum(raw.EXPIR_DATE) - QUOTE_DATE) / 365;

[~, ia] = unique([T_all, K_all, double(is_call)], 'rows', 'stable');
T_all = T_all(ia); K_all = K_all(ia);
mid_all = mid_all(ia); is_call = is_call(ia);

% F e D per scadenza via regressione (C-P) = -D*K + D*F
T_uniq_raw = unique(T_all);
F_per_T = nan(length(T_uniq_raw), 1);
D_per_T = nan(length(T_uniq_raw), 1);

for i = 1:length(T_uniq_raw)
    Ti = T_uniq_raw(i);
    mC = (T_all == Ti) & is_call;
    mP = (T_all == Ti) & ~is_call;
    [Kcommon, iC, iP] = intersect(K_all(mC), K_all(mP));
    if length(Kcommon) < 3, continue; end
    cIdx = find(mC); pIdx = find(mP);
    diff_CP = mid_all(cIdx(iC)) - mid_all(pIdx(iP));

    near = abs(Kcommon/S0_feed - 1) <= 0.25;
    if sum(near) >= 5
        Kreg = Kcommon(near); dreg = diff_CP(near);
    else
        Kreg = Kcommon; dreg = diff_CP;
    end

    coef = [ones(size(Kreg)), Kreg] \ dreg;
    D_per_T(i) = -coef(2);
    F_per_T(i) = coef(1) / D_per_T(i);
end

valid_T = ~isnan(F_per_T) & D_per_T > 0 & D_per_T < 1.05 & F_per_T > 0;
T_uniq_raw = T_uniq_raw(valid_T);
F_per_T    = F_per_T(valid_T);
D_per_T    = D_per_T(valid_T);
fprintf('Scadenze con forward estratto: %d / %d\n', sum(valid_T), length(valid_T));

% Put-side in forward space (convenzione OTM)
F_surf = []; T_surf = []; K_surf = []; Mid_surf = [];
for i = 1:length(T_uniq_raw)
    Ti = T_uniq_raw(i); Fi = F_per_T(i); Di = D_per_T(i);
    mTi  = (T_all == Ti);
    Ki   = K_all(mTi); midi = mid_all(mTi); calli = is_call(mTi);

    P_fwd_i = nan(size(Ki));
    for j = 1:length(Ki)
        if Ki(j) < Fi && ~calli(j)
            P_fwd_i(j) = midi(j) / Di;
        elseif Ki(j) >= Fi && calli(j)
            P_fwd_i(j) = midi(j)/Di - Fi + Ki(j);
        end
    end
    okj = ~isnan(P_fwd_i) & P_fwd_i > 1e-4;
    F_surf   = [F_surf;   Fi*ones(sum(okj),1)];
    T_surf   = [T_surf;   Ti*ones(sum(okj),1)];
    K_surf   = [K_surf;   Ki(okj)];
    Mid_surf = [Mid_surf; P_fwd_i(okj)];
end

IV_Mkt = InvertiIV_OTM(F_surf, K_surf, T_surf, Mid_surf);
N_orig = length(F_surf);

% STEP 1: FILTRI STANDARD + (opzionale) GRIGLIA RIDOTTA
LOGM_CUT  = 0.5;
MIN_K_MAT = 5;
T_MIN     = 0.03;
T_MAX     = 3.0;

logM = log(K_surf ./ F_surf);
idx  = ~isnan(IV_Mkt) & IV_Mkt > 0.01 & IV_Mkt < 2 & ...
       Mid_surf > 1e-4 & abs(logM) < LOGM_CUT & ...
       T_surf >= T_MIN & T_surf <= T_MAX;
F_surf = F_surf(idx); T_surf = T_surf(idx); K_surf = K_surf(idx);
Mid_surf = Mid_surf(idx); IV_Mkt = IV_Mkt(idx);

T_uniq_all = unique(T_surf);
mask_all = false(length(T_surf),1);
for i = 1:length(T_uniq_all)
    m = T_surf == T_uniq_all(i);
    if sum(m) >= MIN_K_MAT, mask_all = mask_all | m; end
end
F_surf=F_surf(mask_all); T_surf=T_surf(mask_all); K_surf=K_surf(mask_all);
Mid_surf=Mid_surf(mask_all); IV_Mkt=IV_Mkt(mask_all);

if REDUCED_GRID
        % ~12 scadenze target log-spaziate; per ognuna la slice disponibile
    T_avail = unique(T_surf);
    T_target = [0.0384 0.0833 0.25 0.5 0.75 1.0 1.5];  % 7 scadenze target
    % per ognuna si prende la slice disponibile piu' vicina, con dedup.
    T_sel = [];
    for tt = T_target
        [~, ii] = min(abs(T_avail - tt));
        T_sel = [T_sel; T_avail(ii)];
    end
    T_sel = unique(T_sel);
    keepm = ismember(T_surf, T_sel);
    F_surf=F_surf(keepm); T_surf=T_surf(keepm); K_surf=K_surf(keepm);
    Mid_surf=Mid_surf(keepm); IV_Mkt=IV_Mkt(keepm);
    fprintf('Griglia ridotta: %d scadenze selezionate\n', length(T_sel));
end

T_uniq = unique(T_surf);
fprintf('Opzioni disponibili (tutti gli strike): %d | Maturities: %d\n', ...
    length(F_surf), length(T_uniq));

% IV ATM DI MERCATO PER SCADENZA (serve alla selezione per moneyness std)
IV_ATM_mkt = zeros(length(T_uniq),1);
for i = 1:length(T_uniq)
    m = (T_surf == T_uniq(i));
    K_i = K_surf(m); IV_i = IV_Mkt(m); F_i = F_surf(find(m,1));
    [~, atm_idx] = min(abs(K_i - F_i));
    IV_ATM_mkt(i) = IV_i(atm_idx);
end

% SELEZIONE 15 STRIKE PER SCADENZA (moneyness standardizzata)
% z = ln(K/F) / (IV_ATM*sqrt(T)).
N_PER_MAT = 15;
Z_LO = -3.0;  Z_HI = 2.0;

keep15 = false(size(T_surf));
for i = 1:length(T_uniq)
    m = find(T_surf == T_uniq(i));
    z = log(K_surf(m) ./ F_surf(m)) ./ (IV_ATM_mkt(i) * sqrt(T_uniq(i)));
    if numel(m) <= N_PER_MAT
        keep15(m) = true; continue
    end
    z_targets = linspace(max(Z_LO, min(z)), min(Z_HI, max(z)), N_PER_MAT);
    for tt = z_targets
        [~, jj] = min(abs(z - tt));
        keep15(m(jj)) = true;
    end
    fprintf('  T=%.4f: z in [%.2f, %.2f], %d strike tenuti\n', ...
        T_uniq(i), min(z(keep15(m))), max(z(keep15(m))), sum(keep15(m)));
end
F_surf=F_surf(keep15); T_surf=T_surf(keep15); K_surf=K_surf(keep15);
Mid_surf=Mid_surf(keep15); IV_Mkt=IV_Mkt(keep15);

logM_surf = log(K_surf ./ F_surf);
R_surf = zeros(size(F_surf));
Q_surf = zeros(size(F_surf));
CP_flag = -1 * ones(size(F_surf));

fprintf('Selezione %d strike/scadenza: %d opzioni in calibrazione\n\n', ...
    N_PER_MAT, length(F_surf));

% STEP 2: PESI - RELATIVE IV (paper eq. 12-13, winsorizzati)
N_valid = length(F_surf);
Vega_Mkt = zeros(N_valid,1);
for i = 1:N_valid
    d1 = (log(F_surf(i)/K_surf(i)) + 0.5*IV_Mkt(i)^2*T_surf(i)) / (IV_Mkt(i)*sqrt(T_surf(i)));
    Vega_Mkt(i) = F_surf(i) * sqrt(T_surf(i)) * normpdf(d1);
end
Vega_Mkt = max(Vega_Mkt, 1e-4);
IV_Mkt_w = max(IV_Mkt, 0.01);
sqrt_w   = 1 ./ (Vega_Mkt .* IV_Mkt_w);
w_cap    = prctile(sqrt_w, 99);
n_capped = sum(sqrt_w > w_cap);
sqrt_w   = min(sqrt_w, w_cap);
fprintf('Pesi winsorizzati al 99mo percentile: %d punti cappati\n', n_capped);

% STEP 3: LOSS (rough Hawkes Heston, 7 params)
% Firma:
% PrezzaSuperficieSINC_yield_new(F,T,K,R,Q,model,params,xi,j,Xc,Xm,pow,CP)
% Per 'rough Hawkes Heston': xi=[] (non usata), j = NSTEPS Eulero.
NSTEPS_FAST = 250;
NSTEPS_EXC  = 400;
[XC_FAST, POW_FAST] = XcPow('rHawkesH', 'fast');
[XC_EXC,  POW_EXC ] = XcPow('rHawkesH', 'exact');

% Guardia anti-NaN: se un individuo GA produce prezzi non finiti, residuo
% grande invece di crash della spline in SINC_fastFT.
Residui_Fast = @(p) ResiduiRHH_safe(p, F_surf, T_surf, K_surf, R_surf, Q_surf, ...
    NSTEPS_FAST, XC_FAST, POW_FAST, CP_flag, Mid_surf, sqrt_w);

Residui_Exact = @(p) ResiduiRHH_safe(p, F_surf, T_surf, K_surf, R_surf, Q_surf, ...
    NSTEPS_EXC, XC_EXC, POW_EXC, CP_flag, Mid_surf, sqrt_w);

Loss_Fast_GA = @(p) sum( Residui_Fast(p).^2 );

% STEP 4: CALIBRAZIONE 3-STEP (7 parametri, GA con warm start)
fprintf('\n=== Calibrazione rough Hawkes Heston (7 params) ===\n');

%        [alpha,  b,      c,     rho,   Lambda, sig0^2, beta ]
lb     = [0.505, -20.0,   1e-4, -0.99,  0.0,    1e-5,   1e-4];
ub     = [0.995, -0.01,   2.00,  0.10,  5.0,    0.10,   0.30];
nvars  = 7;

% Seeds:
%  1) Tabella 1 del paper (SPX/VIX 2017): punto di partenza informato
%  2) Tabella 1 con Lambda~0: test nested (salti scollegati dai returns)
%  3) esplorativo interno
% La superficie obiettivo del rHH e' fortemente multimodale.
pop0 = [ 0.8949, -3.4964, 1.4030, -0.6280, 0.6725, 0.026294, 0.12397;  % ottimo noto
         0.506,  -2.008,  0.156,  -0.737,  0.242,  0.007,    0.048;
         0.506,  -2.008,  0.156,  -0.737,  0.001,  0.007,    0.048;
         0.600,  -1.000,  0.300,  -0.650,  0.500,  0.020,    0.060 ];

fprintf('\n[Step 1] GA (8 gen, pop 40, popolazione seedata)...\n');
% Pre-flight: valutazione NUDA del pricer al seed Tabella 1. Se c'e' un
% errore di path/porting emerge QUI con il messaggio vero, invece di essere
% assorbito dalla guardia e produrre una loss piatta a 1e3^2*N.
fprintf('Pre-flight pricer (senza guardia)...\n');
Pr0 = PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, R_surf, Q_surf, ...
    'rough Hawkes Heston', pop0(1,:), [], NSTEPS_FAST, XC_FAST, 0, POW_FAST, CP_flag);
assert(all(isfinite(Pr0)), 'Pre-flight: prezzi non finiti al seed Tabella 1.');
fprintf('Pre-flight OK: loss fast al seed Tab.1 = %.4e\n', ...
    sum((sqrt_w.*(Pr0 - Mid_surf)).^2));

options_ga = optimoptions('ga', ...
    'MaxGenerations', 8, 'PopulationSize', 40, ...
    'InitialPopulationMatrix', pop0, ...
    'Display', 'iter', 'FunctionTolerance', 1e-9);
t_ga = tic;
theta_GA = ga(Loss_Fast_GA, nvars, [], [], [], [], lb, ub, [], options_ga);
fprintf('Step 1 in %.1f s\n', toc(t_ga));

fprintf('\n[Step 2] lsqnonlin fast (max 100 iter)...\n');
options_lsq = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunctionEvaluations', 2500, 'MaxIterations', 100, 'Display', 'iter', ...
    'StepTolerance', 1e-6, 'FunctionTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-8);
t_lsq1 = tic;
theta_LSQ1 = lsqnonlin(Residui_Fast, theta_GA, lb, ub, options_lsq);
fprintf('Step 2 in %.1f s\n', toc(t_lsq1));

fprintf('\n[Step 3] lsqnonlin exact (max 50 iter)...\n');
options_lsq_exc = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunctionEvaluations', 800, 'MaxIterations', 15, 'Display', 'iter', ...
    'StepTolerance', 1e-6, 'FunctionTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-8);
t_lsq2 = tic;
theta_ottimi = lsqnonlin(Residui_Exact, theta_LSQ1, lb, ub, options_lsq_exc);
fprintf('Step 3 in %.1f s\n', toc(t_lsq2));

al_o  = theta_ottimi(1); b_o   = theta_ottimi(2); c_o  = theta_ottimi(3);
rho_o = theta_ottimi(4); Lam_o = theta_ottimi(5);
s0_o  = theta_ottimi(6); be_o  = theta_ottimi(7);

fprintf('\n=== PARAMETRI ROUGH HAWKES HESTON ===\n');
fprintf('alpha (kernel)     : %.4f   (H_eff = alpha-1/2 = %.4f)\n', al_o, al_o-0.5);
fprintf('b (mean reversion) : %.4f\n', b_o);
fprintf('c (vol-of-vol^2)   : %.4f   (sqrt(c) = %.4f)\n', c_o, sqrt(c_o));
fprintf('rho                : %.4f\n', rho_o);
fprintf('Lambda (jump lev.) : %.4f\n', Lam_o);
fprintf('sigma0^2           : %.5f  (vol spot = %.2f%%)\n', s0_o, sqrt(s0_o)*100);
fprintf('beta (curva g0)    : %.5f\n', be_o);

% Diagnostiche
% Remark 3 paper: QV dei returns = (1 + Lambda^2 * int z^2 nu(dz)) * int
% sigma^2 con nu esponenziale rate 1: int z^2 nu = 2 -> quota salti =
% 2L^2/(1+2L^2)
qv_jump_share = 2*Lam_o^2 / (1 + 2*Lam_o^2);
fprintf('\n--- DIAGNOSTICHE ---\n');
fprintf('Quota QV da salti  : %.1f%%   (2*Lambda^2/(1+2*Lambda^2))\n', 100*qv_jump_share);
fprintf('g0(1y)             : %.5f  (sqrt = %.2f%%)\n', ...
    s0_o + be_o/gamma(al_o+1), sqrt(s0_o + be_o/gamma(al_o+1))*100);

% Boundary check
bnd = {};
pn = {'alpha','b','c','rho','Lambda','sig0','beta'};
for k = 1:nvars
    rng_k = ub(k) - lb(k);
    if theta_ottimi(k) < lb(k) + 1e-3*rng_k, bnd{end+1} = [pn{k} '@lb']; end
    if theta_ottimi(k) > ub(k) - 1e-3*rng_k, bnd{end+1} = [pn{k} '@ub']; end
end
if isempty(bnd)
    fprintf('Boundary check     : tutti i parametri interni\n');
else
    fprintf('Boundary check     : %s\n', strjoin(bnd, ', '));
end
if Lam_o < 0.01
    fprintf('>>> Lambda ~ 0: salti scollegati dai returns (nested verso rH-like)\n');
else
    fprintf('>>> Lambda = %.3f attivo: salti comuni prezzati dal mercato\n', Lam_o);
end

lf = sum(Residui_Fast(theta_ottimi).^2);
le = sum(Residui_Exact(theta_ottimi).^2);
fprintf('Sanity Fast/Exact  : %.4e / %.4e  Ratio=%.2f\n', lf, le, lf/le);

% salvataggio automatico per la PricingSuite salvataggio automatico robusto
cartella_corrente = fileparts(mfilename('fullpath'));
cartella_results = fullfile(cartella_corrente, '..', 'Results');
if ~exist(cartella_results, 'dir'), mkdir(cartella_results); end
save(fullfile(cartella_results, 'calib_rhh_7x15z.mat'), 'theta_ottimi', 'REDUCED_GRID');
fprintf('\nParametri salvati in: %s\n', fullfile(cartella_results, 'calib_rhh_7x15z.mat'));
fprintf('Parametri salvati in calib_rhh_7x15z.mat (REDUCED_GRID=%d)\n', REDUCED_GRID);

% STEP 5: PREZZI E IV MODELLO
fprintf('\nRicalcolo prezzi e inversione IV...\n');
PrezziModello = PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, R_surf, Q_surf, ...
    'rough Hawkes Heston', theta_ottimi, [], NSTEPS_EXC, XC_EXC, 0, POW_EXC, CP_flag);

IV_Mod = InvertiIV_OTM(F_surf, K_surf, T_surf, max(PrezziModello, 1e-6));
valid = ~isnan(IV_Mod) & ~isnan(IV_Mkt);

% STEP 6: ERRORI AGGREGATI + BREAKDOWN (paper eq. 12-14)
err_IV    = IV_Mod(valid) - IV_Mkt(valid);
err_price = PrezziModello(valid) - Mid_surf(valid);

eps_rel   = -err_IV ./ IV_Mkt(valid);

RMSE_paper   = sqrt(mean(eps_rel.^2)) * 100;
MAE_paper    = mean(abs(eps_rel)) * 100;
MaxErr_paper = max(abs(eps_rel)) * 100;

RMSE_Vol_abs = sqrt(mean(err_IV.^2)) * 100;
MAPE_Prezzi  = mean(abs(err_price ./ Mid_surf(valid))) * 100;

if REDUCED_GRID
    fprintf('\n--- RISULTATI (GRIGLIA RIDOTTA: non confrontare 1:1 col full surface) ---\n');
else
    fprintf('\n--- RISULTATI AGGREGATI ---\n');
end
fprintf('N opzioni valide         : %d / %d\n', sum(valid), length(valid));
fprintf('\n[Paper eq. 13-14 - RMSE relativo IV]\n');
fprintf('RMSE_rel IV (paper)      : %.4f%%\n', RMSE_paper);
fprintf('MAE_rel  IV (paper)      : %.4f%%\n', MAE_paper);
fprintf('Max |eps| IV (paper)     : %.4f%%\n', MaxErr_paper);
fprintf('\n[Metriche assolute - confronto]\n');
fprintf('RMSE_abs IV              : %.4f%%\n', RMSE_Vol_abs);
fprintf('MAPE Prezzi fwd          : %.4f%%\n', MAPE_Prezzi);

fprintf('\n--- ERRORE IV PER SCADENZA (eq. 14 paper) ---\n');
fprintf('%-9s  %-6s  %-12s  %-12s  %-12s\n', ...
    'T(yr)','N','RMSE_rel%','MAE_rel%','RMSE_abs%');
T_sorted = sort(unique(T_surf(valid)));
rmse_by_T = zeros(length(T_sorted),1);
for i = 1:length(T_sorted)
    m = (T_surf == T_sorted(i)) & valid;
    e     = IV_Mod(m) - IV_Mkt(m);
    e_rel = -e ./ IV_Mkt(m);
    rmse_by_T(i) = sqrt(mean(e_rel.^2))*100;
    rmse_abs_T   = sqrt(mean(e.^2))*100;
    fprintf('%-9.4f  %-6d  %-12.4f  %-12.4f  %-12.4f\n', ...
        T_sorted(i), sum(m), rmse_by_T(i), mean(abs(e_rel))*100, rmse_abs_T);
end

fprintf('\n--- ERRORE IV PER MONEYNESS BUCKET ---\n');
buckets = [-inf,-0.30,-0.15,-0.05,0.05,0.15,0.30,inf];
labels = {'<-0.30 (deep OTM P)','-0.30..-0.15','-0.15..-0.05', ...
          '-0.05..0.05 (ATM)','0.05..0.15','0.15..0.30','>0.30 (deep OTM C)'};
fprintf('%-22s  %-6s  %-12s  %-12s  %-12s\n', 'Bucket','N','RMSE_rel%','MAE_rel%','Bias_rel%');
for i = 1:length(buckets)-1
    m = (logM_surf >= buckets(i)) & (logM_surf < buckets(i+1)) & valid;
    if sum(m) == 0, continue; end
    e_rel = -(IV_Mod(m) - IV_Mkt(m)) ./ IV_Mkt(m);
    fprintf('%-22s  %-6d  %-12.4f  %-12.4f  %+.4f\n', ...
        labels{i}, sum(m), sqrt(mean(e_rel.^2))*100, mean(abs(e_rel))*100, mean(e_rel)*100);
end

% ERRORE PER MONEYNESS STANDARDIZZATA  z = ln(K/F)/(IV_ATM*sqrt(T))
z_surf = zeros(size(logM_surf));
for i = 1:length(T_uniq)
    m = (T_surf == T_uniq(i));
    z_surf(m) = logM_surf(m) ./ (IV_ATM_mkt(i) * sqrt(T_uniq(i)));
end

fprintf('\n--- ERRORE IV PER MONEYNESS STANDARDIZZATA (z) ---\n');
z_buckets = [-inf, -2, -1, -0.5, 0.5, 1, inf];
z_labels  = {'z < -2 (deep OTM P)','-2..-1','-1..-0.5','-0.5..0.5 (ATM)', ...
             '0.5..1','z > 1 (OTM C)'};
fprintf('%-22s  %-6s  %-12s  %-12s  %-12s  %-12s\n', ...
    'Bucket z','N','RMSE_rel%','MAE_rel%','Bias_rel%','RMSE vol pt');
for i = 1:length(z_buckets)-1
    m = (z_surf >= z_buckets(i)) & (z_surf < z_buckets(i+1)) & valid;
    if sum(m) == 0, continue; end
    e     = IV_Mod(m) - IV_Mkt(m);
    e_rel = -e ./ IV_Mkt(m);
    fprintf('%-22s  %-6d  %-12.4f  %-12.4f  %+-12.4f  %-12.3f\n', ...
        z_labels{i}, sum(m), sqrt(mean(e_rel.^2))*100, mean(abs(e_rel))*100, ...
        mean(e_rel)*100, sqrt(mean(e.^2))*100);
end

fprintf('\n--- composizione bucket z: N punti per (scadenza x bucket) ---\n');
fprintf('%-9s', 'T(yr)');
for i = 1:length(z_labels), fprintf('%12s', z_labels{i}(1:min(11,end))); end
fprintf('\n');
for k = 1:length(T_uniq)
    fprintf('%-9.4f', T_uniq(k));
    for i = 1:length(z_buckets)-1
        m = (T_surf == T_uniq(k)) & (z_surf >= z_buckets(i)) & ...
            (z_surf < z_buckets(i+1)) & valid;
        fprintf('%12d', sum(m));
    end
    fprintf('\n');
end

% RMSE IN VOL POINTS PER SMILE
fprintf('\n--- RMSE IN VOL POINTS PER SMILE ---\n');
fprintf('%-9s  %-5s  %-14s  %-12s  %-12s\n', ...
    'T(yr)','N','RMSE (vol pt)','MAE (vol pt)','Max |err|');
for i = 1:length(T_sorted)
    m = (T_surf == T_sorted(i)) & valid;
    e = (IV_Mod(m) - IV_Mkt(m))*100;
    fprintf('%-9.4f  %-5d  %-14.3f  %-12.3f  %-12.3f\n', ...
        T_sorted(i), sum(m), sqrt(mean(e.^2)), mean(abs(e)), max(abs(e)));
end
e_all = (IV_Mod(valid) - IV_Mkt(valid))*100;
fprintf('%-9s  %-5d  %-14.3f  %-12.3f  %-12.3f\n', ...
    'TUTTE', sum(valid), sqrt(mean(e_all.^2)), mean(abs(e_all)), max(abs(e_all)));

% PLOT
figure('Name','Superficie IV 3D - rHawkes Heston','Color','w','Position',[100 100 900 600]);
scatter3(logM_surf(valid), T_surf(valid), IV_Mkt(valid)*100, 25, 'r', 'o', 'DisplayName','Mercato');
hold on; grid on;
scatter3(logM_surf(valid), T_surf(valid), IV_Mod(valid)*100, 25, [0.85 0.5 0.1], 'filled', 'DisplayName','rHawkes Heston');
xlabel('ln(K/F)'); ylabel('T (yr)'); zlabel('IV (%)');
title(sprintf('rHawkes Heston: \\alpha=%.3f b=%.2f c=%.3f \\rho=%.3f \\Lambda=%.3f', ...
    al_o, b_o, c_o, rho_o, Lam_o));
legend('Location','best'); view(-45,20);

ScadenzePlt = unique(T_surf(valid));
num_plots = length(ScadenzePlt);
cols = min(5, num_plots); rows = ceil(num_plots/cols);
figure('Name','Smiles IV rHawkes Heston','Color','w','Position',[50 50 1800 1000]);
for i = 1:num_plots
    T_target = ScadenzePlt(i);
    mask = (T_surf == T_target) & valid;
    [logM_sorted, sidx] = sort(logM_surf(mask));
    IVm = IV_Mkt(mask); IVm = IVm(sidx);
    IVh = IV_Mod(mask); IVh = IVh(sidx);
    subplot(rows, cols, i);
    plot(logM_sorted, IVm*100, 'ro','MarkerSize',4,'LineWidth',1, 'DisplayName','Mercato'); hold on;
    plot(logM_sorted, IVh*100, '-','Color',[0.85 0.5 0.1],'LineWidth',1.5, 'DisplayName','rHH');
    xline(0, 'k--','F','HandleVisibility','off');
    title(sprintf('T=%.4f (N=%d)', T_target, sum(mask)), 'FontSize', 9);
    ylabel('IV (%)'); xlabel('ln(K/F)'); grid on;
    if i==1, legend('Location','best','FontSize',8); end
end

figure('Name','RMSE IV per T - rHawkes Heston','Color','w','Position',[200 200 900 450]);
bar(T_sorted, rmse_by_T, 'FaceColor', [0.85 0.5 0.1]); grid on;
xlabel('Maturity (years)'); ylabel('RMSE IV (%)');
title('RMSE IV per scadenza - rough Hawkes Heston - SPX 2026-05-15');
yline(RMSE_paper, 'r--', sprintf('RMSE rel (paper): %.4f%%', RMSE_paper), 'LineWidth',1.5);

toc

% FUNZIONE LOCALE: inversione IV lato OTM
function IV = InvertiIV_OTM(F_vec, K_vec, T_vec, P_fwd_vec)
    IV = nan(size(F_vec));
    for i = 1:length(F_vec)
        try
            if K_vec(i) <= F_vec(i)
                IV(i) = blkimpv(F_vec(i), K_vec(i), 0, T_vec(i), ...
                                P_fwd_vec(i), 'Class', 'put');
            else
                C_fwd = max(P_fwd_vec(i) + F_vec(i) - K_vec(i), 1e-6);
                IV(i) = blkimpv(F_vec(i), K_vec(i), 0, T_vec(i), ...
                                C_fwd, 'Class', 'call');
            end
        catch
            IV(i) = NaN;
        end
    end
end

% FUNZIONE LOCALE: residui con guardia anti-NaN (per GA su bounds larghi)
function res = ResiduiRHH_safe(p, F, T, K, R, Q, nst, xc, pw, CP, Mid, sw)
    try
        Pr = PrezzaSuperficieSINC_yield_new(F, T, K, R, Q, ...
            'rough Hawkes Heston', p, [], nst, xc, 0, pw, CP);
        res = sw .* (Pr - Mid);
        res(~isfinite(res)) = 1e3;
    catch
        res = 1e3 * ones(size(Mid));
    end
end
