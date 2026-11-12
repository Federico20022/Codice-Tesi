clearvars; clc;
tic
% ROUGH HESTON-KAPPA + CGMY - SPX 2026-05-15.
% Come Rough_Heston_CMGY_def, con kappa libero: separa l'effetto della mean
% reversion da quello del tipo di salti.
%   params = [H, nu, rho, q, M, Y, kappa],  G = 0.05 fisso,  M >= 5
% CF: ramo 'rough Heston kappa CGMY' di charfun.
% Output: Results/calib_rhcgmykappa_7x15z.mat.
XI_MODE = 'pchip';

REDUCED_GRID = true;

fprintf('\nCaricamento dataset Refinitiv SPX 2026-05-15 (merged)...\n');

CSV_PATH   = '20260515SPX_merged.csv';
QUOTE_DATE = datenum('2026-05-15');

raw = readtable(CSV_PATH, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
rng(42);

% Livello dell'indice dal feed (CF_LAST della riga Root = SPOT)
% NON e' sincrono con lo snapshot delle catene: sta ~1.1% SOPRA lo spot
% implicito dai forward (7501.24 contro ~7418, q implicito 0.85%).
is_spot = strcmp(raw.Root, 'SPOT');
S0_feed = raw.CF_LAST(is_spot);
assert(isscalar(S0_feed), ...
    'Attesa una sola riga SPOT nel CSV, trovate %d.', numel(S0_feed));
fprintf('Livello indice da feed (non sincrono): %.2f | XI_MODE = %s\n', S0_feed, XI_MODE);

raw     = raw(~is_spot, :);

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

% STEP 1: FILTRI
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
    T_avail  = unique(T_surf);
    T_target = [0.0384 0.0833 0.25 0.5 0.75 1.0 1.5];  % 7 scadenze
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
fprintf('Opzioni post-filtro: %d (su %d) | Maturities: %d\n\n', ...
    length(F_surf), N_orig, length(T_uniq));

% STEP 2: ESTRAZIONE Kvar(T) - DEMETERFI
fprintf('Estrazione Kvar(T) via variance swap replication...\n');

K_TAIL_LOW  = 0.05;
K_TAIL_HIGH = 2.50;
N_GRID      = 4000;

Kvar_mkt   = zeros(length(T_uniq),1);
IV_ATM_mkt = zeros(length(T_uniq),1);

for i = 1:length(T_uniq)
    T_i = T_uniq(i);
    m   = (T_surf == T_i);
    K_i = K_surf(m); IV_i = IV_Mkt(m); F_i = F_surf(find(m,1));
    [K_s, srt] = sort(K_i); IV_s = IV_i(srt);
    [~, atm_idx] = min(abs(K_s - F_i));
    IV_ATM_mkt(i) = IV_s(atm_idx);
    K_grid = logspace(log10(K_TAIL_LOW*F_i), log10(K_TAIL_HIGH*F_i), N_GRID)';
    IV_grid = interp1(K_s, IV_s, K_grid, 'linear', NaN);
    IV_grid(K_grid < K_s(1))   = IV_s(1);
    IV_grid(K_grid > K_s(end)) = IV_s(end);
    OTM_prices = zeros(size(K_grid));
    for j = 1:length(K_grid)
        if K_grid(j) < F_i
            [~, p] = blkprice(F_i, K_grid(j), 0, T_i, IV_grid(j));
            OTM_prices(j) = p;
        else
            [c, ~] = blkprice(F_i, K_grid(j), 0, T_i, IV_grid(j));
            OTM_prices(j) = c;
        end
    end
    Kvar_mkt(i) = (2/T_i) * trapz(K_grid, OTM_prices ./ K_grid.^2);
end

% STEP 3: COSTRUZIONE xi SECONDO XI_MODE
if strcmp(XI_MODE, 'pchip')
    fprintf('Costruzione curva xi0(t) PCHIP su W(T)=T*Kvar(T)...\n');
    T_anchor = [0; T_uniq(:)];
    W_anchor = [0; T_uniq(:) .* Kvar_mkt(:)];
    W_pp     = pchip(T_anchor, W_anchor);
    xi0_eval = @(tvec) eval_xi0_pchip(tvec, W_pp);
    fprintf('%-10s  %-10s  %-10s  %-10s\n', 'T(yr)', 'Kvar', 'xi0(T)', 'sqrt(xi0)%');
    for i = 1:length(T_uniq)
        xi_i = xi0_eval(T_uniq(i));
        fprintf('%-10.4f  %-10.5f  %-10.5f  %-10.3f\n', ...
            T_uniq(i), Kvar_mkt(i), xi_i, sqrt(max(xi_i,0))*100);
    end
else
    fprintf('xi flat per maturity: xi(t) = Kvar(T_i) - L_cgmy\n');
    xi0_eval = [];
    fprintf('%-10s  %-10s  %-12s  %-10s\n', 'T(yr)', 'Kvar', 'sqrt(Kvar)%', 'IV_ATM%');
    for i = 1:length(T_uniq)
        fprintf('%-10.4f  %-10.5f  %-12.3f  %-10.3f\n', ...
            T_uniq(i), Kvar_mkt(i), sqrt(max(Kvar_mkt(i),0))*100, IV_ATM_mkt(i)*100);
    end
end

% SELEZIONE 15 STRIKE PER SCADENZA
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
    n_kept = sum(keep15(m));
    if n_kept < N_PER_MAT
        warning('T=%.4f: solo %d strike distinti su %d target', T_uniq(i), n_kept, N_PER_MAT);
    end
    fprintf('  T=%.4f: z in [%.2f, %.2f], %d strike tenuti\n', ...
        T_uniq(i), min(z(keep15(m))), max(z(keep15(m))), sum(keep15(m)));
end
F_surf=F_surf(keep15); T_surf=T_surf(keep15); K_surf=K_surf(keep15);
Mid_surf=Mid_surf(keep15); IV_Mkt=IV_Mkt(keep15);

logM_surf = log(K_surf ./ F_surf);
R_surf  = zeros(size(F_surf));
Q_surf  = zeros(size(F_surf));
CP_flag = -1 * ones(size(F_surf));
fprintf('Selezione %d strike/scadenza: %d opzioni in calibrazione\n\n', ...
    N_PER_MAT, length(F_surf));

% STEP 4: PESI
N_valid = length(F_surf);
Vega_Mkt = zeros(N_valid,1);
for i = 1:N_valid
    d1 = (log(F_surf(i)/K_surf(i)) + 0.5*IV_Mkt(i)^2*T_surf(i)) / ...
         (IV_Mkt(i)*sqrt(T_surf(i)));
    Vega_Mkt(i) = F_surf(i) * sqrt(T_surf(i)) * normpdf(d1);
end
Vega_Mkt = max(Vega_Mkt, 1e-4);
IV_Mkt_w = max(IV_Mkt, 0.01);
sqrt_w   = 1 ./ (Vega_Mkt .* IV_Mkt_w);
w_cap    = prctile(sqrt_w, 99);
n_capped = sum(sqrt_w > w_cap);
sqrt_w   = min(sqrt_w, w_cap);
fprintf('Pesi winsorizzati al 99mo percentile: %d punti cappati\n', n_capped);

% STEP 5: G FISSO + LOSS
% G non identificato (stesso profile del CGMY puro).
G_FIX = 0.05;
g_fix = G_FIX;
Q_MAX = 0.98;  % quota massima della varianza forward minima presa dai salti

% xi_min: minimo della curva forward variance sull'unione delle griglie
% Riccati effettivamente usate.
if strcmp(XI_MODE, 'pchip')
    tt_all = [];
    for k = 1:length(T_uniq)
        tt_all = [tt_all, ((0:399)+0.5)/400*T_uniq(k), ((0:99)+0.5)/100*T_uniq(k)];
    end
    xi_min = min(xi0_eval(tt_all));
else
    xi_min = min(Kvar_mkt);
end
fprintf('min xi0 sulle griglie Riccati : %.6f  (sqrt = %.2f%%)\n', ...
    xi_min, 100*sqrt(xi_min));
fprintf('Vincolo di ammissibilita'': Lcg <= %.2f * xi_min = %.6f\n', ...
    Q_MAX, Q_MAX*xi_min);

% Vedi CMGY_def: nessuna griglia SINC coarse sicura sul range di Y.
% Lo stadio coarse degrada i NODI RICCATI, non la griglia SINC.
[XC, POW_FAST] = XcPow('rH+CGMY', 'fast');
[~,  POW_EXC ] = XcPow('rH+CGMY', 'exact');
% Errore di quadratura sull'integrale di xi, misurato in bp di IV contro una
% soluzione a 6400 nodi sul rough Heston calibrato: 0.25 bp a n=100, 0.09 bp
% a n=200, 0.03 bp a n=400.
N_riccati_fast = 100;  % ~0.25 bp: sufficiente per il warm start
N_riccati_exc  = 400;  % ~0.03 bp: stadio finale

% p = [H, nu, rho, q, M, Y]   (G fisso, C ricavato da q)
Residui_Fast = @(p) sqrt_w .* ( ...
    PrezziRoughHestonCGMY(F_surf, T_surf, K_surf, p, g_fix, xi_min, T_uniq, Kvar_mkt, ...
        xi0_eval, XI_MODE, N_riccati_fast, XC, POW_FAST, CP_flag) - Mid_surf);

Residui_Exact = @(p) sqrt_w .* ( ...
    PrezziRoughHestonCGMY(F_surf, T_surf, K_surf, p, g_fix, xi_min, T_uniq, Kvar_mkt, ...
        xi0_eval, XI_MODE, N_riccati_exc, XC, POW_EXC, CP_flag) - Mid_surf);

Loss_Fast_GA = @(p) sum(Residui_Fast(p).^2);

% STEP 6: CALIBRAZIONE 3-STEP (6 parametri: H, nu, rho, C, M, Y)
fprintf('\n=== Calibrazione rough Heston + CGMY (G=%.2f fisso, 6 params, xi %s) ===\n', ...
    G_FIX, XI_MODE);

%        [H,     nu,    rho,   q,     M,     Y   ]
% q in (0, Q_MAX]: vincolo di ammissibilita'.
% M >= 5: restrizione economica sulla scala del salto positivo (vedi
% header).
lb = [0.005, 0.05, -0.99, 1e-4, 5.0,   0.9,  0.0];
ub = [0.50,  5.00,  0.10, Q_MAX, 200.0, 1.95, 50.0];
nvars = 7;  % ultimo parametro: kappa

% Seed espressi in (C, M, Y) e convertiti in q, cosi' restano leggibili.
%            H       nu      rho      C         M      Y      kappa
seed_CMY = [ 0.0878, 0.3676, -0.7182, 5.681e-4, 21.71, 1.832, 0.0  % ottimo con kappa = 0
             0.2000, 0.6000, -0.6900, 5.681e-4, 21.71, 1.832, 1.0
             0.4000, 1.1000, -0.6900, 5.681e-4, 21.71, 1.832, 3.0 ];  % come rough Heston-kappa
pop0 = zeros(size(seed_CMY));
for k = 1:size(seed_CMY,1)
    qk = seed_CMY(k,4) * LambdaCGMY(G_FIX, seed_CMY(k,5), seed_CMY(k,6)) / xi_min;
    pop0(k,:) = [seed_CMY(k,1:3), min(max(qk,lb(4)),ub(4)), ...
                 min(max(seed_CMY(k,5),lb(5)),ub(5)), seed_CMY(k,6), seed_CMY(k,7)];
end

fprintf('\n[Step 1] GA (3 gen, pop 20, 7 vars)...\n');
options_ga = optimoptions('ga', ...
    'MaxGenerations', 3, 'PopulationSize', 20, ...
    'InitialPopulationMatrix', pop0, ...
    'Display', 'iter', 'FunctionTolerance', 1e-9);
t_ga = tic;
theta_GA = ga(Loss_Fast_GA, nvars, [], [], [], [], lb, ub, [], options_ga);
fprintf('Step 1 in %.1f s\n', toc(t_ga));

fprintf('\n[Step 2] lsqnonlin fast (max 100 iter)...\n');
options_lsq = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunctionEvaluations', 5000, 'MaxIterations', 300, 'Display', 'iter', ...
    'StepTolerance', 1e-8, 'FunctionTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-8);
t_lsq1 = tic;
theta_LSQ1 = lsqnonlin(Residui_Fast, theta_GA, lb, ub, options_lsq);
fprintf('Step 2 in %.1f s\n', toc(t_lsq1));

fprintf('\n[Step 3] lsqnonlin exact (max 50 iter)...\n');
options_lsq_exc = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunctionEvaluations', 2000, 'MaxIterations', 100, 'Display', 'iter', ...
    'StepTolerance', 1e-6, 'FunctionTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-8);
t_lsq2 = tic;
theta_ottimi_6 = lsqnonlin(Residui_Exact, theta_LSQ1, lb, ub, options_lsq_exc);
fprintf('Step 3 in %.1f s\n', toc(t_lsq2));

% Ricostruisci il vettore completo [H, nu, rho, C, G, M, Y] per
% compatibilita' con il pricer e con la PricingSuite: q torna a essere C.
q_opt = theta_ottimi_6(4);
C_da_q = q_opt * xi_min / LambdaCGMY(G_FIX, theta_ottimi_6(5), theta_ottimi_6(6));
theta_ottimi = [theta_ottimi_6(1:3), C_da_q, G_FIX, ...
                theta_ottimi_6(5), theta_ottimi_6(6), theta_ottimi_6(7)];

H_opt  = theta_ottimi(1); nu_opt  = theta_ottimi(2); rho_opt = theta_ottimi(3);
C_opt  = theta_ottimi(4); G_opt   = theta_ottimi(5); M_opt   = theta_ottimi(6);
Y_opt  = theta_ottimi(7);
kap_opt = theta_ottimi(8);

fprintf('\n=== PARAMETRI ROUGH HESTON + CGMY (G=%.2f fisso, xi %s) ===\n', G_FIX, XI_MODE);
fprintf('--- Rough Heston ---\n');
fprintf('H (Hurst)        : %.4f\n', H_opt);
fprintf('nu (vol-of-vol)  : %.4f\n', nu_opt);
fprintf('rho              : %.4f\n', rho_opt);
fprintf('kappa            : %.4f   (scala 1/kappa = %.3f anni)\n', kap_opt, 1/max(kap_opt,1e-12));
fprintf('--- CGMY (G fisso, M vincolato >= %.0f) ---\n', lb(5));
fprintf('q (quota xi_min) : %.4f   (tetto %.2f)\n', q_opt, Q_MAX);
fprintf('C                : %.6g   (= q * xi_min / lambda)\n', C_opt);
fprintf('G                : %.5f  (FISSO)\n', G_opt);
fprintf('M                : %.4f\n', M_opt);
fprintf('Y                : %.4f\n', Y_opt);

% Boundary check su 6 parametri liberi
bnd = {};
pn6 = {'H','nu','rho','q','M','Y','kappa'};
% Soglia RELATIVA per i parametri che coprono ordini di grandezza (C, M):
% con lb=1e-4 e ub=10 una soglia di 1e-3*range vale 0.0101 e segnala "@lb"
% su quasi tutto il dominio plausibile, cioe' non discrimina nulla.
log_scale = [false false false true true false false];  % [H nu rho q M Y kappa]
for k = 1:7
    if log_scale(k)
        at_lb = theta_ottimi_6(k) < lb(k)*1.05;
        at_ub = theta_ottimi_6(k) > ub(k)*0.95;
    else
        rng_k = ub(k) - lb(k);
        at_lb = theta_ottimi_6(k) < lb(k) + 1e-3*rng_k;
        at_ub = theta_ottimi_6(k) > ub(k) - 1e-3*rng_k;
    end
    if at_lb, bnd{end+1} = [pn6{k} '@lb']; end
    if at_ub, bnd{end+1} = [pn6{k} '@ub']; end
end
if isempty(bnd)
    fprintf('Boundary check: tutti interni (G fissato a %.3f)\n', G_FIX);
else
    fprintf('Boundary check: %s (G fissato a %.3f)\n', strjoin(bnd,', '), G_FIX);
end

% DIAGNOSTICHE
% Le due quantita' del secondo ordine della componente CGMY hanno ruoli
% diversi: lc entra nella costruzione di xi0 (log-contract), vcg entra nella
% quadratic variation.
[lc_cgmy_1y, var_cgmy_1y] = ContributiCGMY(C_opt, G_opt, M_opt, Y_opt);

skew_cgmy = C_opt * gamma(3-Y_opt) * (M_opt^(Y_opt-3) - G_opt^(Y_opt-3)) / ...
            max(var_cgmy_1y,1e-12)^1.5;

[~, i1y] = min(abs(T_uniq - 1));
T_1y = T_uniq(i1y);
if strcmp(XI_MODE, 'pchip')
    ti_1y  = ((0:N_riccati_exc-1) + 0.5) / N_riccati_exc * T_1y;
    xie_1y = max(xi0_eval(ti_1y) - lc_cgmy_1y, 1e-6);
else
    xie_1y = max(Kvar_mkt(i1y) - lc_cgmy_1y, 1e-6) * ones(1, N_riccati_exc);
end
% Quota della QUADRATIC VARIATION del modello attribuibile ai salti:
% numeratore e denominatore sono ora entrambi varianze.
QV_mod_1y = mean(xie_1y) + var_cgmy_1y;
quota_j   = 100 * var_cgmy_1y / max(QV_mod_1y, 1e-12);
VSL_cgmy  = var_cgmy_1y / max(lc_cgmy_1y, 1e-12);

fprintf('\n--- DIAGNOSTICHE ---\n');
fprintf('CGMY var contrib 1y  : %.5f -> sqrt = %.2f%%\n', var_cgmy_1y, sqrt(max(var_cgmy_1y,0))*100);
fprintf('CGMY log-contract 1y : %.5f  (sottratto a Kvar per costruire xi0)\n', lc_cgmy_1y);
fprintf('VS/L della componente: %.3f  (~1 -> quasi-diffusione, >>1 -> salti veri)\n', VSL_cgmy);
fprintf('CGMY skewness 1y     : %+.4f\n', skew_cgmy);
fprintf('CGMY G/M             : %.4f  (<1 = piu salti negativi)\n', G_opt/M_opt);
fprintf('rH alpha = H+1/2     : %.4f\n', H_opt+0.5);
fprintf('Quota salti su QV    : %.1f%% della QV del modello @T=%.2f\n', quota_j, T_1y);
fprintf('Lcg / xi_min         : %.1f%%  (tetto di ammissibilita'' %.0f%%)\n', ...
    100*lc_cgmy_1y/xi_min, 100*Q_MAX);
fprintf('Scala salto positivo : 1/M = %.4f in log  (%.0f%% in prezzo)\n', ...
    1/M_opt, 100*(exp(1/M_opt)-1));
fprintf('Scala salto negativo : 1/G = %.4f in log  (%.0f%% in prezzo)\n', ...
    1/G_opt, 100*(1-exp(-1/G_opt)));
if VSL_cgmy < 1
    fprintf(2, ['ATTENZIONE: VS/L della componente < 1. La misura di Levy e''\n' ...
                '  dominata dai salti POSITIVI: componente di rialzo, non di\n' ...
                '  crash. Su un indice azionario e'' il segno sbagliato.\n']);
end

% Due regimi degeneri da riconoscere, entrambi descritti dalla Remark 6.2:
%  - C al lower bound: l'attivita' dei salti e' spenta e i restanti
%    parametri CGMY non sono piu' identificati;
%  - Y al upper bound: la componente non e' spenta ma degenera al limite
%    DIFFUSIVO, perche' per Y -> 2 il processo di Levy converge a un moto
%    browniano. In quel regime gamma(2-Y) diverge e C, M contano solo
%    attraverso la varianza della componente, non separatamente.
% La condizione sui salti innocui usa AND, non OR: una quota di varianza
% alta con skew basso NON e' un salto innocuo.
C_al_lb = theta_ottimi_6(4) < lb(4)*1.05;
Y_al_ub = theta_ottimi_6(6) > ub(6)*0.99;
if Y_al_ub
    fprintf(['>>> VERDETTO: componente CGMY spinta al limite DIFFUSIVO (Y=%.3f al bound).\n' ...
             '    Per Y->2 il processo di Levy degenera in un moto browniano: la\n' ...
             '    componente "salti" sta agendo da seconda diffusione, non da salti.\n' ...
             '    Conferma indipendente dai bound: VS/L della componente = %.3f\n' ...
             '    (un salto vero dara'' un valore ben maggiore di 1).\n' ...
             '    Quota (%.1f%%) e skew (%+.2f) vanno letti con cautela: a Y->2 il\n' ...
             '    fattore gamma(2-Y) diverge e C, M non sono piu'' separatamente\n' ...
             '    identificati (contano solo attraverso la varianza della componente).\n'], ...
             theta_ottimi_6(6), VSL_cgmy, quota_j, skew_cgmy);
elseif C_al_lb
    fprintf(['>>> VERDETTO: attivita' ' dei salti SPENTA (C=%.2e al lower bound).\n' ...
             '    I restanti parametri CGMY non sono identificati: non riportarli\n' ...
             '    come stime. Vedi Remark 6.2.\n'], theta_ottimi_6(4));
elseif quota_j < 3 && abs(skew_cgmy) < 0.5
    fprintf('>>> VERDETTO: salti innocui (quasi-diffusione rough)\n');
else
    fprintf('>>> VERDETTO: salti attivi (quota %.1f%% della QV, skew %+.2f, VS/L %.3f)\n', ...
        quota_j, skew_cgmy, VSL_cgmy);
end

% Coerenza del log-contract: con la sottrazione corretta l'identita'
% mean(xi_eff) + L_cgmy = Kvar(T) deve valere a meno del clip.
fprintf('\n--- COERENZA LOG-CONTRACT (modello vs mercato) ---\n');
fprintf('%-9s  %-12s  %-12s  %-12s  %-10s\n', 'T(yr)','Kvar mkt','Kvar mod','diff','vol pt');
for k = 1:length(T_uniq)
    Tk = T_uniq(k);
    if strcmp(XI_MODE, 'pchip')
        tik = ((0:N_riccati_exc-1) + 0.5) / N_riccati_exc * Tk;
        xie = max(xi0_eval(tik) - lc_cgmy_1y, 1e-6);
    else
        xie = max(Kvar_mkt(k) - lc_cgmy_1y, 1e-6) * ones(1, N_riccati_exc);
    end
    Kmod = mean(xie) + lc_cgmy_1y;
    fprintf('%-9.4f  %-12.5f  %-12.5f  %-12.2e  %+10.3f\n', ...
        Tk, Kvar_mkt(k), Kmod, Kmod - Kvar_mkt(k), ...
        100*(sqrt(max(Kmod,0)) - sqrt(max(Kvar_mkt(k),0))));
end

lf = sum(Residui_Fast(theta_ottimi_6).^2);
le = sum(Residui_Exact(theta_ottimi_6).^2);
fprintf('\nSanity Fast/Exact    : %.4e / %.4e  Ratio=%.2f\n', lf, le, lf/le);

% salvataggio
cartella_corrente = fileparts(mfilename('fullpath'));
cartella_results  = fullfile(cartella_corrente, '..', 'Results');
if ~exist(cartella_results, 'dir'), mkdir(cartella_results); end
save(fullfile(cartella_results, 'calib_rhcgmykappa_7x15z.mat'), ...
     'theta_ottimi', 'REDUCED_GRID', 'T_uniq', 'Kvar_mkt', 'XI_MODE', 'G_FIX', ...
     'xi_min', 'Q_MAX', 'q_opt');
fprintf('\nParametri salvati in: %s (REDUCED_GRID=%d)\n', ...
    fullfile(cartella_results, 'calib_rhcgmykappa_7x15z.mat'), REDUCED_GRID);

% STEP 7: PREZZI E IV MODELLO
fprintf('\nRicalcolo prezzi e inversione IV...\n');
PrezziModello = PrezziRoughHestonCGMY(F_surf, T_surf, K_surf, theta_ottimi_6, g_fix, ...
    xi_min, T_uniq, Kvar_mkt, xi0_eval, XI_MODE, N_riccati_exc, XC, POW_EXC, CP_flag);

IV_Mod = InvertiIV_OTM(F_surf, K_surf, T_surf, max(PrezziModello, 1e-6));
valid  = ~isnan(IV_Mod) & ~isnan(IV_Mkt);

% STEP 8: ERRORI
% NB sul segno: eps_rel = (IV_mkt - IV_mod)/IV_mkt, quindi un Bias_rel%
err_IV    = IV_Mod(valid) - IV_Mkt(valid);
err_price = PrezziModello(valid) - Mid_surf(valid);
eps_rel   = -err_IV ./ IV_Mkt(valid);

RMSE_paper   = sqrt(mean(eps_rel.^2)) * 100;
MAE_paper    = mean(abs(eps_rel)) * 100;
MaxErr_paper = max(abs(eps_rel)) * 100;
RMSE_Vol_abs = sqrt(mean(err_IV.^2)) * 100;
MAPE_Prezzi  = mean(abs(err_price ./ Mid_surf(valid))) * 100;

fprintf('\n--- RISULTATI AGGREGATI ---\n');
fprintf('N opzioni valide         : %d / %d\n', sum(valid), length(valid));
fprintf('\n[Paper eq. 13-14 - RMSE relativo IV]\n');
fprintf('RMSE_rel IV (paper)      : %.4f%%\n', RMSE_paper);
fprintf('MAE_rel  IV (paper)      : %.4f%%\n', MAE_paper);
fprintf('Max |eps| IV (paper)     : %.4f%%\n', MaxErr_paper);
fprintf('\n[Metriche assolute]\n');
fprintf('RMSE_abs IV              : %.4f%%\n', RMSE_Vol_abs);
fprintf('MAPE Prezzi fwd          : %.4f%%\n', MAPE_Prezzi);

fprintf('\n--- ERRORE IV PER SCADENZA ---\n');
fprintf('%-9s  %-6s  %-12s  %-12s  %-12s\n', 'T(yr)','N','RMSE_rel%','MAE_rel%','RMSE_abs%');
T_sorted  = sort(unique(T_surf(valid)));
rmse_by_T = zeros(length(T_sorted),1);
for i = 1:length(T_sorted)
    m = (T_surf == T_sorted(i)) & valid;
    e     = IV_Mod(m) - IV_Mkt(m);
    e_rel = -e ./ IV_Mkt(m);
    rmse_by_T(i) = sqrt(mean(e_rel.^2))*100;
    fprintf('%-9.4f  %-6d  %-12.4f  %-12.4f  %-12.4f\n', ...
        T_sorted(i), sum(m), rmse_by_T(i), mean(abs(e_rel))*100, sqrt(mean(e.^2))*100);
end

fprintf('\n--- ERRORE IV PER MONEYNESS BUCKET ---\n');
buckets = [-inf,-0.30,-0.15,-0.05,0.05,0.15,0.30,inf];
labels  = {'<-0.30 (deep OTM P)','-0.30..-0.15','-0.15..-0.05', ...
           '-0.05..0.05 (ATM)','0.05..0.15','0.15..0.30','>0.30 (deep OTM C)'};
fprintf('%-22s  %-6s  %-12s  %-12s  %-12s\n','Bucket','N','RMSE_rel%','MAE_rel%','Bias_rel%');
for i = 1:length(buckets)-1
    m = (logM_surf >= buckets(i)) & (logM_surf < buckets(i+1)) & valid;
    if sum(m) == 0, continue; end
    e_rel = -(IV_Mod(m) - IV_Mkt(m)) ./ IV_Mkt(m);
    fprintf('%-22s  %-6d  %-12.4f  %-12.4f  %+.4f\n', ...
        labels{i}, sum(m), sqrt(mean(e_rel.^2))*100, mean(abs(e_rel))*100, mean(e_rel)*100);
end

fprintf('\n--- ERRORE IV PER MONEYNESS STANDARDIZZATA (z) ---\n');
z_surf = zeros(size(logM_surf));
for i = 1:length(T_uniq)
    m = (T_surf == T_uniq(i));
    z_surf(m) = logM_surf(m) ./ (IV_ATM_mkt(i) * sqrt(T_uniq(i)));
end
z_buckets = [-inf,-2,-1,-0.5,0.5,1,inf];
z_labels  = {'z<-2 (deep OTM P)','-2..-1','-1..-0.5','-0.5..0.5 (ATM)','0.5..1','z>1 (OTM C)'};
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

fprintf('\n--- RMSE IN VOL POINTS PER SMILE ---\n');
fprintf('%-9s  %-5s  %-14s  %-12s  %-12s\n', 'T(yr)','N','RMSE (vol pt)','MAE (vol pt)','Max |err|');
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
figure('Name','Superficie IV 3D - rH-kappa+CGMY','Color','w','Position',[100 100 900 600]);
scatter3(logM_surf(valid), T_surf(valid), IV_Mkt(valid)*100, 25,'r','o','DisplayName','Mercato');
hold on; grid on;
scatter3(logM_surf(valid), T_surf(valid), IV_Mod(valid)*100, 25,[0.4 0.2 0.6],'filled','DisplayName','rH+CGMY');
xlabel('ln(K/F)'); ylabel('T (yr)'); zlabel('IV (%)');
title(sprintf('rH+CGMY (G=%.2f fisso, xi %s): H=%.3f nu=%.3f rho=%.3f | C=%.4f M=%.2f Y=%.3f', ...
    G_FIX, XI_MODE, H_opt, nu_opt, rho_opt, C_opt, M_opt, Y_opt));
legend('Location','best'); view(-45,20);

ScadenzePlt = unique(T_surf(valid));
num_plots = length(ScadenzePlt);
cols = min(4, num_plots); rows = ceil(num_plots/cols);
figure('Name','Smiles IV rH-kappa+CGMY','Color','w','Position',[50 50 1400 800]);
for i = 1:num_plots
    T_tgt = ScadenzePlt(i);
    mask  = (T_surf == T_tgt) & valid;
    [logM_sorted, sidx] = sort(logM_surf(mask));
    IVm = IV_Mkt(mask); IVm = IVm(sidx);
    IVh = IV_Mod(mask); IVh = IVh(sidx);
    subplot(rows, cols, i);
    plot(logM_sorted, IVm*100,'ro','MarkerSize',4,'LineWidth',1,'DisplayName','Mercato'); hold on;
    plot(logM_sorted, IVh*100,'-','Color',[0.4 0.2 0.6],'LineWidth',1.5,'DisplayName','rH+CGMY');
    xline(0,'k--','F','HandleVisibility','off');
    title(sprintf('T=%.4f (N=%d)', T_tgt, sum(mask)),'FontSize',9);
    ylabel('IV (%)'); xlabel('ln(K/F)'); grid on;
    if i==1, legend('Location','best','FontSize',8); end
end

figure('Name','RMSE IV per T - rH-kappa+CGMY','Color','w','Position',[200 200 900 450]);
bar(T_sorted, rmse_by_T,'FaceColor',[0.4 0.2 0.6]); grid on;
xlabel('Maturity (years)'); ylabel('RMSE IV (%)');
title(sprintf('RMSE IV - rH+CGMY (G=%.2f fisso, xi %s) - SPX 2026-05-15', G_FIX, XI_MODE));
yline(RMSE_paper,'r--',sprintf('RMSE rel: %.4f%%', RMSE_paper),'LineWidth',1.5);

toc

% FUNZIONI LOCALI
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

function xi = eval_xi0_pchip(tvec, W_pp)
    h      = 1e-5;
    tvec   = tvec(:);
    tplus  = max(tvec + h, 0);
    tminus = max(tvec - h, 0);
    denom  = tplus - tminus;
    denom(denom < 1e-12) = 1e-12;
    xi = (ppval(W_pp, tplus) - ppval(W_pp, tminus)) ./ denom;
    xi = max(xi, 1e-8);
end

function lam = LambdaCGMY(G, M, Y)
% Lcg = C * lam(G,M,Y): il log-contract e' LINEARE nell'attivita' C.
    lam = 2*( gamma(-Y) * ((M-1)^Y - M^Y + (G+1)^Y - G^Y) ...
            - gamma(1-Y) * (M^(Y-1) - G^(Y-1)) );
end

function [Lcg, Vcg] = ContributiCGMY(C, G, M, Y)
% Le due quantita' del secondo ordine di una componente CGMY.
%   Lcg = 2 * int (e^x - 1 - x) nu(dx)   contributo al LOG-CONTRACT
%   Vcg =     int x^2 nu(dx)             contributo alla QUADRATIC VARIATION
% La curva xi0 si ottiene derivando W(T) = T*Kvar(T), e Kvar(T) e' lo strike
% del variance swap replicato alla Demeterfi, cioe' il log-contract: (2/T)
% E[-log(F_T/F_0)] = (1/T) int_0^T xi0 + Lcg.
    Lcg = C * LambdaCGMY(G, M, Y);
    Vcg = C * gamma(2-Y) * (M^(Y-2) + G^(Y-2));
end

function P_fwd = PrezziRoughHestonCGMY(F_vec, T_vec, K_vec, params6, G_fix, xi_min, ...
                    T_uniq, Kvar_mkt, xi0_fn, xi_mode, N_riccati, Xc, pow, CP_vec)
% params6 = [H, nu, rho, q, M, Y, kappa] q e' la quota della varianza
% forward minima assorbita dal log-contract dei salti; C si ricava da q, e
% il vincolo q <= Q_MAX garantisce per costruzione che xi0 - Lcg resti
% positiva (vedi header).
% Da Kvar si sottrae il LOG-CONTRACT dei salti, non la varianza.
    H   = params6(1); nu  = params6(2); rho = params6(3);
    q   = params6(4); M   = params6(5); Y   = params6(6);
    G   = G_fix;
    C   = q * xi_min / LambdaCGMY(G, M, Y);
    kap = params6(7);
    params_combined = [H, nu, rho, kap, C, G, M, Y];

    L_cgmy = ContributiCGMY(C, G, M, Y);

    P_fwd = zeros(length(K_vec),1);

    for k = 1:length(T_uniq)
        T   = T_uniq(k);
        idx = (T_vec == T);
        if ~any(idx), continue; end
        F_t  = F_vec(find(idx,1));
        K_t  = K_vec(idx);
        CP_t = CP_vec(idx);

        if strcmp(xi_mode, 'pchip')
            ti = ((0:N_riccati-1) + 0.5) / N_riccati * T;
            xi_grid = max(xi0_fn(ti) - L_cgmy, 1e-6);
            xi_grid = xi_grid(:)';
        else
            xi_grid = max(Kvar_mkt(k) - L_cgmy, 1e-6) * ones(1, N_riccati);
        end

        prezzi_t = PrezzaSuperficieSINC_yield(F_t, T*ones(size(K_t)), K_t, ...
            zeros(size(K_t)), zeros(size(K_t)), ...
            'rough Heston kappa CGMY', params_combined, xi_grid, N_riccati, Xc, 0, pow, CP_t(1));

        P_fwd(idx) = prezzi_t;
    end
end
