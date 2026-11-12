clearvars; clc;

tic
% MULTISTART_RHH  Multistart ridotto del rough Hawkes Heston, sulla stessa
% funzione obiettivo di RoughHawkesHeston_SPX_20260515.
% N_START partenze Sobol + l'ottimo attuale, raffinamento exact dei
% N_REFINE migliori. Tempo: 2-4 ore (USE_PARFOR = true se disponibile).
% Output: Results/multistart_rhh.mat, Figure/tab_app_multistart_rhh.tex.
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

% MULTISTART
N_START    = 10;  % partenze Sobol (oltre all'ottimo attuale)
N_REFINE   = 3;  % quante delle migliori raffinare con la loss exact
TOL_BACINO = 0.05;  % distanza max-norm nello spazio normalizzato
USE_PARFOR = false;  % true se e' installato il Parallel Computing Toolbox

%        [alpha,  b,      c,     rho,   Lambda, sig0^2, beta ]  (= calibrazione)
lb     = [0.505, -20.0,   1e-4, -0.99,  0.0,    1e-5,   1e-4];
ub     = [0.995, -0.01,   2.00,  0.10,  5.0,    0.10,   0.30];
% scala logaritmica per i parametri che coprono ordini di grandezza; per b
% (negativo) si usa la scala log di |b|
SCALA  = {'lin','neglog','log','lin','lin','log','log'};
nomi   = {'alpha','b','c','rho','Lambda','sig0sq','beta'};
nvars  = 7;

cartella_results = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
cartella_figure  = fullfile(fileparts(mfilename('fullpath')), '..', 'Figure');
S = load(fullfile(cartella_results,'calib_rhh_7x15z.mat'),'theta_ottimi');
seed = min(max(S.theta_ottimi(:)', lb), ub);

% punti di partenza
P = sobolset(nvars, 'Skip', 1);
P = scramble(P, 'MatousekAffineOwen');
rng(42);
U = net(P, N_START);
starts = zeros(N_START+1, nvars);
starts(1,:) = seed;
for k = 1:N_START
    starts(k+1,:) = DaUnitario(U(k,:), lb, ub, SCALA);
end

opt_fast = optimoptions('lsqnonlin', ...
    'Algorithm','trust-region-reflective', ...
    'MaxFunctionEvaluations',2500,'MaxIterations',100,'Display','off', ...
    'StepTolerance',1e-6,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);
opt_exc = optimoptions('lsqnonlin', ...
    'Algorithm','trust-region-reflective', ...
    'MaxFunctionEvaluations',800,'MaxIterations',15,'Display','off', ...
    'StepTolerance',1e-6,'FunctionTolerance',1e-6,'OptimalityTolerance',1e-8);

% ottimizzazione fast da ogni partenza
nS = size(starts,1);
theta = nan(nS, nvars);  loss = inf(nS,1);  tempo = nan(nS,1);
fprintf('\n=== MULTISTART rough Hawkes Heston: %d partenze ===\n', nS);
if USE_PARFOR
    parfor k = 1:nS
        t0 = tic;
        try
            th = lsqnonlin(Residui_Fast, starts(k,:), lb, ub, opt_fast);
            theta(k,:) = th;  loss(k) = sum(Residui_Fast(th).^2);
        catch
        end
        tempo(k) = toc(t0);
    end
    for k = 1:nS
        fprintf('%3d  loss fast = %.6e  (%.0f s)  %s\n', k, loss(k), tempo(k), ...
            strjoin(compose('%.4g', theta(k,:)), ' '));
    end
else
    for k = 1:nS
        t0 = tic;
        try
            th = lsqnonlin(Residui_Fast, starts(k,:), lb, ub, opt_fast);
            theta(k,:) = th;  loss(k) = sum(Residui_Fast(th).^2);
        catch ME
            fprintf('%3d  FALLITA (%s)\n', k, ME.message);
        end
        tempo(k) = toc(t0);
        eti = ''; if k == 1, eti = '   <- ottimo attuale'; end
        fprintf('%3d  loss fast = %.6e  (%.0f s)  %s%s\n', k, loss(k), tempo(k), ...
            strjoin(compose('%.4g', theta(k,:)), ' '), eti);
    end
end

ok = isfinite(loss);
assert(any(ok), 'Tutte le partenze sono fallite.');
[loss_ord, ord] = sort(loss(ok));
th_ok = theta(ok,:);  theta_ord = th_ok(ord,:);

% bacini
Z = zeros(size(theta_ord));
for k = 1:size(theta_ord,1), Z(k,:) = AUnitario(theta_ord(k,:), lb, ub, SCALA); end
bacino = zeros(size(Z,1),1);  nb = 0;
for k = 1:size(Z,1)
    for bb = 1:nb
        rap = find(bacino == bb, 1);
        if max(abs(Z(k,:) - Z(rap,:))) < TOL_BACINO, bacino(k) = bb; break; end
    end
    if bacino(k) == 0, nb = nb + 1; bacino(k) = nb; end
end

% raffinamento exact e RMSE relativo in IV
% Si raffina il migliore di ciascun bacino, cosi' il confronto fra bacini e'
% fatto con la loss exact e con la stessa metrica riportata in Table 6.2.
[~, primo] = unique(bacino, 'stable');  % primo = migliore del bacino
primo = primo(1:min(N_REFINE, numel(primo)));
ref_theta = nan(numel(primo), nvars); ref_loss = nan(numel(primo),1);
ref_rmse  = nan(numel(primo),1);
fprintf('\n--- RAFFINAMENTO EXACT (migliore di ciascun bacino) ---\n');
for j = 1:numel(primo)
    k = primo(j);
    try
        th = lsqnonlin(Residui_Exact, theta_ord(k,:), lb, ub, opt_exc);
    catch
        th = theta_ord(k,:);
    end
    ref_theta(j,:) = th;  ref_loss(j) = sum(Residui_Exact(th).^2);
    Pm = PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, R_surf, Q_surf, ...
        'rough Hawkes Heston', th, [], NSTEPS_EXC, XC_EXC, 0, POW_EXC, CP_flag);
    IVm = InvertiIV_OTM(F_surf, K_surf, T_surf, max(Pm,1e-6));
    v = ~isnan(IVm) & ~isnan(IV_Mkt);
    ref_rmse(j) = 100*sqrt(mean(((IV_Mkt(v)-IVm(v))./IV_Mkt(v)).^2));
    fprintf('bacino %d: loss exact %.6e | RMSE rel %.3f%% | alpha %.4f (H_eff %.3f) | Lambda %.4f\n', ...
        bacino(k), ref_loss(j), ref_rmse(j), th(1), th(1)-0.5, th(5));
end

% riepilogo
fprintf('\n=== RIEPILOGO ===\n');
fprintf('%-7s %-9s %-14s %-10s\n','bacino','partenze','loss fast min','x globale');
for bb = 1:nb
    mb = (bacino == bb);
    fprintf('%-7d %-9d %-14.6e %-10.3f\n', bb, sum(mb), min(loss_ord(mb)), ...
        min(loss_ord(mb))/loss_ord(1));
end
seed_loss = sum(Residui_Exact(seed).^2);
[best_loss, ib] = min(ref_loss);
fprintf('\nOttimo attuale (tesi): loss exact %.6e\n', seed_loss);
fprintf('Migliore trovato     : loss exact %.6e (RMSE rel %.3f%%)\n', best_loss, ref_rmse(ib));
if best_loss < 0.999*seed_loss
    fprintf(2, ['\nATTENZIONE: trovato un ottimo MIGLIORE di quello in tesi.\n' ...
                'Non sovrascrivere nulla: mandami l''output prima di aggiornare le tabelle.\n']);
else
    fprintf('\nNessun bacino migliore dell''ottimo attuale: l''ottimo della tesi e'' confermato.\n');
end
if isfinite(loss(1))
    i_seed = find(all(abs(theta_ord - theta(1,:)) < 1e-12, 2), 1);
    n_stesso = sum(bacino == bacino(i_seed));
    fprintf('Partenze finite nel bacino dell''ottimo attuale: %d su %d\n', n_stesso, sum(ok));
end

% salvataggio + tabella LaTeX
save(fullfile(cartella_results,'multistart_rhh.mat'), 'starts','theta','loss','tempo', ...
    'theta_ord','loss_ord','bacino','nb','ref_theta','ref_loss','ref_rmse', ...
    'seed','seed_loss','N_START','TOL_BACINO','lb','ub');
fid = fopen(fullfile(cartella_figure,'tab_app_multistart_rhh.tex'),'w');
fprintf(fid,'%% Generato da Multistart_RHH.m il %s -- non modificare a mano\n', datestr(now));
fprintf(fid,'\\begin{tabular}{rrrrr}\n\\toprule\n');
fprintf(fid,'Basin & Restarts & Relative RMSE (\\%%) & $\\alpha$ & $\\Lambda$\\\\\n\\midrule\n');
for j = 1:numel(primo)
    fprintf(fid,'%d & %d & %.2f & %.3f & %.3f\\\\\n', bacino(primo(j)), ...
        sum(bacino == bacino(primo(j))), ref_rmse(j), ref_theta(j,1), ref_theta(j,5));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);
fprintf('\nSalvati multistart_rhh.mat e tab_app_multistart_rhh.tex\n');
toc

% FUNZIONI LOCALI (le prime due copiate dallo script di calibrazione)
function x = DaUnitario(u, lb, ub, scala)
    x = zeros(size(u));
    for j = 1:numel(u)
        switch scala{j}
            case 'log',    x(j) = exp(log(lb(j)) + u(j)*(log(ub(j))-log(lb(j))));
            case 'neglog', x(j) = -exp(log(-ub(j)) + u(j)*(log(-lb(j))-log(-ub(j))));
            otherwise,     x(j) = lb(j) + u(j)*(ub(j)-lb(j));
        end
    end
end

function u = AUnitario(x, lb, ub, scala)
    u = zeros(size(x));
    for j = 1:numel(x)
        switch scala{j}
            case 'log',    u(j) = (log(x(j))-log(lb(j))) / (log(ub(j))-log(lb(j)));
            case 'neglog', u(j) = (log(-x(j))-log(-ub(j))) / (log(-lb(j))-log(-ub(j)));
            otherwise,     u(j) = (x(j)-lb(j)) / (ub(j)-lb(j));
        end
    end
end

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
