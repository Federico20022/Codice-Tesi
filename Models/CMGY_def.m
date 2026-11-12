clearvars; clc;
tic
% CGMY PURO - SPX 2026-05-15. Calibra C, M, Y con G = 0.05 fisso.
% Qui il profilo in G non e' piatto (Profile_G): G = 0.05 e' una
% convenzione dentro il tratto indistinguibile, non una stima.
REDUCED_GRID = true;  % true = 7 scadenze x 15 strike | false = full

fprintf('\nCaricamento dataset Refinitiv SPX 2026-05-15 (merged)...\n');

CSV_PATH   = '20260515SPX_merged.csv';
QUOTE_DATE = datenum('2026-05-15');

raw = readtable(CSV_PATH, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
rng(42);

% Livello dell'indice dal feed (CF_LAST della riga Root = SPOT)
% Non sincrono con lo snapshot delle catene: sta ~1.1% SOPRA lo spot
% implicito dai forward.
is_spot = strcmp(raw.Root, 'SPOT');
S0_feed = raw.CF_LAST(is_spot);
assert(isscalar(S0_feed), ...
    'Attesa una sola riga SPOT nel CSV, trovate %d.', numel(S0_feed));
fprintf('Livello indice da feed (non sincrono): %.2f\n', S0_feed);
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
    T_avail   = unique(T_surf);
    T_targets = [0.0384 0.0833 0.25 0.5 0.75 1.0 1.5];  % 7 scadenze
    T_sel = [];
    for tt = T_targets
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

% IV ATM + SELEZIONE 15 STRIKE PER SCADENZA
IV_ATM_mkt = zeros(length(T_uniq),1);
for i = 1:length(T_uniq)
    m = (T_surf == T_uniq(i));
    K_i = K_surf(m); IV_i = IV_Mkt(m); F_i = F_surf(find(m,1));
    [~, atm_idx] = min(abs(K_i - F_i));
    IV_ATM_mkt(i) = IV_i(atm_idx);
end

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

% STEP 2: PESI
N_valid = length(F_surf);
Vega_Mkt = zeros(N_valid, 1);
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

% STEP 3: G FISSO + LOSS
% G fissato per convenzione, non stimato.
G_FIX = 0.05;
g_fix = G_FIX;  % copia catturata dalle lambda

% Per il CGMY NON esiste una coppia (Xc,pow) coarse sicura sull'intero range
% di Y che il GA esplora: con Y~1.5 e G=0.05 la coda sinistra decade come
% e^{-0.05|x|} e serve Xc grande; con Y~1.95 a T corta la densita' e' quasi
% gaussiana e strettissima e serve dx = 2*Xc/N <= 0.008. I due vincoli
% tirano in direzioni opposte.
[XC,       POW_FAST] = XcPow('CGMY', 'fast');
[~,        POW_EXC ] = XcPow('CGMY', 'exact');
% Conseguenza da tenere presente: con XC, POW_FAST e POW_EXC identici i due
% residui coincidono, quindi lo Step 3 riparte dall'ottimo dello Step 2 ed
% e' un no-op, e un confronto "fast vs exact" darebbe 1.00 per costruzione.
% Il controllo numerico che serve davvero qui e' un altro: il troncamento
% del dominio, che per il CGMY e' l'errore dominante perche' C ~ 0.017 e la
% coda sinistra decade come e^{-G|x|} con G = 0.05. Si confronta quindi la
% loss contro una configurazione con Xc molto piu' largo.
[XC_GUARD, POW_GUARD] = XcPow('CGMY', 'guard');

Residui_Fast = @(p) sqrt_w .* ( ...
    PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, R_surf, Q_surf, ...
        'CGMY', [p(1), g_fix, p(2), p(3)], [], [], XC, 0, POW_FAST, CP_flag) - Mid_surf );

Residui_Exact = @(p) sqrt_w .* ( ...
    PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, R_surf, Q_surf, ...
        'CGMY', [p(1), g_fix, p(2), p(3)], [], [], XC, 0, POW_EXC, CP_flag) - Mid_surf );

Residui_Guard = @(p) sqrt_w .* ( ...
    PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, R_surf, Q_surf, ...
        'CGMY', [p(1), g_fix, p(2), p(3)], [], [], XC_GUARD, 0, POW_GUARD, CP_flag) - Mid_surf );

Loss_Fast_GA = @(p) sum( Residui_Fast(p).^2 );

% STEP 4: CALIBRAZIONE 3-STEP (C, M, Y)
fprintf('\n=== Calibrazione CGMY (G=%.2f fisso, 3 params: C, M, Y) ===\n', G_FIX);
% M >= 5: stessa restrizione economica del modello composito (salto positivo
% di scala massima 20%).
lb = [1e-4, 5.0, 0.05];
ub = [50.0, 2000, 1.95];
nvars = 3;
guess_iniziale = [0.017, 33.0, 1.51];  % warm start dal profile G=0.05

fprintf('\n[Step 1] GA su SINC fast (10 gen, pop 60)...\n');
options_ga = optimoptions('ga', ...
    'MaxGenerations', 10, 'PopulationSize', 60, ...
    'InitialPopulationMatrix', guess_iniziale, ...
    'Display', 'iter', 'FunctionTolerance', 1e-9);
t_ga = tic;
theta_GA = ga(Loss_Fast_GA, nvars, [], [], [], [], lb, ub, [], options_ga);
fprintf('Step 1 in %.2f s\n', toc(t_ga));

fprintf('\n[Step 2] lsqnonlin su SINC fast...\n');
options_lsq = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunctionEvaluations', 1500, 'Display', 'iter', ...
    'StepTolerance', 1e-8, 'FunctionTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-8);
t_lsq1 = tic;
theta_LSQ1 = lsqnonlin(Residui_Fast, theta_GA, lb, ub, options_lsq);
fprintf('Step 2 in %.2f s\n', toc(t_lsq1));

fprintf('\n[Step 3] lsqnonlin su SINC exact (stessa griglia: rifinitura no-op)...\n');
t_lsq2 = tic;
theta_ottimi_CMY = lsqnonlin(Residui_Exact, theta_LSQ1, lb, ub, options_lsq);
fprintf('Step 3 in %.2f s\n', toc(t_lsq2));

% Ricostruisci vettore completo [C, G, M, Y] per compatibilita' col pricer
theta_ottimi = [theta_ottimi_CMY(1), G_FIX, theta_ottimi_CMY(2), theta_ottimi_CMY(3)];

C_opt = theta_ottimi(1);
G_opt = theta_ottimi(2);  % = G_FIX
M_opt = theta_ottimi(3);
Y_opt = theta_ottimi(4);

fprintf('\n=== PARAMETRI CGMY CALIBRATI ===\n');
fprintf('C  : %.5f\n', C_opt);
fprintf('G  : %.5f  (FISSO per convenzione, non stimato)\n', G_opt);
fprintf('     NB: nel CGMY puro G NON e'' non-identificato. Profile_G.m mostra\n');
fprintf('     che la loss cresce di un fattore 8.2 fra G=0.01 e G=50, con\n');
fprintf('     minimo sul bordo inferiore della griglia: il modello a salti\n');
fprintf('     puri, da solo, chiede una coda sinistra non temperata. La\n');
fprintf('     de-identificazione vale solo nel modello composito.\n');
fprintf('M  : %.5f\n', M_opt);
fprintf('Y  : %.5f\n', Y_opt);

% Boundary check su C, M, Y (G e' fisso)
bnd = {};
pn  = {'C','M','Y'};
th3 = theta_ottimi_CMY;
for k = 1:3
    if th3(k) < lb(k)*(1+1e-3)+1e-6, bnd{end+1} = [pn{k} '@lb']; end
    if th3(k) > ub(k)*(1-1e-3),      bnd{end+1} = [pn{k} '@ub']; end
end
if isempty(bnd)
    fprintf('Boundary check: tutti interni (G fissato a %.3f)\n', G_FIX);
else
    fprintf('Boundary check: %s (G fissato a %.3f)\n', strjoin(bnd,', '), G_FIX);
end

% DIAGNOSTICHE
var_1y  = C_opt * gamma(2-Y_opt) * (M_opt^(Y_opt-2) + G_opt^(Y_opt-2));
skew_1y = C_opt * gamma(3-Y_opt) * (M_opt^(Y_opt-3) - G_opt^(Y_opt-3)) / ...
          max(var_1y,1e-12)^1.5;
fprintf('\n--- DIAGNOSTICHE CGMY ---\n');
fprintf('Var(X_1)       : %.5f  -> sqrt = %.2f%%\n', var_1y, sqrt(max(var_1y,0))*100);
fprintf('Skewness(X_1)  : %+.4f  [solo segno: ~G^(Y-3), G e'' una convenzione]\n', skew_1y);
fprintf('G/M (asymmetry): %.4f  (<1 = piu salti negativi, atteso per SPX)\n', G_opt/M_opt);
fprintf('Regime         : ');
if Y_opt < 1
    fprintf('attivita infinita, variazione finita (0 <= Y < 1)\n');
else
    fprintf('variazione infinita (1 <= Y < 2)\n');
end
% Il ramo "attivita finita" (Y < 0) non e' raggiungibile con lb(3)=0.05.

lg = sum(Residui_Guard(theta_ottimi_CMY).^2);
le = sum(Residui_Exact(theta_ottimi_CMY).^2);
fprintf(['Troncamento    : loss(Xc=%d) %.4e / loss(Xc=%d) %.4e  Ratio=%.4f\n' ...
         '                 (>1 = il dominio a Xc=%d sta tagliando la coda sinistra)\n'], ...
         XC, le, XC_GUARD, lg, le/lg, XC);

% salvataggio
cartella_corrente = fileparts(mfilename('fullpath'));
cartella_results  = fullfile(cartella_corrente, '..', 'Results');
if ~exist(cartella_results, 'dir'), mkdir(cartella_results); end
save(fullfile(cartella_results, 'calib_cgmy_7x15z.mat'), ...
     'theta_ottimi', 'REDUCED_GRID', 'G_FIX');
fprintf('\nParametri salvati in: %s (REDUCED_GRID=%d)\n', ...
    fullfile(cartella_results, 'calib_cgmy_7x15z.mat'), REDUCED_GRID);

% STEP 5: PREZZI E IV MODELLO
fprintf('\nRicalcolo prezzi e inversione IV...\n');
PrezziModello = PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, R_surf, Q_surf, ...
    'CGMY', theta_ottimi, [], [], XC, 0, POW_EXC, CP_flag);

IV_Mod = InvertiIV_OTM(F_surf, K_surf, T_surf, max(PrezziModello, 1e-6));
valid  = ~isnan(IV_Mod) & ~isnan(IV_Mkt);

% STEP 6: ERRORI AGGREGATI
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
fprintf('%-9s  %-6s  %-12s  %-12s  %-12s\n', ...
    'T(yr)','N','RMSE_rel%','MAE_rel%','RMSE_abs%');
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
figure('Name','Superficie IV 3D - CGMY','Color','w','Position',[100 100 900 600]);
scatter3(logM_surf(valid), T_surf(valid), IV_Mkt(valid)*100, 25,'r','o','DisplayName','Mercato');
hold on; grid on;
scatter3(logM_surf(valid), T_surf(valid), IV_Mod(valid)*100, 25,'m','filled','DisplayName','CGMY');
xlabel('ln(K/F)'); ylabel('T (yr)'); zlabel('IV (%)');
title(sprintf('CGMY (G=%.2f fisso): C=%.4f M=%.2f Y=%.3f - SPX 2026-05-15', ...
    G_FIX, C_opt, M_opt, Y_opt));
legend('Location','best'); view(-45,20);

ScadenzePlt = unique(T_surf(valid));
num_plots = length(ScadenzePlt);
cols = min(4, num_plots); rows = ceil(num_plots/cols);
figure('Name','Smiles IV - CGMY','Color','w','Position',[50 50 1400 800]);
for i = 1:num_plots
    T_plot = ScadenzePlt(i);
    mask = (T_surf == T_plot) & valid;
    [logM_sorted, sidx] = sort(logM_surf(mask));
    IVm = IV_Mkt(mask); IVm = IVm(sidx);
    IVh = IV_Mod(mask); IVh = IVh(sidx);
    subplot(rows, cols, i);
    plot(logM_sorted, IVm*100,'ro','MarkerSize',4,'LineWidth',1,'DisplayName','Mercato'); hold on;
    plot(logM_sorted, IVh*100,'m-','LineWidth',1.5,'DisplayName','CGMY');
    xline(0,'k--','F','HandleVisibility','off');
    title(sprintf('T=%.4f (N=%d)', T_plot, sum(mask)),'FontSize',9);
    ylabel('IV (%)'); xlabel('ln(K/F)'); grid on;
    if i==1, legend('Location','best','FontSize',8); end
end

figure('Name','RMSE IV per scadenza - CGMY','Color','w','Position',[200 200 900 450]);
bar(T_sorted, rmse_by_T,'FaceColor',[0.7 0.3 0.5]); grid on;
xlabel('Maturity (years)'); ylabel('RMSE IV (%)');
title(sprintf('RMSE IV - CGMY (G=%.2f fisso) - SPX 2026-05-15', G_FIX));
yline(RMSE_paper,'r--',sprintf('RMSE rel: %.4f%%', RMSE_paper),'LineWidth',1.5);

toc

% FUNZIONE LOCALE
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
