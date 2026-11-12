clearvars; clc;
tic
% HESTON - SPX 2026-05-15 (Refinitiv, SPX + SPXW). Calibra i 5 parametri.
% REDUCED_GRID = true  -> le 7 scadenze comuni a tutti i modelli
% REDUCED_GRID = false -> superficie intera (non mischiare i due risultati)
REDUCED_GRID = true;

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
fprintf('Livello indice da feed (non sincrono): %.2f\n', S0_feed);
raw     = raw(~is_spot, :);

% Flag call/put dalla colonna pre-calcolata
is_call = strcmp(raw.PutCall, 'C');
keep    = ~isnan(raw.STRIKE_PRC) & ...
          raw.CF_BID > 0 & raw.CF_ASK > raw.CF_BID & ...
          raw.CF_BID >= 0.375;  % paper: escludi bid < 3/8 (liquidita')
raw     = raw(keep, :);
is_call = is_call(keep);

K_all   = raw.STRIKE_PRC;
mid_all = 0.5*(raw.CF_BID + raw.CF_ASK);
T_all   = (datenum(raw.EXPIR_DATE) - QUOTE_DATE) / 365;

% Dedupe di sicurezza (T,K,flag)
[~, ia] = unique([T_all, K_all, double(is_call)], 'rows', 'stable');
T_all = T_all(ia); K_all = K_all(ia);
mid_all = mid_all(ia); is_call = is_call(ia);

% Estrazione F e D per scadenza via regressione (C-P) = -D*K + D*F
% La regressione usa solo coppie con strike vicino allo spot (|K/S0_feed -
% 1| <= 25%), dove entrambe le gambe sono liquide.
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
        Kreg = Kcommon; dreg = diff_CP;  % fallback: tutte le coppie
    end

    % Regressione: diff_CP = a + b*K, con a = D*F, b = -D
    coef = [ones(size(Kreg)), Kreg] \ dreg;
    D_per_T(i) = -coef(2);
    F_per_T(i) = coef(1) / D_per_T(i);
end

valid_T = ~isnan(F_per_T) & D_per_T > 0 & D_per_T < 1.05 & F_per_T > 0;
T_uniq_raw = T_uniq_raw(valid_T);
F_per_T    = F_per_T(valid_T);
D_per_T    = D_per_T(valid_T);
fprintf('Scadenze con forward estratto: %d / %d\n', sum(valid_T), length(valid_T));

% Build put-side in forward space (convenzione OTM)
% Per K<F uso put; per K>=F converto call: P = C - D*(F-K) -> P_fwd = P/D
F_surf = []; T_surf = []; K_surf = []; Mid_surf = [];
for i = 1:length(T_uniq_raw)
    Ti = T_uniq_raw(i); Fi = F_per_T(i); Di = D_per_T(i);
    mTi = (T_all == Ti);
    Ki  = K_all(mTi); midi = mid_all(mTi); calli = is_call(mTi);

    P_fwd_i = nan(size(Ki));
    for j = 1:length(Ki)
        if Ki(j) < Fi && ~calli(j)  % put OTM
            P_fwd_i(j) = midi(j) / Di;
        elseif Ki(j) >= Fi && calli(j)  % call OTM -> put via parity
            P_fwd_i(j) = midi(j)/Di - Fi + Ki(j);
        end
    end
    okj = ~isnan(P_fwd_i) & P_fwd_i > 1e-4;
    F_surf   = [F_surf;   Fi*ones(sum(okj),1)];
    T_surf   = [T_surf;   Ti*ones(sum(okj),1)];
    K_surf   = [K_surf;   Ki(okj)];
    Mid_surf = [Mid_surf; P_fwd_i(okj)];
end

% IV di mercato da mid
% Inversione sempre sul lato OTM (put per K<=F, call via parity per
IV_Mkt = InvertiIV_OTM(F_surf, K_surf, T_surf, Mid_surf);

N_orig = length(F_surf);

% STEP 1: FILTRI
LOGM_CUT  = 0.5;
MIN_K_MAT = 5;
T_MIN     = 0.03;  % escludo T<11gg (oscillazioni SINC su CF piccata)
T_MAX     = 3.0;  % cap long-end

logM = log(K_surf ./ F_surf);
idx  = ~isnan(IV_Mkt) & IV_Mkt > 0.01 & IV_Mkt < 2 & ...
       Mid_surf > 1e-4 & abs(logM) < LOGM_CUT & ...
       T_surf >= T_MIN & T_surf <= T_MAX;
F_surf = F_surf(idx); T_surf = T_surf(idx); K_surf = K_surf(idx);
Mid_surf = Mid_surf(idx); IV_Mkt = IV_Mkt(idx);

T_uniq = unique(T_surf);
mask_all = false(length(T_surf),1);
for i = 1:length(T_uniq)
    m = T_surf == T_uniq(i);
    if sum(m) >= MIN_K_MAT, mask_all = mask_all | m; end
end
F_surf = F_surf(mask_all); T_surf = T_surf(mask_all); K_surf = K_surf(mask_all);
Mid_surf = Mid_surf(mask_all); IV_Mkt = IV_Mkt(mask_all);

if REDUCED_GRID
    % Stessa selezione del run rough Hawkes Heston (12 target log-spaziati)
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
fprintf('Opzioni disponibili (tutti gli strike): %d | Maturities: %d\n', ...
    length(F_surf), length(T_uniq));

% IV ATM DI MERCATO PER SCADENZA (per la selezione a moneyness std)
IV_ATM_mkt = zeros(length(T_uniq),1);
for i = 1:length(T_uniq)
    m = (T_surf == T_uniq(i));
    K_i = K_surf(m); IV_i = IV_Mkt(m); F_i = F_surf(find(m,1));
    [~, atm_idx] = min(abs(K_i - F_i));
    IV_ATM_mkt(i) = IV_i(atm_idx);
end

% SELEZIONE 15 STRIKE PER SCADENZA (moneyness standardizzata)
% z = ln(K/F)/(IV_ATM*sqrt(T)): in log-moneyness grezza le scadenze corte
% arriverebbero a -8/-10 deviazioni standard.
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
R_surf  = zeros(size(F_surf));
Q_surf  = zeros(size(F_surf));
CP_flag = -1 * ones(size(F_surf));
fprintf('Selezione %d strike/scadenza: %d opzioni in calibrazione\n\n', ...
    N_PER_MAT, length(F_surf));

fprintf('Opzioni in calibrazione: %d (su %d originali)\n', length(F_surf), N_orig);
fprintf('Maturities: %d | Range T: [%.4f, %.3f] yr\n', ...
    length(unique(T_surf)), min(T_surf), max(T_surf));
fprintf('Range |logM|: [%.3f, %.3f]\n\n', min(abs(logM_surf)), max(abs(logM_surf)));

% STEP 2: PESI - RELATIVE IV (Baschetti et al. 2021, eq. 12-13)
% peso = 1 / (vega * IV_mkt) ~ errore relativo in IV via prezzi.
N_valid = length(F_surf);
Vega_Mkt = zeros(N_valid, 1);
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

% STEP 3: LOSS FUNCTIONS (Mrazek-Pospisil 2017, Sez. 3.5.7)
XC_FAST = 6;   POW_FAST = 10;
XC_EXC  = 12;  POW_EXC  = 12;

Residui_Fast  = @(p) sqrt_w .* ( ...
    PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, R_surf, Q_surf, ...
        'Heston', p, [], [], XC_FAST, 0, POW_FAST, CP_flag) - Mid_surf );

Residui_Exact = @(p) sqrt_w .* ( ...
    PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, R_surf, Q_surf, ...
        'Heston', p, [], [], XC_EXC, 0, POW_EXC, CP_flag) - Mid_surf );

Loss_Fast_GA = @(p) sum( Residui_Fast(p).^2 );

% STEP 4: CALIBRAZIONE 3-STEP
fprintf('=== Calibrazione 3-Step (Mrazek-Pospisil 2017) ===\n');

lb = [0.01, 0.005, 0.01, 0.005, -0.99];
ub = [15.0, 0.30,  3.00, 0.20,   0.10];
nvars = 5;
guess_iniziale = [1.5, 0.025, 0.5, 0.028, -0.7];

fprintf('\n[Step 1] GA su SINC coarse (10 gen, pop 50)...\n');
options_ga = optimoptions('ga', ...
    'MaxGenerations', 10, 'PopulationSize', 50, ...
    'InitialPopulationMatrix', guess_iniziale, ...
    'Display', 'iter', 'FunctionTolerance', 1e-9);
t_ga = tic;
theta_GA = ga(Loss_Fast_GA, nvars, [], [], [], [], lb, ub, [], options_ga);
fprintf('Step 1 in %.2f s\n', toc(t_ga));

fprintf('\n[Step 2] lsqnonlin su SINC coarse...\n');options_lsq = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunctionEvaluations', 3000, 'Display', 'iter', ...
    'StepTolerance', 1e-9, 'FunctionTolerance', 1e-9, ...
    'OptimalityTolerance', 1e-9);
t_lsq1 = tic;
theta_LSQ1 = lsqnonlin(Residui_Fast, theta_GA, lb, ub, options_lsq);
fprintf('Step 2 in %.2f s\n', toc(t_lsq1));

fprintf('\n[Step 3] lsqnonlin su SINC esatto...\n');
t_lsq2 = tic;
theta_ottimi = lsqnonlin(Residui_Exact, theta_LSQ1, lb, ub, options_lsq);
fprintf('Step 3 in %.2f s\n', toc(t_lsq2));

fprintf('\n=== PARAMETRI CALIBRATI ===\n');
fprintf('kappa  : %.5f\ntheta  : %.5f\nsigma  : %.5f\nv0     : %.5f\nrho    : %.5f\n', ...
    theta_ottimi(1), theta_ottimi(2), theta_ottimi(3), theta_ottimi(4), theta_ottimi(5));
feller = 2*theta_ottimi(1)*theta_ottimi(2) - theta_ottimi(3)^2;
if feller > 0, fs = 'SODDISFATTO'; else, fs = 'VIOLATO'; end
fprintf('Feller 2*k*theta - sigma^2 = %.4f  (%s)\n', feller, fs);

lf = sum(Residui_Fast(theta_ottimi).^2);
le = sum(Residui_Exact(theta_ottimi).^2);
fprintf('\nSanity Fast/Exact: Loss_Fast=%.4e  Loss_Exact=%.4e  Ratio=%.1f\n', lf, le, lf/le);
% salvataggio automatico per la PricingSuite (cartella ../Results)
cartella_corrente = fileparts(mfilename('fullpath'));
cartella_results  = fullfile(cartella_corrente, '..', 'Results');
if ~exist(cartella_results, 'dir'), mkdir(cartella_results); end
save(fullfile(cartella_results, 'calib_heston_7x15z.mat'), 'theta_ottimi', 'REDUCED_GRID');
fprintf('\nParametri salvati in: %s (REDUCED_GRID=%d)\n', ...
    fullfile(cartella_results, 'calib_heston_7x15z.mat'), REDUCED_GRID);

% STEP 5: PREZZI E IV MODELLO
fprintf('\nRicalcolo prezzi e inversione IV...\n');
PrezziModello = PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, R_surf, Q_surf, ...
    'Heston', theta_ottimi, [], [], 12, 0, 12, CP_flag);

IV_Mod = InvertiIV_OTM(F_surf, K_surf, T_surf, max(PrezziModello, 1e-6));
valid = ~isnan(IV_Mod) & ~isnan(IV_Mkt);

% STEP 6: ERRORI AGGREGATI - formula paper (eq. 12-14)
err_IV    = IV_Mod(valid) - IV_Mkt(valid);
err_price = PrezziModello(valid) - Mid_surf(valid);

eps_rel   = -err_IV ./ IV_Mkt(valid);  % paper eq. 12

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
fprintf('\n[Metriche assolute - confronto]\n');
fprintf('RMSE_abs IV              : %.4f%%\n', RMSE_Vol_abs);
fprintf('MAPE Prezzi fwd          : %.4f%%\n', MAPE_Prezzi);

% STEP 7: DIAGNOSTIC BREAKDOWN (eq. 14 - per maturity slice)
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
buckets = [-inf, -0.30, -0.15, -0.05, 0.05, 0.15, 0.30, inf];
labels  = {'<-0.30 (deep OTM P)', '-0.30..-0.15', '-0.15..-0.05', ...
           '-0.05..0.05 (ATM)', '0.05..0.15', '0.15..0.30', '>0.30 (deep OTM C)'};
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
figure('Name', 'Superficie IV 3D', 'Color', 'w', 'Position', [100, 100, 900, 600]);
scatter3(logM_surf(valid), T_surf(valid), IV_Mkt(valid)*100, 25, 'r', 'o', 'DisplayName', 'Mercato');
hold on; grid on;
scatter3(logM_surf(valid), T_surf(valid), IV_Mod(valid)*100, 25, 'b', 'filled', 'DisplayName', 'Heston');
xlabel('Log-Moneyness ln(K/F)'); ylabel('Maturity (Years)'); zlabel('IV (%)');
title('Superficie IV - SPX 2026-05-15'); legend('Location', 'best'); view(-45, 20);

ScadenzePlt = unique(T_surf(valid));
num_plots = length(ScadenzePlt);
cols = min(5, num_plots);
rows = ceil(num_plots / cols);
figure('Name', 'Smiles IV', 'Color', 'w', 'Position', [50, 50, 1800, 1000]);
for i = 1:num_plots
    T_target = ScadenzePlt(i);
    mask = (T_surf == T_target) & valid;
    [logM_sorted, sidx] = sort(logM_surf(mask));
    IVm = IV_Mkt(mask); IVm = IVm(sidx);
    IVh = IV_Mod(mask); IVh = IVh(sidx);
    subplot(rows, cols, i);
    plot(logM_sorted, IVm*100, 'ro', 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'Mercato'); hold on;
    plot(logM_sorted, IVh*100, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Heston');
    xline(0, 'k--', 'F', 'HandleVisibility', 'off');
    title(sprintf('T=%.4f yr (N=%d)', T_target, sum(mask)), 'FontSize', 9);
    ylabel('IV (%)'); xlabel('Log-Moneyness'); grid on;
    if i == 1, legend('Location', 'best', 'FontSize', 8); end
end

figure('Name', 'RMSE IV per scadenza', 'Color', 'w', 'Position', [200, 200, 900, 450]);
bar(T_sorted, rmse_by_T, 'FaceColor', [0.2 0.4 0.7]); grid on;
xlabel('Maturity (years)'); ylabel('RMSE IV (%)');
title('RMSE IV per scadenza - Heston (SPX 2026-05-15)');
yline(RMSE_paper, 'r--', sprintf('RMSE rel (paper): %.4f%%', RMSE_paper), 'LineWidth', 1.5);

toc

% FUNZIONE LOCALE: inversione IV lato OTM (coerente mercato/modello)
function IV = InvertiIV_OTM(F_vec, K_vec, T_vec, P_fwd_vec)
    % Input: prezzi PUT in forward space.
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

% FUNZIONE LOCALE: Alos-de Santiago-Vives (2015)
function P_fwd = PrezzaSuperficieAlos(F_vec, T_vec, K_vec, params, CP_vec)
    kappa = params(1); theta = params(2); sigma = params(3);
    v0 = params(4);    rho = params(5);
    N = length(K_vec);
    P_fwd = zeros(N,1);
    inv_sqrt2pi = 1/sqrt(2*pi);
    for i = 1:N
        F = F_vec(i); K = K_vec(i); tau = T_vec(i); CP = CP_vec(i);
        D  = exp(-kappa*tau);
        kt = kappa*tau;
        if kt > 1e-8
            vtilde_sq = theta + (v0 - theta)/kt * (1 - D);
        else
            vtilde_sq = v0 + (theta - v0)*kt*0.5;
        end
        vtilde_sq = max(vtilde_sq, 1e-12);
        vtilde = sqrt(vtilde_sq);
        sqrtT  = sqrt(tau);
        vsT    = vtilde * sqrtT;
        X  = log(F/K);
        dp = (X + 0.5*vtilde_sq*tau) / vsT;
        dm = dp - vsT;
        BS_call = F*normcdf(dp) - K*normcdf(dm);
        U_t = (rho*sigma)/(2*kappa^2) * ( ...
              theta*kt - 2*theta + v0 + D*(2*theta - v0) - kt*D*(v0 - theta) );
        R_t = (sigma^2)/(8*kappa^2) * ( ...
              theta*tau + (v0 - theta)/kappa*(1 - D) - 2*theta/kappa*(1 - D) ...
              - 2*(v0 - theta)*tau*D + theta/(2*kappa)*(1 - D^2) ...
              + (v0 - theta)/kappa*(D - D^2) );
        exp_x = F/K;
        common = exp_x * inv_sqrt2pi / vsT * exp(-0.5*dp^2);
        H_fun = common * (1 - dp/vsT);
        K_fun = common * (dp^2 - dp*vsT - 1) / (vtilde_sq * tau);
        C_fwd = BS_call + K * (H_fun*U_t + K_fun*R_t);
        if CP == 1
            P_fwd(i) = C_fwd;
        else
            P_fwd(i) = C_fwd - (F - K);
        end
    end
end
