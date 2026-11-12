clearvars; clc;
tic
% ROUGH HESTON - SPX 2026-05-15. Calibra H, nu, rho (kappa = 0, Pade 3,3).
% xi0(t) dal PCHIP di W(T) = T*V_fair(T), con V_fair su tutti gli strike.
REDUCED_GRID = true;

fprintf('\nCaricamento dataset Refinitiv SPX 2026-05-15 (merged)...\n');

CSV_PATH   = '20260515SPX_merged.csv';
QUOTE_DATE = datenum('2026-05-15');

raw = readtable(CSV_PATH, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
rng(42);

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
    T_avail  = unique(T_surf);
    T_target = [0.0384 0.0833 0.25 0.5 0.75 1.0 1.5];
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
fprintf('Opzioni disponibili (tutti gli strike): %d | Maturities: %d\n\n', ...
    length(F_surf), length(T_uniq));

% STEP 2: ESTRAZIONE V_fair(T_i) - DEMETERFI su TUTTI gli strike
% V_fair calcolata PRIMA della selezione a 15 strike, usando tutto il
% campione filtrato (logM_cut, T_min/max).
fprintf('Estrazione V_fair(T) via Demeterfi (su tutti gli strike disponibili)...\n');

K_TAIL_LOW  = 0.05;
K_TAIL_HIGH = 2.50;
N_GRID      = 4000;

V_fair     = zeros(length(T_uniq),1);
IV_ATM_mkt = zeros(length(T_uniq),1);

for i = 1:length(T_uniq)
    T_i = T_uniq(i);
    m   = (T_surf == T_i);  % TUTTI gli strike per questa maturity
    K_i = K_surf(m);
    IV_i = IV_Mkt(m);
    F_i = F_surf(find(m,1));

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
    V_fair(i) = (2/T_i) * trapz(K_grid, OTM_prices ./ K_grid.^2);
end

fprintf('V_fair per scadenza (Demeterfi, tutti gli strike):\n');
fprintf('%-10s  %-10s  %-12s  %-10s\n', 'T(yr)', 'V_fair', 'sqrt(Vf)%', 'IV_ATM%');
for i = 1:length(T_uniq)
    fprintf('%-10.4f  %-10.5f  %-12.3f  %-10.3f\n', ...
        T_uniq(i), V_fair(i), sqrt(max(V_fair(i),0))*100, IV_ATM_mkt(i)*100);
end

% CURVA FORWARD VARIANCE xi0(t) = dW/dt
T_anchor = [0; T_uniq(:)];
W_anchor = [0; T_uniq(:) .* V_fair(:)];
if any(diff(W_anchor) <= 0)
    warning('W(T) non monotona: la forward variance potrebbe essere negativa.');
end
W_pp = pchip(T_anchor, W_anchor);

fprintf('\nCurva xi0(t) = dW/dt (PCHIP su W = T*V_fair):\n');
fprintf('%-10s  %-12s  %-12s  %-12s  %-10s\n', ...
    'T(yr)','W=T*Vfair','xi0(T)','sqrt(xi0)%','V_fair');
for i = 1:length(T_uniq)
    xi_i = eval_xi0_pchip(T_uniq(i), W_pp);
    fprintf('%-10.4f  %-12.5f  %-12.5f  %-12.3f  %-10.5f\n', ...
        T_uniq(i), T_uniq(i)*V_fair(i), xi_i, sqrt(max(xi_i,0))*100, V_fair(i));
end

% SELEZIONE 15 STRIKE PER SCADENZA (moneyness standardizzata)
% V_fair e la curva xi0 sono gia' state calcolate sopra su tutti gli strike.
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
fprintf('Selezione %d strike/scadenza: %d opzioni totali\n\n', ...
    N_PER_MAT, length(F_surf));

% STEP 3: PESI
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
fprintf('Pesi winsorizzati al 99mo percentile: %d punti cappati\n\n', n_capped);

% STEP 4: LOSS
% Xc=12 fisso, pow separato per fast/exact
XC_FAST        = 6;    POW_FAST = 10;  % N = 2^10 = 1024, dx = 2*Xc/N = 0.0117
XC_EXC         = 12;   POW_EXC  = 12;  % N = 2^12 = 4096
N_riccati_fast = 100;
N_riccati_exc  = 400;  % 0.03 bp di errore di quadratura contro n=6400
                       % (a n=200 sarebbero 0.09 bp, a n=100 0.25 bp)
Residui_Fast = @(p) sqrt_w .* ( ...
    PrezziRoughHestonXi0(F_surf, T_surf, K_surf, p, T_uniq, W_pp, ...
                         N_riccati_fast, XC_FAST, POW_FAST, CP_flag) - Mid_surf);

Residui_Exact = @(p) sqrt_w .* ( ...
    PrezziRoughHestonXi0(F_surf, T_surf, K_surf, p, T_uniq, W_pp, ...
                         N_riccati_exc, XC_EXC, POW_EXC, CP_flag) - Mid_surf);

Loss_Fast_GA = @(p) sum(Residui_Fast(p).^2);

% STEP 5: CALIBRAZIONE 3-STEP (H, nu, rho)
fprintf('=== Calibrazione rough Heston (3 params, xi0(t) curva PCHIP) ===\n');
lb    = [0.005, 0.05, -0.99];
ub    = [0.49,  5.00,  0.10];
guess = [0.05,  0.40, -0.65];
nvars = 3;

fprintf('\n[Step 1] GA (10 gen, pop 40)...\n');
options_ga = optimoptions('ga', ...
    'MaxGenerations', 10, 'PopulationSize', 40, ...
    'InitialPopulationMatrix', guess, ...
    'Display', 'iter', 'FunctionTolerance', 1e-9);
t_ga = tic;
theta_GA = ga(Loss_Fast_GA, nvars, [], [], [], [], lb, ub, [], options_ga);
fprintf('Step 1 in %.1f s\n', toc(t_ga));

fprintf('\n[Step 2] lsqnonlin fast...\n');
options_lsq = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunctionEvaluations', 1500, 'Display', 'iter', ...
    'StepTolerance', 1e-8, 'FunctionTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-8);
t_lsq1 = tic;
theta_LSQ1 = lsqnonlin(Residui_Fast, theta_GA, lb, ub, options_lsq);
fprintf('Step 2 in %.1f s\n', toc(t_lsq1));

fprintf('\n[Step 3] lsqnonlin exact...\n');
t_lsq2 = tic;
theta_ottimi = lsqnonlin(Residui_Exact, theta_LSQ1, lb, ub, options_lsq);
fprintf('Step 3 in %.1f s\n', toc(t_lsq2));

H_opt   = theta_ottimi(1);
nu_opt  = theta_ottimi(2);
rho_opt = theta_ottimi(3);

fprintf('\n=== PARAMETRI ROUGH HESTON (xi0 curva PCHIP) ===\n');
fprintf('H (Hurst)        : %.4f\n', H_opt);
fprintf('nu (vol-of-vol)  : %.4f\n', nu_opt);
fprintf('rho              : %.4f\n', rho_opt);
fprintf('alpha = H+1/2    : %.4f\n', H_opt + 0.5);

bnd = {};
if H_opt   < lb(1)+1e-3, bnd{end+1}='H@lb';   end
if H_opt   > ub(1)-1e-3, bnd{end+1}='H@ub';   end
if nu_opt  < lb(2)+1e-3, bnd{end+1}='nu@lb';  end
if nu_opt  > ub(2)-1e-3, bnd{end+1}='nu@ub';  end
if rho_opt < lb(3)+1e-3, bnd{end+1}='rho@lb'; end
if rho_opt > ub(3)-1e-3, bnd{end+1}='rho@ub'; end
if isempty(bnd)
    fprintf('Boundary check   : tutti interni\n');
else
    fprintf('Boundary check   : %s\n', strjoin(bnd,', '));
end

lf = sum(Residui_Fast(theta_ottimi).^2);
le = sum(Residui_Exact(theta_ottimi).^2);
fprintf('Sanity Fast/Exact: %.4e / %.4e  Ratio=%.2f\n', lf, le, lf/le);

% salvataggio
cartella_corrente = fileparts(mfilename('fullpath'));
cartella_results  = fullfile(cartella_corrente, '..', 'Results');
if ~exist(cartella_results, 'dir'), mkdir(cartella_results); end
save(fullfile(cartella_results, 'calib_rheston_7x15z_curva.mat'), ...
     'theta_ottimi', 'REDUCED_GRID', 'T_uniq', 'V_fair');
fprintf('\nParametri salvati in: %s (REDUCED_GRID=%d)\n', ...
    fullfile(cartella_results, 'calib_rheston_7x15z_curva.mat'), REDUCED_GRID);

% STEP 6: PREZZI E IV MODELLO
fprintf('\nRicalcolo prezzi e inversione IV...\n');
PrezziModello = PrezziRoughHestonXi0(F_surf, T_surf, K_surf, theta_ottimi, ...
    T_uniq, W_pp, N_riccati_exc, XC_EXC, POW_EXC, CP_flag);

IV_Mod = InvertiIV_OTM(F_surf, K_surf, T_surf, max(PrezziModello, 1e-6));
valid  = ~isnan(IV_Mod) & ~isnan(IV_Mkt);

% STEP 7: ERRORI
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
fprintf('%-9s  %-6s  %-12s  %-12s  %-12s  %-10s\n', ...
    'T(yr)','N','RMSE_rel%','MAE_rel%','RMSE_abs%','V_fair');
T_sorted  = sort(unique(T_surf(valid)));
rmse_by_T = zeros(length(T_sorted),1);
for i = 1:length(T_sorted)
    m = (T_surf == T_sorted(i)) & valid;
    e     = IV_Mod(m) - IV_Mkt(m);
    e_rel = -e ./ IV_Mkt(m);
    rmse_by_T(i) = sqrt(mean(e_rel.^2))*100;
    Vf_here = V_fair(T_uniq == T_sorted(i));
    fprintf('%-9.4f  %-6d  %-12.4f  %-12.4f  %-12.4f  %-10.5f\n', ...
        T_sorted(i), sum(m), rmse_by_T(i), mean(abs(e_rel))*100, sqrt(mean(e.^2))*100, Vf_here);
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
figure('Name','V_fair e curva xi0','Color','w','Position',[100 100 900 450]);
stairs([T_uniq; max(T_uniq)*1.05],[V_fair; V_fair(end)],'b-','LineWidth',2,...
    'DisplayName','V_{fair}(T) (variance swap)'); hold on; grid on;
plot(T_uniq, V_fair,'ro','MarkerSize',8,'LineWidth',1.5,'DisplayName','V_{fair} market');
plot(T_uniq, IV_ATM_mkt.^2,'g^','MarkerSize',8,'LineWidth',1.5,'DisplayName','IV_{ATM}^2');
xlabel('T (years)'); ylabel('variance'); legend('Location','best');
title('V_{fair}(T) via Demeterfi (su tutti gli strike disponibili)');

figure('Name','Superficie IV 3D','Color','w','Position',[100 100 900 600]);
scatter3(logM_surf(valid), T_surf(valid), IV_Mkt(valid)*100, 25,'r','o','DisplayName','Mercato');
hold on; grid on;
scatter3(logM_surf(valid), T_surf(valid), IV_Mod(valid)*100, 25,'b','filled','DisplayName','rHeston');
xlabel('ln(K/F)'); ylabel('T (yr)'); zlabel('IV (%)');
title(sprintf('Rough Heston (xi0 PCHIP): H=%.3f nu=%.3f rho=%.3f', H_opt, nu_opt, rho_opt));
legend('Location','best'); view(-45,20);

ScadenzePlt = unique(T_surf(valid));
num_plots = length(ScadenzePlt);
cols = min(4, num_plots); rows = ceil(num_plots/cols);
figure('Name','Smiles IV rough Heston','Color','w','Position',[50 50 1400 800]);
for i = 1:num_plots
    T_tgt = ScadenzePlt(i);
    mask  = (T_surf == T_tgt) & valid;
    [logM_sorted, sidx] = sort(logM_surf(mask));
    IVm = IV_Mkt(mask); IVm = IVm(sidx);
    IVh = IV_Mod(mask); IVh = IVh(sidx);
    subplot(rows, cols, i);
    plot(logM_sorted, IVm*100,'ro','MarkerSize',4,'LineWidth',1,'DisplayName','Mercato'); hold on;
    plot(logM_sorted, IVh*100,'b-','LineWidth',1.5,'DisplayName','rHeston');
    xline(0,'k--','F','HandleVisibility','off');
    title(sprintf('T=%.4f (N=%d)', T_tgt, sum(mask)),'FontSize',9);
    ylabel('IV (%)'); xlabel('ln(K/F)'); grid on;
    if i==1, legend('Location','best','FontSize',8); end
end

figure('Name','RMSE IV per T','Color','w','Position',[200 200 900 450]);
bar(T_sorted, rmse_by_T,'FaceColor',[0.3 0.6 0.4]); grid on;
xlabel('Maturity (years)'); ylabel('RMSE IV (%)');
title('RMSE IV per scadenza - rough Heston (xi0 PCHIP) - SPX 2026-05-15');
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

function P_fwd = PrezziRoughHestonXi0(F_vec, T_vec, K_vec, params, T_uniq, W_pp, ...
                                        N_riccati, Xc, pow, CP_vec)
% Prezza rough Heston con la curva xi0(t) = dW/dt valutata sui nodi Riccati.
% params = [H, nu, rho]. Xc e pow sono separati.
    H = params(1); nu = params(2); rho = params(3);
    params_cf = [H, nu, rho, 0];

    P_fwd = zeros(length(K_vec),1);

    for k = 1:length(T_uniq)
        T   = T_uniq(k);
        idx = (T_vec == T);
        if ~any(idx), continue; end

        F_t  = F_vec(find(idx,1));
        K_t  = K_vec(idx);
        CP_t = CP_vec(idx);

        ti = ((0:N_riccati-1) + 0.5) / N_riccati * T;
        xi_grid = eval_xi0_pchip(ti, W_pp);
        xi_grid = xi_grid(:)';

        prezzi_t = PrezzaSuperficieSINC_yield(F_t, T*ones(size(K_t)), K_t, ...
            zeros(size(K_t)), zeros(size(K_t)), ...
            'rough Heston', params_cf, xi_grid, N_riccati, Xc, 0, pow, CP_t(1));

        P_fwd(idx) = prezzi_t;
    end
end
