clearvars; clc;
tic
% MULTISTART  Verifica dei minimi di rough Heston, CGMY e rough Heston + CGMY.
% Partenze Sobol nei bound piu' l'ottimo corrente, lsqnonlin fast, poi
% raffinamento exact dei migliori. I risultati sono raggruppati in bacini.
% Output: Results/multistart.mat. Tempo: circa 13 minuti.
% Il rough Hawkes Heston e' in Multistart_RHH.

N_START = 24;  % punti Sobol per modello (oltre al seed)
N_REFINE = 6;  % quanti dei migliori raffinare con la loss exact
TOL_BACINO = 0.05;  % distanza max-norm nello spazio normalizzato

fprintf('\nCaricamento dataset Refinitiv SPX 2026-05-15 (merged)...\n');

CSV_PATH   = '20260515SPX_merged.csv';
QUOTE_DATE = datenum('2026-05-15');

raw = readtable(CSV_PATH, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
rng(42);

is_spot = strcmp(raw.Root, 'SPOT');
S0_feed = raw.CF_LAST(is_spot);
assert(isscalar(S0_feed), 'Attesa una sola riga SPOT nel CSV, trovate %d.', numel(S0_feed));
fprintf('Livello indice da feed (non sincrono): %.2f\n', S0_feed);
raw     = raw(~is_spot, :);

is_call = strcmp(raw.PutCall, 'C');
keep    = ~isnan(raw.STRIKE_PRC) & raw.CF_BID > 0 & ...
          raw.CF_ASK > raw.CF_BID & raw.CF_BID >= 0.375;
raw     = raw(keep, :);
is_call = is_call(keep);

K_all   = raw.STRIKE_PRC;
mid_all = 0.5*(raw.CF_BID + raw.CF_ASK);
T_all   = (datenum(raw.EXPIR_DATE) - QUOTE_DATE) / 365;

[~, ia] = unique([T_all, K_all, double(is_call)], 'rows', 'stable');
T_all = T_all(ia); K_all = K_all(ia);
mid_all = mid_all(ia); is_call = is_call(ia);

T_uniq_raw = unique(T_all);
F_per_T = nan(length(T_uniq_raw),1);
D_per_T = nan(length(T_uniq_raw),1);
for i = 1:length(T_uniq_raw)
    Ti = T_uniq_raw(i);
    mC = (T_all == Ti) & is_call;
    mP = (T_all == Ti) & ~is_call;
    [Kcommon, iC, iP] = intersect(K_all(mC), K_all(mP));
    if length(Kcommon) < 3, continue; end
    cIdx = find(mC); pIdx = find(mP);
    diff_CP = mid_all(cIdx(iC)) - mid_all(pIdx(iP));
    near = abs(Kcommon/S0_feed - 1) <= 0.25;
    if sum(near) >= 5, Kreg = Kcommon(near); dreg = diff_CP(near);
    else,              Kreg = Kcommon;       dreg = diff_CP; end
    coef = [ones(size(Kreg)), Kreg] \ dreg;
    D_per_T(i) = -coef(2);
    F_per_T(i) = coef(1) / D_per_T(i);
end
valid_T = ~isnan(F_per_T) & D_per_T > 0 & D_per_T < 1.05 & F_per_T > 0;
T_uniq_raw = T_uniq_raw(valid_T);
F_per_T = F_per_T(valid_T);  D_per_T = D_per_T(valid_T);

F_surf = []; T_surf = []; K_surf = []; Mid_surf = [];
for i = 1:length(T_uniq_raw)
    Ti = T_uniq_raw(i); Fi = F_per_T(i); Di = D_per_T(i);
    mTi = (T_all == Ti);
    Ki = K_all(mTi); midi = mid_all(mTi); calli = is_call(mTi);
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

logM = log(K_surf ./ F_surf);
idx  = ~isnan(IV_Mkt) & IV_Mkt > 0.01 & IV_Mkt < 2 & Mid_surf > 1e-4 & ...
       abs(logM) < 0.5 & T_surf >= 0.03 & T_surf <= 3.0;
F_surf=F_surf(idx); T_surf=T_surf(idx); K_surf=K_surf(idx);
Mid_surf=Mid_surf(idx); IV_Mkt=IV_Mkt(idx);

T_uniq_all = unique(T_surf);
mask_all = false(length(T_surf),1);
for i = 1:length(T_uniq_all)
    m = T_surf == T_uniq_all(i);
    if sum(m) >= 5, mask_all = mask_all | m; end
end
F_surf=F_surf(mask_all); T_surf=T_surf(mask_all); K_surf=K_surf(mask_all);
Mid_surf=Mid_surf(mask_all); IV_Mkt=IV_Mkt(mask_all);

T_avail = unique(T_surf);
T_sel = [];
for tt = [0.0384 0.0833 0.25 0.5 0.75 1.0 1.5]
    [~, ii] = min(abs(T_avail - tt));
    T_sel = [T_sel; T_avail(ii)];
end
T_sel = unique(T_sel);
keepm = ismember(T_surf, T_sel);
F_surf=F_surf(keepm); T_surf=T_surf(keepm); K_surf=K_surf(keepm);
Mid_surf=Mid_surf(keepm); IV_Mkt=IV_Mkt(keepm);
T_uniq = unique(T_surf);

% Demeterfi + curva xi0
Kvar_mkt   = zeros(length(T_uniq),1);
IV_ATM_mkt = zeros(length(T_uniq),1);
for i = 1:length(T_uniq)
    T_i = T_uniq(i);
    m   = (T_surf == T_i);
    K_i = K_surf(m); IV_i = IV_Mkt(m); F_i = F_surf(find(m,1));
    [K_s, srt] = sort(K_i); IV_s = IV_i(srt);
    [~, atm_idx] = min(abs(K_s - F_i));
    IV_ATM_mkt(i) = IV_s(atm_idx);
    K_grid  = logspace(log10(0.05*F_i), log10(2.50*F_i), 4000)';
    IV_grid = interp1(K_s, IV_s, K_grid, 'linear', NaN);
    IV_grid(K_grid < K_s(1))   = IV_s(1);
    IV_grid(K_grid > K_s(end)) = IV_s(end);
    OTM = zeros(size(K_grid));
    for j = 1:length(K_grid)
        if K_grid(j) < F_i
            [~, p] = blkprice(F_i, K_grid(j), 0, T_i, IV_grid(j)); OTM(j) = p;
        else
            [c, ~] = blkprice(F_i, K_grid(j), 0, T_i, IV_grid(j)); OTM(j) = c;
        end
    end
    Kvar_mkt(i) = (2/T_i) * trapz(K_grid, OTM ./ K_grid.^2);
end
W_pp     = pchip([0; T_uniq(:)], [0; T_uniq(:) .* Kvar_mkt(:)]);
xi0_eval = @(tvec) eval_xi0_pchip(tvec, W_pp);

% selezione 15 strike/scadenza
keep15 = false(size(T_surf));
for i = 1:length(T_uniq)
    m = find(T_surf == T_uniq(i));
    z = log(K_surf(m) ./ F_surf(m)) ./ (IV_ATM_mkt(i) * sqrt(T_uniq(i)));
    if numel(m) <= 15, keep15(m) = true; continue; end
    for tt = linspace(max(-3.0, min(z)), min(2.0, max(z)), 15)
        [~, jj] = min(abs(z - tt));
        keep15(m(jj)) = true;
    end
end
F_surf=F_surf(keep15); T_surf=T_surf(keep15); K_surf=K_surf(keep15);
Mid_surf=Mid_surf(keep15); IV_Mkt=IV_Mkt(keep15);
CP_flag = -1 * ones(size(F_surf));
fprintf('Campione: %d opzioni su %d scadenze\n', length(F_surf), length(T_uniq));

% pesi
Vega_Mkt = zeros(length(F_surf),1);
for i = 1:length(F_surf)
    d1 = (log(F_surf(i)/K_surf(i)) + 0.5*IV_Mkt(i)^2*T_surf(i)) / ...
         (IV_Mkt(i)*sqrt(T_surf(i)));
    Vega_Mkt(i) = F_surf(i) * sqrt(T_surf(i)) * normpdf(d1);
end
sqrt_w = 1 ./ (max(Vega_Mkt,1e-4) .* max(IV_Mkt,0.01));
sqrt_w = min(sqrt_w, prctile(sqrt_w, 99));

resdir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
if isempty(fileparts(mfilename('fullpath'))), resdir = fullfile(pwd,'..','Results'); end

opt_lsq = optimoptions('lsqnonlin', ...
    'Algorithm','trust-region-reflective', ...
    'MaxFunctionEvaluations',2000,'MaxIterations',150,'Display','off', ...
    'StepTolerance',1e-8,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);
[XC_RH,  POW_RH ] = XcPow('rHeston', 'exact');
[XC_RHC, POW_RHC] = XcPow('rH+CGMY', 'exact');
[XC_CG,  POW_CG ] = XcPow('CGMY',    'exact');
% xi_min e vincolo di ammissibilita'
% Deve coincidere con quello calcolato in Rough_Heston_CMGY_def, altrimenti
% il multistart esplora una specificazione diversa da quella calibrata.
Q_MAX = 0.98;
tt_all = [];
for k = 1:length(T_uniq)
    tt_all = [tt_all, ((0:399)+0.5)/400*T_uniq(k), ((0:99)+0.5)/100*T_uniq(k)];
end
xi_min = min(xi0_eval(tt_all));
fprintf('min xi0 sulle griglie Riccati = %.6f | vincolo Lcg <= %.6f\n', ...
    xi_min, Q_MAX*xi_min);

OUT = struct();

% MODELLO 1 - ROUGH HESTON (H, nu, rho)
fprintf('\n========== MULTISTART: rough Heston (3 parametri) ==========\n');
lb1 = [0.005, 0.05, -0.99];
ub1 = [0.49,  5.00,  0.10];
log1 = [true, true, false];  % H e nu coprono ordini di grandezza
nomi1 = {'H','nu','rho'};

c1 = load(fullfile(resdir,'calib_rheston_7x15z_curva.mat'), 'theta_ottimi');
seed1 = min(max(c1.theta_ottimi(:)', lb1), ub1);

r1F = @(p) sqrt_w .* (PrezziRoughHeston(F_surf,T_surf,K_surf,p,T_uniq, ...
            xi0_eval,100,XC_RH,POW_RH,CP_flag) - Mid_surf);
r1E = @(p) sqrt_w .* (PrezziRoughHeston(F_surf,T_surf,K_surf,p,T_uniq, ...
            xi0_eval,400,XC_RH,POW_RH,CP_flag) - Mid_surf);

OUT.rheston = EseguiMultistart(r1F, r1E, lb1, ub1, log1, seed1, ...
                               N_START, N_REFINE, TOL_BACINO, opt_lsq, nomi1);

% MODELLO 2 - ROUGH HESTON + CGMY (H, nu, rho, C, M, Y), G = 0.05
fprintf('\n========== MULTISTART: rough Heston + CGMY (6 parametri) ==========\n');
G_FIX = 0.05;
%       [H,     nu,   rho,   q,     M,     Y   ]
lb2 = [0.005, 0.05, -0.99, 1e-4,  5.0,   0.9 ];
ub2 = [0.49,  5.00,  0.10, Q_MAX, 200.0, 1.95];
log2 = [true, true, false, true, true, false];
nomi2 = {'H','nu','rho','q','M','Y'};

% Il .mat salva C; qui serve q.
c2 = load(fullfile(resdir,'calib_rhcgmy_7x15z.mat'), 'theta_ottimi');
t7 = c2.theta_ottimi;  % [H,nu,rho,C,G,M,Y]
q2 = t7(4) * LambdaCGMY(t7(5), t7(6), t7(7)) / xi_min;
seed2 = min(max([t7(1:3), q2, t7(6), t7(7)], lb2), ub2);
fprintf('seed: q = %.4f (C = %.6g), M = %.3f\n', seed2(4), t7(4), seed2(5));

r2F = @(p) sqrt_w .* (PrezziRoughHestonCGMY(F_surf,T_surf,K_surf,p,G_FIX,xi_min, ...
            T_uniq,xi0_eval,100,XC_RHC,POW_RHC,CP_flag) - Mid_surf);
r2E = @(p) sqrt_w .* (PrezziRoughHestonCGMY(F_surf,T_surf,K_surf,p,G_FIX,xi_min, ...
            T_uniq,xi0_eval,400,XC_RHC,POW_RHC,CP_flag) - Mid_surf);

OUT.rhcgmy = EseguiMultistart(r2F, r2E, lb2, ub2, log2, seed2, ...
                              N_START, N_REFINE, TOL_BACINO, opt_lsq, nomi2);

% MODELLO 3 - CGMY PURO (C, M, Y), G = 0.05
% Bound allineati a CMGY_def, con la stessa restrizione economica M >= 5
% adottata nel modello composito.
fprintf('\n========== MULTISTART: CGMY puro (3 parametri) ==========\n');
lb3 = [1e-4,  5.0,  0.05];
ub3 = [50.0, 2000,  1.95];
log3 = [true, true, false];
nomi3 = {'C','M','Y'};

c3 = load(fullfile(resdir,'calib_cgmy_7x15z.mat'), 'theta_ottimi');
seed3 = min(max(c3.theta_ottimi([1 3 4]), lb3), ub3);  % [C,G,M,Y] -> [C,M,Y]

R_surf = zeros(size(F_surf));
r3 = @(p) sqrt_w .* (PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, ...
        R_surf, R_surf, 'CGMY', [p(1), G_FIX, p(2), p(3)], [], [], ...
        XC_CG, 0, POW_CG, CP_flag) - Mid_surf);

% Per il CGMY la CF e' in forma chiusa e non c'e' stadio fast/exact
% distinto: si passa la stessa funzione due volte.
OUT.cgmy = EseguiMultistart(r3, r3, lb3, ub3, log3, seed3, ...
                            N_START, N_REFINE, TOL_BACINO, opt_lsq, nomi3);

% FIGURE + SALVATAGGIO
nomi_mod = {'rough Heston','rH+CGMY','CGMY puro'};
campi    = {'rheston','rhcgmy','cgmy'};
colori   = {[0 0.45 0.74], [0.4 0.2 0.6], [0.7 0.3 0.5]};

figure('Name','Multistart','Color','w','Position',[100 100 1300 400]);
for i = 1:3
    subplot(1,3,i);
    R = OUT.(campi{i});
    semilogy(1:numel(R.loss_ord), R.loss_ord/min(R.loss_ord), 'o', ...
        'Color', colori{i}); grid on;
    xlabel('partenza (ordinata per loss)'); ylabel('loss / loss_{min}');
    title(sprintf('%s: %d bacini', nomi_mod{i}, R.n_bacini));
end

figure('Name','Primo parametro contro loss','Color','w','Position',[150 150 1300 400]);
for i = 1:3
    subplot(1,3,i);
    R = OUT.(campi{i});
    if i == 3, semilogx(R.theta_ord(:,1), R.loss_ord/min(R.loss_ord), 'o', ...
            'Color', colori{i});
    else,      semilogy(R.theta_ord(:,1), R.loss_ord/min(R.loss_ord), 'o', ...
            'Color', colori{i}); end
    grid on; xlabel(R.nomi{1}); ylabel('loss / loss_{min}'); title(nomi_mod{i});
end

if ~exist(resdir,'dir'), mkdir(resdir); end
save(fullfile(resdir,'multistart.mat'), 'OUT', 'N_START', 'TOL_BACINO', ...
     'xi_min', 'Q_MAX');
fprintf('\nRisultati salvati in %s\n', fullfile(resdir,'multistart.mat'));

toc

% FUNZIONI LOCALI
function R = EseguiMultistart(rFast, rExact, lb, ub, log_scale, seed, ...
                              N_START, N_REFINE, tol, opt_lsq, nomi)
    n = numel(lb);
    % punti di partenza: seed + Sobol
    P = sobolset(n, 'Skip', 1);
    P = scramble(P, 'MatousekAffineOwen');
    U = net(P, N_START);
    starts = zeros(N_START+1, n);
    starts(1,:) = seed;
    for k = 1:N_START
        for j = 1:n
            if log_scale(j)
                starts(k+1,j) = exp(log(lb(j)) + U(k,j)*(log(ub(j))-log(lb(j))));
            else
                starts(k+1,j) = lb(j) + U(k,j)*(ub(j)-lb(j));
            end
        end
    end

    % ottimizzazione FAST da ogni partenza
    nS    = size(starts,1);
    theta = nan(nS, n);
    loss  = inf(nS, 1);
    fprintf('%-5s  %-12s  %s\n', 'i', 'loss (fast)', 'parametri');
    for k = 1:nS
        try
            th = lsqnonlin(rFast, starts(k,:), lb, ub, opt_lsq);
            theta(k,:) = th;
            loss(k)    = sum(rFast(th).^2);
        catch ME
            fprintf('%-5d  FALLITA (%s)\n', k, ME.identifier);
            continue
        end
        eti = '';
        if k == 1, eti = '   <- seed (ottimo corrente)'; end
        fprintf('%-5d  %-12.6e  %s%s\n', k, loss(k), ...
            strjoin(compose('%.4g', theta(k,:)), ' '), eti);
    end

    ok = isfinite(loss);
    assert(any(ok), 'Tutte le partenze sono fallite.');
    [loss_ord, ord] = sort(loss(ok));
    th_ok = theta(ok,:);
    theta_ord = th_ok(ord,:);

    % clustering in bacini (max-norm nello spazio normalizzato)
    Z = theta_ord;
    for j = 1:n
        if log_scale(j)
            Z(:,j) = (log(Z(:,j))-log(lb(j))) / (log(ub(j))-log(lb(j)));
        else
            Z(:,j) = (Z(:,j)-lb(j)) / (ub(j)-lb(j));
        end
    end
    bacino = zeros(size(Z,1),1);
    nb = 0;
    for k = 1:size(Z,1)
        for b = 1:nb
            rap = find(bacino == b, 1);
            if max(abs(Z(k,:) - Z(rap,:))) < tol, bacino(k) = b; break; end
        end
        if bacino(k) == 0, nb = nb + 1; bacino(k) = nb; end
    end

    fprintf('\n--- BACINI DISTINTI (tol max-norm %.2f) ---\n', tol);
    fprintf('%-8s  %-6s  %-14s  %-10s  %s\n', ...
        'bacino','conteggio','loss migliore','vs globale','parametri');
    for b = 1:nb
        m = (bacino == b);
        [lb_best, ib] = min(loss_ord(m));
        idxs = find(m);
        th_b = theta_ord(idxs(ib),:);
        fprintf('%-8d  %-6d  %-14.6e  %-10.4f  %s\n', b, sum(m), lb_best, ...
            lb_best/min(loss_ord), strjoin(compose('%.4g', th_b), ' '));
    end

    % raffinamento EXACT sui migliori
    nref = min(N_REFINE, size(theta_ord,1));
    fprintf('\n--- RAFFINAMENTO EXACT sui %d migliori ---\n', nref);
    fprintf('%-5s  %-14s  %-14s  %s\n', 'i', 'loss fast', 'loss exact', 'parametri');
    th_ref = nan(nref, n);
    ls_ref = inf(nref, 1);
    for k = 1:nref
        try
            th = lsqnonlin(rExact, theta_ord(k,:), lb, ub, opt_lsq);
            th_ref(k,:) = th;
            ls_ref(k)   = sum(rExact(th).^2);
            fprintf('%-5d  %-14.6e  %-14.6e  %s\n', k, loss_ord(k), ls_ref(k), ...
                strjoin(compose('%.4g', th), ' '));
        catch
            fprintf('%-5d  raffinamento fallito\n', k);
        end
    end

    [best_loss, ibest] = min(ls_ref);
    fprintf('\nMIGLIORE (exact): loss = %.6e\n', best_loss);
    for j = 1:n
        fprintf('  %-4s = %.6g\n', nomi{j}, th_ref(ibest,j));
    end

    % confronto col seed
    seed_loss = sum(rExact(seed).^2);
    fprintf('Seed (ottimo corrente), loss exact = %.6e  ->  rapporto %.4f\n', ...
        seed_loss, seed_loss/best_loss);
    if seed_loss > best_loss*1.001
        fprintf(2, ['ATTENZIONE: il multistart ha trovato un ottimo MIGLIORE del\n' ...
                    '  valore attualmente salvato. La calibrazione va rifatta a\n' ...
                    '  partire da questi parametri prima di usarli in tesi.\n']);
    else
        fprintf('L''ottimo corrente e'' confermato come il migliore trovato.\n');
    end

    % escursione dei parametri fra soluzioni entro l'1%
    vicini = ls_ref <= 1.01*best_loss;
    if sum(vicini) > 1
        fprintf('\nParametri fra le soluzioni entro l''1%% dal minimo (%d soluzioni):\n', ...
            sum(vicini));
        for j = 1:n
            v = th_ref(vicini,j);
            fprintf('  %-4s in [%.6g, %.6g]  (fattore %.2f)\n', ...
                nomi{j}, min(v), max(v), max(v)/min(v));
        end
    end

    R = struct('starts',starts,'theta',theta,'loss',loss, ...
               'theta_ord',theta_ord,'loss_ord',loss_ord,'bacino',bacino, ...
               'n_bacini',nb,'theta_ref',th_ref,'loss_ref',ls_ref, ...
               'best',th_ref(ibest,:),'best_loss',best_loss, ...
               'seed',seed,'seed_loss',seed_loss,'nomi',{nomi});
end

function IV = InvertiIV_OTM(F_vec, K_vec, T_vec, P_fwd_vec)
    IV = nan(size(F_vec));
    for i = 1:length(F_vec)
        try
            if K_vec(i) <= F_vec(i)
                IV(i) = blkimpv(F_vec(i),K_vec(i),0,T_vec(i),P_fwd_vec(i),'Class','put');
            else
                C_fwd = max(P_fwd_vec(i) + F_vec(i) - K_vec(i), 1e-6);
                IV(i) = blkimpv(F_vec(i),K_vec(i),0,T_vec(i),C_fwd,'Class','call');
            end
        catch
            IV(i) = NaN;
        end
    end
end

function xi = eval_xi0_pchip(tvec, W_pp)
    h = 1e-5;
    tvec = tvec(:);
    tplus = max(tvec+h,0); tminus = max(tvec-h,0);
    denom = tplus - tminus; denom(denom < 1e-12) = 1e-12;
    xi = max((ppval(W_pp,tplus) - ppval(W_pp,tminus)) ./ denom, 1e-8);
end

function lam = LambdaCGMY(G, M, Y)
% Lcg = C * lam(G,M,Y): il log-contract e' lineare in C.
    lam = 2*( gamma(-Y) * ((M-1)^Y - M^Y + (G+1)^Y - G^Y) ...
            - gamma(1-Y) * (M^(Y-1) - G^(Y-1)) );
end

function [Lcg, Vcg] = ContributiCGMY(C, G, M, Y)
    Lcg = C * LambdaCGMY(G, M, Y);
    Vcg = C * gamma(2-Y) * (M^(Y-2) + G^(Y-2));
end

function P_fwd = PrezziRoughHeston(F_vec, T_vec, K_vec, params3, T_uniq, ...
                                   xi0_fn, N_riccati, Xc, pow, CP_vec)
    params_cf = [params3(1), params3(2), params3(3), 0];
    P_fwd = zeros(length(K_vec),1);
    for k = 1:length(T_uniq)
        T = T_uniq(k);
        idx = (T_vec == T);
        if ~any(idx), continue; end
        F_t = F_vec(find(idx,1)); K_t = K_vec(idx); CP_t = CP_vec(idx);
        ti = ((0:N_riccati-1) + 0.5) / N_riccati * T;
        xi_grid = xi0_fn(ti); xi_grid = xi_grid(:)';
        P_fwd(idx) = PrezzaSuperficieSINC_yield(F_t, T*ones(size(K_t)), K_t, ...
            zeros(size(K_t)), zeros(size(K_t)), ...
            'rough Heston', params_cf, xi_grid, N_riccati, Xc, 0, pow, CP_t(1));
    end
end

function P_fwd = PrezziRoughHestonCGMY(F_vec, T_vec, K_vec, params6, G_fix, xi_min, ...
                    T_uniq, xi0_fn, N_riccati, Xc, pow, CP_vec)
% params6 = [H, nu, rho, q, M, Y]: q e' la quota di xi_min assorbita dal
% log-contract dei salti, C si ricava da q.
    H = params6(1); nu = params6(2); rho = params6(3);
    q = params6(4); M  = params6(5); Y   = params6(6);
    G = G_fix;
    C = q * xi_min / LambdaCGMY(G, M, Y);
    params_combined = [H, nu, rho, 0, C, G, M, Y];
    L_cgmy = ContributiCGMY(C, G, M, Y);
    P_fwd = zeros(length(K_vec),1);
    for k = 1:length(T_uniq)
        T = T_uniq(k);
        idx = (T_vec == T);
        if ~any(idx), continue; end
        F_t = F_vec(find(idx,1)); K_t = K_vec(idx); CP_t = CP_vec(idx);
        ti = ((0:N_riccati-1) + 0.5) / N_riccati * T;
        xi_grid = max(xi0_fn(ti) - L_cgmy, 1e-6); xi_grid = xi_grid(:)';
        P_fwd(idx) = PrezzaSuperficieSINC_yield(F_t, T*ones(size(K_t)), K_t, ...
            zeros(size(K_t)), zeros(size(K_t)), ...
            'rough Heston CGMY', params_combined, xi_grid, N_riccati, Xc, 0, pow, CP_t(1));
    end
end
