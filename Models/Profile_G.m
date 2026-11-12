clearvars; clc;
tic
% PROFILE_G  Profilo della loss in G per CGMY puro e rough Heston + CGMY.
% Per ogni G della griglia si riottimizzano gli altri parametri.
% Profilo piatto: G non identificato, fissarlo e' una convenzione.
% Ogni punto e' ricontrollato con Xc = 30 (colonna Xc_ratio).
% Output: Results/profile_G.mat e due figure.

fprintf('\nCaricamento dataset Refinitiv SPX 2026-05-15 (merged)...\n');

CSV_PATH   = '20260515SPX_merged.csv';
QUOTE_DATE = datenum('2026-05-15');

raw = readtable(CSV_PATH, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
rng(42);

% Livello dell'indice dal feed (CF_LAST della riga Root = SPOT)
% Non sincrono con lo snapshot delle catene; serve solo a centrare la
% finestra +/-25% della regressione put-call parity.
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

LOGM_CUT = 0.5;  MIN_K_MAT = 5;  T_MIN = 0.03;  T_MAX = 3.0;
logM = log(K_surf ./ F_surf);
idx  = ~isnan(IV_Mkt) & IV_Mkt > 0.01 & IV_Mkt < 2 & ...
       Mid_surf > 1e-4 & abs(logM) < LOGM_CUT & ...
       T_surf >= T_MIN & T_surf <= T_MAX;
F_surf=F_surf(idx); T_surf=T_surf(idx); K_surf=K_surf(idx);
Mid_surf=Mid_surf(idx); IV_Mkt=IV_Mkt(idx);

T_uniq_all = unique(T_surf);
mask_all = false(length(T_surf),1);
for i = 1:length(T_uniq_all)
    m = T_surf == T_uniq_all(i);
    if sum(m) >= MIN_K_MAT, mask_all = mask_all | m; end
end
F_surf=F_surf(mask_all); T_surf=T_surf(mask_all); K_surf=K_surf(mask_all);
Mid_surf=Mid_surf(mask_all); IV_Mkt=IV_Mkt(mask_all);

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
T_uniq = unique(T_surf);

% Demeterfi + curva xi0 (identici a Rough_Heston_CMGY_def)
K_TAIL_LOW = 0.05;  K_TAIL_HIGH = 2.50;  N_GRID = 4000;
Kvar_mkt   = zeros(length(T_uniq),1);
IV_ATM_mkt = zeros(length(T_uniq),1);
for i = 1:length(T_uniq)
    T_i = T_uniq(i);
    m   = (T_surf == T_i);
    K_i = K_surf(m); IV_i = IV_Mkt(m); F_i = F_surf(find(m,1));
    [K_s, srt] = sort(K_i); IV_s = IV_i(srt);
    [~, atm_idx] = min(abs(K_s - F_i));
    IV_ATM_mkt(i) = IV_s(atm_idx);
    K_grid  = logspace(log10(K_TAIL_LOW*F_i), log10(K_TAIL_HIGH*F_i), N_GRID)';
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
T_anchor = [0; T_uniq(:)];
W_anchor = [0; T_uniq(:) .* Kvar_mkt(:)];
W_pp     = pchip(T_anchor, W_anchor);
xi0_eval = @(tvec) eval_xi0_pchip(tvec, W_pp);

% selezione 15 strike/scadenza
N_PER_MAT = 15;  Z_LO = -3.0;  Z_HI = 2.0;
keep15 = false(size(T_surf));
for i = 1:length(T_uniq)
    m = find(T_surf == T_uniq(i));
    z = log(K_surf(m) ./ F_surf(m)) ./ (IV_ATM_mkt(i) * sqrt(T_uniq(i)));
    if numel(m) <= N_PER_MAT, keep15(m) = true; continue; end
    z_targets = linspace(max(Z_LO, min(z)), min(Z_HI, max(z)), N_PER_MAT);
    for tt = z_targets
        [~, jj] = min(abs(z - tt));
        keep15(m(jj)) = true;
    end
end
F_surf=F_surf(keep15); T_surf=T_surf(keep15); K_surf=K_surf(keep15);
Mid_surf=Mid_surf(keep15); IV_Mkt=IV_Mkt(keep15);
R_surf  = zeros(size(F_surf));
Q_surf  = zeros(size(F_surf));
CP_flag = -1 * ones(size(F_surf));
fprintf('Campione del profile: %d opzioni su %d scadenze\n', ...
    length(F_surf), length(T_uniq));

% pesi
N_valid  = length(F_surf);
Vega_Mkt = zeros(N_valid,1);
for i = 1:N_valid
    d1 = (log(F_surf(i)/K_surf(i)) + 0.5*IV_Mkt(i)^2*T_surf(i)) / ...
         (IV_Mkt(i)*sqrt(T_surf(i)));
    Vega_Mkt(i) = F_surf(i) * sqrt(T_surf(i)) * normpdf(d1);
end
Vega_Mkt = max(Vega_Mkt, 1e-4);
sqrt_w   = 1 ./ (Vega_Mkt .* max(IV_Mkt, 0.01));
sqrt_w   = min(sqrt_w, prctile(sqrt_w, 99));

% GRIGLIA DI G E OPZIONI
% Log-spaziata su tre ordini di grandezza e contenente 0.05, il valore
% attualmente fissato.
G_grid = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50];
iG0    = find(abs(G_grid - 0.05) < 1e-12, 1);
assert(~isempty(iG0), 'La griglia deve contenere G = 0.05.');

% Profile B (rH+CGMY): Xc = 12 basta, Xc_ratio misurato 1.0007 su tutta la
% griglia perche' C ~ 5.7e-4.
% Profile B: rH+CGMY, stessa configurazione della calibrazione.
[XC,       POW      ] = XcPow('rH+CGMY', 'exact');
[XC_GUARD, POW_GUARD] = XcPow('rH+CGMY', 'guard');
% Profile A: CGMY puro.
[XC_A,       POW_A      ] = XcPow('CGMY', 'exact');
[XC_A_GUARD, POW_A_GUARD] = XcPow('CGMY', 'guard');

N_RIC_F = 100;
N_RIC_E = 400;

opt_lsq = optimoptions('lsqnonlin', ...
    'Algorithm','trust-region-reflective', ...
    'MaxFunctionEvaluations',3000,'MaxIterations',200,'Display','off', ...
    'StepTolerance',1e-8,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);

% PROFILE A - CGMY PURO  (liberi: C, M, Y)
fprintf('\n=== PROFILE A: CGMY puro, liberi C M Y  (Xc=%d) ===\n', XC_A);
lbA = [1e-4, 1.01, 0.05];
ubA = [10.0, 200.0, 1.95];

% Punti di partenza letti dai .mat delle calibrazioni.
resdir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
if isempty(fileparts(mfilename('fullpath'))), resdir = fullfile(pwd,'..','Results'); end
fA  = fullfile(resdir, 'calib_cgmy_7x15z.mat');
assert(isfile(fA), 'Manca %s: lancia prima CMGY_def.', fA);
cA  = load(fA, 'theta_ottimi');
thA = cA.theta_ottimi([1 3 4]);  % [C,G,M,Y] -> [C,M,Y]
thA = min(max(thA(:)', lbA), ubA);
fprintf('Partenza da calib_cgmy: C=%.5f M=%.3f Y=%.4f\n', thA(1), thA(2), thA(3));

ResA = struct('G',[],'loss',[],'theta',[],'ratio',[]);
thA_start = thA;
for dir = [0 1]  % 0 = da G0 verso il basso, 1 = da G0 verso l'alto
    if dir == 0, order = iG0:-1:1; else, order = iG0+1:numel(G_grid); end
    th = thA_start;
    for k = order
        Gk = G_grid(k);
        rA = @(p) sqrt_w .* (PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, ...
                R_surf, Q_surf, 'CGMY', [p(1), Gk, p(2), p(3)], [], [], ...
                XC_A, 0, POW_A, CP_flag) - Mid_surf);
        th = lsqnonlin(rA, th, lbA, ubA, opt_lsq);
        rG = @(p) sqrt_w .* (PrezzaSuperficieSINC_yield_new(F_surf, T_surf, K_surf, ...
                R_surf, Q_surf, 'CGMY', [p(1), Gk, p(2), p(3)], [], [], ...
                XC_A_GUARD, 0, POW_A_GUARD, CP_flag) - Mid_surf);
        ResA.G(k)       = Gk;
        ResA.loss(k)    = sum(rA(th).^2);
        ResA.theta(k,:) = th;
        ResA.ratio(k)   = ResA.loss(k) / sum(rG(th).^2);
        fprintf('  G=%7.3f  loss=%.6e  C=%.5f M=%8.3f Y=%.4f  Xc_ratio=%.4f\n', ...
            Gk, ResA.loss(k), th(1), th(2), th(3), ResA.ratio(k));
    end
end

% PROFILE B - ROUGH HESTON + CGMY  (liberi: H, nu, rho, C, M, Y)
fprintf('\n=== PROFILE B: rough Heston + CGMY, liberi H nu rho C M Y  (Xc=%d) ===\n', XC);
lbB = [0.005, 0.05, -0.99, 1e-4, 1.01, 0.9];
ubB = [0.49,  5.00,  0.10, 10.0, 200.0, 1.95];

fB  = fullfile(resdir, 'calib_rhcgmy_7x15z.mat');
assert(isfile(fB), 'Manca %s: lancia prima Rough_Heston_CMGY_def.', fB);
cB  = load(fB, 'theta_ottimi');
thB = cB.theta_ottimi([1 2 3 4 6 7]);  % [H,nu,rho,C,G,M,Y] -> togli G
thB = min(max(thB(:)', lbB), ubB);
fprintf('Partenza da calib_rhcgmy: H=%.4f nu=%.4f rho=%+.4f C=%.2e M=%.3f Y=%.4f\n', ...
    thB(1), thB(2), thB(3), thB(4), thB(5), thB(6));

ResB = struct('G',[],'loss',[],'theta',[],'ratio',[],'VSL',[],'quota',[]);
thB_start = thB;
for dir = [0 1]
    if dir == 0, order = iG0:-1:1; else, order = iG0+1:numel(G_grid); end
    th = thB_start;
    for k = order
        Gk = G_grid(k);
        rF = @(p) sqrt_w .* (PrezziRoughHestonCGMY(F_surf, T_surf, K_surf, p, Gk, ...
                T_uniq, Kvar_mkt, xi0_eval, N_RIC_F, XC, POW, CP_flag) - Mid_surf);
        rE = @(p) sqrt_w .* (PrezziRoughHestonCGMY(F_surf, T_surf, K_surf, p, Gk, ...
                T_uniq, Kvar_mkt, xi0_eval, N_RIC_E, XC, POW, CP_flag) - Mid_surf);
        rG = @(p) sqrt_w .* (PrezziRoughHestonCGMY(F_surf, T_surf, K_surf, p, Gk, ...
                T_uniq, Kvar_mkt, xi0_eval, N_RIC_E, XC_GUARD, POW_GUARD, CP_flag) - Mid_surf);
        th = lsqnonlin(rF, th, lbB, ubB, opt_lsq);
        th = lsqnonlin(rE, th, lbB, ubB, opt_lsq);

        [Lcg, Vcg] = ContributiCGMY(th(4), Gk, th(5), th(6));
        [~, i1y] = min(abs(T_uniq - 1));
        ti1  = ((0:N_RIC_E-1) + 0.5) / N_RIC_E * T_uniq(i1y);
        xie1 = max(xi0_eval(ti1) - Lcg, 1e-6);

        ResB.G(k)       = Gk;
        ResB.loss(k)    = sum(rE(th).^2);
        ResB.theta(k,:) = th;
        ResB.ratio(k)   = ResB.loss(k) / sum(rG(th).^2);
        ResB.VSL(k)     = Vcg / max(Lcg, 1e-12);
        ResB.quota(k)   = 100 * Vcg / (mean(xie1) + Vcg);
        fprintf(['  G=%7.3f  loss=%.6e  H=%.4f nu=%.4f rho=%+.4f  ' ...
                 'C=%.2e M=%8.3f Y=%.4f  VS/L=%.3f  quota=%4.1f%%  Xc_ratio=%.4f\n'], ...
            Gk, ResB.loss(k), th(1), th(2), th(3), th(4), th(5), th(6), ...
            ResB.VSL(k), ResB.quota(k), ResB.ratio(k));
    end
end

% LETTURA
fprintf('\n=== PROFILE NORMALIZZATI (loss / loss_min) ===\n');
fprintf('%-9s  %-14s  %-12s  %-14s  %-12s\n', ...
    'G', 'CGMY puro', 'Xc_ratio A', 'rH+CGMY', 'Xc_ratio B');
for k = 1:numel(G_grid)
    fprintf('%-9.3f  %-14.4f  %-12.4f  %-14.4f  %-12.4f\n', G_grid(k), ...
        ResA.loss(k)/min(ResA.loss), ResA.ratio(k), ...
        ResB.loss(k)/min(ResB.loss), ResB.ratio(k));
end

[~, kA] = min(ResA.loss);
[~, kB] = min(ResB.loss);
spreadA = max(ResA.loss)/min(ResA.loss) - 1;
spreadB = max(ResB.loss)/min(ResB.loss) - 1;
fprintf('\nCGMY puro : minimo a G = %.3f, escursione del profile = %.2f%%\n', ...
    G_grid(kA), 100*spreadA);
fprintf('rH+CGMY   : minimo a G = %.3f, escursione del profile = %.2f%%\n', ...
    G_grid(kB), 100*spreadB);
fprintf('Loss a G = 0.05 (valore fissato): CGMY %.6e | rH+CGMY %.6e\n', ...
    ResA.loss(iG0), ResB.loss(iG0));
fprintf('Guadagno di RMSE passando da G=0.05 al minimo del profile (rH+CGMY): %.2f%%\n', ...
    100*(1 - sqrt(ResB.loss(kB)/ResB.loss(iG0))));

% Ampiezza del tratto indistinguibile: tutti i G entro l'1% dal minimo.
flatB = ResB.loss <= 1.01*min(ResB.loss);
fprintf('G entro l''1%% dal minimo (rH+CGMY): %s\n', ...
    strjoin(compose('%.3g', G_grid(flatB)), ', '));
fprintf('  H su quei punti: [%.4f, %.4f]  (fattore %.2f)\n', ...
    min(ResB.theta(flatB,1)), max(ResB.theta(flatB,1)), ...
    max(ResB.theta(flatB,1))/min(ResB.theta(flatB,1)));
% L'insieme non e' un intervallo.
mainB = flatB & (G_grid <= 0.5);
fprintf('  solo bacino principale (G<=0.5): H in [%.4f, %.4f]  (fattore %.2f)\n', ...
    min(ResB.theta(mainB,1)), max(ResB.theta(mainB,1)), ...
    max(ResB.theta(mainB,1))/min(ResB.theta(mainB,1)));
fprintf(['\nLettura. Se l''escursione e'' di pochi punti percentuali su tre ordini\n' ...
         'di grandezza di G, il dato non distingue fra quei valori: G non e''\n' ...
         'identificato, fissarlo e'' legittimo e il valore va dichiarato come\n' ...
         'convenzione. Se invece c''e'' un minimo netto, G va liberato e la\n' ...
         'calibrazione a 7 parametri va rifatta. In quel caso verificare che il\n' ...
         'guadagno non sia sovradattamento: su 105 punti un settimo parametro\n' ...
         'deve pagare piu'' della penalita'' di AIC/BIC, oppure reggere un test\n' ...
         'out-of-sample sulle scadenze escluse dalla griglia ridotta.\n' ...
         'Se Xc_ratio si allontana da 1 in qualche punto, li'' il pricer sta\n' ...
         'troncando: quel tratto del profile non e'' affidabile.\n']);

if any(abs(ResA.ratio - 1) > 0.02) || any(abs(ResB.ratio - 1) > 0.02)
    warning(['Troncamento: Xc_ratio fuori dall''1%% in %d punti del profile A ' ...
             'e %d del profile B. Rialza Xc prima di interpretarli.'], ...
             sum(abs(ResA.ratio - 1) > 0.02), sum(abs(ResB.ratio - 1) > 0.02));
end

% FIGURE + SALVATAGGIO
figure('Name','Profile likelihood su G','Color','w','Position',[100 100 1000 430]);
subplot(1,2,1);
semilogx(ResA.G, ResA.loss/min(ResA.loss), 'o-','LineWidth',1.6); grid on; hold on;
xline(0.05,'k--','G fissato','HandleVisibility','off');
xlabel('G'); ylabel('loss / loss_{min}'); title('CGMY puro (C, M, Y riottimizzati)');
subplot(1,2,2);
semilogx(ResB.G, ResB.loss/min(ResB.loss), 'o-','LineWidth',1.6,'Color',[0.4 0.2 0.6]);
grid on; hold on;
xline(0.05,'k--','G fissato','HandleVisibility','off');
xlabel('G'); ylabel('loss / loss_{min}'); title('rough Heston + CGMY (6 par. riottimizzati)');

figure('Name','Parametri lungo il profile (rH+CGMY)','Color','w','Position',[150 150 1100 640]);
nomiB = {'H','nu','rho','C','M','Y'};
for j = 1:6
    subplot(2,3,j);
    if any(j == [4 5]), loglog(ResB.G, ResB.theta(:,j),'o-','LineWidth',1.4);
    else,               semilogx(ResB.G, ResB.theta(:,j),'o-','LineWidth',1.4); end
    grid on; xline(0.05,'k--','HandleVisibility','off');
    xlabel('G'); ylabel(nomiB{j}); title(nomiB{j});
end

if ~exist(resdir, 'dir'), mkdir(resdir); end
save(fullfile(resdir, 'profile_G.mat'), 'G_grid', 'ResA', 'ResB');
fprintf('\nProfile salvato in %s\n', fullfile(resdir, 'profile_G.mat'));

toc

% FUNZIONI LOCALI
function IV = InvertiIV_OTM(F_vec, K_vec, T_vec, P_fwd_vec)
    IV = nan(size(F_vec));
    for i = 1:length(F_vec)
        try
            if K_vec(i) <= F_vec(i)
                IV(i) = blkimpv(F_vec(i), K_vec(i), 0, T_vec(i), P_fwd_vec(i), 'Class','put');
            else
                C_fwd = max(P_fwd_vec(i) + F_vec(i) - K_vec(i), 1e-6);
                IV(i) = blkimpv(F_vec(i), K_vec(i), 0, T_vec(i), C_fwd, 'Class','call');
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

function [Lcg, Vcg] = ContributiCGMY(C, G, M, Y)
% Identica alla funzione omonima in Rough_Heston_CMGY_def.m e nella
% PricingSuite.
    g_m1 = gamma(-Y) * ((M-1)^Y - M^Y + (G+1)^Y - G^Y);
    mu1  = gamma(1-Y) * (M^(Y-1) - G^(Y-1));
    Lcg  = 2*C*(g_m1 - mu1);
    Vcg  = C * gamma(2-Y) * (M^(Y-2) + G^(Y-2));
end

function P_fwd = PrezziRoughHestonCGMY(F_vec, T_vec, K_vec, params6, G_fix, ...
                    T_uniq, Kvar_mkt, xi0_fn, N_riccati, Xc, pow, CP_vec)  %#ok<INUSD>
% params6 = [H, nu, rho, C, M, Y].
    H = params6(1); nu = params6(2); rho = params6(3);
    C = params6(4); M  = params6(5); Y   = params6(6);
    G = G_fix;
    params_combined = [H, nu, rho, 0, C, G, M, Y];

    L_cgmy = ContributiCGMY(C, G, M, Y);
    P_fwd  = zeros(length(K_vec),1);

    for k = 1:length(T_uniq)
        T   = T_uniq(k);
        idx = (T_vec == T);
        if ~any(idx), continue; end
        F_t  = F_vec(find(idx,1));
        K_t  = K_vec(idx);
        CP_t = CP_vec(idx);

        ti = ((0:N_riccati-1) + 0.5) / N_riccati * T;
        xi_grid = max(xi0_fn(ti) - L_cgmy, 1e-6);
        xi_grid = xi_grid(:)';

        P_fwd(idx) = PrezzaSuperficieSINC_yield(F_t, T*ones(size(K_t)), K_t, ...
            zeros(size(K_t)), zeros(size(K_t)), ...
            'rough Heston CGMY', params_combined, xi_grid, N_riccati, Xc, 0, pow, CP_t(1));
    end
end
