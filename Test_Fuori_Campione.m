clearvars; clc;
tic
% TEST_FUORI_CAMPIONE  Modelli calibrati contro dati non usati in stima.
% (A) Skew ATM di mercato e dei modelli su tutte le scadenze, anche sotto le
%     due settimane (CURVA_SKEW = 'calibrata' di default).
% (B) RMSE relativo di IV su tutta la superficie filtrata (0.03-1.6 anni).
% Output: Figure/fig_skew_mkt.pdf, Figure/tab_skew_mkt.tex,
% Figure/tab_fuori_campione.tex, Results/test_fuori_campione.mat.
% Tempo: 20-40 minuti.
FAI_SKEW       = true;  % parte (A)
FAI_SUPERFICIE = true;  % parte (B)
CURVA_SKEW     = 'calibrata';  % 'calibrata' (default) oppure 'mercato'

root = fileparts(mfilename('fullpath')); if isempty(root), root = pwd; end
addpath(fullfile(root,'Functions')); addpath(fullfile(root,'Models'));
addpath(fullfile(root,'Data'));      addpath(fullfile(root,'Results'));
resdir = fullfile(root,'Results');   figdir = fullfile(root,'Figure');
if ~isfolder(figdir), mkdir(figdir); end
NRIC = 400;  % nodi Riccati, come in calibrazione

% PARAMETRI CALIBRATI
fl = {'hes','calib_heston_7x15z.mat'; 'hcur','calib_hestoncurve_7x15z.mat'; ...
      'hck','calib_hestoncurvekappa_7x15z.mat'; 'rh','calib_rheston_7x15z_curva.mat'; ...
      'rk','calib_rhestonkappa_7x15z.mat'; 'cg','calib_cgmy_7x15z.mat'; ...
      'rhc','calib_rhcgmy_7x15z.mat'; 'ck','calib_rhcgmykappa_7x15z.mat'; ...
      'rhh','calib_rhh_7x15z.mat'};
P = struct();
for i = 1:size(fl,1)
    f = fullfile(resdir, fl{i,2});
    assert(isfile(f), 'Manca %s: serve il Run_All gia'' fatto.', fl{i,2});
    S = load(f, 'theta_ottimi'); P.(fl{i,1}) = S.theta_ottimi(:).';
end
t7 = P.rhc;  [Lcg7, ~] = ContributiCGMY(t7(4), t7(5), t7(6), t7(7));
t8 = P.ck;   [Lcg8, ~] = ContributiCGMY(t8(4), t8(5), t8(6), t8(7));  % [H nu rho C G M Y kappa]

% Elenco dei modelli. gruppo: 0 = senza mean reversion, 1 = con. curva: true
% se il modello prende xi0 dal mercato.
MOD = struct('nome',{}, 'gruppo',{}, 'curva',{}, 'dom',{});
MOD(end+1) = struct('nome','Heston',                   'gruppo',1,'curva',false,'dom','Heston');
MOD(end+1) = struct('nome','Heston-curve',             'gruppo',0,'curva',true, 'dom','rHeston');
MOD(end+1) = struct('nome','Heston-curve-\kappa',      'gruppo',1,'curva',true, 'dom','rHeston');
MOD(end+1) = struct('nome','rough Heston',             'gruppo',0,'curva',true, 'dom','rHeston');
MOD(end+1) = struct('nome','rough Heston-\kappa',      'gruppo',1,'curva',true, 'dom','rHeston');
MOD(end+1) = struct('nome','CGMY',                     'gruppo',0,'curva',false,'dom','CGMY');
MOD(end+1) = struct('nome','rough Heston + CGMY',      'gruppo',0,'curva',true, 'dom','rH+CGMY');
MOD(end+1) = struct('nome','rough Heston-\kappa + CGMY','gruppo',1,'curva',true, 'dom','rH+CGMY');
MOD(end+1) = struct('nome','rough Hawkes Heston',      'gruppo',1,'curva',false,'dom','rHawkesH');
NM = numel(MOD);
NOMETEX = strrep({MOD.nome}, '\kappa', '$\kappa$');

% DATI: stesso blocco degli script di calibrazione, ma senza il taglio
% T >= 0.03 (lo si applica dopo, solo dove serve)
fprintf('Caricamento dati...\n');
raw = readtable('20260515SPX_merged.csv','Delimiter',',','VariableNamingRule','preserve');
is_spot = strcmp(raw.Root,'SPOT');
S0_feed = raw.CF_LAST(is_spot);  raw = raw(~is_spot,:);
is_call = strcmp(raw.PutCall,'C');
kp = ~isnan(raw.STRIKE_PRC) & raw.CF_BID>0 & raw.CF_ASK>raw.CF_BID & raw.CF_BID>=0.375;
raw = raw(kp,:); is_call = is_call(kp);
K_all = raw.STRIKE_PRC; mid_all = 0.5*(raw.CF_BID+raw.CF_ASK);
T_all = (datenum(raw.EXPIR_DATE)-datenum('2026-05-15'))/365;
[~,ia] = unique([T_all,K_all,double(is_call)],'rows','stable');
T_all = T_all(ia); K_all = K_all(ia); mid_all = mid_all(ia); is_call = is_call(ia);

Tr = unique(T_all); Fp = nan(numel(Tr),1); Dp = nan(numel(Tr),1);
for i = 1:numel(Tr)
    Ti = Tr(i); mC = (T_all==Ti)&is_call; mP = (T_all==Ti)&~is_call;
    [Kc,iC,iP] = intersect(K_all(mC),K_all(mP));
    if numel(Kc) < 3, continue; end
    cI = find(mC); pI = find(mP); d = mid_all(cI(iC)) - mid_all(pI(iP));
    nr = abs(Kc/S0_feed-1) <= 0.25;
    if sum(nr) >= 5, Kr = Kc(nr); dr = d(nr); else, Kr = Kc; dr = d; end
    co = [ones(size(Kr)),Kr]\dr; Dp(i) = -co(2); Fp(i) = co(1)/Dp(i);
end
v = ~isnan(Fp) & Dp>0 & Dp<1.05 & Fp>0; Tr = Tr(v); Fp = Fp(v); Dp = Dp(v);

F_s=[]; T_s=[]; K_s=[]; Mid_s=[];
for i = 1:numel(Tr)
    Ti = Tr(i); Fi = Fp(i); Di = Dp(i); m = (T_all==Ti);
    Ki = K_all(m); mi = mid_all(m); ci = is_call(m); Pf = nan(size(Ki));
    for j = 1:numel(Ki)
        if Ki(j) < Fi && ~ci(j),     Pf(j) = mi(j)/Di;
        elseif Ki(j) >= Fi && ci(j), Pf(j) = mi(j)/Di - Fi + Ki(j); end
    end
    o = ~isnan(Pf) & Pf > 1e-4;
    F_s=[F_s;Fi*ones(sum(o),1)]; T_s=[T_s;Ti*ones(sum(o),1)];  %#ok<AGROW>
    K_s=[K_s;Ki(o)]; Mid_s=[Mid_s;Pf(o)];  %#ok<AGROW>
end
IV_M = InvertiIV_OTM(F_s,K_s,T_s,Mid_s);
lm = log(K_s./F_s);
id = ~isnan(IV_M) & IV_M>0.01 & IV_M<2 & Mid_s>1e-4 & abs(lm)<0.5 & ...
     T_s >= 2/365 & T_s <= 3.0;
F_s=F_s(id); T_s=T_s(id); K_s=K_s(id); Mid_s=Mid_s(id); IV_M=IV_M(id);
Ta = unique(T_s); ma = false(numel(T_s),1);
for i = 1:numel(Ta), m = T_s==Ta(i); if sum(m) >= 5, ma = ma | m; end, end
F_s=F_s(ma); T_s=T_s(ma); K_s=K_s(ma); Mid_s=Mid_s(ma); IV_M=IV_M(ma);
T_tutte = unique(T_s);
fprintf('Scadenze filtrate da %.0f giorni a %.2f anni: %d, quote: %d\n', ...
    365*min(T_tutte), max(T_tutte), numel(T_tutte), numel(T_s));

% Le 7 scadenze di calibrazione, scelte come negli script di calibrazione
Tc = T_tutte(T_tutte >= 0.03); T_cal = [];
for tt = [0.0384 0.0833 0.25 0.5 0.75 1.0 1.5]
    [~,ii] = min(abs(Tc-tt)); T_cal = [T_cal; Tc(ii)];  %#ok<AGROW>
end
T_cal = unique(T_cal);

% Demeterfi su ogni scadenza (stessa ricetta delle calibrazioni)
Kvar = zeros(numel(T_tutte),1); IVatm = zeros(numel(T_tutte),1);
for i = 1:numel(T_tutte)
    [Kvar(i), IVatm(i)] = Demeterfi(T_tutte(i), T_s, K_s, F_s, IV_M);
end
% curva calibrata: 7 ancore
[~, ic] = ismember(T_cal, T_tutte);
Wcal = pchip([0; T_cal], [0; T_cal.*Kvar(ic)]);
xi_cal = @(tv) eval_xi0_pchip(tv, Wcal);
% curva di mercato: tutte le scadenze fino a 1.6 anni, W reso crescente
Tm = T_tutte(T_tutte <= 1.6); Wm = Tm .* Kvar(T_tutte <= 1.6);
ok = true(size(Tm)); wmax = 0;
for i = 1:numel(Tm)
    if Wm(i) <= wmax, ok(i) = false; else, wmax = Wm(i); end
end
if any(~ok)
    fprintf('Curva di mercato: %d scadenze scartate perche'' W non cresce.\n', sum(~ok));
end
Wmkt = pchip([0; Tm(ok)], [0; Wm(ok)]);
xi_mkt = @(tv) eval_xi0_pchip(tv, Wmkt);
fprintf('Minimo di xi0 sulla curva di mercato: %.5f  (Lcg composto %.5f, composto-kappa %.5f)\n', ...
    min(xi_mkt(linspace(1e-4,1.6,2000))), Lcg7, Lcg8);

% (A) SKEW ATM
if FAI_SKEW
fprintf('\n=== (A) SKEW ATM: mercato contro modelli (curva %s) ===\n', CURVA_SKEW);
if strcmp(CURVA_SKEW,'mercato'), xi_sk = xi_mkt; else, xi_sk = xi_cal; end
Tsk = T_tutte(T_tutte <= 1.6);
ns  = numel(Tsk);
psiM = nan(ns,1); psiM_lo = nan(ns,1); psiM_hi = nan(ns,1); nq = zeros(ns,1); ivM = nan(ns,1);
for i = 1:ns
    m = (T_s == Tsk(i));
    k = log(K_s(m)./F_s(m)); iv = IV_M(m);
    [~, a] = min(abs(k)); sa = iv(a);
    kw = max(0.75*sa*sqrt(Tsk(i)), 0.006);
    [psiM(i), nq(i), ivM(i)] = SkewMercato(k, iv, kw);
    psiM_lo(i) = SkewMercato(k, iv, 0.6*kw);
    psiM_hi(i) = SkewMercato(k, iv, 1.4*kw);
end

F0 = 100;
psi = nan(ns, NM); conv = nan(ns, NM); ivA = nan(ns, NM);
for j = 1:NM
    fprintf('  %-28s', MOD(j).nome); tj = tic;
    for i = 1:ns
        T = Tsk(i);
        dk = min(0.005, 0.25*0.15*sqrt(T));
        Kpm = F0*exp([-dk; dk]);
        if T < 0.1
            cfgA = [2 14]; cfgB = [4 15];  % stesso passo, dominio doppio
        else
            [xc,pw] = XcPow(MOD(j).dom,'exact'); cfgA = [xc pw]; cfgB = [];
        end
        [psi(i,j), ivA(i,j)] = SkewModello(j, MOD, P, T, Kpm, dk, cfgA, xi_sk, Lcg7, Lcg8, NRIC);
        if ~isempty(cfgB) && (i == 1 || Tsk(i-1) >= 0.1 || mod(i,3) == 0)
            pB = SkewModello(j, MOD, P, T, Kpm, dk, cfgB, xi_sk, Lcg7, Lcg8, NRIC);
            conv(i,j) = abs(pB/psi(i,j) - 1);
        end
    end
    fprintf(' %6.1f s   (convergenza del dominio: max %.1e)\n', toc(tj), max(conv(:,j)));
end

% metriche
okM   = ~isnan(psiM) & psiM < 0;
corto = okM & Tsk <  0.0384 - 1e-6;
dentro= okM & Tsk >= 0.0384 - 1e-6;
fin25 = okM & Tsk <= 0.25;
cfM   = polyfit(log(Tsk(fin25)), log(abs(psiM(fin25))), 1);
cfMc  = polyfit(log(Tsk(corto)), log(abs(psiM(corto))), 1);
fprintf('\nMercato: %d scadenze con skew stimato (%d sotto le due settimane)\n', sum(okM), sum(corto));
fprintf('Mercato: pendenza log-log su T<=0.25 = %+.3f, su T<2 settimane = %+.3f\n', cfM(1), cfMc(1));

E_corto = nan(NM,1); E_dentro = nan(NM,1); SL = nan(NM,1); SLc = nan(NM,1); B_corto = nan(NM,1);
for j = 1:NM
    r = log(abs(psi(:,j))) - log(abs(psiM));
    E_corto(j)  = sqrt(mean(r(corto).^2));
    E_dentro(j) = sqrt(mean(r(dentro).^2));
    B_corto(j)  = mean(r(corto));
    cf = polyfit(log(Tsk(fin25)), log(abs(psi(fin25,j))), 1);  SL(j)  = cf(1);
    cf = polyfit(log(Tsk(corto)), log(abs(psi(corto,j))), 1);  SLc(j) = cf(1);
end
fprintf('\n%-28s %9s %9s %12s %12s %12s\n','Modello','pend<=3m','pend<2w', ...
    'logRMSE<2w','bias<2w','logRMSE 2w+');
fprintf('%-28s %+9.3f %+9.3f %12s %12s %12s\n','Mercato',cfM(1),cfMc(1),'-','-','-');
for j = 1:NM
    fprintf('%-28s %+9.3f %+9.3f %12.3f %+12.3f %12.3f\n', MOD(j).nome, SL(j), SLc(j), ...
        E_corto(j), B_corto(j), E_dentro(j));
end
fprintf(['Lettura: logRMSE e'' l''errore medio in log dello skew, 0.10 = circa 10%%.\n' ...
         'bias<2w > 0 vuol dire che il modello ha uno skew piu'' ripido del mercato\n' ...
         'sotto le due settimane, < 0 piu'' piatto.\n']);

fprintf('\nDettaglio per scadenza (|psi|):\n%-8s %8s', 'T', 'mercato');
fprintf(' %8s', 'H','Hc','Hck','rH','rHk','CGMY','rHC','rHkC','rHH'); fprintf('\n');
for i = 1:ns
    fprintf('%-8.4f %8.3f', Tsk(i), abs(psiM(i)));
    fprintf(' %8.3f', abs(psi(i,:))); fprintf('\n');
end
fprintf('\nIV ATM (%%) per scadenza, per controllare il livello:\n%-8s %8s', 'T', 'mercato');
fprintf(' %8s', 'H','Hc','Hck','rH','rHk','CGMY','rHC','rHkC','rHH'); fprintf('\n');
for i = 1:ns
    fprintf('%-8.4f %8.2f', Tsk(i), 100*ivM(i));
    fprintf(' %8.2f', 100*ivA(i,:)); fprintf('\n');
end

% figura
CC = { [0.00 0.45 0.74], [0.30 0.60 0.85], [0.10 0.60 0.55], [0.85 0.33 0.10], ...
       [0.64 0.08 0.18], [0.47 0.67 0.19], [0.49 0.18 0.56], [0.75 0.40 0.75], ...
       [0.93 0.69 0.13] };
f = figure('Color','w','Position',[80 80 1200 470]);
tit = {'Without mean reversion','With mean reversion'};
for g = 0:1
    subplot(1,2,g+1); hold on; box on; grid on
    set(gca,'XScale','log','YScale','log');
    yl = [0.5*min(abs(psiM(okM))) 2*max(abs(psiM(okM)))];
    patch([min(Tsk)*0.9 0.0384 0.0384 min(Tsk)*0.9], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.92 0.92 0.92], 'EdgeColor','none','HandleVisibility','off');
    e1 = abs(psiM) - min(abs(psiM_lo),abs(psiM_hi)); e2 = max(abs(psiM_lo),abs(psiM_hi)) - abs(psiM);
    e1(isnan(e1)) = 0; e2(isnan(e2)) = 0; e1 = max(e1,0); e2 = max(e2,0);
    errorbar(Tsk(okM), abs(psiM(okM)), e1(okM), e2(okM), 'ko', 'MarkerSize',4, ...
        'MarkerFaceColor','k', 'LineStyle','none', 'DisplayName','Market');
    for j = find([MOD.gruppo] == g)
        plot(Tsk, abs(psi(:,j)), '-', 'Color', CC{j}, 'LineWidth', 1.5, 'DisplayName', MOD(j).nome);
    end
    xlim([min(Tsk)*0.9 1.7]); ylim(yl);
    xlabel('$T$ (years)','Interpreter','latex');
    ylabel('$|\partial\sigma_{BS}/\partial k|$ at the money','Interpreter','latex');
    title(tit{g+1});
    legend('Location','southwest','Interpreter','tex','Box','off');
end
Esporta(f, figdir, 'fig_skew_mkt');

% tabella
fid = fopen(fullfile(figdir,'tab_skew_mkt.tex'),'w');
fprintf(fid,'%% Generato da Test_Fuori_Campione.m il %s\n', datestr(now));
fprintf(fid,'\\begin{tabular}{lrrrr}\n\\toprule\n');
fprintf(fid,' & \\multicolumn{2}{c}{Slope} & \\multicolumn{2}{c}{Error in $\\log|\\psi|$}\\\\\n');
fprintf(fid,'\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n');
fprintf(fid,'Model & $T\\le 3$m & $T<2$w & $T<2$w & $T\\ge 2$w\\\\\n\\midrule\n');
fprintf(fid,'Market & %.3f & %.3f & --- & ---\\\\\n\\addlinespace\n', cfM(1), cfMc(1));
for g = 0:1
    for j = find([MOD.gruppo] == g)
        fprintf(fid,'%s & %.3f & %.3f & %.3f & %.3f\\\\\n', NOMETEX{j}, SL(j), SLc(j), E_corto(j), E_dentro(j));
    end
    if g == 0, fprintf(fid,'\\addlinespace\n'); end
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);
fprintf('Scritti fig_skew_mkt.pdf e tab_skew_mkt.tex\n');
else
    Tsk=[]; psiM=[]; psi=[]; E_corto=[]; E_dentro=[]; SL=[]; SLc=[]; cfM=[]; cfMc=[]; B_corto=[]; conv=[]; nq=[];
end

% (B) ERRORE SU TUTTA LA SUPERFICIE FILTRATA
if FAI_SUPERFICIE
fprintf('\n=== (B) RMSE relativo di IV su tutta la superficie filtrata ===\n');
sel = T_s >= 0.03 & T_s <= 1.6;
Fb = F_s(sel); Tb = T_s(sel); Kb = K_s(sel); IVb = IV_M(sel);
Tu = unique(Tb); is_cal = ismember(Tb, T_cal);
fprintf('Quote: %d su %d scadenze (%d sulle 7 scadenze di calibrazione)\n', ...
    numel(Tb), numel(Tu), sum(is_cal));
IVmod = nan(numel(Tb), NM);
for j = 1:NM
    fprintf('  %-28s', MOD(j).nome); tj = tic;
    [xc,pw] = XcPow(MOD(j).dom,'exact');
    for i = 1:numel(Tu)
        m = (Tb == Tu(i));
        Fi = Fb(find(m,1)); Ki = Kb(m);
        C  = PrezziModello(j, MOD, P, Tu(i), Fi, Ki, [xc pw], xi_cal, Lcg7, Lcg8, NRIC);
        Pp = C(:) - Fi + Ki(:);  % put forward
        IVmod(m,j) = InvertiIV_OTM(Fi*ones(sum(m),1), Ki, Tu(i)*ones(sum(m),1), max(Pp,1e-8));
    end
    fprintf(' %6.1f s\n', toc(tj));
end
EPS = (IVb - IVmod) ./ IVb;
rm  = @(x) 100*sqrt(mean(x.^2,'omitnan'));
R_all = zeros(NM,1); R_cal = zeros(NM,1); R_oth = zeros(NM,1);
bande = [0.03 0.1; 0.1 0.5; 0.5 1.61]; R_b = zeros(NM,3);
for j = 1:NM
    R_all(j) = rm(EPS(:,j)); R_cal(j) = rm(EPS(is_cal,j)); R_oth(j) = rm(EPS(~is_cal,j));
    for b = 1:3
        mb = Tb >= bande(b,1) & Tb < bande(b,2); R_b(j,b) = rm(EPS(mb,j));
    end
end
fprintf('\n%-28s %9s %11s %11s %10s %10s %10s\n','Modello','tutte','scad. cal.', ...
    'altre scad.','T<0.1','0.1-0.5','0.5-1.6');
for j = 1:NM
    fprintf('%-28s %8.2f%% %10.2f%% %10.2f%% %9.2f%% %9.2f%% %9.2f%%\n', MOD(j).nome, ...
        R_all(j), R_cal(j), R_oth(j), R_b(j,:));
end
fid = fopen(fullfile(figdir,'tab_fuori_campione.tex'),'w');
fprintf(fid,'%% Generato da Test_Fuori_Campione.m il %s\n', datestr(now));
fprintf(fid,'\\begin{tabular}{lrrr}\n\\toprule\n');
fprintf(fid,'Model & All quotes & Calibration maturities & Other maturities\\\\\n\\midrule\n');
for g = 0:1
    for j = find([MOD.gruppo] == g)
        fprintf(fid,'%s & %.2f & %.2f & %.2f\\\\\n', NOMETEX{j}, R_all(j), R_cal(j), R_oth(j));
    end
    if g == 0, fprintf(fid,'\\addlinespace\n'); end
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);
fprintf('Quote: %d in tutto, %d sulle scadenze di calibrazione, %d sulle altre\n', ...
    numel(Tb), sum(is_cal), sum(~is_cal));
fprintf('Scritto tab_fuori_campione.tex\n');
else
    R_all=[]; R_cal=[]; R_oth=[]; R_b=[]; Tb=[];
end

if ~exist('ivM','var'), ivM = []; ivA = []; end
save(fullfile(resdir,'test_fuori_campione.mat'), 'CURVA_SKEW','ivM','ivA','Tsk','psiM','psi','E_corto','E_dentro', ...
    'B_corto','SL','SLc','cfM','cfMc','conv','nq','R_all','R_cal','R_oth','R_b','T_cal','MOD');
fprintf('\nSalvato Results/test_fuori_campione.mat\n');
toc

% FUNZIONI LOCALI
function C = PrezziModello(j, MOD, P, T, F0, K, cfg, xifn, Lcg7, Lcg8, NRIC)
% Call forward (r = q = 0) del modello j a scadenza T, sugli strike K.
    K = K(:); n = numel(K); Z = zeros(n,1); Tv = T*ones(n,1);
    xc = cfg(1); pw = cfg(2);
    xg = RiccatiXi(xifn, NRIC, T);
    switch MOD(j).nome
        case 'Heston'
            C = PrezzaSuperficieSINC_yield_new(F0, Tv, K, Z, Z, 'Heston', P.hes, [], [], xc, 0, pw, 1);
        case 'Heston-curve'
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, 'Heston curve', P.hcur, xg, NRIC, xc, 0, pw, 1);
        case 'Heston-curve-\kappa'
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, 'Heston curve kappa', P.hck, xg, NRIC, xc, 0, pw, 1);
        case 'rough Heston'
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, 'rough Heston', [P.rh 0], xg, NRIC, xc, 0, pw, 1);
        case 'rough Heston-\kappa'
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, 'rough Heston kappa', P.rk, xg, NRIC, xc, 0, pw, 1);
        case 'CGMY'
            C = PrezzaSuperficieSINC_yield_new(F0, Tv, K, Z, Z, 'CGMY', P.cg, [], [], xc, 0, pw, 1);
        case 'rough Heston + CGMY'
            t = P.rhc;
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, 'rough Heston CGMY', [t(1:3) 0 t(4:7)], ...
                max(xg - Lcg7, 1e-6), NRIC, xc, 0, pw, 1);
        case 'rough Heston-\kappa + CGMY'
            t = P.ck;
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, 'rough Heston kappa CGMY', [t(1:3) t(8) t(4:7)], ...
                max(xg - Lcg8, 1e-6), NRIC, xc, 0, pw, 1);
        case 'rough Hawkes Heston'
            C = PrezzaSuperficieSINC_yield_new(F0, Tv, K, Z, Z, 'rough Hawkes Heston', P.rhh, [], NRIC, xc, 0, pw, 1);
        otherwise
            error('Modello non previsto: %s', MOD(j).nome);
    end
    C = C(:);
end

function [s, ivatm] = SkewModello(j, MOD, P, T, Kpm, dk, cfg, xifn, Lcg7, Lcg8, NRIC)
% Skew ATM del modello: differenza centrata della IV in k = log(K/F).
    F0 = 100;
    C  = PrezziModello(j, MOD, P, T, F0, Kpm, cfg, xifn, Lcg7, Lcg8, NRIC);
    Pp = C - F0 + Kpm(:);
    iv = InvertiIV_OTM(F0*[1;1], Kpm(:), T*[1;1], max(Pp,1e-10));
    s  = (iv(2) - iv(1)) / (2*dk);
    ivatm = 0.5*(iv(1) + iv(2));
end

function [s, n, ivatm] = SkewMercato(k, iv, kw)
% Skew ATM di mercato: pendenza in k = 0 di una parabola su |k| <= kw.
    m = abs(k) <= kw; n = sum(m); ivatm = NaN;
    if n < 5 || sum(k(m) < 0) < 2 || sum(k(m) > 0) < 2, s = NaN; return; end
    p = polyfit(k(m), iv(m), 2);
    s = p(2); ivatm = p(3);
end

function [Kv, ivatm] = Demeterfi(Ti, T_s, K_s, F_s, IV_M)
    m = (T_s == Ti); Fi = F_s(find(m,1));
    [Ks, s] = sort(K_s(m)); IVs = IV_M(m); IVs = IVs(s);
    [~, a] = min(abs(Ks - Fi)); ivatm = IVs(a);
    Kg  = logspace(log10(0.05*Fi), log10(2.50*Fi), 4000)';
    IVg = interp1(Ks, IVs, Kg, 'linear', NaN);
    IVg(Kg < Ks(1)) = IVs(1); IVg(Kg > Ks(end)) = IVs(end);
    O = zeros(size(Kg));
    for j = 1:numel(Kg)
        if Kg(j) < Fi, [~,p] = blkprice(Fi,Kg(j),0,Ti,IVg(j)); O(j) = p;
        else,          [c,~] = blkprice(Fi,Kg(j),0,Ti,IVg(j)); O(j) = c; end
    end
    Kv = (2/Ti) * trapz(Kg, O./Kg.^2);
end

function xg = RiccatiXi(xifn, N, T)
    ti = ((0:N-1)+0.5)/N*T;  xg = xifn(ti);  xg = xg(:)';
end

function xi = eval_xi0_pchip(tvec, W_pp)
    h = 1e-5; tvec = tvec(:);
    tp = max(tvec+h,0); tm = max(tvec-h,0);
    dn = tp - tm; dn(dn < 1e-12) = 1e-12;
    xi = max((ppval(W_pp,tp) - ppval(W_pp,tm))./dn, 1e-8);
end

function [Lcg, Vcg] = ContributiCGMY(C, G, M, Y)
    Lcg = 2*C*( gamma(-Y)*((M-1)^Y - M^Y + (G+1)^Y - G^Y) ...
              - gamma(1-Y)*(M^(Y-1) - G^(Y-1)) );
    Vcg = C * gamma(2-Y) * (M^(Y-2) + G^(Y-2));
end

function IV = InvertiIV_OTM(F,K,T,P)
    IV = nan(size(F));
    for i = 1:numel(F)
        try
            if K(i) <= F(i)
                IV(i) = blkimpv(F(i),K(i),0,T(i),P(i),'Class','put');
            else
                IV(i) = blkimpv(F(i),K(i),0,T(i),max(P(i)+F(i)-K(i),1e-8),'Class','call');
            end
        catch
            IV(i) = NaN;
        end
    end
end

function Esporta(f, cartella, nome)
    set(findall(f,'-property','FontSize'),'FontSize',9);
    try
        exportgraphics(f, fullfile(cartella,[nome '.pdf']),'ContentType','vector');
    catch
        set(f,'PaperPositionMode','auto');
        print(f, fullfile(cartella,[nome '.pdf']), '-dpdf', '-painters', '-bestfit');
    end
    fprintf('  %s.pdf\n', nome);
    close(f);
end
