clearvars; clc;
tic
% ESTENSIONE  Calibrazione con le scadenze sotto le due settimane (Sez. 7.8).
% Campione: le 7 scadenze della tesi piu' 3 scadenze corte (5, 7, 11 gg).
% Ricalibra Heston-curve-kappa, rough Heston, rough Heston-kappa e rough
% Hawkes Heston, poi fa il profilo in H con kappa libero.
% Output: Results/estensione.mat, Figure/fig_estensione.pdf,
%         Figure/tab_estensione.tex. Tempo: circa 1.5 ore.
FAI_RHH        = true;  % il rough Hawkes e' il piu' lento
FAI_PROFILO    = true;  % profilo in H con kappa libero
GIORNI_CORTI   = [5 7 11];  % scadenze corte da aggiungere
SOLO_PROFILO   = false;  % true: rilegge le calibrazioni da Results/estensione.mat
                                % e rifa' solo il profilo in H (circa 20-30
                                % minuti)

root = fileparts(mfilename('fullpath')); if isempty(root), root = pwd; end
addpath(fullfile(root,'Functions')); addpath(fullfile(root,'Data')); addpath(fullfile(root,'Results'));
resdir = fullfile(root,'Results'); figdir = fullfile(root,'Figure');
if ~isfolder(figdir), mkdir(figdir); end

% parametri della tesi (punti di partenza e confronto)
P0 = struct();
P0.hck = CaricaTheta(resdir,'calib_hestoncurvekappa_7x15z.mat', [4.9242 1.4828 -0.6857]);
P0.rh  = CaricaTheta(resdir,'calib_rheston_7x15z_curva.mat',    [0.086568 0.33987 -0.67611]);
P0.rhk = CaricaTheta(resdir,'calib_rhestonkappa_7x15z.mat',     [0.4176 1.1262 -0.6828 3.3175]);
P0.rhh = CaricaTheta(resdir,'calib_rhh_7x15z.mat',              [0.8568 -1.6171 0.45765 -0.63089 0.65392 0.026375 0.10894]);

% DATI (blocco standard), senza il taglio T >= 0.03
fprintf('Caricamento dati...\n');
raw = readtable('20260515SPX_merged.csv','Delimiter',',','VariableNamingRule','preserve');
is_spot = strcmp(raw.Root,'SPOT'); S0_feed = raw.CF_LAST(is_spot); raw = raw(~is_spot,:);
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
id = ~isnan(IV_M) & IV_M>0.01 & IV_M<2 & Mid_s>1e-4 & abs(lm)<0.5 & T_s>=2/365 & T_s<=3.0;
F_s=F_s(id); T_s=T_s(id); K_s=K_s(id); Mid_s=Mid_s(id); IV_M=IV_M(id);
Ta = unique(T_s); ma = false(numel(T_s),1);
for i = 1:numel(Ta), m = T_s==Ta(i); if sum(m) >= 5, ma = ma | m; end, end
F_s=F_s(ma); T_s=T_s(ma); K_s=K_s(ma); Mid_s=Mid_s(ma); IV_M=IV_M(ma);
Ta = unique(T_s);

% scadenze: 7 della tesi + corte
Tc = Ta(Ta >= 0.03); T_std = [];
for tt = [0.0384 0.0833 0.25 0.5 0.75 1.0 1.5]
    [~,ii] = min(abs(Tc-tt)); T_std = [T_std; Tc(ii)];  %#ok<AGROW>
end
T_std = unique(T_std);
Tb = Ta(Ta < T_std(1) - 1e-6); T_cor = [];  % strettamente sotto la prima scadenza della tesi
for gg = GIORNI_CORTI
    [~,ii] = min(abs(Tb - gg/365)); T_cor = [T_cor; Tb(ii)];  %#ok<AGROW>
end
T_cor = unique(T_cor);
T_uniq = sort([T_cor; T_std]);
fprintf('Scadenze corte aggiunte: %s giorni\n', mat2str(round(365*T_cor')));
fprintf('Scadenze del campione esteso: %s\n', mat2str(T_uniq', 4));

% Demeterfi e curva xi0 su tutte le 10 ancore
Kvar = zeros(numel(T_uniq),1); IVatm = zeros(numel(T_uniq),1);
for i = 1:numel(T_uniq), [Kvar(i), IVatm(i)] = Demeterfi(T_uniq(i), T_s, K_s, F_s, IV_M); end
W = T_uniq.*Kvar;
if any(diff([0; W]) <= 0), warning('W(T) non crescente sulle ancore: xi0 potrebbe annullarsi.'); end
W_pp = pchip([0; T_uniq], [0; W]);
fprintf('\n%-8s %-6s %-10s %-10s\n','T','giorni','sqrt(Kvar)','IV ATM');
for i = 1:numel(T_uniq)
    fprintf('%-8.4f %-6.0f %-10.2f %-10.2f\n', T_uniq(i), 365*T_uniq(i), 100*sqrt(Kvar(i)), 100*IVatm(i));
end

% 15 strike in z per scadenza e pesi (come nella tesi)
keep = false(size(T_s));
for i = 1:numel(T_uniq)
    m = find(T_s == T_uniq(i));
    z = log(K_s(m)./F_s(m)) ./ (IVatm(i)*sqrt(T_uniq(i)));
    if numel(m) <= 15, keep(m) = true; continue; end
    for tt = linspace(max(-3.0,min(z)), min(2.0,max(z)), 15)
        [~,jj] = min(abs(z-tt)); keep(m(jj)) = true;
    end
end
F = F_s(keep); T = T_s(keep); K = K_s(keep); Mid = Mid_s(keep); IVm = IV_M(keep);
is_cor = ismember(T, T_cor);
d1 = (log(F./K) + 0.5*IVm.^2.*T) ./ (IVm.*sqrt(T));
vega = max(F.*sqrt(T).*exp(-0.5*d1.^2)/sqrt(2*pi), 1e-4);
sw = 1./(vega.*max(IVm,0.01)); sw = min(sw, prctile(sw,99));
fprintf('\nCampione esteso: %d quote (%d sulle scadenze corte)\n', numel(F), sum(is_cor));

Pr = @(nome,p,stadio) PrezziPut(nome, p, F, T, K, T_uniq, W_pp, stadio);
Res = @(nome,p,stadio) sw .* (Pr(nome,p,stadio) - Mid);
Met = @(nome,p) Metriche(Pr(nome,p,'exact'), F, T, K, IVm, is_cor);

% CALIBRAZIONI
MOD = {'hck','rh','rhk'}; if FAI_RHH, MOD{end+1} = 'rhh'; end
NOME = struct('hck','Heston-curve-kappa','rh','rough Heston','rhk','rough Heston-kappa','rhh','rough Hawkes Heston');
LB = struct('hck',[0 0.05 -0.99], 'rh',[0.005 0.05 -0.99], 'rhk',[0.005 0.05 -0.99 0], ...
            'rhh',[0.505 -20 1e-4 -0.99 0 1e-5 1e-4]);
UB = struct('hck',[50 5 0.10], 'rh',[0.49 5 0.10], 'rhk',[0.50 5 0.10 50], ...
            'rhh',[0.995 -0.01 2 0.10 5 0.10 0.30]);
ST = struct();
ST.hck = [P0.hck; 2 1.0 -0.70; 10 2.2 -0.68];
ST.rh  = [P0.rh; 0.03 0.30 -0.70; 0.20 0.50 -0.68];
ST.rhk = [P0.rhk; 0.0866 0.34 -0.676 0.05; 0.20 0.55 -0.68 0.85];
ST.rhh = P0.rhh;
opF = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','Display','off', ...
      'MaxFunctionEvaluations',600,'StepTolerance',1e-8,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);
opE = optimoptions(opF,'MaxFunctionEvaluations',200,'MaxIterations',25);

if SOLO_PROFILO
    OLD = load(fullfile(resdir,'estensione.mat'));
    R = OLD.R; MOD = OLD.MOD;
    fprintf('Calibrazioni rilette da Results/estensione.mat\n');
else
R = struct();
for im = 1:numel(MOD)
    nm = MOD{im}; tm = tic;
    fprintf('\n=== %s ===\n', NOME.(nm));
    M0 = Met(nm, P0.(nm));
    fprintf('Parametri della tesi sul campione esteso: RMSE %.3f%% (corte %.3f%%, standard %.3f%%)\n', ...
        M0.tot, M0.cor, M0.std);
    best = []; bl = Inf;
    for s = 1:size(ST.(nm),1)
        x0 = min(max(ST.(nm)(s,:), LB.(nm)+1e-6), UB.(nm)-1e-6);
        [x, l] = lsqnonlin(@(p) Res(nm,p,'fast'), x0, LB.(nm), UB.(nm), opF);
        fprintf('  partenza %d: loss fast %.5e  %s\n', s, l, mat2str(x,4));
        if l < bl, bl = l; best = x; end
    end
    [th, le] = lsqnonlin(@(p) Res(nm,p,'exact'), best, LB.(nm), UB.(nm), opE);
    M1 = Met(nm, th);
    fprintf('Ricalibrato: %s\n', mat2str(th,4));
    fprintf('RMSE %.3f%% (corte %.3f%%, standard %.3f%%)   loss exact %.5e   [%.0f s]\n', ...
        M1.tot, M1.cor, M1.std, le, toc(tm));
    bnd = (th <= LB.(nm) + 1e-3*(UB.(nm)-LB.(nm))) | (th >= UB.(nm) - 1e-3*(UB.(nm)-LB.(nm)));
    if any(bnd), fprintf('  ATTENZIONE: parametri al bound: %s\n', mat2str(find(bnd))); end
    R.(nm) = struct('theta0',P0.(nm),'theta',th,'M0',M0,'M1',M1,'loss',le,'bound',bnd);
end
end

% PROFILO IN H CON KAPPA LIBERO (campione esteso)
Hg = [0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.50]; PRO = nan(numel(Hg), 5);
if FAI_PROFILO
    fprintf('\n=== Profilo in H, kappa libero, campione esteso ===\n');
    % Si parte dall'ottimo e ci si allontana un nodo alla volta, nelle due
    % direzioni.
    lbq = [0.05 -0.99 0]; ubq = [5 0.10 50];
    cl  = @(q) min(max(q, lbq+1e-6), ubq-1e-6);
    H0  = R.rhk.theta(1);
    giu = sort(Hg(Hg <  H0), 'descend');
    su  = sort(Hg(Hg >= H0), 'ascend');
    for ramo = {giu, su}
        xprev = R.rhk.theta(2:4);
        for H = ramo{1}
            rr = @(q) Res('rhk', [H q], 'fast');
            starts = [xprev; R.rh.theta(2:3) 0.05; R.rhk.theta(2:4)];
            bl = Inf; xb = xprev;
            for s0 = 1:size(starts,1)
                [xq, lq] = lsqnonlin(rr, cl(starts(s0,:)), lbq, ubq, opF);
                if lq < bl, bl = lq; xb = xq; end
            end
            xprev = xb;
            Mh = Met('rhk', [H xb]);
            ih = find(Hg == H, 1);
            PRO(ih,:) = [H xb Mh.tot];
            fprintf('H=%.2f  nu=%.4f  rho=%.4f  kappa=%.4f  RMSE=%.3f%%  (corte %.3f%%)\n', H, xb, Mh.tot, Mh.cor);
        end
    end
end

% SKEW SOTTO LE DUE SETTIMANE: mercato contro modelli ricalibrati
if SOLO_PROFILO
    Tk = OLD.Tk; psiM = OLD.psiM; PSI = OLD.PSI; SL = OLD.SL; RA = OLD.RA; LE = OLD.LE; cf = OLD.cf;
    ok = ~isnan(psiM) & psiM < 0;
else
fprintf('\n=== Skew ATM sotto le due settimane (modelli ricalibrati) ===\n');
Tk = Ta(Ta <= 0.0385 & Ta >= 2/365);
psiM = nan(numel(Tk),1);
for i = 1:numel(Tk)
    m = (T_s == Tk(i)); k = log(K_s(m)./F_s(m)); iv = IV_M(m);
    [~,a] = min(abs(k)); psiM(i) = SkewMercato(k, iv, max(0.75*iv(a)*sqrt(Tk(i)),0.006));
end
ok = ~isnan(psiM) & psiM < 0;
cf = polyfit(log(Tk(ok)), log(abs(psiM(ok))), 1);
fprintf('Mercato: pendenza %.3f  rapporto psi(%.0fg)/psi(%.0fg) = %.2f\n', cf(1), ...
    365*Tk(find(ok,1)), 365*Tk(find(ok,1,'last')), psiM(find(ok,1))/psiM(find(ok,1,'last')));
PSI = nan(numel(Tk), numel(MOD)); SL = nan(1,numel(MOD)); RA = nan(1,numel(MOD)); LE = nan(1,numel(MOD));
for im = 1:numel(MOD)
    nm = MOD{im};
    for i = 1:numel(Tk)
        PSI(i,im) = SkewModelloPut(nm, R.(nm).theta, Tk(i), W_pp);
    end
    c = polyfit(log(Tk(ok)), log(abs(PSI(ok,im))), 1); SL(im) = c(1);
    RA(im) = PSI(find(ok,1),im)/PSI(find(ok,1,'last'),im);
    LE(im) = sqrt(mean((log(abs(PSI(ok,im))) - log(abs(psiM(ok)))).^2));
    fprintf('%-22s pendenza %.3f  rapporto %.2f  errore log %.3f\n', NOME.(nm), SL(im), RA(im), LE(im));
end

end

% RIEPILOGO, TABELLA, FIGURA
fprintf('\n================ RIEPILOGO ================\n');
fprintf('%-22s %-34s %-34s\n','Modello','tesi (7 scadenze)','esteso (con le corte)');
for im = 1:numel(MOD)
    nm = MOD{im}; r = R.(nm);
    fprintf('%-22s %-34s %-34s\n', NOME.(nm), mat2str(r.theta0,3), mat2str(r.theta,3));
    fprintf('%-22s RMSE campione esteso: %.3f%% -> %.3f%%   (corte %.3f%% -> %.3f%%)\n', '', ...
        r.M0.tot, r.M1.tot, r.M0.cor, r.M1.cor);
end
fprintf('Riferimento Cap. 6 (7 scadenze): rough Heston-kappa H = %.3f, kappa = %.3f\n', P0.rhk(1), P0.rhk(4));
fprintf('Esteso:                         rough Heston-kappa H = %.3f, kappa = %.3f\n', R.rhk.theta(1), R.rhk.theta(4));

fid = fopen(fullfile(figdir,'tab_estensione.tex'),'w');
fprintf(fid,'%% Generato da Estensione.m il %s\n', datestr(now));
fprintf(fid,'\\begin{tabular}{lrrrrrr}\n\\toprule\n');
fprintf(fid,' & & & \\multicolumn{2}{c}{RMSE (\\%%)} & \\multicolumn{2}{c}{Skew $<2$w}\\\\\n');
fprintf(fid,'\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n');
fprintf(fid,'Model & $H$ & $\\kappa$ & All & Short & Slope & Ratio\\\\\n\\midrule\n');
fprintf(fid,'Market & & & & & %.3f & %.2f\\\\\n\\addlinespace\n', cf(1), psiM(find(ok,1))/psiM(find(ok,1,'last')));
for im = 1:numel(MOD)
    nm = MOD{im}; th = R.(nm).theta;
    switch nm
        case 'hck', Hs = '0.5'; ks = sprintf('%.3f', th(1));
        case 'rh',  Hs = sprintf('%.3f', th(1)); ks = '0';
        case 'rhk', Hs = sprintf('%.3f', th(1)); ks = sprintf('%.3f', th(4));
        case 'rhh', Hs = sprintf('%.3f', th(1)-0.5); ks = sprintf('%.3f', -th(2));
    end
    fprintf(fid,'%s & %s & %s & %.2f & %.2f & %.3f & %.2f\\\\\n', strrep(NOME.(nm),'-kappa','-$\kappa$'), ...
        Hs, ks, R.(nm).M1.tot, R.(nm).M1.cor, SL(im), RA(im));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

f = figure('Color','w','Units','centimeters','Position',[2 2 16 6.5]);
subplot(1,2,1); hold on; grid on; box on
H6 = [0.05 0.10 0.20 0.30 0.40 0.50]; R6 = [4.034 3.859 3.553 3.334 3.240 3.294];  % Tab. 6.4
plot(H6, R6, 'o--', 'Color',[0.5 0.5 0.5], 'LineWidth',1.2, 'DisplayName','7 maturities (Table 6.4)');
if FAI_PROFILO
    plot(PRO(:,1), PRO(:,5), 's-', 'Color',[0.64 0.08 0.18], 'LineWidth',1.4, 'DisplayName','with short maturities');
end
xlabel('$H$','Interpreter','latex'); ylabel('Relative RMSE (\%)','Interpreter','latex');
title('Profile in $H$, $\kappa$ free','Interpreter','latex','FontWeight','normal');
legend('Location','best','Box','off','FontSize',7);
subplot(1,2,2); hold on; grid on; box on; set(gca,'XScale','log','YScale','log');
plot(Tk(ok), abs(psiM(ok)), 'ko', 'MarkerFaceColor','k', 'MarkerSize',3, 'DisplayName','Market');
CC = {[0.10 0.60 0.55],[0.85 0.33 0.10],[0.64 0.08 0.18],[0.93 0.69 0.13]};
for im = 1:numel(MOD)
    plot(Tk, abs(PSI(:,im)), '-', 'Color', CC{im}, 'LineWidth',1.3, 'DisplayName', strrep(NOME.(MOD{im}),'kappa','\kappa'));
end
xlabel('$T$ (years)','Interpreter','latex'); ylabel('$|\psi(T)|$','Interpreter','latex');
title('Skew below two weeks, recalibrated','Interpreter','latex','FontWeight','normal');
legend('Location','southwest','Box','off','FontSize',7,'Interpreter','tex');
set(findall(f,'-property','FontSize'),'FontSize',8);
try, exportgraphics(f, fullfile(figdir,'fig_estensione.pdf'),'ContentType','vector');
catch, set(f,'PaperPositionMode','auto'); print(f, fullfile(figdir,'fig_estensione.pdf'),'-dpdf','-painters'); end
close(f);

save(fullfile(resdir,'estensione.mat'),'R','PRO','Hg','Tk','psiM','PSI','SL','RA','LE','cf', ...
    'T_uniq','T_cor','T_std','Kvar','MOD','GIORNI_CORTI');
fprintf('\nScritti Results/estensione.mat, Figure/fig_estensione.pdf, Figure/tab_estensione.tex\n');
toc

% FUNZIONI LOCALI
function th = CaricaTheta(resdir, nome, default)
    f = fullfile(resdir, nome);
    if isfile(f), S = load(f,'theta_ottimi'); th = S.theta_ottimi(:).';
    else, warning('%s non trovato: uso i valori della tesi.', nome); th = default; end
end

function cfg = Config(nm, T, stadio)
% [Xc, pow, nodi Riccati].
    fast = strcmp(stadio,'fast');
    if T < 0.03
        if fast, cfg = [2 12]; else, cfg = [2 14]; end
    elseif strcmp(nm,'rhh')
        if fast, cfg = [12 12]; else, cfg = [18 13]; end
    else
        if fast, cfg = [6 10]; else, cfg = [12 12]; end
    end
    if strcmp(nm,'rhh'), nr = 250 + 150*(~fast); else, nr = 100 + 300*(~fast); end
    cfg = [cfg nr];
end

function P = PrezziPut(nm, p, F, T, K, T_uniq, W_pp, stadio)
% Put forward (CP = -1) del modello nm sulle quote del campione.
    P = zeros(numel(K),1);
    for k = 1:numel(T_uniq)
        Tk = T_uniq(k); idx = (T == Tk);
        if ~any(idx), continue; end
        P(idx) = PutScadenza(nm, p, F(find(idx,1)), Tk, K(idx), W_pp, Config(nm, Tk, stadio));
    end
end

function P = PutScadenza(nm, p, Ft, Tk, Kt, W_pp, cfg)
    Kt = Kt(:); n = numel(Kt); Z = zeros(n,1); Tv = Tk*ones(n,1);
    xc = cfg(1); pw = cfg(2); nr = cfg(3);
    ti = ((0:nr-1)+0.5)/nr*Tk; xg = eval_xi0_pchip(ti, W_pp); xg = xg(:)';
    switch nm
        case 'hck', P = PrezzaSuperficieSINC_yield(Ft, Tv, Kt, Z, Z, 'Heston curve kappa', p, xg, nr, xc, 0, pw, -1);
        case 'rh',  P = PrezzaSuperficieSINC_yield(Ft, Tv, Kt, Z, Z, 'rough Heston', [p 0], xg, nr, xc, 0, pw, -1);
        case 'rhk', P = PrezzaSuperficieSINC_yield(Ft, Tv, Kt, Z, Z, 'rough Heston kappa', p, xg, nr, xc, 0, pw, -1);
        case 'rhh', P = PrezzaSuperficieSINC_yield_new(Ft, Tv, Kt, Z, Z, 'rough Hawkes Heston', p, [], nr, xc, 0, pw, -1);
    end
    P = P(:);
end

function M = Metriche(P, F, T, K, IVm, is_cor)
    IV = InvertiIV_OTM(F, K, T, max(P,1e-8));
    e  = (IVm - IV) ./ IVm;
    r  = @(x) 100*sqrt(mean(x.^2,'omitnan'));
    M  = struct('tot', r(e), 'cor', r(e(is_cor)), 'std', r(e(~is_cor)));
end

function s = SkewModelloPut(nm, p, Tk, W_pp)
    dk = min(0.005, 0.25*0.15*sqrt(Tk)); F0 = 100; Kp = F0*exp([-dk; dk]);
    P = PutScadenza(nm, p, F0, Tk, Kp, W_pp, [2 14 400]);
    iv = InvertiIV_OTM([F0;F0], Kp, [Tk;Tk], max(P,1e-10));
    s = (iv(2) - iv(1))/(2*dk);
end

function s = SkewMercato(k, iv, kw)
    m = abs(k) <= kw;
    if sum(m) < 5 || sum(k(m) < 0) < 2 || sum(k(m) > 0) < 2, s = NaN; return; end
    p = polyfit(k(m), iv(m), 2); s = p(2);
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

function xi = eval_xi0_pchip(tvec, W_pp)
    h = 1e-5; tvec = tvec(:);
    tp = max(tvec+h,0); tm = max(tvec-h,0);
    dn = tp - tm; dn(dn < 1e-12) = 1e-12;
    xi = max((ppval(W_pp,tp) - ppval(W_pp,tm))./dn, 1e-8);
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
