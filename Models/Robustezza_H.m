clearvars; clc;
tic
% ROBUSTEZZA_H  Due controlli sulla ruvidezza.
% Parte 1: profilo della loss in H per rough Heston e rough Heston + CGMY.
% Parte 2: multistart del composito sulla curva SSVI.
% Richiede Results/sensibilita_modelli.mat (da Sensibilita_SVI).
% Output: Results/profile_H.mat, Results/multistart_rc_ssvi.mat, figure e
% tabelle in Figure/. Tempo: parte 1 3-5 ore, parte 2 2-4 ore.
FAI_PROFILO    = false;  % parte 1 (true = ricalcolo completo, 3-5 ore)
RIDISEGNA_PROFILO = true;  % true con FAI_PROFILO = false
                           % e rifa' solo figura e tabella (pochi secondi),
                           % con banda e RMSE del rough Hawkes aggiornati
FAI_MULTISTART = false;  % parte 2 (true = 2-4 ore)

% parte 1
CURVE   = [1 3];  % 1 = piatta (base), 2 = lineare, 3 = SSVI
H_GRID  = [0.02 0.05 0.087 0.12 0.16 0.20 0.25 0.30 0.35 0.385 0.42 0.46];
RAFFINA_EXACT = true;  % false = solo pricer fast
% RMSE_RHH e HEFF_RHH sono letti da multistart_rhh.mat piu' sotto (dopo
% resdir): prima erano scritti a mano e non si aggiornavano quando il rough
% Hawkes Heston veniva ricalibrato.

% parte 2
CURVA      = 3;  % curva del multistart (3 = SSVI)
N_START    = 10;  % partenze Sobol (oltre all'ottimo attuale)
N_REFINE   = 3;  % bacini da raffinare con il pricer exact
TOL_BACINO = 0.05;  % distanza max-norm nello spazio normalizzato
USE_PARFOR = false;  % true se e' installato il Parallel Computing Toolbox

ETI = {'Flat','Linear','SSVI'};  % etichette brevi delle tre curve

root = fileparts(mfilename('fullpath'));
if isempty(root), root = pwd; end
resdir = fullfile(root, '..', 'Results');
figdir = fullfile(root, '..', 'Figure');
MS_RHH   = load(fullfile(resdir,'multistart_rhh.mat'), 'theta', 'ref_rmse', 'loss');
RMSE_RHH = min(MS_RHH.ref_rmse);  % rough Hawkes Heston
% banda di H_eff: solo le partenze nella valle dell'ottimo, cosi' un bacino
% degenere (es. alpha ~ 0.51 con Lambda al bound) non la allarga
ok_ms    = MS_RHH.loss <= 1.01*min(MS_RHH.loss);
HEFF_RHH = [min(MS_RHH.theta(ok_ms,1)), max(MS_RHH.theta(ok_ms,1))] - 0.5;
fprintf('Banda H_eff del rough Hawkes: [%.3f, %.3f] (%d partenze su %d)\n', ...
    HEFF_RHH, sum(ok_ms), numel(ok_ms));

% DATI: blocco copiato da Figure_Tesi.m (stesso campione di 105)
fprintf('Caricamento dati...\n');
raw = readtable('20260515SPX_merged.csv','Delimiter',',','VariableNamingRule','preserve');
is_spot = strcmp(raw.Root,'SPOT');
S0_feed = raw.CF_LAST(is_spot);  raw = raw(~is_spot,:);
is_call = strcmp(raw.PutCall,'C');
kp = ~isnan(raw.STRIKE_PRC) & raw.CF_BID>0 & raw.CF_ASK>raw.CF_BID & raw.CF_BID>=0.375;
raw = raw(kp,:); is_call = is_call(kp);
K_all=raw.STRIKE_PRC; mid_all=0.5*(raw.CF_BID+raw.CF_ASK); bid_all=raw.CF_BID;
T_all=(datenum(raw.EXPIR_DATE)-datenum('2026-05-15'))/365;
[~,ia]=unique([T_all,K_all,double(is_call)],'rows','stable');
T_all=T_all(ia); K_all=K_all(ia); mid_all=mid_all(ia);
bid_all=bid_all(ia); is_call=is_call(ia);

Tr=unique(T_all); Fp=nan(numel(Tr),1); Dp=nan(numel(Tr),1);
for i=1:numel(Tr)
    Ti=Tr(i); mC=(T_all==Ti)&is_call; mP=(T_all==Ti)&~is_call;
    [Kc,iC,iP]=intersect(K_all(mC),K_all(mP));
    if numel(Kc)<3, continue; end
    cI=find(mC); pI=find(mP); d=mid_all(cI(iC))-mid_all(pI(iP));
    nr=abs(Kc/S0_feed-1)<=0.25;
    if sum(nr)>=5, Kr=Kc(nr); dr=d(nr); else, Kr=Kc; dr=d; end
    co=[ones(size(Kr)),Kr]\dr; Dp(i)=-co(2); Fp(i)=co(1)/Dp(i);
end
v=~isnan(Fp)&Dp>0&Dp<1.05&Fp>0; Tr=Tr(v); Fp=Fp(v); Dp=Dp(v);

F_s=[];T_s=[];K_s=[];Mid_s=[];Bid_s=[];
for i=1:numel(Tr)
    Ti=Tr(i); Fi=Fp(i); Di=Dp(i); m=(T_all==Ti);
    Ki=K_all(m); mi=mid_all(m); bi=bid_all(m); ci=is_call(m); P=nan(size(Ki));
    for j=1:numel(Ki)
        if Ki(j)<Fi && ~ci(j),     P(j)=mi(j)/Di;
        elseif Ki(j)>=Fi && ci(j), P(j)=mi(j)/Di-Fi+Ki(j); end
    end
    o=~isnan(P)&P>1e-4;
    F_s=[F_s;Fi*ones(sum(o),1)]; T_s=[T_s;Ti*ones(sum(o),1)];
    K_s=[K_s;Ki(o)]; Mid_s=[Mid_s;P(o)]; Bid_s=[Bid_s;bi(o)];
end
IV_M=InvertiIV_OTM(F_s,K_s,T_s,Mid_s);
lm=log(K_s./F_s);
id=~isnan(IV_M)&IV_M>0.01&IV_M<2&Mid_s>1e-4&abs(lm)<0.5&T_s>=0.03&T_s<=3.0;
F_s=F_s(id);T_s=T_s(id);K_s=K_s(id);Mid_s=Mid_s(id);Bid_s=Bid_s(id);IV_M=IV_M(id);
Ta=unique(T_s); ma=false(numel(T_s),1);
for i=1:numel(Ta), m=T_s==Ta(i); if sum(m)>=5, ma=ma|m; end, end
F_s=F_s(ma);T_s=T_s(ma);K_s=K_s(ma);Mid_s=Mid_s(ma);Bid_s=Bid_s(ma);IV_M=IV_M(ma);
Ta=unique(T_s); Ts=[];
for tt=[0.0384 0.0833 0.25 0.5 0.75 1.0 1.5]
    [~,ii]=min(abs(Ta-tt)); Ts=[Ts;Ta(ii)];
end
Ts=unique(Ts); km=ismember(T_s,Ts);
F_s=F_s(km);T_s=T_s(km);K_s=K_s(km);Mid_s=Mid_s(km);Bid_s=Bid_s(km);IV_M=IV_M(km);
T_uniq=unique(T_s);

% Demeterfi -> V_fair e curva xi0
Kvar=zeros(numel(T_uniq),1); IV_ATM=zeros(numel(T_uniq),1);
for i=1:numel(T_uniq)
    Ti=T_uniq(i); m=(T_s==Ti);
    Ki=K_s(m); IVi=IV_M(m); Fi=F_s(find(m,1));
    [Ks,s]=sort(Ki); IVs=IVi(s);
    [~,a]=min(abs(Ks-Fi)); IV_ATM(i)=IVs(a);
    Kg=logspace(log10(0.05*Fi),log10(2.50*Fi),4000)';
    IVg=interp1(Ks,IVs,Kg,'linear',NaN);
    IVg(Kg<Ks(1))=IVs(1); IVg(Kg>Ks(end))=IVs(end);
    O=zeros(size(Kg));
    for j=1:numel(Kg)
        if Kg(j)<Fi, [~,p]=blkprice(Fi,Kg(j),0,Ti,IVg(j)); O(j)=p;
        else,        [c,~]=blkprice(Fi,Kg(j),0,Ti,IVg(j)); O(j)=c; end
    end
    Kvar(i)=(2/Ti)*trapz(Kg,O./Kg.^2);
end
W_pp=pchip([0;T_uniq(:)],[0;T_uniq(:).*Kvar(:)]);
xifn=@(tv) eval_xi0_pchip(tv,W_pp);

% 15 strike/scadenza
k15=false(size(T_s));
for i=1:numel(T_uniq)
    m=find(T_s==T_uniq(i));
    z=log(K_s(m)./F_s(m))./(IV_ATM(i)*sqrt(T_uniq(i)));
    if numel(m)<=15, k15(m)=true; continue; end
    for tt=linspace(max(-3.0,min(z)),min(2.0,max(z)),15)
        [~,jj]=min(abs(z-tt)); k15(m(jj))=true;
    end
end
F=F_s(k15); T=T_s(k15); K=K_s(k15); Mid=Mid_s(k15); IVmkt=IV_M(k15);
CP=-1*ones(size(F)); logM=log(K./F);
zz=zeros(size(logM));
for i=1:numel(T_uniq)
    m=(T==T_uniq(i)); zz(m)=logM(m)./(IV_ATM(i)*sqrt(T_uniq(i)));
end

% pesi del campione fisso (identici alla calibrazione)
Vega = zeros(numel(F),1);
for i = 1:numel(F)
    d1 = (log(F(i)/K(i))+0.5*IVmkt(i)^2*T(i))/(IVmkt(i)*sqrt(T(i)));
    Vega(i) = F(i)*sqrt(T(i))*normpdf(d1);
end
sqrt_w = 1./(max(Vega,1e-4).*max(IVmkt,0.01));
sqrt_w = min(sqrt_w, prctile(sqrt_w,99));

% curve xi0 da Sensibilita_modelli
SM = load(fullfile(resdir,'sensibilita_modelli.mat'),'RES','V','KVAR','T_uniq');
assert(numel(SM.T_uniq)==numel(T_uniq) && max(abs(SM.T_uniq(:)-T_uniq(:)))<1e-10, ...
    'Le scadenze non coincidono con sensibilita_modelli.mat.');
fprintf('Controllo curva piatta: max |Kvar - KVAR(:,1)|/Kvar = %.2e (atteso ~0)\n', ...
    max(abs(Kvar(:)-SM.KVAR(:,1))./Kvar(:)));

% pricer e bound (identici a Sensibilita_modelli.m)
[XC_RH, POW_RH] = XcPow('rHeston', 'exact');
G_FIX = 0.05;  Q_MAX = 0.98;
[XC_RC, POW_RC_F] = XcPow('rH+CGMY', 'fast');
[~,     POW_RC_E] = XcPow('rH+CGMY', 'exact');
lb_rh = [0.05, -0.99];  % [nu, rho]
ub_rh = [5.00,  0.10];
lb_rc = [0.05, -0.99, 1e-4,  5.0,   0.9 ];  % [nu, rho, q, M, Y]
ub_rc = [5.00,  0.10, Q_MAX, 200.0, 1.95];

opt_f = optimoptions('lsqnonlin','Algorithm','trust-region-reflective', ...
    'MaxFunctionEvaluations',1500,'MaxIterations',100,'Display','off', ...
    'StepTolerance',1e-8,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);
opt_e = optimoptions('lsqnonlin','Algorithm','trust-region-reflective', ...
    'MaxFunctionEvaluations',600,'MaxIterations',30,'Display','off', ...
    'StepTolerance',1e-8,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);

if FAI_PROFILO || RIDISEGNA_PROFILO
fprintf('\n\n############ PARTE 1: PROFILO IN H ############\n');
if FAI_PROFILO
nH = numel(H_GRID);  nC = numel(CURVE);
[~, i1y] = min(abs(T_uniq - 1));
PR = struct();

for ic = 1:nC
    iv  = CURVE(ic);
    eti = ETI{iv};
    Wv  = pchip([0; T_uniq(:)], [0; T_uniq(:).*SM.KVAR(:,iv)]);
    xifn_v = @(tv) eval_xi0_pchip(tv, Wv);
    tt_all = [];
    for k = 1:numel(T_uniq)
        tt_all = [tt_all, ((0:399)+0.5)/400*T_uniq(k), ((0:99)+0.5)/100*T_uniq(k)];
    end
    xi_min = min(xifn_v(tt_all));
    fprintf('\n================ Curva: %s ================\n', eti);

    % rough Heston
    th0 = SM.RES(iv).rh.theta;  % [H, nu, rho]
    A = struct('theta',nan(nH,3),'rmse',nan(nH,1),'loss',nan(nH,1),'bias',nan(nH,6), ...
               'H_free',th0(1),'rmse_free',SM.RES(iv).rh.rmse);
    [~, i0] = min(abs(H_GRID - th0(1)));
    fprintf('\n--- rough Heston (ottimo libero: H = %.4f, RMSE %.2f%%) ---\n', th0(1), A.rmse_free);
    for dir = [0 1]
        if dir == 0, order = i0:-1:1; else, order = i0+1:nH; end
        p = min(max(th0(2:3), lb_rh), ub_rh);
        for k = order
            Hk = H_GRID(k);
            rF = @(x) sqrt_w .* (PrezziRoughHeston(F,T,K,[Hk x],T_uniq,xifn_v,100,XC_RH,POW_RH,CP) - Mid);
            rE = @(x) sqrt_w .* (PrezziRoughHeston(F,T,K,[Hk x],T_uniq,xifn_v,400,XC_RH,POW_RH,CP) - Mid);
            p = lsqnonlin(rF, p, lb_rh, ub_rh, opt_f);
            if RAFFINA_EXACT, p = lsqnonlin(rE, p, lb_rh, ub_rh, opt_e); end
            Pm = PrezziRoughHeston(F,T,K,[Hk p],T_uniq,xifn_v,400,XC_RH,POW_RH,CP);
            R  = Valuta([Hk p], Pm, F, K, T, IVmkt, zz, T_uniq);
            A.theta(k,:) = [Hk p];  A.rmse(k) = R.rmse;  A.bias(k,:) = R.bias;
            A.loss(k) = sum((sqrt_w.*(Pm - Mid)).^2);
            fprintf('  H = %.3f | nu %.4f rho %+.4f | RMSE %6.2f%%  (%.0f s)\n', Hk, p, R.rmse, toc);
        end
    end
    PR(ic).rh = A;

    % rough Heston + CGMY
    th0 = SM.RES(iv).rc.theta;  % [H, nu, rho, q, M, Y]
    B = struct('theta',nan(nH,6),'rmse',nan(nH,1),'loss',nan(nH,1),'bias',nan(nH,6), ...
               'vsl',nan(nH,1),'quota',nan(nH,1),'al_bound',false(nH,1), ...
               'H_free',th0(1),'rmse_free',SM.RES(iv).rc.rmse);
    [~, i0] = min(abs(H_GRID - th0(1)));
    fprintf('\n--- rough Heston + CGMY (ottimo libero: H = %.4f, RMSE %.2f%%) ---\n', th0(1), B.rmse_free);
    for dir = [0 1]
        if dir == 0, order = i0:-1:1; else, order = i0+1:nH; end
        p = min(max(th0(2:6), lb_rc), ub_rc);
        for k = order
            Hk = H_GRID(k);
            rF = @(x) sqrt_w .* (PrezziRoughHestonCGMY(F,T,K,[Hk x],G_FIX,xi_min,T_uniq,xifn_v, ...
                     100,XC_RC,POW_RC_F,CP) - Mid);
            rE = @(x) sqrt_w .* (PrezziRoughHestonCGMY(F,T,K,[Hk x],G_FIX,xi_min,T_uniq,xifn_v, ...
                     400,XC_RC,POW_RC_E,CP) - Mid);
            p = lsqnonlin(rF, p, lb_rc, ub_rc, opt_f);
            if RAFFINA_EXACT, p = lsqnonlin(rE, p, lb_rc, ub_rc, opt_e); end
            Pm = PrezziRoughHestonCGMY(F,T,K,[Hk p],G_FIX,xi_min,T_uniq,xifn_v,400,XC_RC,POW_RC_E,CP);
            R  = Valuta([Hk p], Pm, F, K, T, IVmkt, zz, T_uniq);
            Ccg = p(3)*xi_min/LambdaCGMY(G_FIX, p(4), p(5));
            [Lcg, Vcg] = ContributiCGMY(Ccg, G_FIX, p(4), p(5));
            B.theta(k,:) = [Hk p];  B.rmse(k) = R.rmse;  B.bias(k,:) = R.bias;
            B.loss(k)  = sum((sqrt_w.*(Pm - Mid)).^2);
            B.vsl(k)   = Vcg / Lcg;
            B.quota(k) = 100 * Vcg / (SM.KVAR(i1y,iv) - Lcg + Vcg);
            B.al_bound(k) = any(abs(p - lb_rc) < 1e-3*(ub_rc-lb_rc) | abs(p - ub_rc) < 1e-3*(ub_rc-lb_rc));
            flag = ''; if B.al_bound(k), flag = '  <- parametro al bound'; end
            fprintf('  H = %.3f | nu %.3f rho %+.3f q %.3f M %6.2f Y %.3f | RMSE %6.2f%% | QV salti %4.1f%%  (%.0f s)%s\n', ...
                Hk, p, R.rmse, B.quota(k), toc, flag);
        end
    end
    PR(ic).rc  = B;
    PR(ic).eti = eti;  PR(ic).iv = iv;
    save(fullfile(resdir,'profile_H.mat'), 'PR','H_GRID','CURVE','RMSE_RHH','HEFF_RHH');
    fprintf('(salvataggio intermedio profile_H.mat)\n');
end
else
    % solo ridisegno: profilo gia' calcolato (non dipende dal rough Hawkes)
    OLD = load(fullfile(resdir,'profile_H.mat'), 'PR', 'H_GRID', 'CURVE');
    PR = OLD.PR; H_GRID = OLD.H_GRID; CURVE = OLD.CURVE;
    nH = numel(H_GRID); nC = numel(CURVE);
    fprintf('Profilo riletto da profile_H.mat: si rifanno solo figura e tabella.\n');
end

% RIEPILOGO
fprintf('\n\n================ RIEPILOGO ================\n');
fprintf('Relative RMSE (%%) con H fissato\n%-8s', 'H');
for ic = 1:nC, fprintf(' %13s %13s', ['rH ' PR(ic).eti], ['rC ' PR(ic).eti]); end
fprintf('\n');
for k = 1:nH
    fprintf('%-8.3f', H_GRID(k));
    for ic = 1:nC, fprintf(' %13.2f %13.2f', PR(ic).rh.rmse(k), PR(ic).rc.rmse(k)); end
    fprintf('\n');
end
fprintf('(rH = rough Heston, rC = rough Heston + CGMY; rough Hawkes Heston = %.2f%%)\n', RMSE_RHH);

inBand = H_GRID >= HEFF_RHH(1) & H_GRID <= HEFF_RHH(2);
fprintf('\nNella banda H in [%.2f, %.2f] del rough Hawkes Heston:\n', HEFF_RHH);
for ic = 1:nC
    a = min(PR(ic).rh.rmse(inBand));  b = min(PR(ic).rc.rmse(inBand));
    fprintf('  %-26s rough Heston %.2f%% (ottimo %.2f%%, x%.2f) | composito %.2f%% (ottimo %.2f%%, x%.2f)\n', ...
        PR(ic).eti, a, PR(ic).rh.rmse_free, a/PR(ic).rh.rmse_free, ...
        b, PR(ic).rc.rmse_free, b/PR(ic).rc.rmse_free);
end
fprintf(['\nLettura: se il composito nella banda e'' molto sopra il suo ottimo e sopra\n' ...
         'il %.2f%% del rough Hawkes Heston, la differenza di ruvidezza fra i due\n' ...
         'meccanismi e'' identificata dai prezzi. Se il profilo e'' piatto, non lo e''.\n'], RMSE_RHH);

% FIGURA
c_rh = [0.000 0.447 0.741];  c_rc = [0.850 0.325 0.098];
f = figure('Color','w','Position',[100 100 430*nC 360]);
ymax = 1.08 * max([arrayfun(@(s) max(s.rh.rmse), PR), arrayfun(@(s) max(s.rc.rmse), PR), RMSE_RHH]);
for ic = 1:nC
    subplot(1, nC, ic); hold on; box on; grid on;
    patch([HEFF_RHH(1) HEFF_RHH(2) HEFF_RHH(2) HEFF_RHH(1)], [0 0 ymax ymax], ...
        [0.92 0.92 0.92], 'EdgeColor','none', 'HandleVisibility','off');
    h1 = plot(H_GRID, PR(ic).rh.rmse, 'o-', 'Color', c_rh, 'LineWidth', 1.4, 'MarkerSize', 4);
    h2 = plot(H_GRID, PR(ic).rc.rmse, 's-', 'Color', c_rc, 'LineWidth', 1.4, 'MarkerSize', 4);
    plot(PR(ic).rh.H_free, PR(ic).rh.rmse_free, 'o', 'MarkerFaceColor', c_rh, ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'HandleVisibility','off');
    plot(PR(ic).rc.H_free, PR(ic).rc.rmse_free, 's', 'MarkerFaceColor', c_rc, ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'HandleVisibility','off');
    h3 = yline(RMSE_RHH, 'k--', 'LineWidth', 1.1);
    xlim([0 0.5]); ylim([0 ymax]);
    xlabel('$H$', 'Interpreter', 'latex');
    ylabel('Relative RMSE (\%)', 'Interpreter', 'latex');
    title([PR(ic).eti ' wing rule'], 'Interpreter', 'none');
    if ic == 1
        legend([h1 h2 h3], {'rough Heston','rough Heston + CGMY','rough Hawkes Heston'}, ...
            'Location','northwest','Box','off','Interpreter','none');
    end
end
set(findall(f,'-property','FontSize'),'FontSize',9);
try
    exportgraphics(f, fullfile(figdir,'fig_profile_H.pdf'), 'ContentType','vector');
catch
    set(f,'PaperPositionMode','auto');
    print(f, fullfile(figdir,'fig_profile_H.pdf'), '-dpdf', '-painters', '-bestfit');
end

% TABELLA LaTeX (appendice)
fid = fopen(fullfile(figdir,'tab_app_profile_H.tex'),'w');
fprintf(fid,'%% Generato da Profilo_H.m il %s -- non modificare a mano\n', datestr(now));
fprintf(fid,'\\begin{tabular}{l%s}\n\\toprule\n', repmat('cc',1,nC));
fprintf(fid,' ');
for ic = 1:nC
    fprintf(fid,' & \\multicolumn{2}{c}{%s}', PR(ic).eti);
end
fprintf(fid,'\\\\\n');
for ic = 1:nC, fprintf(fid,'\\cmidrule(lr){%d-%d}', 2*ic, 2*ic+1); end
fprintf(fid,'\n$H$');
for ic = 1:nC, fprintf(fid,' & rough Heston & + CGMY'); end
fprintf(fid,'\\\\\n\\midrule\n');
for k = 1:nH
    fprintf(fid,'%.3f', H_GRID(k));
    for ic = 1:nC, fprintf(fid,' & %.2f & %.2f', PR(ic).rh.rmse(k), PR(ic).rc.rmse(k)); end
    fprintf(fid,'\\\\\n');
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

save(fullfile(resdir,'profile_H.mat'), 'PR','H_GRID','CURVE','RMSE_RHH','HEFF_RHH');
fprintf('\nSalvati profile_H.mat, fig_profile_H.pdf, tab_app_profile_H.tex\n');
end  % FAI_PROFILO

if FAI_MULTISTART
fprintf('\n\n############ PARTE 2: MULTISTART COMPOSITO (curva %s) ############\n', ETI{CURVA});
% curva xi0 da Sensibilita_modelli
eti = ETI{CURVA};
Wv  = pchip([0; T_uniq(:)], [0; T_uniq(:).*SM.KVAR(:,CURVA)]);
xifn_v = @(tv) eval_xi0_pchip(tv, Wv);
tt_all = [];
for k = 1:numel(T_uniq)
    tt_all = [tt_all, ((0:399)+0.5)/400*T_uniq(k), ((0:99)+0.5)/100*T_uniq(k)];
end
xi_min = min(xifn_v(tt_all));
[~, i1y] = min(abs(T_uniq - 1));

% pricer, bound e scala (identici a Multistart.m)
G_FIX = 0.05;  Q_MAX = 0.98;
[XC_RC, POW_RC_F] = XcPow('rH+CGMY', 'fast');
[~,     POW_RC_E] = XcPow('rH+CGMY', 'exact');
%      [H,     nu,   rho,   q,     M,     Y   ]
lb   = [0.005, 0.05, -0.99, 1e-4,  5.0,   0.9 ];
ub   = [0.49,  5.00,  0.10, Q_MAX, 200.0, 1.95];
LOGS = [true,  true, false, true,  true,  false];
nvars = 6;

rF = @(p) ResiduiRC_safe(p, F,T,K,G_FIX,xi_min,T_uniq,xifn_v,100,XC_RC,POW_RC_F,CP,Mid,sqrt_w);
rE = @(p) ResiduiRC_safe(p, F,T,K,G_FIX,xi_min,T_uniq,xifn_v,400,XC_RC,POW_RC_E,CP,Mid,sqrt_w);

seed = min(max(SM.RES(CURVA).rc.theta(:)', lb), ub);
fprintf('Curva: %s | ottimo locale attuale: H %.4f, RMSE %.2f%%\n', eti, seed(1), SM.RES(CURVA).rc.rmse);

% punti di partenza
rng(42);
P = scramble(sobolset(nvars, 'Skip', 1), 'MatousekAffineOwen');
U = net(P, N_START);
starts = zeros(N_START+1, nvars);
starts(1,:) = seed;
for k = 1:N_START, starts(k+1,:) = DaUnitario(U(k,:), lb, ub, LOGS); end

opt_fm = optimoptions('lsqnonlin','Algorithm','trust-region-reflective', ...
    'MaxFunctionEvaluations',2000,'MaxIterations',150,'Display','off', ...
    'StepTolerance',1e-8,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);
opt_em = optimoptions('lsqnonlin','Algorithm','trust-region-reflective', ...
    'MaxFunctionEvaluations',800,'MaxIterations',40,'Display','off', ...
    'StepTolerance',1e-8,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);

% ottimizzazione fast da ogni partenza
nS = size(starts,1);
theta = nan(nS, nvars);  loss = inf(nS,1);  tempo = nan(nS,1);
fprintf('\n=== MULTISTART rough Heston + CGMY (%s): %d partenze ===\n', eti, nS);
if USE_PARFOR
    parfor k = 1:nS
        t0 = tic;
        try
            th = lsqnonlin(rF, starts(k,:), lb, ub, opt_fm);
            theta(k,:) = th;  loss(k) = sum(rF(th).^2);
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
            th = lsqnonlin(rF, starts(k,:), lb, ub, opt_fm);
            theta(k,:) = th;  loss(k) = sum(rF(th).^2);
        catch ME
            fprintf('%3d  FALLITA (%s)\n', k, ME.message);
        end
        tempo(k) = toc(t0);
        e = ''; if k == 1, e = '   <- ottimo locale attuale'; end
        fprintf('%3d  loss fast = %.6e  (%.0f s)  H = %.4f | %s%s\n', k, loss(k), tempo(k), ...
            theta(k,1), strjoin(compose('%.4g', theta(k,:)), ' '), e);
    end
end

ok = isfinite(loss);
assert(any(ok), 'Tutte le partenze sono fallite.');
th_ok = theta(ok,:);
[loss_ord, ord] = sort(loss(ok));
theta_ord = th_ok(ord,:);

% bacini
Z = zeros(size(theta_ord));
for k = 1:size(theta_ord,1), Z(k,:) = AUnitario(theta_ord(k,:), lb, ub, LOGS); end
bacino = zeros(size(Z,1),1);  nb = 0;
for k = 1:size(Z,1)
    for bb = 1:nb
        rap = find(bacino == bb, 1);
        if max(abs(Z(k,:) - Z(rap,:))) < TOL_BACINO, bacino(k) = bb; break; end
    end
    if bacino(k) == 0, nb = nb + 1; bacino(k) = nb; end
end

% raffinamento exact del migliore di ogni bacino
[~, primo] = unique(bacino, 'stable');
primo = primo(1:min(N_REFINE, numel(primo)));
nR = numel(primo);
ref_theta = nan(nR, nvars);  ref_loss = nan(nR,1);  ref_rmse = nan(nR,1);
ref_bias  = nan(nR, 6);      ref_vsl  = nan(nR,1);  ref_quota = nan(nR,1);
fprintf('\n--- RAFFINAMENTO EXACT (migliore di ciascun bacino) ---\n');
for j = 1:nR
    k = primo(j);
    try
        th = lsqnonlin(rE, theta_ord(k,:), lb, ub, opt_em);
    catch
        th = theta_ord(k,:);
    end
    Pm = PrezziRoughHestonCGMY(F,T,K,th,G_FIX,xi_min,T_uniq,xifn_v,400,XC_RC,POW_RC_E,CP);
    R  = Valuta(th, Pm, F, K, T, IVmkt, zz, T_uniq);
    Ccg = th(4)*xi_min/LambdaCGMY(G_FIX, th(5), th(6));
    [Lcg, Vcg] = ContributiCGMY(Ccg, G_FIX, th(5), th(6));
    ref_theta(j,:) = th;  ref_loss(j) = sum(rE(th).^2);  ref_rmse(j) = R.rmse;
    ref_bias(j,:)  = R.bias;  ref_vsl(j) = Vcg/Lcg;
    ref_quota(j)   = 100 * Vcg / (SM.KVAR(i1y,CURVA) - Lcg + Vcg);
    fprintf(['bacino %d: loss exact %.6e | RMSE %.3f%% | H %.4f nu %.4f rho %+.4f ' ...
             'q %.4f M %.2f Y %.4f | VS/L salti %.3f | QV salti %.1f%%\n'], ...
        bacino(k), ref_loss(j), R.rmse, th, ref_vsl(j), ref_quota(j));
    fprintf('          bias per bucket (%%): %s\n', strjoin(compose('%+.2f', R.bias), '  '));
end

% riepilogo
fprintf('\n=== RIEPILOGO ===\n');
fprintf('%-7s %-9s %-14s %-10s %-8s\n','bacino','partenze','loss fast min','x globale','H');
for bb = 1:nb
    mb = (bacino == bb);  kb = find(mb, 1);
    fprintf('%-7d %-9d %-14.6e %-10.3f %-8.4f\n', bb, sum(mb), loss_ord(kb), ...
        loss_ord(kb)/loss_ord(1), theta_ord(kb,1));
end
seed_loss = sum(rE(seed).^2);
[best_loss, ib] = min(ref_loss);
fprintf('\nOttimo locale attuale: loss exact %.6e (RMSE %.2f%%, H %.4f)\n', ...
    seed_loss, SM.RES(CURVA).rc.rmse, seed(1));
fprintf('Migliore trovato     : loss exact %.6e (RMSE %.2f%%, H %.4f)\n', ...
    best_loss, ref_rmse(ib), ref_theta(ib,1));
if best_loss < 0.999*seed_loss
    fprintf(2, '\nTrovato un ottimo MIGLIORE di quello attuale: mandami l''output prima di aggiornare la tesi.\n');
else
    fprintf('\nNessun bacino migliore: l''ottimo locale e'' confermato come globale nel box.\n');
end
vicini = loss_ord <= 1.10*loss_ord(1);
fprintf('Partenze entro il 10%% della loss migliore: %d su %d, H tra %.4f e %.4f\n', ...
    sum(vicini), numel(loss_ord), min(theta_ord(vicini,1)), max(theta_ord(vicini,1)));
fprintf('Partenze finite con H >= 0.30: %d (loss minima fra queste: x%.2f la migliore)\n', ...
    sum(theta_ord(:,1) >= 0.30), min([inf; loss_ord(theta_ord(:,1) >= 0.30)])/loss_ord(1));

% salvataggio + tabella LaTeX
save(fullfile(resdir,'multistart_rc_ssvi.mat'), 'starts','theta','loss','tempo', ...
    'theta_ord','loss_ord','bacino','nb','ref_theta','ref_loss','ref_rmse','ref_bias', ...
    'ref_vsl','ref_quota','seed','seed_loss','N_START','TOL_BACINO','lb','ub','CURVA');
fid = fopen(fullfile(figdir,'tab_app_multistart_rc_ssvi.tex'),'w');
fprintf(fid,'%% Generato da Multistart_RC_SSVI.m il %s -- non modificare a mano\n', datestr(now));
fprintf(fid,'\\begin{tabular}{rrrrrr}\n\\toprule\n');
fprintf(fid,'Basin & Restarts & Relative RMSE (\\%%) & $H$ & $q$ & $\\mathrm{VS}/L$\\\\\n\\midrule\n');
for j = 1:nR
    fprintf(fid,'%d & %d & %.2f & %.3f & %.3f & %.2f\\\\\n', bacino(primo(j)), ...
        sum(bacino == bacino(primo(j))), ref_rmse(j), ref_theta(j,1), ref_theta(j,4), ref_vsl(j));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);
fprintf('\nSalvati multistart_rc_ssvi.mat e tab_app_multistart_rc_ssvi.tex\n');
end  % FAI_MULTISTART

toc

% FUNZIONI LOCALI
function res = ResiduiRC_safe(p, F,T,K,G,xi_min,T_uniq,xifn,N,Xc,pw,CP,Mid,sw)
    try
        res = sw .* (PrezziRoughHestonCGMY(F,T,K,p,G,xi_min,T_uniq,xifn,N,Xc,pw,CP) - Mid);
        res(~isfinite(res)) = 1e3;
    catch
        res = 1e3 * ones(size(Mid));
    end
end

function x = DaUnitario(u, lb, ub, logs)
    x = zeros(size(u));
    for j = 1:numel(u)
        if logs(j), x(j) = exp(log(lb(j)) + u(j)*(log(ub(j))-log(lb(j))));
        else,       x(j) = lb(j) + u(j)*(ub(j)-lb(j)); end
    end
end

function u = AUnitario(x, lb, ub, logs)
    u = zeros(size(x));
    for j = 1:numel(x)
        if logs(j), u(j) = (log(x(j))-log(lb(j))) / (log(ub(j))-log(lb(j)));
        else,       u(j) = (x(j)-lb(j)) / (ub(j)-lb(j)); end
    end
end

function R = Valuta(th, Pm, F, K, T, IVmkt, zz, T_uniq)
    IVv = InvertiIV_OTM(F, K, T, max(Pm, 1e-6));
    ok  = ~isnan(IVv) & ~isnan(IVmkt);
    e   = (IVmkt - IVv) ./ IVmkt;  % positivo = modello sotto il mercato
    ZB  = [-inf, -2, -1, -0.5, 0.5, 1, inf];
    bias = nan(1, numel(ZB)-1);
    for b = 1:numel(ZB)-1
        mb = ok & zz > ZB(b) & zz <= ZB(b+1);
        bias(b) = 100*mean(e(mb));
    end
    rT = nan(1, numel(T_uniq));
    for it = 1:numel(T_uniq)
        mt = ok & (T == T_uniq(it));
        rT(it) = 100*sqrt(mean(e(mt).^2));
    end
    R = struct('theta', th, 'rmse', 100*sqrt(mean(e(ok).^2)), 'bias', bias, 'rmse_T', rT);
end

function lam = LambdaCGMY(G, M, Y)
% Copiata da Multistart.m / Rough_Heston_CMGY_def.m
    lam = 2*( gamma(-Y) * ((M-1)^Y - M^Y + (G+1)^Y - G^Y) ...
            - gamma(1-Y) * (M^(Y-1) - G^(Y-1)) );
end

function [Lcg, Vcg] = ContributiCGMY(C, G, M, Y)
    Lcg = C * LambdaCGMY(G, M, Y);
    Vcg = C * gamma(2-Y) * (M^(Y-2) + G^(Y-2));
end

function P_fwd = PrezziRoughHestonCGMY(F_vec, T_vec, K_vec, params6, G_fix, xi_min, ...
                    T_uniq, xi0_fn, N_riccati, Xc, pow, CP_vec)
% Copiata da Multistart.m (ramo pchip di Rough_Heston_CMGY_def.m)
    H = params6(1); nu = params6(2); rho = params6(3);
    q = params6(4); M  = params6(5); Y   = params6(6);
    G = G_fix;
    C = q * xi_min / LambdaCGMY(G, M, Y);
    params_combined = [H, nu, rho, 0, C, G, M, Y];
    L_cgmy = ContributiCGMY(C, G, M, Y);
    P_fwd = zeros(length(K_vec),1);
    for k = 1:length(T_uniq)
        T = T_uniq(k); idx = (T_vec == T);
        if ~any(idx), continue; end
        F_t = F_vec(find(idx,1)); K_t = K_vec(idx); CP_t = CP_vec(idx);
        ti = ((0:N_riccati-1) + 0.5) / N_riccati * T;
        xi_grid = max(xi0_fn(ti) - L_cgmy, 1e-6); xi_grid = xi_grid(:)';
        P_fwd(idx) = PrezzaSuperficieSINC_yield(F_t, T*ones(size(K_t)), K_t, ...
            zeros(size(K_t)), zeros(size(K_t)), ...
            'rough Heston CGMY', params_combined, xi_grid, N_riccati, Xc, 0, pow, CP_t(1));
    end
end

function P = PrezziRoughHeston(F_vec,T_vec,K_vec,p3,T_uniq,xifn,N_ric,Xc,pw,CP_vec)
    pc = [p3(1), p3(2), p3(3), 0];
    P = zeros(numel(K_vec),1);
    for k = 1:numel(T_uniq)
        T = T_uniq(k); idx = (T_vec == T);
        if ~any(idx), continue; end
        Ft = F_vec(find(idx,1)); Kt = K_vec(idx); CPt = CP_vec(idx);
        ti = ((0:N_ric-1)+0.5)/N_ric*T;
        xg = xifn(ti); xg = xg(:)';
        P(idx) = PrezzaSuperficieSINC_yield(Ft, T*ones(size(Kt)), Kt, ...
            zeros(size(Kt)), zeros(size(Kt)), ...
            'rough Heston', pc, xg, N_ric, Xc, 0, pw, CPt(1));
    end
end

function xi = eval_xi0_pchip(tvec, W_pp)
    h=1e-5; tvec=tvec(:);
    tp=max(tvec+h,0); tm=max(tvec-h,0);
    dn=tp-tm; dn(dn<1e-12)=1e-12;
    xi=max((ppval(W_pp,tp)-ppval(W_pp,tm))./dn,1e-8);
end

function IV = InvertiIV_OTM(F,K,T,P)
    IV=nan(size(F));
    for i=1:numel(F)
        try
            if K(i)<=F(i)
                IV(i)=blkimpv(F(i),K(i),0,T(i),P(i),'Class','put');
            else
                IV(i)=blkimpv(F(i),K(i),0,T(i),max(P(i)+F(i)-K(i),1e-6),'Class','call');
            end
        catch
            IV(i)=NaN;
        end
    end
end
