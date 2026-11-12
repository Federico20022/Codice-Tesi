clearvars; clc;
tic
% FIGURE_TESI  Figure e tabelle dei Cap. 6 e 7 e dell'Appendice C.
% Uso: Run_All, poi Figure_Tesi.
% Riprezza i modelli calibrati con i domini di XcPow e sovrappone le IV sugli
% stessi assi. Include la sensibilita' di V_fair alle regole sulle ali.
% Output: Figure/fig_*.pdf, Figure/tab_*.tex, Results/sensibilita_Vfair.mat.
% Tempo: circa 9 minuti (4 con RIFAI_SENSIBILITA = false).

RIFAI_SENSIBILITA = true;

root = fileparts(mfilename('fullpath'));
if isempty(root), root = pwd; end
resdir = fullfile(root, '..', 'Results');
figdir = fullfile(root, '..', 'Figure');
if ~exist(figdir,'dir'), mkdir(figdir); end

% Pulizia della cartella: niente backup, in Figure restano solo i file di
% QUESTA corsa.
vecchi = [dir(fullfile(figdir,'*.pdf')); dir(fullfile(figdir,'*.png')); ...
          dir(fullfile(figdir,'tab_*.tex')); dir(fullfile(figdir,'tabelle.tex'))];
for i = 1:numel(vecchi)
    delete(fullfile(figdir, vecchi(i).name));
end
bk = dir(fullfile(figdir,'_backup_*'));
for i = 1:numel(bk)
    if bk(i).isdir, rmdir(fullfile(figdir, bk(i).name), 's'); end
end
if ~isempty(vecchi) || ~isempty(bk)
    fprintf('Cartella Figure ripulita (%d file, %d backup rimossi)\n', ...
        numel(vecchi), numel(bk));
end

COL_MKT = [0.15 0.15 0.15];
CC  = { [0.00 0.45 0.74], ...   % Heston
        [0.30 0.60 0.85], ...   % Heston-curve
        [0.85 0.33 0.10], ...   % rough Heston
        [0.47 0.67 0.19], ...   % CGMY
        [0.49 0.18 0.56], ...   % rH+CGMY
        [0.93 0.69 0.13], ...   % rough Hawkes Heston
        [0.10 0.60 0.55] };  % Heston-curve-kappa
MOD   = {'Heston','Heston-curve','rough Heston','CGMY','rough Heston + CGMY','rough Hawkes Heston','Heston-curve-\kappa'};
% nomi per le etichette con interprete latex e per le tabelle LaTeX
MODTEX = strrep(MOD, '\kappa', '$\kappa$');
SIGLA = {'heston','hestoncurve','rheston','cgmy','rhcgmy','rhh','hestoncurvekappa'};
NPAR  = [5 2 3 3 6 7 3];
NM    = numel(MOD);
% Il modello con kappa e' in coda negli indici, ma nelle tabelle e nelle
% figure compare subito dopo la Heston-curve: ORD e' l'ordine di stampa.
ORD   = [1 2 7 3 4 5 6];

% DATI (blocco standard, identico agli script di calibrazione)
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
fprintf('Campione: %d opzioni su %d scadenze\n', numel(F), numel(T_uniq));

% PREZZI E IV DEI SEI MODELLI
C.hes  = load(fullfile(resdir,'calib_heston_7x15z.mat'),'theta_ottimi');
C.hcur = load(fullfile(resdir,'calib_hestoncurve_7x15z.mat'),'theta_ottimi');
C.rh   = load(fullfile(resdir,'calib_rheston_7x15z_curva.mat'),'theta_ottimi');
C.cg   = load(fullfile(resdir,'calib_cgmy_7x15z.mat'),'theta_ottimi');
C.rhc  = load(fullfile(resdir,'calib_rhcgmy_7x15z.mat'),'theta_ottimi');
C.rhh  = load(fullfile(resdir,'calib_rhh_7x15z.mat'),'theta_ottimi');
C.hck  = load(fullfile(resdir,'calib_hestoncurvekappa_7x15z.mat'),'theta_ottimi');

Z=zeros(size(F));
IVmod=nan(numel(F),NM);
fprintf('Prezzatura dei sette modelli...\n');

% 1) Heston
[xc,pw]=XcPow('Heston','exact');
P=PrezzaSuperficieSINC_yield_new(F,T,K,Z,Z,'Heston',C.hes.theta_ottimi,[],[],xc,0,pw,CP);
IVmod(:,1)=InvertiIV_OTM(F,K,T,max(P,1e-6));

% 2) Heston-curve (H = 1/2 esatto, forma chiusa, stessa curva xi0)
th=C.hcur.theta_ottimi;
[xc,pw]=XcPow('rHeston','exact');
IVmod(:,2)=IVdaModello(@(Ft,Tt,Kt,CPt) PrezzaSuperficieSINC_yield(Ft,Tt,Kt, ...
    zeros(size(Kt)),zeros(size(Kt)),'Heston curve',th, ...
    RiccatiXi(xifn,400,Tt(1)),400,xc,0,pw,CPt), F,T,K,CP,T_uniq);

% 3) rough Heston
[xc,pw]=XcPow('rHeston','exact');
IVmod(:,3)=IVdaModello(@(Ft,Tt,Kt,CPt) PrezzaSuperficieSINC_yield(Ft,Tt,Kt, ...
    zeros(size(Kt)),zeros(size(Kt)),'rough Heston',[C.rh.theta_ottimi,0], ...
    RiccatiXi(xifn,400,Tt(1)),400,xc,0,pw,CPt), F,T,K,CP,T_uniq);

% 4) CGMY
[xc,pw]=XcPow('CGMY','exact');
P=PrezzaSuperficieSINC_yield_new(F,T,K,Z,Z,'CGMY',C.cg.theta_ottimi,[],[],xc,0,pw,CP);
IVmod(:,4)=InvertiIV_OTM(F,K,T,max(P,1e-6));

% 5) rough Heston + CGMY
t7=C.rhc.theta_ottimi; Lcg=ContributiCGMY(t7(4),t7(5),t7(6),t7(7));
pc=[t7(1:3),0,t7(4:7)];
[xc,pw]=XcPow('rH+CGMY','exact');
IVmod(:,5)=IVdaModello(@(Ft,Tt,Kt,CPt) PrezzaSuperficieSINC_yield(Ft,Tt,Kt, ...
    zeros(size(Kt)),zeros(size(Kt)),'rough Heston CGMY',pc, ...
    max(RiccatiXi(xifn,400,Tt(1))-Lcg,1e-6),400,xc,0,pw,CPt), F,T,K,CP,T_uniq);

% 6) rough Hawkes Heston
[xc,pw]=XcPow('rHawkesH','exact');
P=PrezzaSuperficieSINC_yield_new(F,T,K,Z,Z,'rough Hawkes Heston', ...
    C.rhh.theta_ottimi,[],400,xc,0,pw,CP);
IVmod(:,6)=InvertiIV_OTM(F,K,T,max(P,1e-6));

% 7) Heston-curve-kappa (stessa curva xi0, kappa libero, forma chiusa)
[xc,pw]=XcPow('rHeston','exact');
IVmod(:,7)=IVdaModello(@(Ft,Tt,Kt,CPt) PrezzaSuperficieSINC_yield(Ft,Tt,Kt, ...
    zeros(size(Kt)),zeros(size(Kt)),'Heston curve kappa',C.hck.theta_ottimi, ...
    RiccatiXi(xifn,400,Tt(1)),400,xc,0,pw,CPt), F,T,K,CP,T_uniq);

% eq. 12 del paper: eps = (IV_mkt - IV_mod)/IV_mkt, positivo = modello sotto
EPS = -(IVmod - IVmkt) ./ IVmkt;
RMSE_tot = 100*sqrt(mean(EPS.^2,1,'omitnan'));
fprintf('RMSE_rel ricalcolati:\n');
for j=1:NM, fprintf('  %-22s %.4f%%\n', MOD{j}, RMSE_tot(j)); end
fprintf(['(devono coincidere con runall.txt; attesi circa 2.899 / 10.346 / 4.012 / ' ...
         '5.592 / 3.760 / 2.22 / 3.29)\n\n']);

% SENSIBILITA' DI V_fair ALLE IPOTESI SULLE ALI
% La curva xi0 non e' calibrata: e' imposta model-free dalla replica di
if RIFAI_SENSIBILITA
    fprintf('Sensibilita di V_fair: 8 varianti, 8 ricalibrazioni...\n');
    rng(42);

    V = {};
    V{end+1} = struct('eti','BASE (codice attuale)','etiEN','BASE (current code)','lo',0.05,'hi',2.50,'estr','flat','ngrid',4000,'logmcut',0.50,'bidmin',0);
    V{end+1} = struct('eti','troncamento stretto','etiEN','narrow truncation','lo',0.30,'hi',2.00,'estr','flat','ngrid',4000,'logmcut',0.50,'bidmin',0);
    V{end+1} = struct('eti','troncamento largo','etiEN','wide truncation','lo',0.02,'hi',4.00,'estr','flat','ngrid',4000,'logmcut',0.50,'bidmin',0);
    V{end+1} = struct('eti','estrap. lineare in var','etiEN','linear-in-variance extrap.','lo',0.05,'hi',2.50,'estr','linvar','ngrid',4000,'logmcut',0.50,'bidmin',0);
    V{end+1} = struct('eti','nessuna estrapolazione','etiEN','no extrapolation','lo',0.05,'hi',2.50,'estr','nessuna','ngrid',4000,'logmcut',0.50,'bidmin',0);
    V{end+1} = struct('eti','griglia 16000 punti','etiEN','16000-point grid','lo',0.05,'hi',2.50,'estr','flat','ngrid',16000,'logmcut',0.50,'bidmin',0);
    V{end+1} = struct('eti','smile |k|<0.35','etiEN','smile $|k|<0.35$','lo',0.05,'hi',2.50,'estr','flat','ngrid',4000,'logmcut',0.35,'bidmin',0);
    V{end+1} = struct('eti','bid >= 1.00','etiEN','bid $\geq 1.00$','lo',0.05,'hi',2.50,'estr','flat','ngrid',4000,'logmcut',0.50,'bidmin',1.00);

    nV = numel(V);
    KVAR = nan(numel(T_uniq), nV);
    for iv = 1:nV
        KVAR(:,iv) = VFairVariante(V{iv}, T_s, K_s, F_s, IV_M, Bid_s, T_uniq);
    end

    % sqrt(V_fair) in punti di volatilita', per variante
    fprintf('%-9s', 'T(yr)');
    for iv=1:nV, fprintf('%12s', sprintf('var%d', iv)); end
    fprintf('\n');
    for it = 1:numel(T_uniq)
        fprintf('%-9.4f', T_uniq(it));
        for iv=1:nV, fprintf('%12.3f', 100*sqrt(max(KVAR(it,iv),0))); end
        fprintf('\n');
    end
    fprintf('%-9s', 'scarto');
    for iv=1:nV
        d = 100*(sqrt(max(KVAR(:,iv),0)) - sqrt(max(KVAR(:,1),0)));
        fprintf('%12.3f', max(abs(d)));
    end
    fprintf('   <- max |scarto| vs BASE, punti di vol\n');
    for iv=1:nV, fprintf('  var%d = %s\n', iv, V{iv}.eti); end

    % pesi vega del campione fisso (identici alla calibrazione)
    Vega = zeros(numel(F),1);
    for i = 1:numel(F)
        d1=(log(F(i)/K(i))+0.5*IVmkt(i)^2*T(i))/(IVmkt(i)*sqrt(T(i)));
        Vega(i)=F(i)*sqrt(T(i))*normpdf(d1);
    end
    sqrt_w = 1./(max(Vega,1e-4).*max(IVmkt,0.01));
    sqrt_w = min(sqrt_w, prctile(sqrt_w,99));

    seed = C.rh.theta_ottimi(:)';
    lb = [0.005, 0.05, -0.99];
    ub = [0.49,  5.00,  0.10];
    [XC_RH, POW_RH] = XcPow('rHeston', 'exact');
    opt = optimoptions('lsqnonlin','Algorithm','trust-region-reflective', ...
        'MaxFunctionEvaluations',1500,'MaxIterations',100,'Display','off', ...
        'StepTolerance',1e-8,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);

    TH = nan(nV,3); RMSE = nan(nV,1); XIMIN = nan(nV,1); XIMAX = nan(nV,1);
    fprintf('%-24s %8s %8s %9s %10s %10s\n', ...
        'variante','H','nu','rho','RMSE_rel%','d vs BASE');
    for iv = 1:nV
        Wv = pchip([0; T_uniq(:)], [0; T_uniq(:).*KVAR(:,iv)]);
        xifn_v = @(tv) eval_xi0_pchip(tv, Wv);
        tt = [];
        for k=1:numel(T_uniq), tt=[tt, ((0:399)+0.5)/400*T_uniq(k)]; end
        XIMIN(iv)=min(xifn_v(tt)); XIMAX(iv)=max(xifn_v(tt));

        rF = @(p) sqrt_w .* (PrezziRoughHeston(F,T,K,p,T_uniq,xifn_v,100,XC_RH,POW_RH,CP) - Mid);
        rE = @(p) sqrt_w .* (PrezziRoughHeston(F,T,K,p,T_uniq,xifn_v,400,XC_RH,POW_RH,CP) - Mid);
        thv = lsqnonlin(rF, seed, lb, ub, opt);
        thv = lsqnonlin(rE, thv,  lb, ub, opt);
        TH(iv,:) = thv;

        Pm = PrezziRoughHeston(F,T,K,thv,T_uniq,xifn_v,400,XC_RH,POW_RH,CP);
        IVv = InvertiIV_OTM(F,K,T,max(Pm,1e-6));
        ok = ~isnan(IVv) & ~isnan(IVmkt);
        RMSE(iv) = 100*sqrt(mean(((IVmkt(ok)-IVv(ok))./IVmkt(ok)).^2));

        if iv==1, ds='   (riferimento)';
        else,     ds=sprintf('%+9.4f', RMSE(iv)-RMSE(1)); end
        fprintf('%-24s %8.4f %8.4f %+9.4f %10.4f %s\n', ...
            V{iv}.eti, thv(1), thv(2), thv(3), RMSE(iv), ds);
    end

    fprintf('Escursione dei parametri sulle varianti:\n');
    nomi={'H','nu','rho'};
    for j=1:3
        fprintf('  %-4s in [%.4f, %.4f]  (fattore %.2f)\n', ...
            nomi{j}, min(TH(:,j)), max(TH(:,j)), max(abs(TH(:,j)))/min(abs(TH(:,j))));
    end
    fprintf('  RMSE_rel in [%.4f%%, %.4f%%]\n\n', min(RMSE), max(RMSE));

    if ~exist(resdir,'dir'), mkdir(resdir); end
    save(fullfile(resdir,'sensibilita_Vfair.mat'), ...
        'V','KVAR','TH','RMSE','T_uniq','XIMIN','XIMAX');
    fprintf('Salvato %s\n\n', fullfile(resdir,'sensibilita_Vfair.mat'));
else
    fprintf('Sensibilita saltata (RIFAI_SENSIBILITA = false)\n\n');
end

% FIGURA 1 - SMILE: mercato contro i sei modelli
f=figure('Color','w','Position',[50 50 1300 700]);
for i=1:numel(T_uniq)
    subplot(2,4,i); hold on; grid on; box on
    m=(T==T_uniq(i)); [ks,si]=sort(logM(m));
    im=IVmkt(m); im=im(si);
    plot(ks,100*im,'o','Color',COL_MKT,'MarkerSize',4,'LineWidth',0.8);
    for j=ORD
        ij=IVmod(m,j); ij=ij(si);
        plot(ks,100*ij,'-','Color',CC{j},'LineWidth',1.3);
    end
    xline(0,'k:');
    title(sprintf('$T = %.4f$ (%d days)', T_uniq(i), round(365*T_uniq(i))), ...
        'Interpreter','latex');
    xlabel('$\log(K/F)$','Interpreter','latex');
    ylabel('IV (\%)','Interpreter','latex');
end
% Pannello 8: solo legenda.
ax=subplot(2,4,8); hold(ax,'on'); axis(ax,'off');
lh=gobjects(NM+1,1);
lh(1)=plot(ax,nan,nan,'o','Color',COL_MKT,'MarkerSize',5);
for j=1:NM, lh(j+1)=plot(ax,nan,nan,'-','Color',CC{ORD(j)},'LineWidth',2); end
legend(ax,lh,[{'Market'},MOD(ORD)],'Location','west','Interpreter','tex','Box','off');
Esporta(f, figdir, 'fig_smiles');

% FIGURA 1-bis - uno smile per modello, assi comuni
% Nella sovrapposta sei curve sugli stessi assi nascondono gli scarti
% piccoli.
ylim_com = zeros(numel(T_uniq),2);
for i=1:numel(T_uniq)
    m=(T==T_uniq(i));
    tutte=100*[IVmkt(m); reshape(IVmod(m,:),[],1)];
    tutte=tutte(isfinite(tutte));
    ylim_com(i,:)=[min(tutte)-0.5, max(tutte)+0.5];
end
for j=1:NM
    f=figure('Color','w','Position',[50 50 1300 620]);
    for i=1:numel(T_uniq)
        subplot(2,4,i); hold on; grid on; box on
        m=(T==T_uniq(i)); [ks,si]=sort(logM(m));
        im=IVmkt(m); im=im(si); ij=IVmod(m,j); ij=ij(si);
        plot(ks,100*im,'o','Color',COL_MKT,'MarkerSize',4,'LineWidth',0.8);
        plot(ks,100*ij,'-','Color',CC{j},'LineWidth',1.6);
        xline(0,'k:'); ylim(ylim_com(i,:));
        title(sprintf('$T = %.4f$ (%d days)', T_uniq(i), round(365*T_uniq(i))), ...
            'Interpreter','latex');
        xlabel('$\log(K/F)$','Interpreter','latex');
        ylabel('IV (\%)','Interpreter','latex');
    end
    % Ottavo pannello: residuo relativo contro z.
    subplot(2,4,8); hold on; grid on; box on
    plot(zz,100*EPS(:,j),'.','Color',CC{j},'MarkerSize',9);
    yline(0,'k-');
    xlabel('$z$','Interpreter','latex');
    ylabel('relative $\varepsilon$ (\%)','Interpreter','latex');
    title(sprintf('%s: RMSE $%.3f\\%%$', MODTEX{j}, RMSE_tot(j)),'Interpreter','latex');
    Esporta(f, figdir, ['fig_smiles_' SIGLA{j}]);
end

% FIGURA 2 - confronto nested: Heston-curve contro rough Heston
% Stessa curva xi0, stesso campione, stessi pesi, stesso pricer: l'unica
% differenza e' H.
f=figure('Color','w','Position',[50 50 1300 620]);
for i=1:numel(T_uniq)
    subplot(2,4,i); hold on; grid on; box on
    m=(T==T_uniq(i)); [ks,si]=sort(logM(m));
    im=IVmkt(m); im=im(si);
    i2=IVmod(m,2); i2=i2(si); i3=IVmod(m,3); i3=i3(si); i7=IVmod(m,7); i7=i7(si);
    plot(ks,100*im,'o','Color',COL_MKT,'MarkerSize',4,'LineWidth',0.8);
    plot(ks,100*i2,'-','Color',CC{2},'LineWidth',1.6);
    plot(ks,100*i7,'--','Color',CC{7},'LineWidth',1.6);
    plot(ks,100*i3,'-','Color',CC{3},'LineWidth',1.6);
    xline(0,'k:'); ylim(ylim_com(i,:));
    title(sprintf('$T = %.4f$', T_uniq(i)),'Interpreter','latex');
    xlabel('$\log(K/F)$','Interpreter','latex');
    ylabel('IV (\%)','Interpreter','latex');
end
ax=subplot(2,4,8); hold(ax,'on'); axis(ax,'off');
lh=gobjects(4,1);
lh(1)=plot(ax,nan,nan,'o','Color',COL_MKT,'MarkerSize',5);
lh(2)=plot(ax,nan,nan,'-','Color',CC{2},'LineWidth',2);
lh(3)=plot(ax,nan,nan,'--','Color',CC{7},'LineWidth',2);
lh(4)=plot(ax,nan,nan,'-','Color',CC{3},'LineWidth',2);
legend(ax,lh,{'Market', ...
    sprintf('H = 1/2, \\kappa = 0 (RMSE %.2f%%)',RMSE_tot(2)), ...
    sprintf('H = 1/2, \\kappa = %.2f (RMSE %.2f%%)',C.hck.theta_ottimi(1),RMSE_tot(7)), ...
    sprintf('H = %.4f, \\kappa = 0 (RMSE %.2f%%)',C.rh.theta_ottimi(1),RMSE_tot(3))}, ...
    'Location','west','Interpreter','tex','Box','off');
Esporta(f, figdir, 'fig_nested');

% FIGURA 3 - RMSE relativo per scadenza
RM=zeros(numel(T_uniq),NM);
for i=1:numel(T_uniq)
    m=(T==T_uniq(i));
    RM(i,:)=100*sqrt(mean(EPS(m,:).^2,1,'omitnan'));
end
f=figure('Color','w','Position',[100 100 900 400]);
b=bar(RM(:,ORD)); grid on; box on
for j=1:NM, b(j).FaceColor=CC{ORD(j)}; b(j).EdgeColor='none'; end
set(gca,'XTickLabel',arrayfun(@(t) sprintf('%.3f',t),T_uniq,'UniformOutput',false));
xlabel('Maturity (years)');
ylabel('Relative IV RMSE (\%)','Interpreter','latex');
legend(MOD(ORD),'Location','northwest','Interpreter','tex','Box','off');
Esporta(f, figdir, 'fig_rmse_scadenza');

% FIGURA 4 - RMSE relativo per moneyness standardizzata
zb=[-inf -2 -1 -0.5 0.5 1 inf];
zl={'$z<-2$','$-2..-1$','$-1..-0.5$','ATM','$0.5..1$','$z>1$'};
RZ=zeros(numel(zb)-1,NM);
for i=1:numel(zb)-1
    m=(zz>=zb(i))&(zz<zb(i+1));
    RZ(i,:)=100*sqrt(mean(EPS(m,:).^2,1,'omitnan'));
end
f=figure('Color','w','Position',[100 100 900 400]);
b=bar(RZ(:,ORD)); grid on; box on
for j=1:NM, b(j).FaceColor=CC{ORD(j)}; b(j).EdgeColor='none'; end
set(gca,'XTickLabel',zl,'TickLabelInterpreter','latex');
xlabel('Standardised moneyness $z=\log(K/F)/(\sigma_{ATM}\sqrt{T})$','Interpreter','latex');
ylabel('Relative IV RMSE (\%)','Interpreter','latex');
legend(MOD(ORD),'Location','northwest','Interpreter','tex','Box','off');
Esporta(f, figdir, 'fig_rmse_moneyness');

% FIGURA 5 - BIAS relativo per moneyness standardizzata
% Il bias, a differenza dell'RMSE, ha un segno: dice se il modello sta sopra
% o sotto il mercato.
BZ=zeros(numel(zb)-1,NM);
for i=1:numel(zb)-1
    m=(zz>=zb(i))&(zz<zb(i+1));
    BZ(i,:)=100*mean(EPS(m,:),1,'omitnan');
end
f=figure('Color','w','Position',[100 100 900 400]);
b=bar(BZ(:,ORD)); grid on; box on
for j=1:NM, b(j).FaceColor=CC{ORD(j)}; b(j).EdgeColor='none'; end
yline(0,'k-','LineWidth',1);
set(gca,'XTickLabel',zl,'TickLabelInterpreter','latex');
xlabel('$z$','Interpreter','latex');
ylabel({'Relative bias (\%)','positive: model below market'}, ...
    'Interpreter','latex');
legend(MOD(ORD),'Location','northeast','Interpreter','tex','Box','off');
Esporta(f, figdir, 'fig_bias_moneyness');

% FIGURA 6 - curva forward variance e ancore del variance swap
tg=linspace(1e-4,max(T_uniq),600);
f=figure('Color','w','Position',[100 100 700 400]);
plot(tg,100*sqrt(xifn(tg)),'-','Color',CC{3},'LineWidth',1.6); hold on; grid on; box on
plot(T_uniq,100*sqrt(Kvar),'ks','MarkerSize',7,'MarkerFaceColor','w','LineWidth',1.1);
plot(T_uniq,100*IV_ATM,'^','Color',COL_MKT,'MarkerSize',6,'MarkerFaceColor','w');
xlabel('$t$ (years)','Interpreter','latex');
ylabel('Volatility (\%)','Interpreter','latex');
legend({'$\sqrt{\xi_0(t)}$ (PCHIP on $W=T\,K_{var}$)', ...
        '$\sqrt{K_{var}(T)}$ (Demeterfi)','Market ATM IV'}, ...
       'Interpreter','latex','Location','southeast','Box','off');
Esporta(f, figdir, 'fig_xi0');

% FIGURE DAI .mat GIA' PRODOTTI
% La PricingSuite lavora su cinque modelli: i colori vanno rimappati,
% altrimenti le curve avrebbero tinte diverse fra figure.
CC5 = {CC{1}, CC{3}, CC{4}, CC{5}, CC{6}};
MOD5 = {'Heston','rough Heston','CGMY','rough Heston + CGMY','rough Hawkes Heston'};
MK5  = {'o','s','^','d','v'};
NM5  = 5;

S=load(fullfile(resdir,'PricingSuite_results.mat'));
f=figure('Color','w','Position',[100 100 700 420]);
for j=1:NM5
    loglog(S.T_skew, abs(S.SkewATM(:,j)),'o-','Color',CC5{j},'LineWidth',1.4); hold on
end
grid on; box on
xlabel('$T$ (years)','Interpreter','latex');
ylabel('$|\partial \sigma_{BS}/\partial k|$ ATM','Interpreter','latex');
legend(MOD5,'Location','southwest','Interpreter','none','Box','off');
Esporta(f, figdir, 'fig_skew_term');

f=figure('Color','w','Position',[100 100 700 400]);
b=bar(squeeze(S.ResVS(:,:,3))); grid on; box on
for j=1:NM5, b(j).FaceColor=CC5{j}; b(j).EdgeColor='none'; end
set(gca,'XTickLabel',arrayfun(@(t) sprintf('$T=%.3f$',t),S.T_vs,'UniformOutput',false), ...
    'TickLabelInterpreter','latex');
ylabel('Variance swap / log-contract','Interpreter','latex');
yline(1,'k--','LineWidth',1);
ylim([0.95 max(1.3,1.05*max(S.ResVS(:,:,3),[],'all'))]);
legend(MOD5,'Location','northwest','Interpreter','none','Box','off');
Esporta(f, figdir, 'fig_vs_logcontract');

if isfile(fullfile(resdir,'profile_G.mat'))
    Pg=load(fullfile(resdir,'profile_G.mat'));
    f=figure('Color','w','Position',[100 100 900 380]);
    subplot(1,2,1);
    semilogx(Pg.ResA.G,Pg.ResA.loss/min(Pg.ResA.loss),'o-','Color',CC{4},'LineWidth',1.4);
    grid on; box on; xline(0.05,'k--');
    xlabel('$G$','Interpreter','latex');
    ylabel('loss / loss$_{\min}$','Interpreter','latex');
    title('Pure CGMY','Interpreter','none');
    subplot(1,2,2);
    semilogx(Pg.ResB.G,Pg.ResB.loss/min(Pg.ResB.loss),'o-','Color',CC{5},'LineWidth',1.4);
    grid on; box on; xline(0.05,'k--');
    xlabel('$G$','Interpreter','latex');
    ylabel('loss / loss$_{\min}$','Interpreter','latex');
    title('rough Heston + CGMY','Interpreter','none');
    Esporta(f, figdir, 'fig_profile_G');
end

if isfile(fullfile(resdir,'sensibilita_Vfair.mat'))
    Vs=load(fullfile(resdir,'sensibilita_Vfair.mat'));
    f=figure('Color','w','Position',[100 100 900 380]);
        nV = numel(Vs.RMSE);
    CV = lines(nV); CV(1,:) = [0 0 0];
    subplot(1,2,1); hold on; grid on; box on
    hV = gobjects(nV,1);
    for k=1:nV
        hV(k) = plot(Vs.T_uniq,100*sqrt(max(Vs.KVAR(:,k),0)),'o-', ...
            'Color',CV(k,:),'MarkerFaceColor',CV(k,:),'MarkerSize',3,'LineWidth',1.2);
    end
    xlabel('$T$ (years)','Interpreter','latex');
    ylabel('$\sqrt{K_{var}}$ (\%)','Interpreter','latex');
    title('Variance swap by variant','Interpreter','none');
    if isfield(Vs,'V')
        eti = cellfun(@(s) s.etiEN, Vs.V, 'UniformOutput', false);
        eti{1} = 'base (flat extrapolation)';
        eti = arrayfun(@(k) sprintf('%d\\ %s',k,eti{k}), 1:nV, 'UniformOutput', false);
        legend(hV, eti, 'Location','southeast','FontSize',6,'Interpreter','latex','Box','off');
    end
    subplot(1,2,2); hold on; grid on; box on
    for k=1:nV
        plot(Vs.TH(k,1),Vs.RMSE(k),'o','MarkerSize',9, ...
            'MarkerFaceColor',CV(k,:),'MarkerEdgeColor','none');
        text(Vs.TH(k,1),Vs.RMSE(k),sprintf('  %d',k),'FontSize',8);
    end
    xlabel('Calibrated $H$','Interpreter','latex');
    ylabel('Relative RMSE (\%)','Interpreter','latex');
    title('rough Heston recalibrated','Interpreter','none');
    Esporta(f, figdir, 'fig_sensibilita_Vfair');
end

% CAPITOLO 7 - PRICING FUORI DAL CAMPIONE
% Tutto da PricingSuite_results.mat: qui non si riprezza nulla.

% dispersione relativa dei prezzi delle call
% In termini assoluti i cinque prezzi sono indistinguibili quasi ovunque: la
% struttura si vede solo come scarto dalla media fra modelli.
REL = nan(numel(S.T_call), numel(S.Kf), NM5);
for i = 1:numel(S.T_call)
    for k = 1:numel(S.Kf)
        p = squeeze(S.CallPrices(i,k,:));
        REL(i,k,:) = 100*(p - mean(p))/mean(p);
    end
end
yl = [min(REL(:)) max(REL(:))]; yl = yl + 0.08*diff(yl)*[-1 1];

f = figure('Color','w','Position',[80 80 1200 380]);
for i = 1:numel(S.T_call)
    subplot(1,3,i); hold on; grid on; box on
    for j = 1:NM5
        plot(S.Kf, squeeze(REL(i,:,j)), ['-' MK5{j}], 'Color', CC5{j}, ...
            'MarkerSize',4,'MarkerFaceColor','w','LineWidth',1.3);
    end
    yline(0,'k-'); ylim(yl);
    title(sprintf('$T = %.2f$', S.T_call(i)),'Interpreter','latex');
    xlabel('$K/F$','Interpreter','latex');
    ylabel('deviation from cross-model mean (\%)','Interpreter','latex');
    if i==1
        legend(MOD5,'Location','northwest','Interpreter','none', ...
            'Box','off','FontSize',7);
    end
end
Esporta(f, figdir, 'fig_calls_rel');

% digitali: livello e dispersione
% A sinistra le probabilita' (le cinque curve sono visivamente sovrapposte:
% e' il punto della figura).
DISP = nan(numel(S.T_dig), numel(S.Kf_d));
for i = 1:numel(S.T_dig)
    for k = 1:numel(S.Kf_d)
        p = squeeze(S.Digitals(i,k,:));
        DISP(i,k) = 100*(max(p)-min(p))/mean(p);
    end
end
LS = {'-','--'};

f = figure('Color','w','Position',[80 80 1100 400]);
subplot(1,2,1); hold on; grid on; box on
for i = 1:numel(S.T_dig)
    for j = 1:NM5
        plot(S.Kf_d, squeeze(S.Digitals(i,:,j)), LS{i}, 'Color', CC5{j}, ...
            'LineWidth',1.3);
    end
end
xlabel('$K/F$','Interpreter','latex');
ylabel('$Q(S_T > K)$','Interpreter','latex');  % \mathbb non e' supportato da MATLAB
title('Digital prices','Interpreter','none');
hh = gobjects(NM5+2,1);
for j = 1:NM5, hh(j) = plot(nan,nan,'-','Color',CC5{j},'LineWidth',2); end
hh(NM5+1) = plot(nan,nan,'-','Color',[0.4 0.4 0.4],'LineWidth',1.3);
hh(NM5+2) = plot(nan,nan,'--','Color',[0.4 0.4 0.4],'LineWidth',1.3);
legend(hh,[MOD5, arrayfun(@(t) sprintf('T = %.2f',t), S.T_dig(:)', ...
    'UniformOutput', false)],'Location','northeast', ...
    'Interpreter','none','Box','off','FontSize',7);

subplot(1,2,2); hold on; grid on; box on
for i = 1:numel(S.T_dig)
    plot(S.Kf_d, DISP(i,:), ['-' MK5{i}], 'Color', [0.15 0.15 0.15], ...
        'LineStyle', LS{i}, 'MarkerSize',5,'MarkerFaceColor','w','LineWidth',1.4);
end
set(gca,'YScale','log');
xlabel('$K/F$','Interpreter','latex');
ylabel('dispersion across models, range/mean (\%)','Interpreter','latex');
title('Dispersion','Interpreter','none');
legend(arrayfun(@(t) sprintf('T = %.2f',t), S.T_dig, 'UniformOutput', false), ...
    'Location','northwest','Interpreter','none','Box','off','FontSize',8);
Esporta(f, figdir, 'fig_digitali');

% struttura a termine dello skew ATM, con i fit a legge di potenza
% Il fit e' ristretto a T <= 0.25 perche' gamma = 1/2 - H e' asintotica per
% T -> 0. La banda oltre T = 1.59 marca la zona in cui xi0 e' la cubica
% terminale dell'interpolante, non un dato.
T_FIT = 0.25;
T_ANC = 1.59;
SL = nan(NM5,1);
mfit = S.T_skew <= T_FIT;
for j = 1:NM5
    cfit = polyfit(log(S.T_skew(mfit)), log(abs(S.SkewATM(mfit,j))), 1);
    SL(j) = cfit(1);
end

f = figure('Color','w','Position',[100 100 780 460]);
hold on; grid on; box on
yy = abs(S.SkewATM(:)); yl2 = [0.8*min(yy) 1.25*max(yy)];
patch([T_ANC max(S.T_skew) max(S.T_skew) T_ANC], ...
      [yl2(1) yl2(1) yl2(2) yl2(2)], [0.85 0.85 0.85], ...
      'EdgeColor','none','FaceAlpha',0.45);
hl = gobjects(NM5,1);
for j = 1:NM5
    hl(j) = plot(S.T_skew, abs(S.SkewATM(:,j)), ['-' MK5{j}], ...
        'Color', CC5{j}, 'MarkerSize',4,'MarkerFaceColor','w','LineWidth',1.3);
end
for j = 1:NM5
    cfit = polyfit(log(S.T_skew(mfit)), log(abs(S.SkewATM(mfit,j))), 1);
    tgd  = linspace(min(S.T_skew), T_FIT, 50);
    plot(tgd, exp(polyval(cfit, log(tgd))), '--', 'Color', CC5{j}, 'LineWidth',1.0);
end
set(gca,'XScale','log','YScale','log'); ylim(yl2);
xlabel('$T$ (years)','Interpreter','latex');
ylabel('$|\partial \sigma_{BS}/\partial k|$ ATM','Interpreter','latex');
legend(hl, arrayfun(@(j) sprintf('%s  (slope %.3f)', MOD5{j}, SL(j)), ...
    1:NM5, 'UniformOutput', false), ...
    'Location','southwest','Interpreter','none','Box','off','FontSize',8);
Esporta(f, figdir, 'fig_skew_term_ann');

% premio del variance swap sul log-contract, con la forma chiusa
RATIO  = squeeze(S.ResVS(:,:,3));
LAMBDA = C.rhh.theta_ottimi(5);
VSL_CF = (1+2*LAMBDA^2)/(1 + 2*LAMBDA^2/(1+LAMBDA));

f = figure('Color','w','Position',[100 100 820 420]);
b = bar(RATIO); grid on; box on
for j = 1:NM5, b(j).FaceColor = CC5{j}; b(j).EdgeColor = 'none'; end
yline(1,'k--','LineWidth',1.2);
yline(VSL_CF,':','Color',CC5{5},'LineWidth',1.4);
set(gca,'XTickLabel',arrayfun(@(t) sprintf('$T=%.3f$',t), S.T_vs, ...
    'UniformOutput', false),'TickLabelInterpreter','latex');
ylabel('$\mathrm{VS}/L$','Interpreter','latex');
ylim([0.95 1.05*max(RATIO(:))]);
text(0.55, VSL_CF, sprintf(' closed form %.6f', VSL_CF), ...
    'FontSize',7,'VerticalAlignment','bottom');
legend(MOD5,'Location','northwest','Interpreter','none','Box','off');
Esporta(f, figdir, 'fig_vs_logcontract_ann');

% TABELLE LATEX
% Un file per tabella: nel documento si scrive
fid = ApriTab(figdir,'tab_parametri');
fprintf(fid,'\\begin{tabular}{llr}\n\\toprule\nModel & Parameter & Estimate\\\\\n\\midrule\n');
Scriv(fid,'Heston',{'$\kappa$','$\bar v$','$\eta$','$v_0$','$\rho$'},C.hes.theta_ottimi);
Scriv(fid,'Heston-curve ($H=1/2$)',{'$\nu$','$\rho$'},C.hcur.theta_ottimi);
Scriv(fid,'Heston-curve-$\kappa$ ($H=1/2$)',{'$\kappa$','$\nu$','$\rho$'},C.hck.theta_ottimi);
Scriv(fid,'rough Heston',{'$H$','$\nu$','$\rho$'},C.rh.theta_ottimi);
Scriv(fid,'CGMY',{'$C$','$G$ (fixed)','$M$','$Y$'},C.cg.theta_ottimi);
Scriv(fid,'rough Heston + CGMY',{'$H$','$\nu$','$\rho$','$C$','$G$ (fixed)','$M$','$Y$'},t7);
Scriv(fid,'rough Hawkes Heston',{'$\alpha$','$b$','$c$','$\rho$','$\Lambda$','$\sigma_0^2$','$\beta$'},C.rhh.theta_ottimi);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

fid = ApriTab(figdir,'tab_errori');
fprintf(fid,'\\begin{tabular}{lrrr}\n\\toprule\n');
fprintf(fid,'Model & Params & RMSE rel. (\\%%) & RMSE abs. (\\%%)\\\\\n\\midrule\n');
for j=ORD
    e=IVmod(:,j)-IVmkt; ok=~isnan(e);
    fprintf(fid,'%s & %d & %.3f & %.3f\\\\\n', MODTEX{j}, NPAR(j), ...
        RMSE_tot(j), 100*sqrt(mean(e(ok).^2)));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

fid = ApriTab(figdir,'tab_rmse_scadenza');
fprintf(fid,'\\begin{tabular}{r%s}\n\\toprule\n$T$ (yr)',repmat('r',1,NM));
for j=ORD, fprintf(fid,' & %s',MODTEX{j}); end
fprintf(fid,'\\\\\n\\midrule\n');
for i=1:numel(T_uniq)
    fprintf(fid,'%.4f',T_uniq(i));
    for j=ORD, fprintf(fid,' & %.3f',RM(i,j)); end
    fprintf(fid,'\\\\\n');
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% tabella dell'effetto della mean reversion (kappa)
% Tre famiglie, ciascuna con kappa = 0 e kappa libero, piu' il rough Hawkes
% (che ha mean reversion b).
f_rk = fullfile(resdir,'calib_rhestonkappa_7x15z.mat');
f_ck = fullfile(resdir,'calib_rhcgmykappa_7x15z.mat');
if isfile(f_rk) && isfile(f_ck)
    RK = load(f_rk,'theta_ottimi','TH_K0','RMSE_K0');
    CK = load(f_ck,'theta_ottimi');
    rmse_rel = @(iv) 100*sqrt(mean(((IVmkt - iv)./IVmkt).^2,'omitnan'));
    [xc,pw] = XcPow('rHeston','exact');
    IVrk = IVdaModello(@(Ft,Tt,Kt,CPt) PrezzaSuperficieSINC_yield(Ft,Tt,Kt, ...
        zeros(size(Kt)),zeros(size(Kt)),'rough Heston kappa',RK.theta_ottimi, ...
        RiccatiXi(xifn,400,Tt(1)),400,xc,0,pw,CPt), F,T,K,CP,T_uniq);
    t8  = CK.theta_ottimi;  % [H,nu,rho,C,G,M,Y,kappa]
    Lck = ContributiCGMY(t8(4),t8(5),t8(6),t8(7));
    pck = [t8(1:3), t8(8), t8(4:7)];
    [xc,pw] = XcPow('rH+CGMY','exact');
    IVck = IVdaModello(@(Ft,Tt,Kt,CPt) PrezzaSuperficieSINC_yield(Ft,Tt,Kt, ...
        zeros(size(Kt)),zeros(size(Kt)),'rough Heston kappa CGMY',pck, ...
        max(RiccatiXi(xifn,400,Tt(1))-Lck,1e-6),400,xc,0,pw,CPt), F,T,K,CP,T_uniq);
    RMSE_rk = rmse_rel(IVrk);  RMSE_ck = rmse_rel(IVck);
    fprintf('\nEffetto di kappa (RMSE rel.):\n');
    fprintf('  Heston-curve      %.3f%%  -> con kappa %.3f%%\n', RMSE_tot(2), RMSE_tot(7));
    fprintf('  rough Heston      %.3f%%  -> con kappa %.3f%%  (kappa = 0 stesso solver: %.3f%%)\n', ...
        RMSE_tot(3), RMSE_rk, RK.RMSE_K0);
    fprintf('  rough Heston+CGMY %.3f%%  -> con kappa %.3f%%\n', RMSE_tot(5), RMSE_ck);

    fid = ApriTab(figdir,'tab_kappa');
    fprintf(fid,'\\begin{tabular}{lrrrr}\n\\toprule\n');
    fprintf(fid,'Model & Params & $H$ & $\\kappa$ & RMSE rel. (\\%%)\\\\\n\\midrule\n');
    fprintf(fid,'Heston-curve & 2 & 0.5 & 0 & %.3f\\\\\n', RMSE_tot(2));
    fprintf(fid,'Heston-curve-$\\kappa$ & 3 & 0.5 & %.3f & %.3f\\\\\n', C.hck.theta_ottimi(1), RMSE_tot(7));
    fprintf(fid,'\\addlinespace\n');
    fprintf(fid,'rough Heston & 3 & %.4f & 0 & %.3f\\\\\n', C.rh.theta_ottimi(1), RMSE_tot(3));
    fprintf(fid,'rough Heston-$\\kappa$ & 4 & %.4f & %.3f & %.3f\\\\\n', RK.theta_ottimi(1), RK.theta_ottimi(4), RMSE_rk);
    fprintf(fid,'\\addlinespace\n');
    fprintf(fid,'rough Heston + CGMY & 6 & %.4f & 0 & %.3f\\\\\n', t7(1), RMSE_tot(5));
    fprintf(fid,'rough Heston-$\\kappa$ + CGMY & 7 & %.4f & %.3f & %.3f\\\\\n', t8(1), t8(8), RMSE_ck);
    fprintf(fid,'\\addlinespace\n');
    fprintf(fid,'rough Hawkes Heston ($\\kappa=-b$) & 7 & %.4f & %.3f & %.3f\\\\\n', ...
        C.rhh.theta_ottimi(1)-0.5, -C.rhh.theta_ottimi(2), RMSE_tot(6));
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);
else
    warning('Mancano calib_rhestonkappa / calib_rhcgmykappa: tab_kappa saltata.');
end

if isfile(fullfile(resdir,'diagnostica_Xc.mat'))
    D=load(fullfile(resdir,'diagnostica_Xc.mat'));
    fid = ApriTab(figdir,'tab_dominio');
    fprintf(fid,'\\begin{tabular}{lrrr}\n\\toprule\n');
    fprintf(fid,'Model & $X_c$ & $N$ & Error at $X_c=12$ (bp)\\\\\n\\midrule\n');
    nn=fieldnames(D.RES);
        NOMI_DOM = struct('Heston','Heston','rHeston','rough Heston','CGMY','CGMY', ...
                      'rHCGMY','rough Heston + CGMY','rHawkesH','rough Hawkes Heston');
    for j=1:numel(nn)
        R=D.RES.(nn{j});
        if isfield(NOMI_DOM,nn{j}), nome=NOMI_DOM.(nn{j}); else, nome=nn{j}; end
        fprintf(fid,'%s & %g & $2^{%d}$ & %.2f\\\\\n', nome, max(R.Xc), ...
            max(R.pow), max(R.dIV_usato_bp));
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);
end

fid = ApriTab(figdir,'tab_incertezza_H');
fprintf(fid,'\\begin{tabular}{lc}\n\\toprule\nSource & Range of $H$\\\\\n\\midrule\n');
fprintf(fid,'Optimiser (25 Sobol restarts) & factor $1.00$\\\\\n');
if isfile(fullfile(resdir,'profile_G.mat'))
    Pg=load(fullfile(resdir,'profile_G.mat'));
    Hg = Pg.ResB.theta(Pg.ResB.G <= 1, 1);
    fprintf(fid,'Tempering convention $G$ & factor $%.2f$\\\\\n', max(Hg)/min(Hg));
end
if isfile(fullfile(resdir,'sensibilita_Vfair.mat'))
    Vs=load(fullfile(resdir,'sensibilita_Vfair.mat'));
    estr = cellfun(@(s) s.estr, Vs.V, 'UniformOutput', false);
    i_flat = find(strcmp(estr,'flat'),1);
    i_lin  = find(strcmp(estr,'linvar'),1);
    fprintf(fid,'Wing extrapolation & $[%.3f,\\,%.3f]$, factor $%.2f$\\\\\n', ...
        Vs.TH(i_flat,1), Vs.TH(i_lin,1), Vs.TH(i_lin,1)/Vs.TH(i_flat,1));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% tabelle del capitolo 7 intestazioni su due righe, cosi' le tabelle larghe
% restano nella pagina
MOD5H = {'Heston','\shortstack[r]{rough\\Heston}','CGMY', ...
         '\shortstack[r]{rough Heston\\+ CGMY}','\shortstack[r]{rough Hawkes\\Heston}'};
fid = ApriTab(figdir,'tab_calls');
fprintf(fid,'\\begin{tabular}{rr%s}\n\\toprule\n$T$ & $K/F$',repmat('r',1,NM5));
for j=1:NM5, fprintf(fid,' & %s',MOD5H{j}); end
fprintf(fid,'\\\\\n\\midrule\n');
for i = 1:numel(S.T_call)
    for k = 1:numel(S.Kf)
        if k==1, et=sprintf('%.2f',S.T_call(i)); else, et=''; end
        fprintf(fid,'%s & %.2f',et,S.Kf(k));
        for j=1:NM5, fprintf(fid,' & %.4f',S.CallPrices(i,k,j)); end
        fprintf(fid,'\\\\\n');
    end
    if i<numel(S.T_call), fprintf(fid,'\\addlinespace\n'); end
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

fid = ApriTab(figdir,'tab_digitali');
fprintf(fid,'\\begin{tabular}{rr%s}\n\\toprule\n$T$ & $K/F$',repmat('r',1,NM5));
for j=1:NM5, fprintf(fid,' & %s',MOD5H{j}); end
fprintf(fid,'\\\\\n\\midrule\n');
for i = 1:numel(S.T_dig)
    for k = 1:numel(S.Kf_d)
        if k==1, et=sprintf('%.2f',S.T_dig(i)); else, et=''; end
        fprintf(fid,'%s & %.2f',et,S.Kf_d(k));
        for j=1:NM5, fprintf(fid,' & %.4f',S.Digitals(i,k,j)); end
        fprintf(fid,'\\\\\n');
    end
    if i<numel(S.T_dig), fprintf(fid,'\\addlinespace\n'); end
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% psi ~ T^(-gamma) con gamma = 1/2 - H: la pendenza di log|psi| su log T e'
% H - 1/2, quindi l'H implicito e' slope + 1/2.
H_cal = [NaN, C.rh.theta_ottimi(1), NaN, C.rhc.theta_ottimi(1), C.rhh.theta_ottimi(1)-0.5];
fid = ApriTab(figdir,'tab_skew');
fprintf(fid,'\\begin{tabular}{lrrl}\n\\toprule\n');
fprintf(fid,'Model & Slope & Implied $H$ & Calibrated\\\\\n\\midrule\n');
for j = 1:NM5
    if isnan(H_cal(j)), cal = '---';
    elseif j==5,        cal = sprintf('$H_{\\mathrm{eff}} = %.3f$', H_cal(j));
    else,               cal = sprintf('$H = %.4f$', H_cal(j));
    end
    fprintf(fid,'%s & %.3f & %.3f & %s\\\\\n', MOD5{j}, SL(j), SL(j)+0.5, cal);
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% ResVS(:,:,1) = L e (:,:,2) = VS sono VARIANZE: in tabella vanno le radici.
fid = ApriTab(figdir,'tab_vsl');
fprintf(fid,'\\begin{tabular}{lrrrr}\n\\toprule\n');
fprintf(fid,'Model & $T$ & $\\sqrt{L}$ & $\\sqrt{\\mathrm{VS}}$ & $\\mathrm{VS}/L$\\\\\n');
fprintf(fid,'\\midrule\n');
for j = 1:NM5
    for i = 1:numel(S.T_vs)
        if i==1, et=MOD5{j}; else, et=''; end
        fprintf(fid,'%s & %.4f & %.4f & %.4f & %.6f\\\\\n', et, S.T_vs(i), ...
            sqrt(max(S.ResVS(i,j,1),0)), sqrt(max(S.ResVS(i,j,2),0)), S.ResVS(i,j,3));
    end
    if j<NM5, fprintf(fid,'\\addlinespace\n'); end
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% TABELLE DELL'APPENDICE C (validazioni numeriche)
% Legge Results/valida_pricer.mat e Results/multistart.mat, gia' prodotti da
% Valida_Pricer e Multistart dentro Run_All: qui non si ricalcola nulla, si
% mettono in LaTeX i numeri dei test.
Wv = load(fullfile(resdir,'valida_pricer.mat'));
VA = Wv.VALIDA;

% C.1.1  SINC contro l'inversione di Lewis 24 combinazioni strike/scadenza
% sotto Heston.
fid = ApriTab(figdir,'tab_app_lewis');
fprintf(fid,'\\begin{tabular}{rrrrr}\n\\toprule\n');
fprintf(fid,'$T$ & $K/F$ & Lewis & SINC & Rel. error\\\\\n\\midrule\n');
for i = 1:numel(VA.A.T)
    fprintf(fid,'%.4f & %.2f & %.8f & %.8f & %.2e\\\\\n', ...
        VA.A.T(i), VA.A.KF(i), VA.A.lewis(i), VA.A.sinc(i), VA.A.err_rel(i));
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{4}{l}{Maximum relative error} & %.2e\\\\\n', VA.A.err_rel_max);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.1.2  convergenza interna in N
% Non dipende da alcun riferimento esterno: se pow=14 e pow=16 coincidono,
% la configurazione accurata e' convergente.
fid = ApriTab(figdir,'tab_app_convergenza');
fprintf(fid,'\\begin{tabular}{rrrrrr}\n\\toprule\n');
fprintf(fid,'$T$ & $K/F$ & $N=2^{12}$ & $N=2^{14}$ & $N=2^{16}$ & $|2^{12}-2^{16}|$\\\\\n');
fprintf(fid,'\\midrule\n');
for i = 1:numel(VA.Abis.T)
    fprintf(fid,'%.4f & %.2f & %.8f & %.8f & %.8f & %.2e\\\\\n', ...
        VA.Abis.T(i), VA.Abis.KF(i), VA.Abis.pow12(i), VA.Abis.pow14(i), ...
        VA.Abis.pow16(i), abs(VA.Abis.pow12(i)-VA.Abis.pow16(i)));
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{5}{l}{Maximum difference} & %.2e\\\\\n', VA.Abis.max_diff_12_16);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.1.3  errore della configurazione coarse
% Una cella "n.a." segnala un prezzo non positivo: quella configurazione non
% e' utilizzabile in una loss pesata, i residui la' sono artefatti.
fid = ApriTab(figdir,'tab_app_coarse');
ncfg = size(VA.B.cfg,1);
fprintf(fid,'\\begin{tabular}{rr%s}\n\\toprule\n',repmat('r',1,ncfg));
fprintf(fid,'$T$ & $K/F$');
for j = 1:ncfg
    fprintf(fid,' & $X_c=%g$, $2^{%g}$', VA.B.cfg(j,1), VA.B.cfg(j,2));
end
fprintf(fid,'\\\\\n\\midrule\n');
for i = 1:numel(VA.B.T)
    fprintf(fid,'%.4f & %.2f', VA.B.T(i), VA.B.KF(i));
    for j = 1:ncfg
        if isnan(VA.B.dIV_bp(i,j)), fprintf(fid,' & n.a.');
        else,                       fprintf(fid,' & %.1f', VA.B.dIV_bp(i,j)); end
    end
    fprintf(fid,'\\\\\n');
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.1.4  Pade al bordo H = 1/2
% Questo NON e' un test di convergenza di griglia: il Pade converge al
% proprio valore, che ad alpha = 1 non e' la soluzione esatta.
fid = ApriTab(figdir,'tab_app_pade');
fprintf(fid,'\\begin{tabular}{rrrrr}\n\\toprule\n');
fprintf(fid,'$T$ & $K/F$ & Closed-form Heston & Pad\\''e $(3,3)$ & $\\Delta$IV (bp)\\\\\n');
fprintf(fid,'\\midrule\n');
for i = 1:numel(VA.C.T)
    fprintf(fid,'%.4f & %.2f & %.6f & %.6f & %+.1f\\\\\n', ...
        VA.C.T(i), VA.C.KF(i), VA.C.heston_chiuso(i), VA.C.pade(i), VA.C.dIV_bp(i));
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{4}{l}{Maximum absolute discrepancy} & %.1f\\\\\n', VA.C.max_abs_bp);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.1.5  sensibilita' del caso nested a H = 1/2 - eps
fid = ApriTab(figdir,'tab_app_eps');
fprintf(fid,'\\begin{tabular}{r%s}\n\\toprule\n$\\varepsilon$',repmat('r',1,numel(VA.F.T)));
for j = 1:numel(VA.F.T), fprintf(fid,' & ATM $T=%.3f$', VA.F.T(j)); end
fprintf(fid,'\\\\\n\\midrule\n');
for i = 1:numel(VA.F.eps)
    fprintf(fid,'%.0e', VA.F.eps(i));
    for j = 1:numel(VA.F.T), fprintf(fid,' & %.6f', VA.F.prezzi(i,j)); end
    fprintf(fid,'\\\\\n');
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.2.1  martingala e put-call parity
% Per Heston la martingala e' esatta; per il ramo Pade il residuo e' un
% limite e scala linearmente con l'offset usato per valutare phi(-i).
fid = ApriTab(figdir,'tab_app_martingala');
fprintf(fid,'\\begin{tabular}{lrrr}\n\\toprule\n');
fprintf(fid,'Model & $T$ & $|\\varphi_T(-i)-1|$ & $\\max|C-P-(F-K)|$\\\\\n\\midrule\n');
nomiE = cellstr(VA.E.modello);
for i = 1:numel(VA.E.T)
    fprintf(fid,'%s & %.4f & %.2e & %.2e\\\\\n', ...
        nomiE{i}, VA.E.T(i), VA.E.martingala(i), VA.E.parity(i));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.2.2  convenzione della griglia di Riccati
% Con i punti medi l'errore deve scalare come 2^(1+alpha): un rapporto
% vicino a 2 segnalerebbe griglie disallineate fra routine diverse.
fid = ApriTab(figdir,'tab_app_griglia');
fprintf(fid,'\\begin{tabular}{rrrr}\n\\toprule\n');
fprintf(fid,'$T$ & $N$ & ATM price & $|p(N)-p(2N)|$\\\\\n\\midrule\n');
for r = 1:size(VA.G2.prezzi,1)
    for j = 1:numel(VA.G2.N)
        if j==1, et=sprintf('%.4f',VA.G2.T(r)); else, et=''; end
        if j < numel(VA.G2.N)
            fprintf(fid,'%s & %d & %.8f & %.3e\\\\\n', et, VA.G2.N(j), ...
                VA.G2.prezzi(r,j), abs(VA.G2.prezzi(r,j)-VA.G2.prezzi(r,j+1)));
        else
            fprintf(fid,'%s & %d & %.8f & ---\\\\\n', et, VA.G2.N(j), VA.G2.prezzi(r,j));
        end
    end
    fprintf(fid,'\\multicolumn{3}{l}{measured ratio} & %.2f\\\\\n', VA.G2.rapporto_misurato(r));
    if r < size(VA.G2.prezzi,1), fprintf(fid,'\\addlinespace\n'); end
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{3}{l}{expected ratio $2^{1+\\alpha}$, $\\alpha=%.4f$} & %.2f\\\\\n', ...
    VA.G2.alpha, VA.G2.rapporto_atteso);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.3.1  variance swap estratte
fid = ApriTab(figdir,'tab_app_kvar');
fprintf(fid,'\\begin{tabular}{rrr}\n\\toprule\n');
fprintf(fid,'$T$ (yr) & $\\mathcal{K}_{\\mathrm{var}}(T)$ & $\\sqrt{\\mathcal{K}_{\\mathrm{var}}}$ (\\%%)\\\\\n');
fprintf(fid,'\\midrule\n');
for i = 1:numel(VA.T_uniq)
    fprintf(fid,'%.4f & %.5f & %.3f\\\\\n', VA.T_uniq(i), VA.V_fair(i), ...
        100*sqrt(VA.V_fair(i)));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.3.2  curva PCHIP contro livelli flat
% Stessa varianza totale a ogni scadenza, forma diversa: la tabella misura
% quanto vale la differenza in punti base di IV.
fid = ApriTab(figdir,'tab_app_pchipflat');
fprintf(fid,'\\begin{tabular}{rrrrr}\n\\toprule\n');
fprintf(fid,'$T$ & $K/F$ & IV pchip (\\%%) & IV flat (\\%%) & $\\Delta$IV (bp)\\\\\n');
fprintf(fid,'\\midrule\n');
for i = 1:numel(VA.D.T)
    fprintf(fid,'%.4f & %.2f & %.3f & %.3f & %+.1f\\\\\n', ...
        VA.D.T(i), VA.D.KF(i), 100*VA.D.iv_pchip(i), 100*VA.D.iv_flat(i), VA.D.dIV_bp(i));
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{4}{l}{Maximum absolute discrepancy} & %.1f\\\\\n', VA.D.max_abs_bp);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.4  multistart
% Una riga per parametro: quante partenze, quanti bacini, e l'escursione del
% parametro fra i restart entro il bacino.
if isfile(fullfile(resdir,'multistart.mat'))
    Ms  = load(fullfile(resdir,'multistart.mat'));
    mm  = {'rheston','rhcgmy','cgmy'};
    eti = {'rough Heston','rough Heston + CGMY','CGMY'};
    fid = ApriTab(figdir,'tab_app_multistart');
    fprintf(fid,'\\begin{tabular}{llrrl}\n\\toprule\n');
    fprintf(fid,'Model & Parameter & Restarts & Basins & Range across restarts\\\\\n');
    fprintf(fid,'\\midrule\n');
    for k = 1:numel(mm)
        Sm = Ms.OUT.(mm{k});
        for j = 1:numel(Sm.nomi)
            if j==1
                et1 = eti{k};
                et2 = sprintf('%d', numel(Sm.loss));
                et3 = sprintf('%d', Sm.n_bacini);
            else
                et1=''; et2=''; et3='';
            end
            fprintf(fid,'%s & $%s$ & %s & %s & $[%.6g,\\; %.6g]$\\\\\n', ...
                et1, Sm.nomi{j}, et2, et3, min(Sm.theta(:,j)), max(Sm.theta(:,j)));
        end
        if k < numel(mm), fprintf(fid,'\\addlinespace\n'); end
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);
end

% controlli dell'appendice
fprintf('\nControlli appendice C (devono coincidere col testo):\n');
fprintf('  err. rel. max SINC vs Lewis %.2e  (atteso 1.25e-04)\n', VA.A.err_rel_max);
fprintf('  Pade al bordo, max |dIV| %.1f bp  (atteso 65.5)\n', VA.C.max_abs_bp);
fprintf('  pchip vs flat, max |dIV| %.1f bp  (atteso 32.0)\n', VA.D.max_abs_bp);
fprintf('  rapporto griglia Riccati %.2f / %.2f  (atteso %.2f)\n', ...
    VA.G2.rapporto_misurato(1), VA.G2.rapporto_misurato(2), VA.G2.rapporto_atteso);

% controlli contro i numeri citati nel capitolo 7
fprintf('\nControlli capitolo 7 (devono coincidere col testo):\n');
fprintf('  call T=0.5 K/F=1.10  Heston %.4f  rH+CGMY %.4f  (atteso 0.8765 / 1.1807)\n', ...
    S.CallPrices(2,6,1), S.CallPrices(2,6,4));
fprintf('  digitale T=0.25 K/F=1.10  Heston %.4f  rH+CGMY %.4f  (atteso 0.0583 / 0.0953)\n', ...
    S.Digitals(2,7,1), S.Digitals(2,7,4));
fprintf('  slope skew rough Heston %.3f  (atteso -0.517)\n', SL(2));
fprintf('  VS/L  CGMY %.3f   rHH %.6f   forma chiusa %.6f\n', ...
    S.ResVS(1,3,3), S.ResVS(1,5,3), VSL_CF);

fprintf('\nFigure (pdf + png) e tabelle in %s\n', figdir);
fprintf('In LaTeX, con \\graphicspath{{./Images/}}:\n');
fprintf('  \\includegraphics{fig_nome}   (senza cartella e senza estensione)\n');
fprintf('  \\input{Images/tab_nome.tex}  (su \\input il graphicspath non vale)\n');
toc

% FUNZIONI LOCALI
function Esporta(f, cartella, nome)
% Salva ogni figura solo in PDF vettoriale: nitido a ogni zoom e in stampa,
% ed e' il formato che pdfLaTeX preferisce.
    set(findall(f,'-property','FontSize'),'FontSize',9);
    try
        exportgraphics(f, fullfile(cartella,[nome '.pdf']),'ContentType','vector');
        fprintf('  %s.pdf\n', nome);
    catch
        % exportgraphics richiede R2020a: ripiego su print.
        set(f,'PaperPositionMode','auto');
        print(f, fullfile(cartella,[nome '.pdf']), '-dpdf', '-painters', '-bestfit');
        fprintf('  %s.pdf  (via print)\n', nome);
    end
    close(f);  % altrimenti restano aperte venti finestre
end

function fid = ApriTab(cartella, nome)
% Apre un file tabella in scrittura (tronca il precedente) e ci mette
% l'intestazione che ricorda di non modificarlo a mano.
    fid = fopen(fullfile(cartella,[nome '.tex']),'w');
    fprintf(fid,'%% Generato da Figure_Tesi.m il %s -- non modificare a mano\n', ...
        datestr(now));
end

function Scriv(fid, modello, nomi, val)
    for i=1:numel(nomi)
        if i==1, et=modello; else, et=''; end
        fprintf(fid,'%s & %s & %.5g\\\\\n', et, nomi{i}, val(i));
    end
    fprintf(fid,'\\addlinespace\n');
end

function xg = RiccatiXi(xifn, N, T)
    ti = ((0:N-1)+0.5)/N*T;  xg = xifn(ti);  xg = xg(:)';
end

function IV = IVdaModello(fn, F, T, K, CP, T_uniq)
    P = zeros(numel(K),1);
    for k = 1:numel(T_uniq)
        idx = (T == T_uniq(k));
        if ~any(idx), continue; end
        P(idx) = fn(F(find(idx,1)), T_uniq(k)*ones(sum(idx),1), K(idx), CP(find(idx,1)));
    end
    IV = InvertiIV_OTM(F, K, T, max(P,1e-6));
end

function [Lcg, Vcg] = ContributiCGMY(C, G, M, Y)
    Lcg = 2*C*( gamma(-Y)*((M-1)^Y - M^Y + (G+1)^Y - G^Y) ...
              - gamma(1-Y)*(M^(Y-1) - G^(Y-1)) );
    Vcg = C * gamma(2-Y) * (M^(Y-2) + G^(Y-2));
end

function xi = eval_xi0_pchip(tvec, W_pp)
    h=1e-5; tvec=tvec(:);
    tp=max(tvec+h,0); tm=max(tvec-h,0);
    dn=tp-tm; dn(dn<1e-12)=1e-12;
    xi=max((ppval(W_pp,tp)-ppval(W_pp,tm))./dn,1e-8);
end

function Kvar = VFairVariante(o, T_s, K_s, F_s, IV_M, Bid_s, T_uniq)
% Replica di Demeterfi con ipotesi parametrizzate su troncamento,
% estrapolazione, risoluzione, ampiezza dello smile e filtro sul bid.
    Kvar = zeros(numel(T_uniq),1);
    for i = 1:numel(T_uniq)
        Ti = T_uniq(i);
        m  = (T_s == Ti) & (abs(log(K_s./F_s)) < o.logmcut) & (Bid_s >= o.bidmin);
        Fi = F_s(find(T_s == Ti, 1));
        [Ks, s] = sort(K_s(m)); IVs = IV_M(m); IVs = IVs(s);
        if numel(Ks) < 5
            warning('T=%.4f: solo %d quote dopo i filtri della variante.', Ti, numel(Ks));
            Kvar(i) = NaN; continue
        end
        if strcmp(o.estr, 'nessuna')
            lo = Ks(1)/Fi;  hi = Ks(end)/Fi;  % integra solo dove si quota
        else
            lo = o.lo;      hi = o.hi;
        end
        Kg  = logspace(log10(lo*Fi), log10(hi*Fi), o.ngrid)';
        IVg = interp1(Ks, IVs, Kg, 'linear', NaN);
        switch o.estr
            case {'flat','nessuna'}
                IVg(Kg < Ks(1))   = IVs(1);
                IVg(Kg > Ks(end)) = IVs(end);
            case 'linvar'
                % varianza totale w = IV^2 estesa linearmente in k con la
                % pendenza ai due bordi dello smile quotato (floor a 1e-6).
                kk = log(Kg/Fi); k0 = log(Ks(1)/Fi); k1 = log(Ks(end)/Fi);
                sl = (IVs(2)^2   - IVs(1)^2)   / (log(Ks(2)/Fi)     - k0);
                sr = (IVs(end)^2 - IVs(end-1)^2)/(k1 - log(Ks(end-1)/Fi));
                L = Kg < Ks(1); R = Kg > Ks(end);
                IVg(L) = sqrt(max(IVs(1)^2   + sl*(kk(L)-k0), 1e-6));
                IVg(R) = sqrt(max(IVs(end)^2 + sr*(kk(R)-k1), 1e-6));
        end
        OTM = zeros(size(Kg));
        for j = 1:numel(Kg)
            if Kg(j) < Fi, [~,p] = blkprice(Fi,Kg(j),0,Ti,IVg(j)); OTM(j)=p;
            else,          [c,~] = blkprice(Fi,Kg(j),0,Ti,IVg(j)); OTM(j)=c; end
        end
        Kvar(i) = (2/Ti) * trapz(Kg, OTM ./ Kg.^2);
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
