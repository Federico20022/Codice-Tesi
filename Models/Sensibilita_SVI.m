clearvars; clc;
tic
% SENSIBILITA_SVI  Ricalibra Heston-curve, rough Heston e rough Heston +
% CGMY con le tre regole sulle ali di xi0 (piatte, lineari, SSVI).
% Si parte dagli ottimi della tesi, lsqnonlin fast e poi exact.
% Output: Results/sensibilita_modelli.mat. Tempo: 30-60 minuti.

root = fileparts(mfilename('fullpath'));
if isempty(root), root = pwd; end
resdir = fullfile(root, '..', 'Results');
figdir = fullfile(root, '..', 'Figure');

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

% FIT SVI PER SCADENZA
nT  = numel(T_uniq);
SVI = nan(nT,5);  SVI_rmse_bp = nan(nT,1);
fprintf('\n%-9s %9s %8s %8s %8s %8s %10s %12s\n','T(yr)','a','b','rho','m','sig', ...
    'b(1+|rho|)','fit IV (bp)');
for i = 1:nT
    Ti = T_uniq(i); m = (T_s == Ti);
    Fi = F_s(find(m,1));
    [Ks, s] = sort(K_s(m)); IVs = IV_M(m); IVs = IVs(s);
    [SVI(i,:), SVI_rmse_bp(i)] = FitSVI(log(Ks/Fi), IVs, Ti);
    fprintf('%-9.4f %9.5f %8.4f %8.4f %8.4f %8.4f %10.3f %12.1f\n', Ti, SVI(i,:), ...
        SVI(i,2)*(1+abs(SVI(i,3))), SVI_rmse_bp(i));
end

fprintf(['\nNB: il fit SVI per singola scadenza e'' solo diagnostico (vedi\n' ...
         'intestazione): se rho > 0 o b(1+|rho|) grande, l''ala e'' non identificata.\n']);

% FIT SSVI CONGIUNTO (Gatheral-Jacquier 2014)
THETA = nan(nT,1); KQ = cell(nT,1); IVQ = cell(nT,1);
for i = 1:nT
    Ti = T_uniq(i); m = (T_s == Ti); Fi = F_s(find(m,1));
    [Ks, s] = sort(K_s(m)); IVs = IV_M(m); IVs = IVs(s);
    kq = log(Ks/Fi);
    THETA(i) = interp1(kq, IVs, 0, 'linear', 'extrap')^2 * Ti;  % w ATM di mercato
    KQ{i} = kq; IVQ{i} = IVs;
end
if any(diff(THETA) <= 0)
    warning('theta ATM non crescente in T: SSVI non e'' priva di calendar arbitrage.');
end
[SSVI, SSVI_rmse_bp] = FitSSVI(KQ, IVQ, THETA, T_uniq);
fprintf('\nSSVI: rho = %.4f  eta = %.4f  gamma = %.4f  | eta(1+|rho|) = %.3f (<= 2)\n', ...
    SSVI, SSVI(2)*(1+abs(SSVI(1))));
fprintf('%-9s %14s\n','T(yr)','fit IV (bp)');
for i = 1:nT, fprintf('%-9.4f %14.1f\n', T_uniq(i), SSVI_rmse_bp(i)); end

% TRE REGOLE DI ESTRAPOLAZIONE -> V_fair -> xi0 -> ricalibrazione rough
% Heston
V = {};
V{end+1} = struct('eti','piatta (base)','etiEN','flat','lo',0.05,'hi',2.50,'estr','flat','ngrid',4000,'logmcut',0.50,'bidmin',0);
V{end+1} = struct('eti','lineare in varianza','etiEN','linear in total variance','lo',0.05,'hi',2.50,'estr','linvar','ngrid',4000,'logmcut',0.50,'bidmin',0);
V{end+1} = struct('eti','SSVI','etiEN','arbitrage-free SSVI','lo',0.05,'hi',2.50,'estr','ssvi','ngrid',4000,'logmcut',0.50,'bidmin',0);
nV = numel(V);

KVAR = nan(nT, nV);
for iv = 1:nV
    KVAR(:,iv) = VFairSVI(V{iv}, T_s, K_s, F_s, IV_M, Bid_s, T_uniq, SSVI, THETA);
end
fprintf('\nsqrt(V_fair) in punti di volatilita''\n%-9s', 'T(yr)');
for iv = 1:nV, fprintf('%14s', V{iv}.eti); end, fprintf('\n');
for it = 1:nT
    fprintf('%-9.4f', T_uniq(it));
    for iv = 1:nV, fprintf('%14.3f', 100*sqrt(max(KVAR(it,iv),0))); end
    fprintf('\n');
end

% pesi del campione fisso (identici alla calibrazione)
Vega = zeros(numel(F),1);
for i = 1:numel(F)
    d1 = (log(F(i)/K(i))+0.5*IVmkt(i)^2*T(i))/(IVmkt(i)*sqrt(T(i)));
    Vega(i) = F(i)*sqrt(T(i))*normpdf(d1);
end
sqrt_w = 1./(max(Vega,1e-4).*max(IVmkt,0.01));
sqrt_w = min(sqrt_w, prctile(sqrt_w,99));

C_rh = load(fullfile(resdir,'calib_rheston_7x15z_curva.mat'),'theta_ottimi');
seed = C_rh.theta_ottimi(:)';
lb = [0.005, 0.05, -0.99];
ub = [0.49,  5.00,  0.10];
[XC_RH, POW_RH] = XcPow('rHeston', 'exact');
opt = optimoptions('lsqnonlin','Algorithm','trust-region-reflective', ...
    'MaxFunctionEvaluations',1500,'MaxIterations',100,'Display','off', ...
    'StepTolerance',1e-8,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);

% configurazioni dei tre modelli (identiche agli script di tesi)
% Heston-curve: Check_Heston_Curve.m usa fast 6/10, exact 12/12
XC_HC_F = 6;  POW_HC_F = 10;  XC_HC_E = 12;  POW_HC_E = 12;
C_hc  = load(fullfile(resdir,'calib_hestoncurve_7x15z.mat'),'theta_ottimi');
seed_hc = C_hc.theta_ottimi(:)';
lb_hc = [0.05, -0.99];  ub_hc = [5.00, 0.10];

% composto: Rough_Heston_CMGY_def.m, p = [H, nu, rho, q, M, Y], G = 0.05
G_FIX = 0.05;  Q_MAX = 0.98;
[XC_RC, POW_RC_F] = XcPow('rH+CGMY', 'fast');
[~,     POW_RC_E] = XcPow('rH+CGMY', 'exact');
lb_rc = [0.005, 0.05, -0.99, 1e-4,  5.0,   0.9 ];
ub_rc = [0.49,  5.00,  0.10, Q_MAX, 200.0, 1.95];
C_rc = load(fullfile(resdir,'calib_rhcgmy_7x15z.mat'),'theta_ottimi');
t7 = C_rc.theta_ottimi(:)';  % [H,nu,rho,C,G,M,Y]

% ciclo sulle tre curve
NOMI_MOD = {'Heston-curve','rough Heston','rH+CGMY'};
RES = struct();
for iv = 1:nV
    Wv = pchip([0; T_uniq(:)], [0; T_uniq(:).*KVAR(:,iv)]);
    xifn_v = @(tv) eval_xi0_pchip(tv, Wv);
    fprintf('\n===== Curva: %s =====\n', V{iv}.eti);

    % 1) Heston-curve
    rF = @(p) sqrt_w .* (PrezziHestonCurve(F,T,K,p,T_uniq,Wv,100,XC_HC_F,POW_HC_F,CP) - Mid);
    rE = @(p) sqrt_w .* (PrezziHestonCurve(F,T,K,p,T_uniq,Wv,400,XC_HC_E,POW_HC_E,CP) - Mid);
    th = lsqnonlin(rF, seed_hc, lb_hc, ub_hc, opt);
    th = lsqnonlin(rE, th,      lb_hc, ub_hc, opt);
    Pm = PrezziHestonCurve(F,T,K,th,T_uniq,Wv,400,XC_HC_E,POW_HC_E,CP);
    RES(iv).hc = Valuta(th, Pm, F, K, T, IVmkt, zz, T_uniq);
    fprintf('Heston-curve : nu %.4f rho %+.4f | RMSE %.3f%%\n', th, RES(iv).hc.rmse);

    % 2) rough Heston
    rF = @(p) sqrt_w .* (PrezziRoughHeston(F,T,K,p,T_uniq,xifn_v,100,XC_RH,POW_RH,CP) - Mid);
    rE = @(p) sqrt_w .* (PrezziRoughHeston(F,T,K,p,T_uniq,xifn_v,400,XC_RH,POW_RH,CP) - Mid);
    th = lsqnonlin(rF, seed, lb, ub, opt);
    th = lsqnonlin(rE, th,   lb, ub, opt);
    Pm = PrezziRoughHeston(F,T,K,th,T_uniq,xifn_v,400,XC_RH,POW_RH,CP);
    RES(iv).rh = Valuta(th, Pm, F, K, T, IVmkt, zz, T_uniq);
    fprintf('rough Heston : H %.4f nu %.4f rho %+.4f | RMSE %.3f%%\n', th, RES(iv).rh.rmse);

    % 3) rough Heston + CGMY xi_min dipende dalla curva: si ricalcola, e il
    % seed in q e' la stessa QUOTA q dell'ottimo di tesi (non lo stesso C).
    tt_all = [];
    for k = 1:nT
        tt_all = [tt_all, ((0:399)+0.5)/400*T_uniq(k), ((0:99)+0.5)/100*T_uniq(k)];
    end
    xi_min = min(xifn_v(tt_all));
    if iv == 1
        q_seed = t7(4) * LambdaCGMY(t7(5), t7(6), t7(7)) / xi_min;
    end
    seed_rc = min(max([t7(1:3), q_seed, t7(6), t7(7)], lb_rc), ub_rc);
    rF = @(p) sqrt_w .* (PrezziRoughHestonCGMY(F,T,K,p,G_FIX,xi_min,T_uniq,xifn_v, ...
                 100,XC_RC,POW_RC_F,CP) - Mid);
    rE = @(p) sqrt_w .* (PrezziRoughHestonCGMY(F,T,K,p,G_FIX,xi_min,T_uniq,xifn_v, ...
                 400,XC_RC,POW_RC_E,CP) - Mid);
    th = lsqnonlin(rF, seed_rc, lb_rc, ub_rc, opt);
    th = lsqnonlin(rE, th,      lb_rc, ub_rc, opt);
    Pm = PrezziRoughHestonCGMY(F,T,K,th,G_FIX,xi_min,T_uniq,xifn_v,400,XC_RC,POW_RC_E,CP);
    RES(iv).rc = Valuta(th, Pm, F, K, T, IVmkt, zz, T_uniq);
    Ccg = th(4)*xi_min/LambdaCGMY(G_FIX, th(5), th(6));
    [Lcg, Vcg] = ContributiCGMY(Ccg, G_FIX, th(5), th(6));
    [~, i1] = min(abs(T_uniq - 1));
    RES(iv).rc.quota_qv_1y = 100 * Vcg / (KVAR(i1,iv) - Lcg + Vcg);  % quota QV dei salti a 1 anno (14.6% in tesi)
    RES(iv).rc.vsl_jump    = Vcg / Lcg;  % VS/L della sola componente di salto
    fprintf('rH+CGMY      : H %.4f nu %.4f rho %+.4f q %.4f M %.2f Y %.4f | RMSE %.3f%% | VS/L salti %.3f\n', ...
        th, RES(iv).rc.rmse, RES(iv).rc.vsl_jump);
end

% tabelle riassuntive
fprintf('\n\n================ RIEPILOGO ================\n');
fprintf('%-24s %14s %14s %14s %12s %12s\n','curva','RMSE H-curve','RMSE rough H', ...
    'RMSE rH+CGMY','rid. rough','rid. salti');
for iv = 1:nV
    a = RES(iv).hc.rmse; b = RES(iv).rh.rmse; c = RES(iv).rc.rmse;
    fprintf('%-24s %13.2f%% %13.2f%% %13.2f%% %11.1f%% %11.1f%%\n', V{iv}.eti, a, b, c, ...
        100*(a-b)/a, 100*(b-c)/b);
end
fprintf('(rid. rough = riduzione dell''errore Heston-curve -> rough Heston;\n');
fprintf(' rid. salti = riduzione rough Heston -> rough Heston + CGMY)\n');

fprintf('\n%-24s %10s %10s %10s %12s %12s\n','curva','H rough H','H rH+CGMY','q','VS/L salti','quota QV 1y');
for iv = 1:nV
    fprintf('%-24s %10.4f %10.4f %10.4f %12.3f %11.1f%%\n', V{iv}.eti, RES(iv).rh.theta(1), ...
        RES(iv).rc.theta(1), RES(iv).rc.theta(4), RES(iv).rc.vsl_jump, RES(iv).rc.quota_qv_1y);
end

fprintf('\nRMSE per scadenza (%%)\n%-34s', 'curva / modello');
fprintf(' %7.4f', T_uniq); fprintf('\n');
campi = {'hc','rh','rc'};
for iv = 1:nV
    for jm = 1:3
        fprintf('%-34s', [V{iv}.eti ' / ' NOMI_MOD{jm}]);
        fprintf(' %7.2f', RES(iv).(campi{jm}).rmse_T); fprintf('\n');
    end
end

fprintf('\nBias relativo per bucket di z (%%, positivo = modello sotto il mercato)\n');
fprintf('%-34s %8s %8s %8s %8s %8s %8s\n','curva / modello','z<-2','-2..-1','-1..-.5','ATM','.5..1','z>1');
for iv = 1:nV
    for jm = 1:3
        fprintf('%-34s', [V{iv}.eti ' / ' NOMI_MOD{jm}]);
        fprintf(' %+8.2f', RES(iv).(campi{jm}).bias); fprintf('\n');
    end
end

save(fullfile(resdir,'sensibilita_modelli.mat'), 'RES','V','KVAR','T_uniq','NOMI_MOD');
fprintf('\nSalvato sensibilita_modelli.mat\n');
toc

% FUNZIONI LOCALI
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

function P_fwd = PrezziHestonCurve(F_vec, T_vec, K_vec, params, T_uniq, W_pp, ...
                                   N_riccati, Xc, pow, CP_vec)
% Copiata da Check_Heston_Curve.m: soluzione esatta della Riccati a alpha=1.
    params_cf = [params(1), params(2)];
    P_fwd = zeros(length(K_vec),1);
    for k = 1:length(T_uniq)
        T = T_uniq(k); idx = (T_vec == T);
        if ~any(idx), continue; end
        F_t = F_vec(find(idx,1)); K_t = K_vec(idx); CP_t = CP_vec(idx);
        ti = ((0:N_riccati-1) + 0.5) / N_riccati * T;
        xi_grid = eval_xi0_pchip(ti, W_pp); xi_grid = xi_grid(:)';
        P_fwd(idx) = PrezzaSuperficieSINC_yield(F_t, T*ones(size(K_t)), K_t, ...
            zeros(size(K_t)), zeros(size(K_t)), ...
            'Heston curve', params_cf, xi_grid, N_riccati, Xc, 0, pow, CP_t(1));
    end
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

function w = SviW(p, k)
    w = p(1) + p(2)*( p(3)*(k - p(4)) + sqrt((k - p(4)).^2 + p(5)^2) );
end

function w = SsviW(p, theta, k)
    rho = p(1); eta = p(2); gam = p(3);
    phi = eta / (theta^gam * (1 + theta)^(1 - gam));
    w = theta/2 * (1 + rho*phi*k + sqrt((phi*k + rho).^2 + 1 - rho^2));
end

function [p, rmse_bp] = FitSSVI(KQ, IVQ, THETA, T_uniq)
% Fit congiunto di (rho, eta, gam) sulle IV quotate di tutte le scadenze.
% Vincoli di assenza di arbitraggio statico: gam in (0, 1/2], eta*(1+|rho|)
% <= 2 (penalita').
    lb = [-0.999, 1e-3, 1e-3];
    ub = [ 0.999, 4.0,  0.5 ];
    x0s = [-0.7 1.0 0.4; -0.5 0.5 0.3; -0.9 1.0 0.5];
    opts = optimoptions('lsqnonlin','Display','off','MaxIterations',2000, ...
        'MaxFunctionEvaluations',2e4,'FunctionTolerance',1e-14,'StepTolerance',1e-12);
    best = inf; p = x0s(1,:);
    for j = 1:size(x0s,1)
        try
            [x, rn] = lsqnonlin(@(x) ResSSVI(x, KQ, IVQ, THETA, T_uniq), x0s(j,:), lb, ub, opts);
        catch
            continue
        end
        if rn < best, best = rn; p = x; end
    end
    rmse_bp = nan(numel(T_uniq),1);
    for i = 1:numel(T_uniq)
        iv = sqrt(max(SsviW(p, THETA(i), KQ{i}),0) / T_uniq(i));
        rmse_bp(i) = 1e4 * sqrt(mean((iv - IVQ{i}(:)).^2));
    end
end

function r = ResSSVI(x, KQ, IVQ, THETA, T_uniq)
    r = [];
    for i = 1:numel(T_uniq)
        iv = sqrt(max(SsviW(x, THETA(i), KQ{i}(:)), 1e-12) / T_uniq(i));
        r = [r; iv - IVQ{i}(:)];
    end
    r = [r; 100*max(0, x(2)*(1 + abs(x(1))) - 2)];
end

function [p, rmse_bp] = FitSVI(k, IVs, T)
% Fit SVI raw sulla varianza totale, con piu' partenze.
    k = k(:); w = IVs(:).^2 * T;
    lb = [-1, 0, -0.999, -1, 1e-4];
    ub = [ 1, 2,  0.999,  1, 2   ];
    x0s = [0.5*min(w), 0.10, -0.50,  0.00, 0.10;
           min(w),     0.05, -0.80,  0.05, 0.20;
           0,          0.20, -0.30, -0.10, 0.05;
           0.8*min(w), 0.30, -0.90,  0.10, 0.30;
           0.5*min(w), 0.50, -0.70,  0.00, 0.50];
    opts = optimoptions('lsqnonlin','Display','off','MaxIterations',3000, ...
        'MaxFunctionEvaluations',2e4,'FunctionTolerance',1e-14,'StepTolerance',1e-12);
    res = @(x) [SviW(x, k) - w; 1e3*max(0, -(x(1) + x(2)*x(5)*sqrt(1 - x(3)^2)))];
    best = inf; p = x0s(1,:);
    for j = 1:size(x0s,1)
        x0 = min(max(x0s(j,:), lb), ub);
        try
            [x, rn] = lsqnonlin(res, x0, lb, ub, opts);
        catch
            continue
        end
        if rn < best, best = rn; p = x; end
    end
    rmse_bp = 1e4 * sqrt(mean((sqrt(max(SviW(p,k),0)/T) - IVs(:)).^2));
end

function Kvar = VFairSVI(o, T_s, K_s, F_s, IV_M, Bid_s, T_uniq, SSVI, THETA)
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
            case 'ssvi'
                % dentro il range quotato: interpolazione delle IV di
                % mercato, IDENTICA alla base.
                kk = log(Kg/Fi);
                wv = SsviW(SSVI, THETA(i), kk);
                LR = (Kg < Ks(1)) | (Kg > Ks(end));
                IVg(LR) = sqrt(max(wv(LR), 1e-8) / Ti);
        end
        OTM = zeros(size(Kg));
        for j = 1:numel(Kg)
            if Kg(j) < Fi, [~,p] = blkprice(Fi,Kg(j),0,Ti,IVg(j)); OTM(j)=p;
            else,          [c,~] = blkprice(Fi,Kg(j),0,Ti,IVg(j)); OTM(j)=c; end
        end
        Kvar(i) = (2/Ti) * trapz(Kg, OTM ./ Kg.^2);
    end
end

function xi = eval_xi0_pchip(tvec, W_pp)
    h=1e-5; tvec=tvec(:);
    tp=max(tvec+h,0); tm=max(tvec-h,0);
    dn=tp-tm; dn(dn<1e-12)=1e-12;
    xi=max((ppval(W_pp,tp)-ppval(W_pp,tm))./dn,1e-8);
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
