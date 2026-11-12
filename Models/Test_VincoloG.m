clearvars; clc;
tic
% TEST_VINCOLOG  Rough Heston + CGMY: effetto del vincolo su q e di G libero.
%   A) G = 0.05 fisso, senza vincolo
%   B) G = 0.05 fisso, con vincolo q = Lcg / min xi0 <= Q_MAX
%   C) G libero, con vincolo
% Output: Results/test_vincolo_G.mat.

Q_MAX = 0.98;  % quota massima della varianza forward minima presa dai salti
G_FIX = 0.05;

fprintf('\nCaricamento dataset Refinitiv SPX 2026-05-15 (merged)...\n');
CSV_PATH   = '20260515SPX_merged.csv';
QUOTE_DATE = datenum('2026-05-15');
raw = readtable(CSV_PATH, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
rng(42);

is_spot = strcmp(raw.Root,'SPOT');
S0_feed = raw.CF_LAST(is_spot);
assert(isscalar(S0_feed), 'Attesa una sola riga SPOT, trovate %d.', numel(S0_feed));
raw = raw(~is_spot,:);
is_call = strcmp(raw.PutCall,'C');
kp = ~isnan(raw.STRIKE_PRC) & raw.CF_BID>0 & raw.CF_ASK>raw.CF_BID & raw.CF_BID>=0.375;
raw = raw(kp,:); is_call = is_call(kp);

K_all = raw.STRIKE_PRC;
mid_all = 0.5*(raw.CF_BID + raw.CF_ASK);
T_all = (datenum(raw.EXPIR_DATE) - QUOTE_DATE)/365;
[~,ia] = unique([T_all,K_all,double(is_call)],'rows','stable');
T_all=T_all(ia); K_all=K_all(ia); mid_all=mid_all(ia); is_call=is_call(ia);

Tr = unique(T_all); Fp = nan(numel(Tr),1); Dp = nan(numel(Tr),1);
for i = 1:numel(Tr)
    Ti=Tr(i); mC=(T_all==Ti)&is_call; mP=(T_all==Ti)&~is_call;
    [Kc,iC,iP] = intersect(K_all(mC),K_all(mP));
    if numel(Kc)<3, continue; end
    cI=find(mC); pI=find(mP); d=mid_all(cI(iC))-mid_all(pI(iP));
    nr = abs(Kc/S0_feed-1)<=0.25;
    if sum(nr)>=5, Kr=Kc(nr); dr=d(nr); else, Kr=Kc; dr=d; end
    co=[ones(size(Kr)),Kr]\dr; Dp(i)=-co(2); Fp(i)=co(1)/Dp(i);
end
v = ~isnan(Fp) & Dp>0 & Dp<1.05 & Fp>0;
Tr=Tr(v); Fp=Fp(v); Dp=Dp(v);

F_surf=[];T_surf=[];K_surf=[];Mid_surf=[];
for i=1:numel(Tr)
    Ti=Tr(i); Fi=Fp(i); Di=Dp(i); m=(T_all==Ti);
    Ki=K_all(m); mi=mid_all(m); ci=is_call(m); P=nan(size(Ki));
    for j=1:numel(Ki)
        if Ki(j)<Fi && ~ci(j),     P(j)=mi(j)/Di;
        elseif Ki(j)>=Fi && ci(j), P(j)=mi(j)/Di-Fi+Ki(j); end
    end
    o = ~isnan(P) & P>1e-4;
    F_surf=[F_surf;Fi*ones(sum(o),1)]; T_surf=[T_surf;Ti*ones(sum(o),1)];
    K_surf=[K_surf;Ki(o)];             Mid_surf=[Mid_surf;P(o)];
end

IV_Mkt = InvertiIV_OTM(F_surf,K_surf,T_surf,Mid_surf);
lm = log(K_surf./F_surf);
id = ~isnan(IV_Mkt)&IV_Mkt>0.01&IV_Mkt<2&Mid_surf>1e-4&abs(lm)<0.5& ...
     T_surf>=0.03&T_surf<=3.0;
F_surf=F_surf(id);T_surf=T_surf(id);K_surf=K_surf(id);
Mid_surf=Mid_surf(id);IV_Mkt=IV_Mkt(id);
Ta=unique(T_surf); ma=false(numel(T_surf),1);
for i=1:numel(Ta), m=T_surf==Ta(i); if sum(m)>=5, ma=ma|m; end, end
F_surf=F_surf(ma);T_surf=T_surf(ma);K_surf=K_surf(ma);
Mid_surf=Mid_surf(ma);IV_Mkt=IV_Mkt(ma);
Ta=unique(T_surf); Ts=[];
for tt=[0.0384 0.0833 0.25 0.5 0.75 1.0 1.5]
    [~,ii]=min(abs(Ta-tt)); Ts=[Ts;Ta(ii)];
end
Ts=unique(Ts); km=ismember(T_surf,Ts);
F_surf=F_surf(km);T_surf=T_surf(km);K_surf=K_surf(km);
Mid_surf=Mid_surf(km);IV_Mkt=IV_Mkt(km);
T_uniq=unique(T_surf);

% Demeterfi + xi0
Kvar_mkt=zeros(numel(T_uniq),1); IV_ATM=zeros(numel(T_uniq),1);
for i=1:numel(T_uniq)
    Ti=T_uniq(i); m=(T_surf==Ti);
    Ki=K_surf(m); IVi=IV_Mkt(m); Fi=F_surf(find(m,1));
    [Ks,s]=sort(Ki); IVs=IVi(s);
    [~,a2]=min(abs(Ks-Fi)); IV_ATM(i)=IVs(a2);
    Kg=logspace(log10(0.05*Fi),log10(2.50*Fi),4000)';
    IVg=interp1(Ks,IVs,Kg,'linear',NaN);
    IVg(Kg<Ks(1))=IVs(1); IVg(Kg>Ks(end))=IVs(end);
    O=zeros(size(Kg));
    for j=1:numel(Kg)
        if Kg(j)<Fi, [~,p]=blkprice(Fi,Kg(j),0,Ti,IVg(j)); O(j)=p;
        else,        [c,~]=blkprice(Fi,Kg(j),0,Ti,IVg(j)); O(j)=c; end
    end
    Kvar_mkt(i)=(2/Ti)*trapz(Kg,O./Kg.^2);
end
W_pp = pchip([0;T_uniq(:)],[0;T_uniq(:).*Kvar_mkt(:)]);
xi0_eval = @(tv) eval_xi0_pchip(tv,W_pp);

k15=false(size(T_surf));
for i=1:numel(T_uniq)
    m=find(T_surf==T_uniq(i));
    z=log(K_surf(m)./F_surf(m))./(IV_ATM(i)*sqrt(T_uniq(i)));
    if numel(m)<=15, k15(m)=true; continue; end
    for tt=linspace(max(-3.0,min(z)),min(2.0,max(z)),15)
        [~,jj]=min(abs(z-tt)); k15(m(jj))=true;
    end
end
F_surf=F_surf(k15);T_surf=T_surf(k15);K_surf=K_surf(k15);
Mid_surf=Mid_surf(k15);IV_Mkt=IV_Mkt(k15);
CP_flag=-1*ones(size(F_surf));

Vega=zeros(numel(F_surf),1);
for i=1:numel(F_surf)
    d1=(log(F_surf(i)/K_surf(i))+0.5*IV_Mkt(i)^2*T_surf(i))/(IV_Mkt(i)*sqrt(T_surf(i)));
    Vega(i)=F_surf(i)*sqrt(T_surf(i))*normpdf(d1);
end
sqrt_w = 1./(max(Vega,1e-4).*max(IV_Mkt,0.01));
sqrt_w = min(sqrt_w, prctile(sqrt_w,99));
fprintf('Campione: %d opzioni su %d scadenze\n', numel(F_surf), numel(T_uniq));

% xi_min: il minimo della curva forward variance sull'unione delle griglie
% Riccati effettivamente usate.
N_RIC_F = 100;  N_RIC_E = 400;
tt_all = [];
for k = 1:numel(T_uniq)
    tt_all = [tt_all, ((0:N_RIC_E-1)+0.5)/N_RIC_E*T_uniq(k)];
    tt_all = [tt_all, ((0:N_RIC_F-1)+0.5)/N_RIC_F*T_uniq(k)];
end
xi_min = min(xi0_eval(tt_all));
fprintf('min xi0 sulle griglie Riccati = %.6f  (sqrt = %.2f%%)\n', ...
    xi_min, 100*sqrt(xi_min));
fprintf('Vincolo: Lcg <= %.3f * xi_min = %.6f\n\n', Q_MAX, Q_MAX*xi_min);

% BOUND
%        [H,     nu,   rho,   q,     G,      M,     Y   ]
lb7 = [0.005, 0.05, -0.99, 1e-4, 0.005, 5,  0.9 ];
ub7 = [0.49,  5.00,  0.10, Q_MAX, 50.0, 200.0, 1.95];
nomi7 = {'H','nu','rho','q','G','M','Y'};
lb6 = lb7([1 2 3 4 6 7]);  ub6 = ub7([1 2 3 4 6 7]);
nomi6 = nomi7([1 2 3 4 6 7]);

% Il lower bound su M resta 1.01, cioe' il vincolo MATEMATICO di esistenza
% di E[exp(X)].

% PUNTI DI PARTENZA
% L'ottimo interno e quello d'angolo trovati dal multistart, espressi in C,
% vengono convertiti in q.
int_C = [0.087785, 0.36758, -0.7182, 0.00056804, 21.7069, 1.8315];  % [H nu rho C M Y]
ang_C = [0.0387062, 0.30203, -0.99,  0.000365274, 1.01,   1.92912];

P0 = {};
P0{end+1} = struct('eti','interno, G=0.05',  'p', Conv(int_C, G_FIX, xi_min), 'G', G_FIX);
P0{end+1} = struct('eti','ANGOLO,  G=0.05',  'p', Conv(ang_C, G_FIX, xi_min), 'G', G_FIX);
P0{end+1} = struct('eti','ANGOLO,  G=0.50',  'p', Conv(ang_C, 0.50,  xi_min), 'G', 0.50);
P0{end+1} = struct('eti','ANGOLO,  G=5.00',  'p', Conv(ang_C, 5.00,  xi_min), 'G', 5.00);
P0{end+1} = struct('eti','interno, G=2.00',  'p', Conv(int_C, 2.00,  xi_min), 'G', 2.00);
P0{end+1} = struct('eti','neutro,  G=1.00',  'p', [0.10 0.40 -0.70 0.15 20 1.60], 'G', 1.00);

opt = optimoptions('lsqnonlin','Algorithm','trust-region-reflective', ...
    'MaxFunctionEvaluations',3000,'MaxIterations',200,'Display','off', ...
    'StepTolerance',1e-8,'FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);

% RESIDUI
% 6 parametri, G fisso, SENZA vincolo: q non esiste, si usa C direttamente.
r6_libero = @(p, Nr) sqrt_w .* (Prezza(F_surf,T_surf,K_surf, ...
        [p(1:3), p(4), G_FIX, p(5), p(6)], T_uniq, xi0_eval, Nr, CP_flag) - Mid_surf);

% 6 parametri, G fisso, CON vincolo: p(4) = q.
r6_vinc = @(p, Nr) sqrt_w .* (Prezza(F_surf,T_surf,K_surf, ...
        DaQ([p(1:3), p(4), G_FIX, p(5), p(6)], xi_min), T_uniq, xi0_eval, Nr, CP_flag) - Mid_surf);

% 7 parametri, G libero, CON vincolo: p(4) = q, p(5) = G.
r7_vinc = @(p, Nr) sqrt_w .* (Prezza(F_surf,T_surf,K_surf, ...
        DaQ(p, xi_min), T_uniq, xi0_eval, Nr, CP_flag) - Mid_surf);

% ESECUZIONE
SPEC = { 'A) 6p, G=0.05, senza vincolo', 6, r6_libero, false
         'B) 6p, G=0.05, con vincolo',   6, r6_vinc,   true
         'C) 7p, G libero, con vincolo', 7, r7_vinc,   true };

RES = struct();
for is = 1:size(SPEC,1)
    eti = SPEC{is,1}; np = SPEC{is,2}; rfun = SPEC{is,3}; usaQ = SPEC{is,4};
    fprintf('\n================ %s ================\n', eti);
    if np == 6, lbx = lb6; ubx = ub6; nomi = nomi6;
    else,       lbx = lb7; ubx = ub7; nomi = nomi7; end
    if ~usaQ
        % specificazione A: il quarto parametro e' C, non q
        lbx(4) = 1e-4; ubx(4) = 10.0; nomi{4} = 'C';
    end

    fprintf('%-20s %-13s %-13s  %s\n','partenza','loss fast','loss exact','parametri');
    best = inf; bth = [];
    TH = nan(numel(P0), np); LS = inf(numel(P0),1);
    for k = 1:numel(P0)
        p0 = P0{k}.p;  % [H nu rho q M Y] con G a parte
        if np == 7, x0 = [p0(1:4), P0{k}.G, p0(5:6)];
        else,       x0 = p0; end
        if ~usaQ  % converti q -> C per la spec A
            Gk = G_FIX;
            x0(4) = x0(4)*xi_min / Lambda(Gk, x0(end-1), x0(end));
        end
        x0 = min(max(x0, lbx), ubx);
        try
            th = lsqnonlin(@(p) rfun(p,N_RIC_F), x0, lbx, ubx, opt);
            lf = sum(rfun(th,N_RIC_F).^2);
            th = lsqnonlin(@(p) rfun(p,N_RIC_E), th, lbx, ubx, opt);
            le = sum(rfun(th,N_RIC_E).^2);
        catch ME
            fprintf('%-20s FALLITA (%s)\n', P0{k}.eti, ME.message); continue
        end
        TH(k,:) = th; LS(k) = le;
        fprintf('%-20s %-13.6e %-13.6e  %s\n', P0{k}.eti, lf, le, ...
            strjoin(compose('%.4g', th), ' '));
        if le < best, best = le; bth = th; end
    end

    % parametri al bound
  tolb = 1e-3;
    log_sc = false(1,np);
    for j = 1:np
        log_sc(j) = ismember(nomi{j}, {'nu','C','q','G','M'});
    end
    alb = {};
    for j = 1:np
        if log_sc(j)
            atlb = bth(j) < lbx(j)*1.05;
            atub = bth(j) > ubx(j)*0.95;
        else
            rj = ubx(j) - lbx(j);
            atlb = bth(j) < lbx(j) + tolb*rj;
            atub = bth(j) > ubx(j) - tolb*rj;
        end
        if atlb, alb{end+1} = [nomi{j} '@lb']; end
        if atub, alb{end+1} = [nomi{j} '@ub']; end
    end
    fprintf('\nMIGLIORE: loss = %.6e\n', best);
    for j = 1:np, fprintf('  %-4s = %.6g\n', nomi{j}, bth(j)); end
    if isempty(alb), fprintf('Bound: tutti interni\n');
    else,            fprintf(2, 'Bound attivi: %s\n', strjoin(alb,', ')); end

    % diagnostiche della componente di salto
    if usaQ
        full = DaQ(SettePar(bth, np, G_FIX), xi_min);
    else
        full = [bth(1:3), bth(4), G_FIX, bth(5), bth(6)];
    end
    [Lcg, Vcg] = ContributiCGMY(full(4), full(5), full(6), full(7));
    fprintf('C = %.6g | G = %.4g | Lcg = %.5f (%.1f%% di xi_min) | VS/L comp = %.3f\n', ...
        full(4), full(5), Lcg, 100*Lcg/xi_min, Vcg/max(Lcg,1e-12));

    RES.(sprintf('spec%c', 'A'+is-1)) = struct('eti',eti,'theta',TH,'loss',LS, ...
        'best',bth,'best_loss',best,'full',full,'bound',{alb},'nomi',{nomi});
end

% CONFRONTO
fprintf('\n\n==================================================\n');
fprintf('CONFRONTO DELLE TRE SPECIFICAZIONI\n');
fprintf('==================================================\n');
fprintf('%-32s %-14s %-10s %s\n','specificazione','loss exact','RMSE rel','bound attivi');
LA = RES.specA.best_loss;
for s = {'specA','specB','specC'}
    R = RES.(s{1});
    if isempty(R.bound), bs = 'nessuno'; else, bs = strjoin(R.bound,', '); end
    fprintf('%-32s %-14.6e %-10.4f %s\n', R.eti, R.best_loss, ...
        sqrt(R.best_loss/LA), bs);
end
fprintf(['\nLettura. La colonna RMSE rel e'' il rapporto di RMSE rispetto alla\n' ...
         'specificazione A. Se in C il bound su M sparisce, liberare G scioglie\n' ...
         'l''angolo e la strada e'' quella: va fatto un multistart completo a 7\n' ...
         'parametri. Se M resta al bound (o ci si aggiunge G), la strada e''\n' ...
         'chiusa e serve il bound economico su M.\n']);

resdir = fullfile(fileparts(mfilename('fullpath')),'..','Results');
if isempty(fileparts(mfilename('fullpath'))), resdir = fullfile(pwd,'..','Results'); end
if ~exist(resdir,'dir'), mkdir(resdir); end
save(fullfile(resdir,'test_vincolo_G.mat'),'RES','xi_min','Q_MAX','P0');
fprintf('\nSalvato in %s\n', fullfile(resdir,'test_vincolo_G.mat'));
toc

% FUNZIONI LOCALI
function lam = Lambda(G, M, Y)
% Lcg = C * lam(G,M,Y).
    lam = 2*( gamma(-Y)*((M-1)^Y - M^Y + (G+1)^Y - G^Y) ...
            - gamma(1-Y)*(M^(Y-1) - G^(Y-1)) );
end

function p7 = DaQ(pq, xi_min)
% [H nu rho q G M Y] -> [H nu rho C G M Y]
    p7 = pq;
    p7(4) = pq(4) * xi_min / Lambda(pq(5), pq(6), pq(7));
end

function pq = Conv(pC6, G, xi_min)
% [H nu rho C M Y] + G -> [H nu rho q M Y] con q = C*lambda/xi_min
    q = pC6(4) * Lambda(G, pC6(5), pC6(6)) / xi_min;
    pq = [pC6(1:3), min(max(q,1e-4),0.98), pC6(5:6)];
end

function p7 = SettePar(th, np, G_FIX)
    if np == 7, p7 = th;
    else,       p7 = [th(1:4), G_FIX, th(5), th(6)]; end
end

function [Lcg, Vcg] = ContributiCGMY(C, G, M, Y)
    Lcg = C * Lambda(G, M, Y);
    Vcg = C * gamma(2-Y) * (M^(Y-2) + G^(Y-2));
end

function P = Prezza(F_vec, T_vec, K_vec, p7, T_uniq, xi0_fn, N_ric, CP_vec)
% p7 = [H, nu, rho, C, G, M, Y]
    pc = [p7(1:3), 0, p7(4), p7(5), p7(6), p7(7)];
    Lcg = ContributiCGMY(p7(4), p7(5), p7(6), p7(7));
    P = zeros(numel(K_vec),1);
    for k = 1:numel(T_uniq)
        T = T_uniq(k); idx = (T_vec == T);
        if ~any(idx), continue; end
        Ft = F_vec(find(idx,1)); Kt = K_vec(idx); CPt = CP_vec(idx);
        ti = ((0:N_ric-1)+0.5)/N_ric*T;
        xg = max(xi0_fn(ti) - Lcg, 1e-6); xg = xg(:)';
        P(idx) = PrezzaSuperficieSINC_yield(Ft, T*ones(size(Kt)), Kt, ...
            zeros(size(Kt)), zeros(size(Kt)), ...
            'rough Heston CGMY', pc, xg, N_ric, 12, 0, 12, CPt(1));
    end
end

function xi = eval_xi0_pchip(tvec, W_pp)
    h = 1e-5; tvec = tvec(:);
    tp = max(tvec+h,0); tm = max(tvec-h,0);
    dn = tp - tm; dn(dn<1e-12) = 1e-12;
    xi = max((ppval(W_pp,tp)-ppval(W_pp,tm))./dn, 1e-8);
end

function IV = InvertiIV_OTM(F,K,T,P)
    IV = nan(size(F));
    for i = 1:numel(F)
        try
            if K(i) <= F(i)
                IV(i) = blkimpv(F(i),K(i),0,T(i),P(i),'Class','put');
            else
                IV(i) = blkimpv(F(i),K(i),0,T(i),max(P(i)+F(i)-K(i),1e-6),'Class','call');
            end
        catch
            IV(i) = NaN;
        end
    end
end
