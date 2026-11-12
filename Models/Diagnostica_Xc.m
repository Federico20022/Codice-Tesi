clearvars; clc;
tic
% DIAGNOSTICA_XC  Applica SelezionaXc ai cinque modelli calibrati, scadenza
% per scadenza, e confronta il dominio richiesto con quello usato.
% Output: Results/diagnostica_Xc.mat. Tempo: 5-15 minuti.

TOL_BP = 0.5;  % tolleranza di convergenza in bp di volatilita' implicita
N_RIC  = 400;  % nodi Riccati, come negli stadi exact

fprintf('\nCaricamento dataset Refinitiv SPX 2026-05-15 (merged)...\n');

CSV_PATH   = '20260515SPX_merged.csv';
QUOTE_DATE = datenum('2026-05-15');
raw = readtable(CSV_PATH, 'Delimiter', ',', 'VariableNamingRule', 'preserve');

is_spot = strcmp(raw.Root, 'SPOT');
S0_feed = raw.CF_LAST(is_spot);
assert(isscalar(S0_feed), 'Attesa una sola riga SPOT nel CSV, trovate %d.', numel(S0_feed));
raw     = raw(~is_spot, :);

is_call = strcmp(raw.PutCall, 'C');
keep    = ~isnan(raw.STRIKE_PRC) & raw.CF_BID > 0 & ...
          raw.CF_ASK > raw.CF_BID & raw.CF_BID >= 0.375;
raw = raw(keep, :);  is_call = is_call(keep);

K_all   = raw.STRIKE_PRC;
mid_all = 0.5*(raw.CF_BID + raw.CF_ASK);
T_all   = (datenum(raw.EXPIR_DATE) - QUOTE_DATE) / 365;
[~, ia] = unique([T_all, K_all, double(is_call)], 'rows', 'stable');
T_all = T_all(ia); K_all = K_all(ia); mid_all = mid_all(ia); is_call = is_call(ia);

T_uniq_raw = unique(T_all);
F_per_T = nan(length(T_uniq_raw),1);  D_per_T = nan(length(T_uniq_raw),1);
for i = 1:length(T_uniq_raw)
    Ti = T_uniq_raw(i);
    mC = (T_all == Ti) & is_call;  mP = (T_all == Ti) & ~is_call;
    [Kc, iC, iP] = intersect(K_all(mC), K_all(mP));
    if length(Kc) < 3, continue; end
    cI = find(mC); pI = find(mP);
    dCP = mid_all(cI(iC)) - mid_all(pI(iP));
    near = abs(Kc/S0_feed - 1) <= 0.25;
    if sum(near) >= 5, Kr = Kc(near); dr = dCP(near); else, Kr = Kc; dr = dCP; end
    co = [ones(size(Kr)), Kr] \ dr;
    D_per_T(i) = -co(2);  F_per_T(i) = co(1)/D_per_T(i);
end
vT = ~isnan(F_per_T) & D_per_T > 0 & D_per_T < 1.05 & F_per_T > 0;
T_uniq_raw = T_uniq_raw(vT); F_per_T = F_per_T(vT); D_per_T = D_per_T(vT);

F_surf = []; T_surf = []; K_surf = []; Mid_surf = [];
for i = 1:length(T_uniq_raw)
    Ti = T_uniq_raw(i); Fi = F_per_T(i); Di = D_per_T(i);
    mTi = (T_all == Ti);
    Ki = K_all(mTi); mi = mid_all(mTi); ci = is_call(mTi);
    Pf = nan(size(Ki));
    for j = 1:length(Ki)
        if Ki(j) < Fi && ~ci(j),      Pf(j) = mi(j)/Di;
        elseif Ki(j) >= Fi && ci(j),  Pf(j) = mi(j)/Di - Fi + Ki(j);
        end
    end
    ok = ~isnan(Pf) & Pf > 1e-4;
    F_surf=[F_surf; Fi*ones(sum(ok),1)]; T_surf=[T_surf; Ti*ones(sum(ok),1)];
    K_surf=[K_surf; Ki(ok)];             Mid_surf=[Mid_surf; Pf(ok)];
end

IV_Mkt = InvertiIV_OTM(F_surf, K_surf, T_surf, Mid_surf);
lm = log(K_surf ./ F_surf);
id = ~isnan(IV_Mkt) & IV_Mkt > 0.01 & IV_Mkt < 2 & Mid_surf > 1e-4 & ...
     abs(lm) < 0.5 & T_surf >= 0.03 & T_surf <= 3.0;
F_surf=F_surf(id); T_surf=T_surf(id); K_surf=K_surf(id);
Mid_surf=Mid_surf(id); IV_Mkt=IV_Mkt(id);

Ta = unique(T_surf); ma = false(length(T_surf),1);
for i = 1:length(Ta)
    m = T_surf == Ta(i);  if sum(m) >= 5, ma = ma | m; end
end
F_surf=F_surf(ma); T_surf=T_surf(ma); K_surf=K_surf(ma);
Mid_surf=Mid_surf(ma); IV_Mkt=IV_Mkt(ma);

Ta = unique(T_surf); T_sel = [];
for tt = [0.0384 0.0833 0.25 0.5 0.75 1.0 1.5]
    [~, ii] = min(abs(Ta - tt));  T_sel = [T_sel; Ta(ii)];
end
T_sel = unique(T_sel);  km = ismember(T_surf, T_sel);
F_surf=F_surf(km); T_surf=T_surf(km); K_surf=K_surf(km);
Mid_surf=Mid_surf(km); IV_Mkt=IV_Mkt(km);
T_uniq = unique(T_surf);

% Demeterfi + xi0
Kvar_mkt = zeros(length(T_uniq),1);  IV_ATM = zeros(length(T_uniq),1);
for i = 1:length(T_uniq)
    Ti = T_uniq(i); m = (T_surf == Ti);
    Ki = K_surf(m); IVi = IV_Mkt(m); Fi = F_surf(find(m,1));
    [Ks, s] = sort(Ki); IVs = IVi(s);
    [~, ia2] = min(abs(Ks - Fi)); IV_ATM(i) = IVs(ia2);
    Kg = logspace(log10(0.05*Fi), log10(2.50*Fi), 4000)';
    IVg = interp1(Ks, IVs, Kg, 'linear', NaN);
    IVg(Kg < Ks(1)) = IVs(1);  IVg(Kg > Ks(end)) = IVs(end);
    O = zeros(size(Kg));
    for j = 1:length(Kg)
        if Kg(j) < Fi, [~,p] = blkprice(Fi,Kg(j),0,Ti,IVg(j)); O(j)=p;
        else,          [c,~] = blkprice(Fi,Kg(j),0,Ti,IVg(j)); O(j)=c; end
    end
    Kvar_mkt(i) = (2/Ti) * trapz(Kg, O ./ Kg.^2);
end
W_pp = pchip([0; T_uniq(:)], [0; T_uniq(:).*Kvar_mkt(:)]);
xi0_eval = @(tv) eval_xi0_pchip(tv, W_pp);

% 15 strike/scadenza
k15 = false(size(T_surf));
for i = 1:length(T_uniq)
    m = find(T_surf == T_uniq(i));
    z = log(K_surf(m)./F_surf(m)) ./ (IV_ATM(i)*sqrt(T_uniq(i)));
    if numel(m) <= 15, k15(m) = true; continue; end
    for tt = linspace(max(-3.0,min(z)), min(2.0,max(z)), 15)
        [~, jj] = min(abs(z - tt));  k15(m(jj)) = true;
    end
end
F_surf=F_surf(k15); T_surf=T_surf(k15); K_surf=K_surf(k15);
fprintf('Campione: %d opzioni su %d scadenze\n', length(F_surf), length(T_uniq));

% calibrazioni
resdir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
if isempty(fileparts(mfilename('fullpath'))), resdir = fullfile(pwd,'..','Results'); end
CAL.hes  = load(fullfile(resdir,'calib_heston_7x15z.mat'),'theta_ottimi');
CAL.rh   = load(fullfile(resdir,'calib_rheston_7x15z_curva.mat'),'theta_ottimi');
CAL.cgmy = load(fullfile(resdir,'calib_cgmy_7x15z.mat'),'theta_ottimi');
CAL.rhc  = load(fullfile(resdir,'calib_rhcgmy_7x15z.mat'),'theta_ottimi');
CAL.rhh  = load(fullfile(resdir,'calib_rhh_7x15z.mat'),'theta_ottimi');

% Configurazione attualmente usata nelle calibrazioni, per il confronto.
USATO = struct('Heston',[12 12],'rHeston',[12 12],'CGMY',[12 12], ...
               'rHCGMY',[12 12],'rHawkesH',[12 12]);

MODELLI = {'Heston','rHeston','CGMY','rHCGMY','rHawkesH'};
RES = struct();

for im = 1:numel(MODELLI)
    nome = MODELLI{im};
    fprintf('\n================ %s ================\n', nome);
    opts = struct('TOL_BP', TOL_BP, 'CP', -1, 'VERBOSE', true);
    % Il rough Hawkes Heston usa uno schema di Eulero il cui numero di passi
    % cresce con la frequenza massima, quindi con la larghezza del dominio:
    % la scala e' accorciata per non farlo esplodere in tempo.
    % Lo schema di Eulero del rough Hawkes Heston adatta il numero di passi
    % alla frequenza massima, che cresce con Xc: a Xc=24 e pow=13 arriva a
    % ~880 passi su 2048 frequenze, e il calcolo va fatto due volte per
    % gradino.
    if strcmp(nome,'rHawkesH'), opts.XC_LADDER = [9 12 18]; end

    Xc_per_T = nan(length(T_uniq),1);  pow_per_T = nan(length(T_uniq),1);
    dIV_usato = nan(length(T_uniq),1); coda_usata = nan(length(T_uniq),1);

    for it = 1:length(T_uniq)
        T  = T_uniq(it);
        m  = (T_surf == T);
        K  = K_surf(m);
        F  = F_surf(find(m,1));
        fprintf('\n-- T = %.4f (%d strike) --\n', T, numel(K));

        [pf, cf] = CostruisciHandle(nome, CAL, T, K, F, xi0_eval, N_RIC);
        cf=[];
        [Xc_per_T(it), pow_per_T(it), INFO] = SelezionaXc(pf, cf, T, K, F, opts);

        % scarto della configurazione oggi in uso, allo stesso riferimento
        ju = find(INFO.Xc == USATO.(nome)(1), 1);
        if ~isempty(ju)
            dIV_usato(it)  = INFO.dIV_bp(ju);
            coda_usata(it) = INFO.coda_sx(ju) + INFO.coda_dx(ju);
        end
    end

    RES.(nome) = struct('T',T_uniq,'Xc',Xc_per_T,'pow',pow_per_T, ...
                        'dIV_usato_bp',dIV_usato,'coda_usata',coda_usata);
    fprintf('\n>> %s: Xc richiesto = %g (max sulle scadenze), pow = %d\n', ...
        nome, max(Xc_per_T), max(pow_per_T));
end

% TABELLA RIASSUNTIVA
fprintf('\n\n========================================================\n');
fprintf('SELEZIONE DEL DOMINIO - tolleranza %.2f bp di IV\n', TOL_BP);
fprintf('========================================================\n');
fprintf('%-10s  %-8s  %-6s  %-9s  %-14s  %-12s\n', ...
    'Modello','Xc rich.','pow','dx','max dIV a Xc=12','coda a Xc=12');
for im = 1:numel(MODELLI)
    nome = MODELLI{im};  R = RES.(nome);
    Xcm = max(R.Xc);  powm = max(R.pow);
    d = max(R.dIV_usato_bp);
    if isinf(d) || isnan(d), ds = '      n.d.'; else, ds = sprintf('%10.2f', d); end
    fprintf('%-10s  %-8g  %-6d  %-9.5f  %s      %-12.3e\n', ...
        nome, Xcm, powm, 2*Xcm/2^powm, ds, max(R.coda_usata));
end
fprintf(['\nLettura. "Xc rich." e'' il minimo dominio che porta la IV entro %.2f bp\n' ...
         'dal riferimento su tutti gli strike e tutte le scadenze. La colonna\n' ...
         '"max dIV a Xc=12" e'' l''errore della configurazione oggi in uso, allo\n' ...
         'stesso riferimento: e'' il numero da riportare come budget di errore\n' ...
         'del pricer per quel modello. La massa di coda e'' un indicatore, non\n' ...
         'il criterio di selezione.\n'], TOL_BP);

if ~exist(resdir,'dir'), mkdir(resdir); end
save(fullfile(resdir,'diagnostica_Xc.mat'), 'RES', 'TOL_BP', 'USATO', 'T_uniq');
fprintf('\nSalvato in %s\n', fullfile(resdir,'diagnostica_Xc.mat'));
toc

% FUNZIONI LOCALI
function [prezzaFn, cfFn] = CostruisciHandle(nome, CAL, T, K, F, xi0_fn, N_RIC)
    n  = numel(K);  Z = zeros(n,1);  Tv = T*ones(n,1);
    ti = ((0:N_RIC-1) + 0.5) / N_RIC * T;
    switch nome
        case 'Heston'
            th = CAL.hes.theta_ottimi;
            prezzaFn = @(Xc,pw) PrezzaSuperficieSINC_yield_new(F*ones(n,1), Tv, K(:), ...
                Z, Z, 'Heston', th, [], [], Xc, 0, pw, -1);
            cfFn = @(u) charfun('Heston', th, [], [], u, T);
        case 'CGMY'
            th = CAL.cgmy.theta_ottimi;
            prezzaFn = @(Xc,pw) PrezzaSuperficieSINC_yield_new(F*ones(n,1), Tv, K(:), ...
                Z, Z, 'CGMY', th, [], [], Xc, 0, pw, -1);
            cfFn = @(u) charfun('CGMY', th, [], [], u, T);
        case 'rHawkesH'
            th = CAL.rhh.theta_ottimi;
            prezzaFn = @(Xc,pw) PrezzaSuperficieSINC_yield_new(F*ones(n,1), Tv, K(:), ...
                Z, Z, 'rough Hawkes Heston', th, [], N_RIC, Xc, 0, pw, -1);
            cfFn = @(u) charfun('rough Hawkes Heston', th, [], N_RIC, u, T);
        case 'rHeston'
            th  = [CAL.rh.theta_ottimi(:)', 0];
            xig = xi0_fn(ti);  xig = xig(:)';
            prezzaFn = @(Xc,pw) PrezzaSuperficieSINC_yield(F, Tv, K(:), Z, Z, ...
                'rough Heston', th, xig, N_RIC, Xc, 0, pw, -1);
            cfFn = @(u) charfun('rough Heston', th, xig, N_RIC, u, T);
        case 'rHCGMY'
            t7  = CAL.rhc.theta_ottimi;  % [H,nu,rho,C,G,M,Y]
            Lcg = ContributiCGMY(t7(4), t7(5), t7(6), t7(7));
            xig = max(xi0_fn(ti) - Lcg, 1e-6);  xig = xig(:)';
            pc  = [t7(1:3), 0, t7(4:7)];
            prezzaFn = @(Xc,pw) PrezzaSuperficieSINC_yield(F, Tv, K(:), Z, Z, ...
                'rough Heston CGMY', pc, xig, N_RIC, Xc, 0, pw, -1);
            cfFn = @(u) charfun('rough Heston CGMY', pc, xig, N_RIC, u, T);
        otherwise
            error('Modello non riconosciuto: %s', nome);
    end
end

function [Lcg, Vcg] = ContributiCGMY(C, G, M, Y)
    g_m1 = gamma(-Y) * ((M-1)^Y - M^Y + (G+1)^Y - G^Y);
    mu1  = gamma(1-Y) * (M^(Y-1) - G^(Y-1));
    Lcg  = 2*C*(g_m1 - mu1);
    Vcg  = C * gamma(2-Y) * (M^(Y-2) + G^(Y-2));
end

function xi = eval_xi0_pchip(tvec, W_pp)
    h = 1e-5;  tvec = tvec(:);
    tp = max(tvec+h,0); tm = max(tvec-h,0);
    dn = tp - tm; dn(dn < 1e-12) = 1e-12;
    xi = max((ppval(W_pp,tp) - ppval(W_pp,tm)) ./ dn, 1e-8);
end

function IV = InvertiIV_OTM(F_vec, K_vec, T_vec, P_fwd_vec)
    IV = nan(size(F_vec));
    for i = 1:length(F_vec)
        try
            if K_vec(i) <= F_vec(i)
                IV(i) = blkimpv(F_vec(i),K_vec(i),0,T_vec(i),P_fwd_vec(i),'Class','put');
            else
                Cf = max(P_fwd_vec(i) + F_vec(i) - K_vec(i), 1e-6);
                IV(i) = blkimpv(F_vec(i),K_vec(i),0,T_vec(i),Cf,'Class','call');
            end
        catch
            IV(i) = NaN;
        end
    end
end
