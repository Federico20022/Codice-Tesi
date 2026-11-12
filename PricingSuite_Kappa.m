tic
% PRICINGSUITE_KAPPA  Pricing suite delle tre varianti con kappa libero,
% confrontate con Heston e rough Hawkes Heston (Sez. 7.6).
%   1) call europee      T = 0.10, 0.50, 1.00   K/F = 0.80 ... 1.10
%   2) digitali Q(S>K)   T = 0.05, 0.25         K/F = 0.90 ... 1.10
%   3) skew ATM          pendenza log-log su T <= 0.25
%   4) VS / log-contract T = 1/12, 0.25, 1
% Input: Results/calib_*.mat. Output: Results/PricingSuite_kappa.mat, il log
% .txt e Figure/tab_kappa_*.tex. Tempo: pochi minuti.
root_path = fileparts(mfilename('fullpath'));
if isempty(root_path), root_path = pwd; end
addpath(fullfile(root_path, 'Functions'));
addpath(fullfile(root_path, 'Models'));
addpath(fullfile(root_path, 'Results'));
addpath(fullfile(root_path, 'Data'));
RES_PATH = fullfile(root_path, 'Results');
FIG_PATH = fullfile(root_path, 'Figure');
if ~isfolder(FIG_PATH), mkdir(FIG_PATH); end
LOG = fullfile(RES_PATH, 'PricingSuite_kappa_log.txt');
if isfile(LOG), delete(LOG); end
diary(LOG);

% 2. CARICAMENTO CALIBRAZIONI
req = { 'calib_heston_7x15z.mat',            'hes';  ...
        'calib_hestoncurvekappa_7x15z.mat',  'hck';  ...
        'calib_rheston_7x15z_curva.mat',     'rh';   ...
        'calib_rhestonkappa_7x15z.mat',      'rhk';  ...
        'calib_cgmy_7x15z.mat',              'cgmy'; ...
        'calib_rhcgmy_7x15z.mat',            'rhc';  ...
        'calib_rhcgmykappa_7x15z.mat',       'rhck'; ...
        'calib_rhh_7x15z.mat',               'rhh' };
CAL = struct();
for i = 1:size(req,1)
    f = fullfile(RES_PATH, req{i,1});
    assert(isfile(f), 'Manca %s in %s: runna la calibrazione corrispondente.', req{i,1}, RES_PATH);
    CAL.(req{i,2}) = load(f);
end

F0 = 100;
% Ordine: prima i cinque modelli CON mean reversion (sono quelli delle
% tabelle nuove), poi i tre senza.
MODELLI = {'Heston','HestonCurveK','rHestonK','rHK+CGMY','rHawkesH', ...
           'rHeston','CGMY','rH+CGMY'};
NOMI_TEX = {'Heston','Heston-curve-$\kappa$','rough Heston-$\kappa$', ...
            'rough Heston-$\kappa$ + CGMY','rough Hawkes Heston', ...
            'rough Heston','CGMY','rough Heston + CGMY'};
CON_K = 1:5;  % gruppo con mean reversion
NM = numel(MODELLI);

fprintf('=== PARAMETRI CARICATI (da %s) ===\n', RES_PATH);
for i = 1:size(req,1)
    fprintf('%-6s: %s\n', req{i,2}, mat2str(CAL.(req{i,2}).theta_ottimi, 5));
end
fprintf('Nodi Riccati: %d\n\n', N_RIC());

% 1) CALL EUROPEE
T_call = [0.10, 0.50, 1.00];
Kf     = [0.80, 0.90, 0.95, 1.00, 1.05, 1.10];
fprintf('=== 1) CALL EUROPEE (F=100, prezzi forward) ===\n');
CallPrices = zeros(numel(T_call), numel(Kf), NM);
CallIV     = nan(numel(T_call), numel(Kf), NM);
for it = 1:numel(T_call)
    T = T_call(it); K = F0*Kf(:);
    for imod = 1:NM
        CallPrices(it,:,imod) = PriceCallsModel(MODELLI{imod}, CAL, F0, T, K, 1);
        for ik = 1:numel(Kf)
            CallIV(it,ik,imod) = IVcall(F0, K(ik), T, CallPrices(it,ik,imod));
        end
    end
    StampaBlocco(sprintf('T = %.2f', T), Kf, squeeze(CallPrices(it,:,:)), MODELLI);
end
fprintf('\nDispersione (max-min)/min nel gruppo CON mean reversion:\n');
StampaDispersione(T_call, Kf, CallPrices(:,:,CON_K), CallIV(:,:,CON_K));

% 2) DIGITALI
T_dig = [0.05, 0.25];
Kf_d  = [0.90, 0.95, 0.98, 1.00, 1.02, 1.05, 1.10];
fprintf('\n=== 2) DIGITALI Q(S_T > K) ===\n');
Digitals = zeros(numel(T_dig), numel(Kf_d), NM);
for it = 1:numel(T_dig)
    T = T_dig(it); K = F0*Kf_d(:);
    for imod = 1:NM
        Digitals(it,:,imod) = PriceCallsModel(MODELLI{imod}, CAL, F0, T, K, 0);
    end
    StampaBlocco(sprintf('T = %.2f', T), Kf_d, squeeze(Digitals(it,:,:)), MODELLI);
end
fprintf('\nDispersione (max-min)/min nel gruppo CON mean reversion:\n');
StampaDispersione(T_dig, Kf_d, Digitals(:,:,CON_K), []);

% 3) SKEW ATM
fprintf('\n=== 3) SKEW ATM: pendenza log-log su T <= 0.25 ===\n');
T_skew = logspace(log10(0.02), log10(2.5), 14);
dk     = 0.005;
SkewATM = zeros(numel(T_skew), NM);
for it = 1:numel(T_skew)
    T = T_skew(it); Kpm = F0*exp([-dk; dk]);
    for imod = 1:NM
        C2 = PriceCallsModel(MODELLI{imod}, CAL, F0, T, Kpm, 1);
        SkewATM(it,imod) = (IVcall(F0,Kpm(2),T,C2(2)) - IVcall(F0,Kpm(1),T,C2(1)))/(2*dk);
    end
end
mshort = T_skew(:) <= 0.25;
Slope = zeros(NM,1);
fprintf('%-10s  %12s  %12s  %10s  %10s\n','Modello','slope','H implicito','|psi| 0.02','|psi| 0.25');
for imod = 1:NM
    cf = polyfit(log(T_skew(mshort)), log(abs(SkewATM(mshort,imod))), 1);
    Slope(imod) = cf(1);
    fprintf('%-10s  %+12.4f  %12.4f  %10.4f  %10.4f\n', MODELLI{imod}, cf(1), cf(1)+0.5, ...
        abs(SkewATM(1,imod)), abs(interp1(log(T_skew), SkewATM(:,imod), log(0.25))));
end

% 4) VARIANCE SWAP vs LOG-CONTRACT
fprintf('\n=== 4) VARIANCE SWAP vs LOG-CONTRACT ===\n');
T_vs = [1/12, 0.25, 1.0];
h_cf = 1e-3;
ResVS = zeros(numel(T_vs), NM, 3);
fprintf('%-8s %-10s %13s %13s %10s\n', 'T', 'Modello', 'sqrt(L)%', 'sqrt(VS)%', 'VS/L');
for it = 1:numel(T_vs)
    T = T_vs(it);
    for imod = 1:NM
        phi_h = CharfunModel(MODELLI{imod}, CAL, h_cf, T);
        L  = -2/T*imag(phi_h)/h_cf;
        VS = VarianceSwapModel(MODELLI{imod}, CAL, T);
        ResVS(it,imod,:) = [L, VS, VS/L];
        fprintf('%-8.4f %-10s %12.3f%% %12.3f%% %10.4f\n', T, MODELLI{imod}, ...
            100*sqrt(max(L,0)), 100*sqrt(max(VS,0)), VS/L);
    end
    fprintf('\n');
end
th = CAL.rhck.theta_ottimi;
[Lj, Vj] = ContributiCGMY(th(4), th(5), th(6), th(7));
fprintf('Componente di salto del composto-kappa da sola: VS/L = %.4f\n', Vj/Lj);
th = CAL.rhc.theta_ottimi;
[Lj, Vj] = ContributiCGMY(th(4), th(5), th(6), th(7));
fprintf('Componente di salto del composto (kappa = 0)  : VS/L = %.4f\n', Vj/Lj);

% 5) LA VALLE DEL PROFILO: fit simile, meccanismo diverso
% Lungo il profilo in H di rough Heston-kappa (Tab. 6.4) l'errore cambia di
% meno del 3% fra H = 0.3 e H = 0.5, mentre kappa passa da circa 1.7 a 4.9.
VALLE = [];
if isfield(CAL.rhk, 'PROF') && isfield(CAL.rhk, 'H_GRID')
    HG = CAL.rhk.H_GRID(:); PR = CAL.rhk.PROF;  % PR = [nu, rho, kappa, RMSE]
    fprintf('\n=== 5) VALLE DEL PROFILO DI rough Heston-kappa ===\n');
    fprintf('%6s %8s %8s %8s %8s | %10s %10s | %10s %10s | %9s %9s\n', 'H', 'nu', 'rho', ...
        'kappa', 'RMSE%', 'C(.5,1.05)', 'C(.5,1.10)', 'D(.25,1.05)', 'D(.25,1.10)', '|psi|.02', 'slope');
    VALLE = nan(numel(HG), 11);
    for ih = 1:numel(HG)
        thv = [HG(ih), PR(ih,1), PR(ih,2), PR(ih,3)];
        Cv  = PrezzoRHK(thv, CAL, F0, 0.50, F0*[1.05; 1.10], 1);
        Dv  = PrezzoRHK(thv, CAL, F0, 0.25, F0*[1.05; 1.10], 0);
        skv = zeros(numel(T_skew),1);
        for it = 1:numel(T_skew)
            T = T_skew(it); Kpm = F0*exp([-dk; dk]);
            C2 = PrezzoRHK(thv, CAL, F0, T, Kpm, 1);
            skv(it) = (IVcall(F0,Kpm(2),T,C2(2)) - IVcall(F0,Kpm(1),T,C2(1)))/(2*dk);
        end
        cf = polyfit(log(T_skew(mshort)), log(abs(skv(mshort))), 1);
        VALLE(ih,:) = [HG(ih), PR(ih,1:4), Cv(1), Cv(2), Dv(1), Dv(2), abs(skv(1)), cf(1)];
        fprintf('%6.2f %8.4f %8.4f %8.4f %8.3f | %10.4f %10.4f | %10.4f %10.4f | %9.4f %+9.4f\n', VALLE(ih,:));
    end
    ScriviTabellaValle(fullfile(FIG_PATH,'tab_kappa_valle.tex'), VALLE);
else
    fprintf('\n[5) saltato: calib_rhestonkappa_7x15z.mat non contiene PROF e H_GRID]\n');
end

% 6) TABELLE LATEX (gruppo con mean reversion)
INTEST = {'',        'Heston-',          'rough',            'rough Heston-',    'rough Hawkes'; ...
          'Heston',  'curve-$\kappa$',   'Heston-$\kappa$', '$\kappa$ + CGMY', 'Heston'};
ScriviTabellaPrezzi(fullfile(FIG_PATH,'tab_kappa_call.tex'), T_call, Kf, ...
    CallPrices(:,:,CON_K), INTEST, '%.4f');
ScriviTabellaPrezzi(fullfile(FIG_PATH,'tab_kappa_dig.tex'), T_dig, Kf_d, ...
    Digitals(:,:,CON_K), INTEST, '%.4f');
ScriviTabellaSkew(fullfile(FIG_PATH,'tab_kappa_skew.tex'), MODELLI, NOMI_TEX, Slope, CAL);
ScriviTabellaVSL(fullfile(FIG_PATH,'tab_kappa_vsl.tex'), T_vs, ResVS, NOMI_TEX, [2 3 4]);
fprintf('\nTabelle scritte in %s: tab_kappa_call/dig/skew/vsl.tex\n', FIG_PATH);

save(fullfile(RES_PATH,'PricingSuite_kappa.mat'), 'MODELLI','T_call','Kf','CallPrices','CallIV', ...
    'T_dig','Kf_d','Digitals','T_skew','SkewATM','Slope','T_vs','ResVS','VALLE');
fprintf('Risultati salvati in Results/PricingSuite_kappa.mat\n');
toc
diary off

% FUNZIONI LOCALI
function n = N_RIC()
    n = 400;
end

function iv = IVcall(F, K, T, C)
% Volatilita' implicita di una call forward (r = q = 0), OTM via parita'.
    try
        if K >= F
            iv = blkimpv(F, K, 0, T, C, 'Class', 'call');
        else
            iv = blkimpv(F, K, 0, T, max(C - F + K, 1e-12), 'Class', 'put');
        end
    catch
        iv = NaN;
    end
end

function C = PriceCallsModel(nome, CAL, F0, T, K, cp)
% cp = 1 call, cp = 0 digitale cash-or-nothing Q(S_T > K).
    K = K(:); n = numel(K); Z = zeros(n,1); Tv = T*ones(n,1);
    switch nome
        case 'Heston'
            [Xc,pw] = XcPow('Heston');
            C = PrezzaSuperficieSINC_yield_new(F0*ones(n,1), Tv, K, Z, Z, ...
                'Heston', CAL.hes.theta_ottimi, [], [], Xc, 0, pw, cp);
        case 'CGMY'
            [Xc,pw] = XcPow('CGMY');
            C = PrezzaSuperficieSINC_yield_new(F0*ones(n,1), Tv, K, Z, Z, ...
                'CGMY', CAL.cgmy.theta_ottimi, [], [], Xc, 0, pw, cp);
        case 'rHawkesH'
            [Xc,pw] = XcPow('rHawkesH');
            C = PrezzaSuperficieSINC_yield_new(F0*ones(n,1), Tv, K, Z, Z, ...
                'rough Hawkes Heston', CAL.rhh.theta_ottimi, [], N_RIC(), Xc, 0, pw, cp);
        case 'rHeston'
            [Xc,pw] = XcPow('rHeston');
            xig = XiPchip(CAL.rh.T_uniq, CAL.rh.V_fair, T, N_RIC());
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, ...
                'rough Heston', [CAL.rh.theta_ottimi(:).', 0], xig, N_RIC(), Xc, 0, pw, cp);
        case 'HestonCurveK'
            [Xc,pw] = XcPow('rHeston');
            xig = XiPchip(CAL.hck.T_uniq, CAL.hck.V_fair, T, N_RIC());
            th  = CAL.hck.theta_ottimi(:).';  % [kappa, nu, rho]
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, ...
                'Heston curve kappa', th(1:3), xig, N_RIC(), Xc, 0, pw, cp);
        case 'rHestonK'
            [Xc,pw] = XcPow('rHeston');
            xig = XiPchip(CAL.rhk.T_uniq, CAL.rhk.V_fair, T, N_RIC());
            th  = CAL.rhk.theta_ottimi(:).';  % [H, nu, rho, kappa]
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, ...
                'rough Heston kappa', th(1:4), xig, N_RIC(), Xc, 0, pw, cp);
        case 'rH+CGMY'
            [Xc,pw] = XcPow('rH+CGMY');
            th  = CAL.rhc.theta_ottimi(:).';  % [H,nu,rho,C,G,M,Y]
            Lcg = ContributiCGMY(th(4), th(5), th(6), th(7));
            xig = max(XiPchip(CAL.rhc.T_uniq, CAL.rhc.Kvar_mkt, T, N_RIC()) - Lcg, 1e-6);
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, ...
                'rough Heston CGMY', [th(1:3), 0, th(4:7)], xig, N_RIC(), Xc, 0, pw, cp);
        case 'rHK+CGMY'
            [Xc,pw] = XcPow('rH+CGMY');
            th  = CAL.rhck.theta_ottimi(:).';  % [H,nu,rho,C,G,M,Y,kappa]
            Lcg = ContributiCGMY(th(4), th(5), th(6), th(7));
            xig = max(XiPchip(CAL.rhck.T_uniq, CAL.rhck.Kvar_mkt, T, N_RIC()) - Lcg, 1e-6);
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, ...
                'rough Heston kappa CGMY', [th(1:3), th(8), th(4:7)], xig, N_RIC(), Xc, 0, pw, cp);
        otherwise
            error('Modello non riconosciuto: %s', nome);
    end
    C = C(:).';
end

function phi = CharfunModel(nome, CAL, u, T)
    switch nome
        case 'Heston'
            phi = charfun('Heston', CAL.hes.theta_ottimi, [], [], u, T);
        case 'CGMY'
            phi = charfun('CGMY', CAL.cgmy.theta_ottimi, [], [], u, T);
        case 'rHawkesH'
            phi = charfun('rough Hawkes Heston', CAL.rhh.theta_ottimi, [], N_RIC(), u, T);
        case 'rHeston'
            xig = XiPchip(CAL.rh.T_uniq, CAL.rh.V_fair, T, N_RIC());
            phi = charfun('rough Heston', [CAL.rh.theta_ottimi(:).', 0], xig, N_RIC(), u, T);
        case 'HestonCurveK'
            xig = XiPchip(CAL.hck.T_uniq, CAL.hck.V_fair, T, N_RIC());
            th  = CAL.hck.theta_ottimi(:).';
            phi = charfun('Heston curve kappa', th(1:3), xig, N_RIC(), u, T);
        case 'rHestonK'
            xig = XiPchip(CAL.rhk.T_uniq, CAL.rhk.V_fair, T, N_RIC());
            th  = CAL.rhk.theta_ottimi(:).';
            phi = charfun('rough Heston kappa', th(1:4), xig, N_RIC(), u, T);
        case 'rH+CGMY'
            th  = CAL.rhc.theta_ottimi(:).';
            Lcg = ContributiCGMY(th(4), th(5), th(6), th(7));
            xig = max(XiPchip(CAL.rhc.T_uniq, CAL.rhc.Kvar_mkt, T, N_RIC()) - Lcg, 1e-6);
            phi = charfun('rough Heston CGMY', [th(1:3), 0, th(4:7)], xig, N_RIC(), u, T);
        case 'rHK+CGMY'
            th  = CAL.rhck.theta_ottimi(:).';
            Lcg = ContributiCGMY(th(4), th(5), th(6), th(7));
            xig = max(XiPchip(CAL.rhck.T_uniq, CAL.rhck.Kvar_mkt, T, N_RIC()) - Lcg, 1e-6);
            phi = charfun('rough Heston kappa CGMY', [th(1:3), th(8), th(4:7)], xig, N_RIC(), u, T);
        otherwise
            error('Modello non riconosciuto: %s', nome);
    end
end

function VS = VarianceSwapModel(nome, CAL, T)
% E[QV]/T in forma chiusa.
    switch nome
        case 'Heston'
            th = CAL.hes.theta_ottimi; k = th(1); vb = th(2); v0 = th(4);
            VS = vb + (v0 - vb)*(1 - exp(-k*T))/(k*T);
        case 'rHeston'
            VS = WsuT(CAL.rh.T_uniq, CAL.rh.V_fair, T);
        case 'HestonCurveK'
            VS = WsuT(CAL.hck.T_uniq, CAL.hck.V_fair, T);
        case 'rHestonK'
            VS = WsuT(CAL.rhk.T_uniq, CAL.rhk.V_fair, T);
        case 'CGMY'
            th = CAL.cgmy.theta_ottimi; [~, VS] = ContributiCGMY(th(1), th(2), th(3), th(4));
        case 'rH+CGMY'
            VS = VScomposto(CAL.rhc.theta_ottimi, CAL.rhc.T_uniq, CAL.rhc.Kvar_mkt, T);
        case 'rHK+CGMY'
            VS = VScomposto(CAL.rhck.theta_ottimi, CAL.rhck.T_uniq, CAL.rhck.Kvar_mkt, T);
        case 'rHawkesH'
            th = CAL.rhh.theta_ottimi;
            a = th(1); b = th(2); Lam = th(5); s0 = th(6); bet = th(7); z = b*T^a;
            I = s0*T*MittagLeffler(a,2,z) + bet*T^(a+1)*MittagLeffler(a,a+2,z);
            VS = (1 + 2*Lam^2) * I / T;
        otherwise
            error('Modello non riconosciuto: %s', nome);
    end
end

function v = WsuT(T_uniq, Kvar, T)
    W_pp = pchip([0; T_uniq(:)], [0; T_uniq(:).*Kvar(:)]);
    v = ppval(W_pp, T) / T;
end

function VS = VScomposto(th, T_uniq, Kvar, T)
    [Lcg, Vcg] = ContributiCGMY(th(4), th(5), th(6), th(7));
    xi_raw = XiPchip(T_uniq, Kvar, T, N_RIC());
    VS = mean(max(xi_raw - Lcg, 1e-6)) + Vcg;
end

function xig = XiPchip(T_uniq, Kvar, T, N)
    W_pp  = pchip([0; T_uniq(:)], [0; T_uniq(:).*Kvar(:)]);
    h     = 1e-5;
    ti    = ((0:N-1) + 0.5)/N * T;
    tp    = max(ti + h, 0); tm = max(ti - h, 0);
    den   = tp - tm; den(den < 1e-12) = 1e-12;
    xig   = max((ppval(W_pp, tp) - ppval(W_pp, tm)) ./ den, 1e-8);
end

function [Lcg, Vcg] = ContributiCGMY(C, G, M, Y)
    g_m1 = gamma(-Y) * ((M-1)^Y - M^Y + (G+1)^Y - G^Y);
    mu1  = gamma(1-Y) * (M^(Y-1) - G^(Y-1));
    Lcg  = 2*C*(g_m1 - mu1);
    Vcg  = C * gamma(2-Y) * (M^(Y-2) + G^(Y-2));
end

function StampaBlocco(titolo, Kf, X, MODELLI)
    fprintf('\n-- %s --\n%-6s', titolo, 'K/F');
    fprintf('%12s', MODELLI{:}); fprintf('\n');
    for ik = 1:numel(Kf)
        fprintf('%-6.2f', Kf(ik)); fprintf('%12.4f', X(ik,:)); fprintf('\n');
    end
end

function StampaDispersione(Tg, Kf, P, IV)
    for it = 1:numel(Tg)
        fprintf('  T = %.2f:', Tg(it));
        for ik = 1:numel(Kf)
            p = squeeze(P(it,ik,:));
            fprintf('  %.2f: %5.1f%%', Kf(ik), 100*(max(p)-min(p))/min(p));
            if ~isempty(IV)
                v = squeeze(IV(it,ik,:));
                fprintf(' (%.2f vol)', 100*(max(v)-min(v)));
            end
        end
        fprintf('\n');
    end
end

function ScriviTabellaPrezzi(file, Tg, Kf, P, nomi, fmt)
% Stesso formato delle Tab. 7.1 e 7.2, piu' una colonna con lo scarto
% relativo (max - min)/min fra i modelli. nomi: cella 2 x nm con
% l'intestazione su due righe, come la Tab. 7.1.
    nm = size(nomi, 2);
    fid = fopen(file, 'w');
    fprintf(fid, '%% Generato da PricingSuite_Kappa.m il %s\n', datestr(now));
    fprintf(fid, '\\begin{tabular}{rr%sr}\n\\toprule\n', repmat('r', 1, nm));
    fprintf(fid, ' & ');
    for j = 1:nm, fprintf(fid, ' & %s', nomi{1,j}); end
    fprintf(fid, ' & Spread\\\\\n');
    fprintf(fid, '$T$ & $K/F$');
    for j = 1:nm, fprintf(fid, ' & %s', nomi{2,j}); end
    fprintf(fid, ' & (\\%%)\\\\\n\\midrule\n');
    for it = 1:numel(Tg)
        for ik = 1:numel(Kf)
            if ik == 1, fprintf(fid, '%.2f', Tg(it)); end
            fprintf(fid, ' & %.2f', Kf(ik));
            p = squeeze(P(it,ik,:));
            for j = 1:nm, fprintf(fid, [' & ' fmt], p(j)); end
            fprintf(fid, ' & %.1f\\\\\n', 100*(max(p)-min(p))/min(p));
        end
        if it < numel(Tg), fprintf(fid, '\\addlinespace\n'); end
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fclose(fid);
end

function ScriviTabellaSkew(file, MODELLI, nomi, Slope, CAL)
    fid = fopen(file, 'w');
    fprintf(fid, '%% Generato da PricingSuite_Kappa.m il %s\n', datestr(now));
    fprintf(fid, '\\begin{tabular}{lrrl}\n\\toprule\n');
    fprintf(fid, 'Model & Slope & Implied $H$ & Calibrated\\\\\n\\midrule\n');
    for j = 1:numel(MODELLI)
        switch MODELLI{j}
            case 'HestonCurveK', c = sprintf('$H = 1/2$, $\\kappa = %.2f$', CAL.hck.theta_ottimi(1));
            case 'rHestonK',     c = sprintf('$H = %.3f$, $\\kappa = %.2f$', CAL.rhk.theta_ottimi(1), CAL.rhk.theta_ottimi(4));
            case 'rHK+CGMY',     c = sprintf('$H = %.3f$, $\\kappa = %.2f$', CAL.rhck.theta_ottimi(1), CAL.rhck.theta_ottimi(8));
            case 'rHeston',      c = sprintf('$H = %.4f$', CAL.rh.theta_ottimi(1));
            case 'rH+CGMY',      c = sprintf('$H = %.4f$', CAL.rhc.theta_ottimi(1));
            case 'rHawkesH',     c = sprintf('$H_{\\mathrm{eff}} = %.3f$', CAL.rhh.theta_ottimi(1) - 0.5);
            otherwise,           c = '---';
        end
        s = strrep(sprintf('%.3f & %.3f', Slope(j), Slope(j)+0.5), '-', '$-$');
        fprintf(fid, '%s & %s & %s\\\\\n', nomi{j}, s, c);
        if j == 5, fprintf(fid, '\\midrule\n'); end
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fclose(fid);
end

function ScriviTabellaVSL(file, T_vs, ResVS, nomi, quali)
    fid = fopen(file, 'w');
    fprintf(fid, '%% Generato da PricingSuite_Kappa.m il %s\n', datestr(now));
    fprintf(fid, '\\begin{tabular}{lrrrr}\n\\toprule\n');
    fprintf(fid, 'Model & $T$ & $\\sqrt{L}$ (\\%%) & $\\sqrt{VS}$ (\\%%) & VS/L\\\\\n\\midrule\n');
    for q = 1:numel(quali)
        j = quali(q);
        for it = 1:numel(T_vs)
            if it == 1, fprintf(fid, '%s', nomi{j}); end
            fprintf(fid, ' & %.3f & %.2f & %.2f & %.4f\\\\\n', T_vs(it), ...
                100*sqrt(ResVS(it,j,1)), 100*sqrt(ResVS(it,j,2)), ResVS(it,j,3));
        end
        if q < numel(quali), fprintf(fid, '\\addlinespace\n'); end
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fclose(fid);
end

function C = PrezzoRHK(th, CAL, F0, T, K, cp)
% rough Heston-kappa con un vettore di parametri qualsiasi [H, nu, rho,
% kappa].
    K = K(:); n = numel(K); Z = zeros(n,1); Tv = T*ones(n,1);
    [Xc,pw] = XcPow('rHeston');
    xig = XiPchip(CAL.rhk.T_uniq, CAL.rhk.V_fair, T, N_RIC());
    C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, 'rough Heston kappa', th(1:4), xig, N_RIC(), Xc, 0, pw, cp);
    C = C(:).';
end

function ScriviTabellaValle(file, V)
    fid = fopen(file, 'w');
    fprintf(fid, '%% Generato da PricingSuite_Kappa.m il %s\n', datestr(now));
    fprintf(fid, '\\begin{tabular}{rrrrrrr}\n\\toprule\n');
    fprintf(fid, '$H$ & $\\kappa$ & RMSE (\\%%) & Call & Digital & $|\\psi(0.02)|$ & Slope\\\\\n\\midrule\n');
    for i = 1:size(V,1)
        s = strrep(sprintf('%.3f', V(i,11)), '-', '$-$');
        fprintf(fid, '%.2f & %.2f & %.2f & %.4f & %.4f & %.3f & %s\\\\\n', ...
            V(i,1), V(i,4), V(i,5), V(i,7), V(i,9), V(i,10), s);
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fclose(fid);
end
