tic

% PRICINGSUITE_SPX_20260515  Pricing suite dei 5 modelli calibrati.
% Legge i parametri da Results/calib_*.mat.
% 400 nodi Riccati per i modelli rough, dominio SINC da XcPow.
% CGMY: da xi0 si toglie il log-contract dei salti, alla QV si somma la
% loro varianza (ContributiCGMY).
root_path = fileparts(mfilename('fullpath'));
if isempty(root_path), root_path = pwd; end
addpath(fullfile(root_path, 'Functions'));
addpath(fullfile(root_path, 'Models'));
addpath(fullfile(root_path, 'Results'));
addpath(fullfile(root_path, 'Data'));
RES_PATH = fullfile(root_path, 'Results');

% 2. CARICAMENTO CALIBRAZIONI
req = { 'calib_heston_7x15z.mat',        'hes'; ...
        'calib_rheston_7x15z_curva.mat', 'rh';  ...
        'calib_cgmy_7x15z.mat',          'cgmy';...
        'calib_rhcgmy_7x15z.mat',        'rhc'; ...
        'calib_rhh_7x15z.mat',           'rhh' };

CAL = struct();
for i = 1:size(req,1)
    f = fullfile(RES_PATH, req{i,1});
    assert(isfile(f), ...
        'Manca %s in %s: runna la calibrazione corrispondente (o Run_All).', ...
        req{i,1}, RES_PATH);
    CAL.(req{i,2}) = load(f);
end

F0 = 100;
MODELLI = {'Heston','rHeston','CGMY','rH+CGMY','rHawkesH'};
NM = numel(MODELLI);

fprintf('=== PARAMETRI CARICATI (da %s) ===\n', RES_PATH);
fprintf('Heston  : %s\n', mat2str(CAL.hes.theta_ottimi, 5));
fprintf('rHeston : %s\n', mat2str(CAL.rh.theta_ottimi, 5));
fprintf('CGMY    : %s\n', mat2str(CAL.cgmy.theta_ottimi, 5));
fprintf('rH+CGMY : %s\n', mat2str(CAL.rhc.theta_ottimi, 5));
fprintf('rHH     : %s\n', mat2str(CAL.rhh.theta_ottimi, 5));
fns = {'hes','rh','cgmy','rhc','rhh'};
for i = 1:numel(fns)
    C = CAL.(fns{i});
    if isfield(C,'REDUCED_GRID') && C.REDUCED_GRID
        fprintf('  [nota] %s: 7 scadenze x 15 strike\n', req{i,1});
    end
end
fprintf('Nodi Riccati usati dalla suite: %d (come gli stadi exact)\n\n', N_RIC());

% 1) CALL EUROPEE
T_call = [0.10, 0.50, 1.00];
% K/F = 1.20 rimosso: a T=0.10 il CGMY dava un prezzo negativo a Xc=12.
% Non era risoluzione della griglia ma troncamento del dominio: con G=0.05
% la coda sinistra decade come e^{-0.05|x|}.
Kf     = [0.80, 0.90, 0.95, 1.00, 1.05, 1.10];

fprintf('=== 1) CALL EUROPEE (F=100, prezzi forward) ===\n');
CallPrices = zeros(numel(T_call), numel(Kf), NM);
for it = 1:numel(T_call)
    T = T_call(it); K = F0*Kf(:);
    for imod = 1:NM
        CallPrices(it,:,imod) = PriceCallsModel(MODELLI{imod}, CAL, F0, T, K);
    end
    fprintf('\n-- T = %.2f --\n', T);
    fprintf('%-6s', 'K/F');
    fprintf('%12s', MODELLI{:}); fprintf('\n');
    for ik = 1:numel(Kf)
        fprintf('%-6.2f', Kf(ik));
        fprintf('%12.4f', squeeze(CallPrices(it,ik,:)));
        fprintf('\n');
    end
end

% 2) DIGITALI Q(S_T > K)
T_dig = [0.05, 0.25];
Kf_d  = [0.90, 0.95, 0.98, 1.00, 1.02, 1.05, 1.10];

fprintf('\n=== 2) DIGITALI Q(S_T > K) ===\n');
Digitals = zeros(numel(T_dig), numel(Kf_d), NM);
% Lette direttamente dalla gamba cash-or-nothing del SINC (CP = 0 in
% SINC_fastFT): esatte, nessuna differenza finita.
for it = 1:numel(T_dig)
    T = T_dig(it); K = F0*Kf_d(:);
    for imod = 1:NM
        Digitals(it,:,imod) = PriceCallsModel(MODELLI{imod}, CAL, F0, T, K, 0);
    end
    fprintf('\n-- T = %.2f --\n', T);
    fprintf('%-6s', 'K/F');
    fprintf('%12s', MODELLI{:}); fprintf('\n');
    for ik = 1:numel(Kf_d)
        fprintf('%-6.2f', Kf_d(ik));
        fprintf('%12.4f', squeeze(Digitals(it,ik,:)));
        fprintf('\n');
    end
end

% 3) TERM STRUCTURE SKEW ATM
fprintf('\n=== 3) SKEW ATM: fit potenza su T<=0.25 (atteso ~T^(H-1/2)) ===\n');
% L'ultimo anchor delle curve xi0 e' a T~1.59, quindi i punti oltre quel
% valore usano l'estrapolazione della cubica PCHIP terminale.
T_skew = logspace(log10(0.02), log10(2.5), 14);
dk     = 0.005;
SkewATM = zeros(numel(T_skew), NM);
for it = 1:numel(T_skew)
    T = T_skew(it);
    Kpm = F0*exp([-dk; dk]);
    for imod = 1:NM
        C2 = PriceCallsModel(MODELLI{imod}, CAL, F0, T, Kpm);
        ivm = blkimpv(F0, Kpm(1), 0, T, C2(1), 'Class', 'call');
        ivp = blkimpv(F0, Kpm(2), 0, T, C2(2), 'Class', 'call');
        SkewATM(it,imod) = (ivp - ivm)/(2*dk);
    end
end
mshort = T_skew(:) <= 0.25;
fprintf('%-10s  %-14s  %-14s\n', 'Modello', 'slope(logT)', 'H implicito');
for imod = 1:NM
    cf = polyfit(log(T_skew(mshort)), log(abs(SkewATM(mshort,imod))), 1);
    fprintf('%-10s  %+14.4f  %14.4f\n', MODELLI{imod}, cf(1), cf(1)+0.5);
end

% Controllo della Sez. 7.4: la stessa pendenza per il rough Heston con una
% curva xi0 PIATTA.
xi_flat = CAL.rh.V_fair(end-1);  % livello a 1 anno
SkewFlat = zeros(numel(T_skew),1);
[Xc_r, pw_r] = XcPow('rHeston');
for it = 1:numel(T_skew)
    T = T_skew(it); Kpm = F0*exp([-dk; dk]);
    C2 = PrezzaSuperficieSINC_yield(F0, T*ones(2,1), Kpm, [0;0], [0;0], ...
        'rough Heston', [CAL.rh.theta_ottimi, 0], xi_flat*ones(1,N_RIC()), N_RIC(), Xc_r, 0, pw_r, 1);
    ivm = blkimpv(F0, Kpm(1), 0, T, C2(1), 'Class', 'call');
    ivp = blkimpv(F0, Kpm(2), 0, T, C2(2), 'Class', 'call');
    SkewFlat(it) = (ivp - ivm)/(2*dk);
end
cf = polyfit(log(T_skew(mshort)), log(abs(SkewFlat(mshort))), 1);
fprintf('rHeston con xi0 piatta (%.4f): slope %+.4f  (asintotico H-1/2 = %+.4f)\n', ...
    xi_flat, cf(1), CAL.rh.theta_ottimi(1) - 0.5);

figure('Name','ATM skew term structure','Color','w','Position',[100 100 800 520]);
loglog(T_skew, abs(SkewATM), 'LineWidth', 1.6); grid on;
xlabel('T (years)'); ylabel('|ATM skew| = |dIV/dk|');
legend(MODELLI, 'Location', 'southwest');
title('Term structure dello skew ATM (log-log): firma T^{H-1/2}');

% 4) VARIANCE SWAP vs LOG-CONTRACT
fprintf('\n=== 4) VARIANCE SWAP vs LOG-CONTRACT ===\n');
T_vs = [1/12, 0.25, 1.0];
h_cf = 1e-3;
fprintf('%-8s %-10s', 'T', 'Modello');
fprintf('%14s %14s %10s\n', 'sqrt(L)%', 'sqrt(VS)%', 'VS/L');
ResVS = zeros(numel(T_vs), NM, 3);
for it = 1:numel(T_vs)
    T = T_vs(it);
    for imod = 1:NM
        phi_h = CharfunModel(MODELLI{imod}, CAL, h_cf, T);
        EX = imag(phi_h)/h_cf;
        L  = -2/T*EX;
        VS = VarianceSwapModel(MODELLI{imod}, CAL, T, L);
        ResVS(it,imod,:) = [L, VS, VS/L];
        fprintf('%-8.4f %-10s %13.3f%% %13.3f%% %10.4f\n', ...
            T, MODELLI{imod}, 100*sqrt(max(L,0)), 100*sqrt(max(VS,0)), VS/L);
        if strcmp(MODELLI{imod}, 'rHawkesH')
            % L numerico (dalla CF) contro L in forma chiusa
            % (Mittag-Leffler): e' il controllo indipendente del solver, il
            % rapporto VS/L no.
            [~, Lcf] = MomentiRHH(CAL.rhh.theta_ottimi, T);
            fprintf('%-8s %-10s sqrt(L) forma chiusa %.3f%%  (diff %.1f bp)\n', '', '', ...
                100*sqrt(Lcf), 1e4*(sqrt(L) - sqrt(Lcf)));
        end
    end
    fprintf('\n');
end
fprintf(['Lettura: VS/L = 1 se e solo se il processo del prezzo e'' continuo,\n' ...
         'quindi il rapporto e'' una misura del contenuto di salto che non\n' ...
         'dipende dalla parametrizzazione ne'' dai bound imposti in calibrazione.\n' ...
         'Heston e rHeston devono dare 1.0000: e'' il check numerico della suite.\n' ...
         'Per rHH il valore e'' in forma chiusa (Remark 3 paper):\n' ...
         '   VS/L = (1 + 2*Lambda^2) / (1 + 2*Lambda^2/(1+Lambda)).\n' ...
         'Confrontare il CGMY da solo con il CGMY dentro rH+CGMY: la distanza\n' ...
         'fra i due rapporti dice se la componente di salto, una volta messa\n' ...
         'sopra la volatilita'' rough, si comporta ancora da salto.\n']);

figure('Name','VS vs log-contract','Color','w','Position',[150 150 800 450]);
bar(squeeze(ResVS(:,:,3))); grid on;
set(gca,'XTickLabel',arrayfun(@(t) sprintf('T=%.3f',t), T_vs, 'UniformOutput', false));
ylabel('VS / log-contract');
ylim([0.95, max(1.3, max(ResVS(:,:,3),[],'all')*1.05)]);
legend(MODELLI, 'Location', 'northwest');
title('Premio del variance swap sullo strip VIX^2: effetto salti');

save(fullfile(RES_PATH,'PricingSuite_results.mat'), ...
    'T_call','Kf','CallPrices','T_dig','Kf_d','Digitals', ...
    'T_skew','SkewATM','T_vs','ResVS','CAL');
fprintf('\nRisultati salvati in PricingSuite_results.mat\n');
toc

% FUNZIONI LOCALI

function [VS, L] = MomentiRHH(th, T)
% Momenti del rough Hawkes Heston in forma chiusa, senza passare per la CF.
%   E[sigma_t^2]       = s0 E_a(b t^a) + beta t^a E_{a,a+1}(b t^a)
%   int_0^T E[sigma^2] = s0 T E_{a,2}(b T^a) + beta T^(a+1) E_{a,a+2}(b T^a)
%   VS = (1 + 2 Lambda^2) * int / T,   L = (1 + 2 Lambda^2/(1+Lambda)) * int / T
    a = th(1); b = th(2); Lam = th(5); s0 = th(6); bet = th(7);
    z = b*T^a;
    I = s0*T*MittagLeffler(a,2,z) + bet*T^(a+1)*MittagLeffler(a,a+2,z);
    VS = (1 + 2*Lam^2) * I / T;
    L  = (1 + 2*Lam^2/(1+Lam)) * I / T;
end

function n = N_RIC()
% Nodi Riccati per i modelli rough.
    n = 400;
end

function C = PriceCallsModel(nome, CAL, F0, T, K, cp)
% Call europee forward (r=q=0), replicando le chiamate exact delle
% calibrazioni. cp = 1 call (default), cp = 0 digitale cash-or-nothing Q(S_T
% > K).
    if nargin < 6, cp = 1; end
    K = K(:); n = numel(K);
    Z = zeros(n,1); Tv = T*ones(n,1);
    [Xc, pw] = XcPow(nome);
    switch nome
        case 'Heston'
            C = PrezzaSuperficieSINC_yield_new(F0*ones(n,1), Tv, K, Z, Z, ...
                'Heston', CAL.hes.theta_ottimi, [], [], Xc, 0, pw, cp);
        case 'CGMY'
            C = PrezzaSuperficieSINC_yield_new(F0*ones(n,1), Tv, K, Z, Z, ...
                'CGMY', CAL.cgmy.theta_ottimi, [], [], Xc, 0, pw, cp);
        case 'rHawkesH'
            C = PrezzaSuperficieSINC_yield_new(F0*ones(n,1), Tv, K, Z, Z, ...
                'rough Hawkes Heston', CAL.rhh.theta_ottimi, [], N_RIC(), Xc, 0, pw, cp);
        case 'rHeston'
            % Usa PCHIP come nella calibrazione (non flat)
            xig = XiPchip(CAL.rh.T_uniq, CAL.rh.V_fair, T, N_RIC());
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, ...
                'rough Heston', [CAL.rh.theta_ottimi, 0], xig, N_RIC(), Xc, 0, pw, cp);
        case 'rH+CGMY'
            th  = CAL.rhc.theta_ottimi;  % [H,nu,rho,C,G,M,Y]
            Lcg = ContributiCGMY(th(4), th(5), th(6), th(7));
            xig = max(XiPchip(CAL.rhc.T_uniq, CAL.rhc.Kvar_mkt, T, N_RIC()) - Lcg, 1e-6);
            pc  = [th(1:3), 0, th(4:7)];
            C = PrezzaSuperficieSINC_yield(F0, Tv, K, Z, Z, ...
                'rough Heston CGMY', pc, xig, N_RIC(), Xc, 0, pw, cp);
        otherwise
            error('Modello non riconosciuto: %s', nome);
    end
    C = C(:).';
end

function xig = XiPchip(T_uniq, Kvar, T, N)
% xi0(t) = dW/dt con W(T)=T*Kvar(T) interpolata PCHIP.
% Identica a eval_xi0_pchip negli script di calibrazione.
    T_anchor = [0; T_uniq(:)];
    W_anchor = [0; T_uniq(:) .* Kvar(:)];
    W_pp     = pchip(T_anchor, W_anchor);
    h        = 1e-5;
    ti = ((0:N-1) + 0.5)/N * T;
    tplus    = max(ti + h, 0);
    tminus   = max(ti - h, 0);
    denom    = tplus - tminus;
    denom(denom < 1e-12) = 1e-12;
    xig      = (ppval(W_pp, tplus) - ppval(W_pp, tminus)) ./ denom;
    xig      = max(xig, 1e-8);
end

function [Lcg, Vcg] = ContributiCGMY(C, G, M, Y)
% Due quantita' diverse e non intercambiabili per una componente CGMY:
% Lcg = 2 int (e^x - 1 - x) nu(dx)  -> contributo al LOG-CONTRACT.
%   Vcg = int x^2 nu(dx)              -> contributo alla QUADRATIC VARIATION.
%         E' cio' che va sommato per ottenere E[QV]/T.
    g_m1 = gamma(-Y) * ((M-1)^Y - M^Y + (G+1)^Y - G^Y);
    mu1  = gamma(1-Y) * (M^(Y-1) - G^(Y-1));
    Lcg  = 2*C*(g_m1 - mu1);
    Vcg  = C * gamma(2-Y) * (M^(Y-2) + G^(Y-2));
end

function phi = CharfunModel(nome, CAL, u, T)
% Funzione caratteristica per E[X_T] via derivata numerica.
    switch nome
        case 'Heston'
            phi = charfun('Heston', CAL.hes.theta_ottimi, [], [], u, T);
        case 'CGMY'
            phi = charfun('CGMY', CAL.cgmy.theta_ottimi, [], [], u, T);
        case 'rHawkesH'
            phi = charfun('rough Hawkes Heston', CAL.rhh.theta_ottimi, [], N_RIC(), u, T);
        case 'rHeston'
            % Usa PCHIP come nella calibrazione
            xig = XiPchip(CAL.rh.T_uniq, CAL.rh.V_fair, T, N_RIC());
            phi = charfun('rough Heston', [CAL.rh.theta_ottimi, 0], xig, N_RIC(), u, T);
        case 'rH+CGMY'
            th  = CAL.rhc.theta_ottimi;
            Lcg = ContributiCGMY(th(4), th(5), th(6), th(7));
            xig = max(XiPchip(CAL.rhc.T_uniq, CAL.rhc.Kvar_mkt, T, N_RIC()) - Lcg, 1e-6);
            phi = charfun('rough Heston CGMY', [th(1:3), 0, th(4:7)], xig, N_RIC(), u, T);
        otherwise
            error('Modello non riconosciuto: %s', nome);
    end
end

function VS = VarianceSwapModel(nome, CAL, T, L)
% E[QV]/T in forma chiusa per modello.
    switch nome
        case 'Heston'
            th = CAL.hes.theta_ottimi;  % [kappa, vbar, eta, v0, rho]
            k = th(1); vb = th(2); v0 = th(4);
            VS = vb + (v0 - vb)*(1 - exp(-k*T))/(k*T);
        case 'rHeston'
            % VS = W(T)/T con W(T)=T*Kvar(T) interpolata PCHIP
            T_anchor = [0; CAL.rh.T_uniq(:)];
            W_anchor = [0; CAL.rh.T_uniq(:) .* CAL.rh.V_fair(:)];
            W_pp = pchip(T_anchor, W_anchor);
            VS   = ppval(W_pp, T) / T;
        case 'CGMY'
            % Per un Levy puro E[QV]/T = int x^2 nu(dx).
            th = CAL.cgmy.theta_ottimi;  % [C,G,M,Y]
            [~, VS] = ContributiCGMY(th(1), th(2), th(3), th(4));
        case 'rH+CGMY'
            % E[QV]/T = media della forward variance EFFETTIVA + VARIANZA
            % dei salti.
            th = CAL.rhc.theta_ottimi;
            [Lcg, Vcg] = ContributiCGMY(th(4), th(5), th(6), th(7));
            xi_raw = XiPchip(CAL.rhc.T_uniq, CAL.rhc.Kvar_mkt, T, N_RIC());
            if any(xi_raw - Lcg < 1e-6)
                warning(['rH+CGMY: il clip su xi morde a T=%.4f. La varianza ' ...
                         'totale del modello non coincide piu'' con Kvar(T).'], T);
            end
            VS = mean(max(xi_raw - Lcg, 1e-6)) + Vcg;
        case 'rHawkesH'
            % X salta di -Lambda*z con z ~ Exp(1) e intensita' sigma^2:
            % E[QV]/T = (1 + 2*Lambda^2) * E[int sigma^2]/T, con E[int
            % sigma^2] in forma chiusa (MomentiRHH).
            VS = MomentiRHH(CAL.rhh.theta_ottimi, T);
        otherwise
            error('Modello non riconosciuto: %s', nome);
    end
end
