clearvars; clc;
tic
% VALIDA_PRICER  Validazione numerica del pricer (Appendice C).
%   A) SINC contro inversione di Lewis
%   B) errore della configurazione coarse al variare di (Xc, pow)
%   C) Pade (3,3) a H = 1/2 contro Heston in forma chiusa
%   D) xi0 PCHIP contro xi0 piatta
%   E) martingala e put-call parity
%   F) stabilita' rispetto a H = 1/2 - eps
%   G) convenzione della griglia di xi0
% Output: Results/valida_pricer.mat.

root = fileparts(mfilename('fullpath'));
if isempty(root), root = pwd; end
addpath(fullfile(root,'Functions'), fullfile(root,'Models'), ...
        fullfile(root,'Results'),   fullfile(root,'Data'));

assert(exist('SINC_fastFT','file')==2, 'SINC_fastFT non trovata sul path.');
assert(exist('charfun','file')==2,     'charfun non trovata sul path.');

F0 = 100;
VALIDA = struct();
VALIDA.data = datestr(now);
VALIDA.F0   = F0;

% parametri di prova
% Heston: se esiste la calibrazione la si usa, altrimenti valori tipici SPX.
fh = fullfile(root,'Results','calib_heston_7x15z.mat');
if isfile(fh)
    S = load(fh); pH = S.theta_ottimi(:).';
    VALIDA.fonte_heston = 'calib_heston_7x15z.mat';
else
    pH = [3.16, 0.0615, 1.275, 0.0281, -0.691];
    VALIDA.fonte_heston = 'default';
    fprintf('[avviso] calib_heston_7x15z.mat assente: uso parametri di default.\n');
end

% rough Heston + curva xi0
fr = fullfile(root,'Results','calib_rheston_7x15z_curva.mat');
if isfile(fr)
    R = load(fr);
    pR = R.theta_ottimi(:).';  % [H, nu, rho]
    T_uniq = R.T_uniq(:); V_fair = R.V_fair(:);
    VALIDA.fonte_rheston = 'calib_rheston_7x15z_curva.mat';
else
    pR = [0.0866, 0.341, -0.677];
    T_uniq = [0.0384;0.0849;0.2685;0.5178;0.7671;1.0164;1.5918];
    V_fair = [0.03026;0.03481;0.04336;0.04690;0.04701;0.04731;0.04670];
    VALIDA.fonte_rheston = 'curva sintetica';
    fprintf('[avviso] calib_rheston_7x15z_curva.mat assente: curva sintetica.\n');
end
W_pp = pchip([0; T_uniq], [0; T_uniq.*V_fair]);

% rough Hawkes Heston: serve solo al test E (martingala e parity)
fk = fullfile(root,'Results','calib_rhh_7x15z.mat');
if isfile(fk)
    SK = load(fk); pK = SK.theta_ottimi(:).';
else
    pK = [];
    fprintf('[avviso] calib_rhh_7x15z.mat assente.\n');
end
NRIC_RHH = 400;
VALIDA.pH = pH;  VALIDA.pR = pR;  VALIDA.pK = pK;
VALIDA.T_uniq = T_uniq;  VALIDA.V_fair = V_fair;

Tset = [0.0384 0.0833 0.25 0.50 1.00 1.50];
kset = [0.90 0.95 1.00 1.05];

% A) SINC accurato (Xc=12, pow=12) contro Lewis con CF indipendente
riga('A) SINC accurato (Xc=12, pow=12) vs inversione di Lewis indipendente');
fprintf('%-8s %-6s %14s %14s %12s %12s\n','T','K/F','Lewis','SINC','diff','diff rel');
nA = numel(Tset)*numel(kset);
A_T = nan(nA,1); A_KF = nan(nA,1); A_lew = nan(nA,1);
A_sinc = nan(nA,1); A_rel = nan(nA,1);
emax = 0; iA = 0;
for T = Tset
    for kf = kset
        iA = iA + 1;
        K = F0*kf;
        ref = CallLewisHeston(F0, K, T, pH);
        sinc = SINC_fastFT(F0, T, K, 0, 0, 'Heston', pH, [], [], 12, 0, 12, 1);
        d = sinc - ref; emax = max(emax, abs(d)/max(ref,1e-12));
        A_T(iA)=T; A_KF(iA)=kf; A_lew(iA)=ref; A_sinc(iA)=sinc;
        A_rel(iA)=d/max(ref,1e-12);
        fprintf('%-8.4f %-6.2f %14.8f %14.8f %12.2e %12.2e\n', ...
                T, kf, ref, sinc, d, d/max(ref,1e-12));
    end
end
fprintf('\n  --> errore relativo massimo: %.2e\n', emax);
fprintf(['  Atteso: <1e-5. Se il valore e'' molto piu'' grande SOLO sulle scadenze\n' ...
         '  corte e deep OTM, sospettare il riferimento prima del pricer: la CF di\n' ...
         '  Heston decade esponenzialmente e il troncamento della quadratura e'' la\n' ...
         '  prima cosa da controllare (vedi CallLewisHeston).\n']);
VALIDA.A = struct('T',A_T,'KF',A_KF,'lewis',A_lew,'sinc',A_sinc, ...
                  'err_rel',A_rel,'err_rel_max',emax);

% controllo di convergenza interno al pricer, indipendente dal riferimento
riga('A-bis) Convergenza interna: pow=12 vs pow=14 vs pow=16 (Xc=12)');
fprintf('%-8s %-6s %14s %14s %14s %12s\n','T','K/F','pow=12','pow=14','pow=16','|12-16|');
Tbis = [Tset(1) Tset(3) Tset(6)];
nB2 = numel(Tbis)*numel(kset);
B2_T = nan(nB2,1); B2_KF = nan(nB2,1);
B2_12 = nan(nB2,1); B2_14 = nan(nB2,1); B2_16 = nan(nB2,1);
iB2 = 0;
for T = Tbis
    for kf = kset
        iB2 = iB2 + 1;
        K = F0*kf;
        c12 = SINC_fastFT(F0,T,K,0,0,'Heston',pH,[],[],12,0,12,1);
        c14 = SINC_fastFT(F0,T,K,0,0,'Heston',pH,[],[],12,0,14,1);
        c16 = SINC_fastFT(F0,T,K,0,0,'Heston',pH,[],[],12,0,16,1);
        B2_T(iB2)=T; B2_KF(iB2)=kf; B2_12(iB2)=c12; B2_14(iB2)=c14; B2_16(iB2)=c16;
        fprintf('%-8.4f %-6.2f %14.8f %14.8f %14.8f %12.2e\n', T, kf, c12, c14, c16, abs(c12-c16));
    end
end
fprintf(['\n  Questo test non dipende da alcun riferimento esterno: se pow=14 e\n' ...
         '  pow=16 coincidono, la configurazione accurata e'' convergente.\n']);
VALIDA.Abis = struct('T',B2_T,'KF',B2_KF,'pow12',B2_12,'pow14',B2_14, ...
                     'pow16',B2_16,'max_diff_12_16',max(abs(B2_12-B2_16)));

% B) Errore del pricer coarse al variare di (Xc, pow)
riga('B) Errore in bp di IV del pricer coarse (riferimento: Lewis)');
cfg = [12 8; 5 8; 5 10; 4 9; 6 10];
fprintf('%-8s %-6s', 'T', 'K/F');
for c = 1:size(cfg,1), fprintf('%14s', sprintf('Xc%d/pow%d', cfg(c,1), cfg(c,2))); end
fprintf('\n');
TB = Tset(1:5); KB = [0.95 1.00 1.05];
B_T = nan(numel(TB)*numel(KB),1); B_KF = B_T;
B_bp = nan(numel(TB)*numel(KB), size(cfg,1));
iB = 0;
for T = TB
    for kf = KB
        iB = iB + 1;
        K = F0*kf;
        ivref = ivcall(F0, K, T, CallLewisHeston(F0, K, T, pH));
        B_T(iB)=T; B_KF(iB)=kf;
        fprintf('%-8.4f %-6.2f', T, kf);
        for c = 1:size(cfg,1)
            p = SINC_fastFT(F0, T, K, 0, 0, 'Heston', pH, [], [], cfg(c,1), 0, cfg(c,2), 1);
            if p <= 0
                fprintf('%14s', 'PREZZO<=0');  % B_bp resta NaN
            else
                B_bp(iB,c) = 1e4*(ivcall(F0,K,T,p) - ivref);
                fprintf('%14.1f', B_bp(iB,c));
            end
        end
        fprintf('\n');
    end
end
fprintf(['\n  Lettura: la configurazione coarse va scelta come la piu'' economica\n' ...
         '  che resta sotto il bp. Una colonna con "PREZZO<=0" non e'' utilizzabile\n' ...
         '  in una loss pesata: i residui su quel punto sono artefatti.\n']);
VALIDA.B = struct('cfg',cfg,'T',B_T,'KF',B_KF,'dIV_bp',B_bp, ...
                  'nota','NaN = prezzo non positivo, configurazione inammissibile');

% C) Pade (3,3) al bordo H=1/2 contro Heston in forma chiusa (kappa=0)
riga('C) Errore del Pade (3,3) al bordo H=1/2 (Appendice C.1)');
H_FIX = 0.5 - 1e-6; nuB = 0.60; rhoB = -0.70; Nric = 400;
fprintf('  parametri: H=%.6f  nu=%.3f  rho=%.3f  (nodi Riccati=%d)\n\n', H_FIX, nuB, rhoB, Nric);
fprintf('%-8s %-6s %10s %14s %14s %12s\n','T','K/F','xi flat','Heston chiuso','Pade (3,3)','diff IV bp');
KC = [0.95 1.00 1.05];
nC = numel(T_uniq)*numel(KC);
C_T = nan(nC,1); C_KF = nan(nC,1); C_ex = nan(nC,1); C_pa = nan(nC,1); C_bp = nan(nC,1);
iC = 0;
for it = 1:numel(T_uniq)
    T = T_uniq(it); xi0 = V_fair(it);
    for kf = KC
        iC = iC + 1;
        K = F0*kf;
        % Heston chiuso equivalente: kappa=0, v0=xi0, eta=nu, rho
        pEq = [0, xi0, nuB, xi0, rhoB];
        cEx = SINC_fastFT(F0, T, K, 0, 0, 'Heston', pEq, [], [], 12, 0, 12, 1);
        % Pade con xi flat allo stesso livello
        cPa = SINC_fastFT(F0, T, K, 0, 0, 'rough Heston', [H_FIX nuB rhoB 0], ...
                          xi0*ones(1,Nric), Nric, 12, 0, 12, 1);
        d = 1e4*(ivcall(F0,K,T,cPa) - ivcall(F0,K,T,cEx));
        C_T(iC)=T; C_KF(iC)=kf; C_ex(iC)=cEx; C_pa(iC)=cPa; C_bp(iC)=d;
        fprintf('%-8.4f %-6.2f %10.5f %14.6f %14.6f %12.1f\n', T, kf, xi0, cEx, cPa, d);
    end
end
fprintf(['\n  Nota: questo NON e'' un test di convergenza di griglia. Il Pade converge\n' ...
         '  al proprio valore, che al bordo alpha=1 non e'' la soluzione esatta.\n' ...
         '  Se l''errore qui e'' dello stesso ordine delle differenze di fit fra i\n' ...
         '  modelli, il confronto nested va prezzato con la forma chiusa.\n']);
VALIDA.C = struct('H_FIX',H_FIX,'nu',nuB,'rho',rhoB,'N_riccati',Nric, ...
                  'T',C_T,'KF',C_KF,'heston_chiuso',C_ex,'pade',C_pa, ...
                  'dIV_bp',C_bp,'max_abs_bp',max(abs(C_bp)));

% D) Trattamento del livello: xi da PCHIP contro xi flat
riga('D) Stesso modello, xi da PCHIP vs xi flat allo stesso livello V_fair');
H = pR(1); nu = pR(2); rho = pR(3); Nric = 200;
fprintf('  parametri: H=%.4f  nu=%.4f  rho=%.4f\n\n', H, nu, rho);
fprintf('%-8s %-6s %11s %11s %12s %10s %10s\n', ...
        'T','K/F','IV pchip %','IV flat %','diff bp','min xi0','max xi0');
KD = [0.95 1.00 1.05];
nD = numel(T_uniq)*numel(KD);
D_T = nan(nD,1); D_KF = nan(nD,1); D_ivp = nan(nD,1); D_ivf = nan(nD,1); D_bp = nan(nD,1);
iD = 0;
for it = 1:numel(T_uniq)
    T  = T_uniq(it);
    ti = ((0:Nric-1) + 0.5)/Nric*T;  % punti medi (vedi nota in charfun)
    xp = eval_xi0(ti, W_pp);
    xf = V_fair(it)*ones(1,Nric);
    for kf = KD
        iD = iD + 1;
        K = F0*kf;
        a = SINC_fastFT(F0, T, K, 0, 0, 'rough Heston', [H nu rho 0], xp, Nric, 12, 0, 12, -1);
        b = SINC_fastFT(F0, T, K, 0, 0, 'rough Heston', [H nu rho 0], xf, Nric, 12, 0, 12, -1);
        ia = ivput(F0,K,T,a); ib = ivput(F0,K,T,b);
        D_T(iD)=T; D_KF(iD)=kf; D_ivp(iD)=ia; D_ivf(iD)=ib; D_bp(iD)=1e4*(ia-ib);
        fprintf('%-8.4f %-6.2f %11.3f %11.3f %12.1f %10.5f %10.5f\n', ...
                T, kf, 100*ia, 100*ib, 1e4*(ia-ib), min(xp), max(xp));
    end
end
fprintf(['\n  I due trattamenti hanno la STESSA varianza totale (int xi = T*V_fair)\n' ...
         '  ma forma diversa. Se questa tabella non e'' trascurabile, la frase\n' ...
         '  "share the curve" del confronto nested va corretta o va allineato il codice.\n']);
VALIDA.D = struct('H',H,'nu',nu,'rho',rho,'N_riccati',Nric, ...
                  'T',D_T,'KF',D_KF,'iv_pchip',D_ivp,'iv_flat',D_ivf, ...
                  'dIV_bp',D_bp,'max_abs_bp',max(abs(D_bp)));

% E) Martingala e put-call parity
riga('E) Condizione di martingala phi(-i)=1 e put-call parity');
% Il ramo 'rough Heston' NON e' valutabile in u = -i esatto.
% In a = -i si ha b1 = b2 = b3 = 0 e g0 = g1 = g2 = 0, quindi il
% denominatore dei coefficienti di Pade in dh_Pade33_coeff e' 0/0 -> NaN.
DELTA = 1e-6;
fprintf('  offset usato per la valutazione: u = -i + %.0e\n\n', DELTA);
fprintf('%-22s %-8s %14s %14s\n','modello','T','|phi(-i)-1|','max |C-P-(F-K)|');
mods = { 'Heston',       pH,                       [], [];
         'rough Heston', [H nu rho 0],             [], Nric };
if ~isempty(pK)
    mods(end+1,:) = {'rough Hawkes Heston', pK, [], NRIC_RHH};
end
TE = [0.0833 1.00];
E_mod = strings(size(mods,1)*numel(TE),1);
E_T = nan(size(E_mod)); E_mart = nan(size(E_mod)); E_par = nan(size(E_mod));
iE = 0;
for im = 1:size(mods,1)
    for T = TE
        iE = iE + 1;
        if isempty(mods{im,4})
                        ph = charfun(mods{im,1}, mods{im,2}, [], [], -1i, T);
        elseif strcmp(mods{im,1},'rough Hawkes Heston')
            ph = charfun(mods{im,1}, mods{im,2}, [], mods{im,4}, -1i + DELTA, T);
        else
            ti = ((0:Nric-1) + 0.5)/Nric*T;  % punti medi (vedi nota in charfun)
            ph = charfun(mods{im,1}, mods{im,2}, eval_xi0(ti,W_pp), Nric, -1i + DELTA, T);
        end
        Kv = F0*kset;
        if isempty(mods{im,4})
            c = SINC_fastFT(F0,T,Kv,0,0,mods{im,1},mods{im,2},[],[],12,0,12, 1);
                        p = SINC_fastFT(F0,T,Kv,0,0,mods{im,1},mods{im,2},[],[],12,0,12,-1);
        elseif strcmp(mods{im,1},'rough Hawkes Heston')
            c = SINC_fastFT(F0,T,Kv,0,0,mods{im,1},mods{im,2},[],mods{im,4},12,0,12, 1);
            p = SINC_fastFT(F0,T,Kv,0,0,mods{im,1},mods{im,2},[],mods{im,4},12,0,12,-1);
        else
            xg = eval_xi0(((0:Nric-1) + 0.5)/Nric*T, W_pp);
            c = SINC_fastFT(F0,T,Kv,0,0,mods{im,1},mods{im,2},xg,Nric,12,0,12, 1);
            p = SINC_fastFT(F0,T,Kv,0,0,mods{im,1},mods{im,2},xg,Nric,12,0,12,-1);
        end
        par = max(abs(c(:) - p(:) - (F0 - Kv(:))));
        E_mod(iE) = string(mods{im,1}); E_T(iE) = T;
        E_mart(iE) = abs(ph-1); E_par(iE) = par;
        fprintf('%-22s %-8.4f %14.2e %14.2e\n', mods{im,1}, T, abs(ph-1), par);
    end
end
fprintf(['\n  Atteso: parity a precisione macchina (~1e-13). Per Heston la martingala\n' ...
         '  e'' esatta; per il ramo Pade il residuo deve scalare LINEARMENTE con\n' ...
         '  l''offset (dimezza l''offset -> dimezza il residuo): e'' un limite, non\n' ...
         '  un''identita'' puntuale. In Appendice C va riportato cosi''.\n']);
VALIDA.E = struct('DELTA',DELTA,'modello',E_mod,'T',E_T, ...
                  'martingala',E_mart,'parity',E_par);

% F) Stabilita' rispetto alla scelta di H_FIX = 1/2 - eps
riga('F) Sensibilita'' del caso nested a H_FIX = 1/2 - eps');
epsv = [1e-3 1e-4 1e-5 1e-6 1e-8 1e-10];
TF   = [0.0384 0.25 1.00 1.50];
F_pr = nan(numel(epsv), numel(TF));
fprintf('%-12s', 'eps');
for T = TF, fprintf('%14s', sprintf('ATM T=%.3f',T)); end
fprintf('\n');
for ie = 1:numel(epsv)
    fprintf('%-12.1e', epsv(ie));
    for jt = 1:numel(TF)
        c = SINC_fastFT(F0, TF(jt), F0, 0, 0, 'rough Heston', [0.5-epsv(ie) nuB rhoB 0], ...
                        0.03*ones(1,400), 400, 12, 0, 12, 1);
        F_pr(ie,jt) = c;
        fprintf('%14.6f', c);
    end
    fprintf('\n');
end
fprintf(['\n  Il valore scelto per eps e'' arbitrario: va documentato che i prezzi\n' ...
         '  si stabilizzano ben prima di 1e-6, cosi'' la scelta non e'' un parametro libero.\n']);
VALIDA.F = struct('eps',epsv(:),'T',TF(:),'prezzi',F_pr, ...
                  'max_spread',max(max(F_pr)-min(F_pr)));

% G) Coerenza della convenzione di griglia xi
riga('G) Coerenza della convenzione di griglia xi');

% G1: controllo DETERMINISTICO sui sorgenti
% Il controllo numerico qui sotto vede solo la griglia costruita DENTRO
% questo script.
if exist('ControllaGrigliaXi','file') == 2
    ControllaGrigliaXi(root);
    VALIDA.G1 = 'ControllaGrigliaXi eseguito (vedi log)';
else
    fprintf(2, '  [ControllaGrigliaXi.m non trovato: salto il controllo sui sorgenti]\n');
    VALIDA.G1 = 'ControllaGrigliaXi non trovato';
end

% G2: controllo NUMERICO: ordine di convergenza in N
% ATTENZIONE ALLA SOGLIA.
xi_smooth = @(tt) 0.030*(1 + 0.6*tt);
alpha_rh  = H + 0.5;
att       = 2^(1 + alpha_rh);
Ns        = [100 200 400 800 1600];
TG        = [T_uniq(1) T_uniq(end)];
G_pr = nan(numel(TG), numel(Ns));
G_rt = nan(numel(TG), 1);
fprintf('  alpha = H+1/2 = %.4f  ->  rapporto atteso 2^(1+alpha) = %.2f\n', alpha_rh, att);
fprintf('  (un rapporto vicino a 2 indica griglie disallineate)\n\n');
fprintf('%-9s %-7s %14s %14s %10s\n', 'T', 'N', 'prezzo ATM', '|p(N)-p(2N)|', 'rapporto');
for ig = 1:numel(TG)
    T = TG(ig);
    pv = zeros(1, numel(Ns));
    for a = 1:numel(Ns)
        Na = Ns(a);
        ti = ((0:Na-1) + 0.5)/Na * T;
        pv(a) = SINC_fastFT(F0, T, F0, 0, 0, 'rough Heston', [H nu rho 0], ...
                            xi_smooth(ti), Na, 12, 0, 12, 1);
    end
    G_pr(ig,:) = pv;
    dv = abs(diff(pv));
    for a = 1:numel(Ns)
        if a <= numel(dv), s1 = sprintf('%14.3e', dv(a)); else, s1 = sprintf('%14s','-'); end
        if a <= numel(dv)-1
            s2 = sprintf('%10.2f', dv(a)/max(dv(a+1), realmin));
        else
            s2 = sprintf('%10s','-');
        end
        fprintf('%-9.4f %-7d %14.8f %s %s\n', T, Ns(a), pv(a), s1, s2);
    end
    rt = median(dv(1:end-1) ./ max(dv(2:end), realmin));
    G_rt(ig) = rt;
    if rt > 2.7
        fprintf('  --> rapporto ~%.2f (atteso %.2f): COERENTE.\n\n', rt, att);
    else
        fprintf(2, '  --> rapporto ~%.2f contro %.2f atteso: possibile disallineamento.\n\n', rt, att);
    end
end
VALIDA.G2 = struct('alpha',alpha_rh,'rapporto_atteso',att,'T',TG(:), ...
                   'N',Ns(:),'prezzi',G_pr,'rapporto_misurato',G_rt);

% SALVATAGGIO
% Run_All verifica che questo file esista e sia stato scritto durante lo
% step: senza, marca "MAT MANCANTE" anche a esecuzione riuscita.
cartella_results = fullfile(root, 'Results');
if ~exist(cartella_results, 'dir'), mkdir(cartella_results); end
save(fullfile(cartella_results, 'valida_pricer.mat'), 'VALIDA');
fprintf('\nRisultati dei test salvati in %s\n', ...
    fullfile(cartella_results, 'valida_pricer.mat'));
fprintf(['Sintesi: err.rel.max test A = %.2e | max |dIV| Pade (C) = %.1f bp | ' ...
         'max |dIV| pchip-flat (D) = %.1f bp\n'], ...
    VALIDA.A.err_rel_max, VALIDA.C.max_abs_bp, VALIDA.D.max_abs_bp);

riga('FINE');
toc

% FUNZIONI LOCALI
function riga(t)
    fprintf('\n=========================================================================\n');
    fprintf('%s\n', t);
    fprintf('=========================================================================\n');
end

function C = CallLewisHeston(F, K, T, p)
% Riferimento indipendente: formula di Lewis con CF di Heston riscritta qui
% dentro (NON usa charfun), quadratura adattiva.
% C = F - sqrt(F K)/pi * int_0^inf Re[e^{iuk} phi(u - i/2)] / (u^2 + 1/4) du
    k = log(F/K);
    fun = @(u) real(exp(1i*u*k) .* cfHestonLocale(u - 0.5i, T, p)) ./ (u.^2 + 0.25);
    I    = integral(fun, 0, Inf,  'RelTol', 1e-12, 'AbsTol', 1e-13);
    Ichk = integral(fun, 0, 5000, 'RelTol', 1e-12, 'AbsTol', 1e-13);
    if abs(I - Ichk) > 1e-9 * max(1, abs(I))
        warning('CallLewisHeston:Troncamento', ...
                ['La quadratura non e'' stabile (T=%.4f, K/F=%.3f): ' ...
                 'Inf vs 5000 differiscono di %.2e. Riferimento inaffidabile.'], ...
                T, K/F, abs(I - Ichk));
    end
    C = F - sqrt(F*K)/pi * I;
end

function phi = cfHestonLocale(u, t, p)
% CF di Heston in parametrizzazione "little trap", scritta indipendentemente
% da charfun per avere un vero controllo incrociato.
    kap = p(1); vb = p(2); eta = p(3); v0 = p(4); rho = p(5);
    b = kap - rho*eta*1i*u;
    a = -0.5*u.^2 - 0.5*1i*u;
    d = sqrt(b.^2 - 2*eta^2*a);
    rm = (b - d)/eta^2;  rp = (b + d)/eta^2;
    g = rm./rp;
    e = exp(-d*t);
    Cc = rm*t - (2/eta^2)*log((1 - g.*e)./(1 - g));
    Dd = rm .* (1 - e) ./ (1 - g.*e);
    phi = exp(kap*vb*Cc + v0*Dd);
end

function xi = eval_xi0(tvec, W_pp)
% Identica alla eval_xi0_pchip degli script di calibrazione.
    h = 1e-5;
    tvec = tvec(:);
    tp = max(tvec + h, 0); tm = max(tvec - h, 0);
    dn = tp - tm; dn(dn < 1e-12) = 1e-12;
    xi = (ppval(W_pp, tp) - ppval(W_pp, tm)) ./ dn;
    xi = max(xi, 1e-8);
    xi = xi(:).';
end

function s = ivcall(F, K, T, price)
    try, s = blkimpv(F, K, 0, T, price, 'Class', 'call'); catch, s = NaN; end
end

function s = ivput(F, K, T, price)
    try, s = blkimpv(F, K, 0, T, price, 'Class', 'put'); catch, s = NaN; end
end
fprintf(['\n  I due trattamenti hanno la STESSA varianza totale\n' ...
         '  (int xi = T*V_fair) ma forma diversa. Tutti gli script usano la\n' ...
         '  PCHIP: il confronto qui misura quanto costerebbe la scelta flat,\n' ...
         '  fino a %.0f bp di IV contro un RMSE assoluto del rough Heston di\n' ...
         '  87 bp. Non e'' quindi una convenzione indifferente, e va dichiarata\n' ...
         '  come scelta di modellazione, non come dettaglio implementativo.\n'], ...
         max(abs(D_bp)));
