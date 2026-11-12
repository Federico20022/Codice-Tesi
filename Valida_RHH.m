clearvars; clc;
% VALIDA_RHH  Validazione della CF del rough Hawkes Heston (Appendice C.2).
%  1. Kernel multi-fattore contro il kernel esatto.
%  2. Log-contract contro la forma chiusa con Mittag-Leffler.
%  3. Limite di Heston (alpha = 1, salti spenti).
%  4. Riccati-Volterra con lo schema di Adams.
% Piu' martingala e convergenza nei passi. Output: Results/valida_rhh.mat.
root = fileparts(mfilename('fullpath')); if isempty(root), root = pwd; end
addpath(fullfile(root,'Functions'));
resdir = fullfile(root,'Results');

f = fullfile(resdir,'calib_rhh_7x15z.mat');
if isfile(f)
    C = load(f); p = C.theta_ottimi(:).';
else
    warning('calib_rhh_7x15z.mat non trovato: uso parametri di riferimento.');
    p = [0.88518 -3.11311 1.18604 -0.62739 0.67747 0.026304 0.11903];
end
alpha = p(1); b = p(2); c = p(3); rho = p(4); Lam = p(5); s0 = p(6); bet = p(7);
fprintf('Parametri: alpha=%.4f b=%.4f c=%.4f rho=%.4f Lambda=%.4f sigma0^2=%.5f beta=%.4f\n\n', p);
V = struct();

% TEST 1
fprintf('=== TEST 1: kernel multi-fattore contro kernel esatto ===\n');
NSTEPS = 400; NN = 40;
V.kernel = [];
for a = unique([0.60 0.75 alpha 0.95])
    for T = [1/12 1.0]
        edges = logspace(log10(1e-4/T), log10(20*NSTEPS/T), NN+1); edges(1) = 0;
        cst = 1/(gamma(a)*gamma(1-a));
        p1 = (edges(2:end).^(1-a) - edges(1:end-1).^(1-a))/(1-a);
        p2 = (edges(2:end).^(2-a) - edges(1:end-1).^(2-a))/(2-a);
        m = cst*p1; x = p2./p1;
        s = logspace(log10(T/NSTEPS), log10(T), 200);
        Kmf = sum(m(:).*exp(-x(:)*s), 1);
        Kex = s.^(a-1)/gamma(a);
        err = max(abs(Kmf./Kex - 1));
        V.kernel(end+1,:) = [a T err];
        fprintf('alpha=%.3f T=%.4f  max errore relativo = %.2e\n', a, T, err);
    end
end

% TEST 2
fprintf('\n=== TEST 2: log-contract, CF numerica contro forma chiusa ===\n');
c1 = 0.5 + Lam^2/(1+Lam);
h = 1e-4;
V.logc = [];
for T = [1/12 0.25 1.0 1.5]
    I  = s0*T*MittagLeffler(alpha,2,b*T^alpha) + bet*T^(alpha+1)*MittagLeffler(alpha,alpha+2,b*T^alpha);
    Lcf = 2*c1*I/T;
    for ns = [400 1600]
        phi = CharFunRHH(p, [h -h], T, ns);
        Lnum = -2/T * (imag(phi(1)) - imag(phi(2)))/(2*h);
        V.logc(end+1,:) = [T ns sqrt(Lnum) sqrt(Lcf)];
        fprintf('T=%.4f nsteps=%4d  sqrt(L) CF = %.5f   forma chiusa = %.5f   diff = %+.1f bp\n', ...
            T, ns, sqrt(Lnum), sqrt(Lcf), 1e4*(sqrt(Lnum)-sqrt(Lcf)));
    end
end

% TEST 3
fprintf('\n=== TEST 3: limite alpha = 1 senza salti contro Heston chiuso ===\n');
ph = [1, -3.0, 1.2, -0.65, 0, 0.03, 0.12];  % b<0: kappa=3, theta=0.04
hes = [-ph(2), -ph(7)/ph(2), sqrt(ph(3)), ph(6), ph(4)];
Xc = 18; pw = 13; N = 2^pw; nn = 1:2:(N/2); wn = nn/(2*Xc); u = 2*pi*wn;
V.heston = [];
for T = [1/12 1.0]
    for ns = [400 1600]
        d1 = max(abs(CharFunRHH(ph, u, T, ns, true) - charfun('Heston', hes, [], [], u, T)));
        d2 = max(abs(CharFunRHH(ph, u-1i, T, ns, true) - charfun('Heston', hes, [], [], u-1i, T)));
        V.heston(end+1,:) = [T ns max(d1,d2)];
        fprintf('T=%.4f nsteps=%4d  max|phi_RHH - phi_Heston| sulla griglia SINC = %.2e\n', T, ns, max(d1,d2));
    end
end

% TEST 4
fprintf('\n=== TEST 4: multi-fattore contro Adams frazionario ===\n');
uu = [0.5 2 5 10 20];
V.adams = [];
for T = [0.25 1.0]
    for wsh = [0 1]  % u reale e u - i
        uc = uu - 1i*wsh;
        lp_mf = log(CharFunRHH(p, uc, T, 1600));
        lp_ad = LogPhiAdams(p, 1i*uc, T, 4000);
        d = max(abs(lp_mf - lp_ad));
        V.adams(end+1,:) = [T wsh d];
        fprintf('T=%.2f  %s  max|log phi_MF - log phi_Adams| = %.2e\n', T, ...
            ifelse(wsh, 'u-i', 'u  '), d);
    end
end

% EXTRA
fprintf('\n=== Martingala e convergenza in nsteps (IV ATM) ===\n');
for T = [1/12 1.0]
    fprintf('T=%.4f  |phi(-i)-1| = %.2e\n', T, abs(CharFunRHH(p, -1i + 1e-6, T, 400) - 1));
    K = 100*[0.95 1 1.05];
    P4 = SINC_fastFT(100, T, K, 0, 0, 'rough Hawkes Heston', p, [], 400,  18, 0, 13, 1);
    P16= SINC_fastFT(100, T, K, 0, 0, 'rough Hawkes Heston', p, [], 1600, 18, 0, 13, 1);
    iv4 = arrayfun(@(i) IVcall(100, K(i), T, P4(i)), 1:3);
    iv16= arrayfun(@(i) IVcall(100, K(i), T, P16(i)), 1:3);
    fprintf('          IV 400 vs 1600 passi: max |diff| = %.2f bp\n', 1e4*max(abs(iv4-iv16)));
end

save(fullfile(resdir,'valida_rhh.mat'), 'V', 'p');
fprintf('\nSalvato in Results/valida_rhh.mat\n');

function lp = LogPhiAdams(p, w, T, N)
% log phi con la Riccati-Volterra psi = K * F(w, psi) risolta con Adams
% predictor-corrector frazionario, poi int F(w,psi(s)) g0(T-s) ds (trapezi).
    a = p(1); b = p(2); c = p(3); rho = p(4); Lam = p(5); s0 = p(6); bet = p(7);
    hst = T/N; w = w(:).';
    Ff = @(v) 0.5*(w.^2 - w) + (b + rho*sqrt(c)*w).*v + 0.5*c*v.^2 ...
            + 1./(1 - v + Lam*w) + w*Lam/(1+Lam) - 1 - v;
    psi = zeros(N+1, numel(w)); Fv = zeros(N+1, numel(w)); Fv(1,:) = Ff(psi(1,:));
    k = (0:N).';
    for n = 0:N-1
        j = (0:n).';
        bw = hst^a/a * ((n+1-j).^a - (n-j).^a);  % predictor
        pred = (bw.' * Fv(1:n+1,:)) / gamma(a);
        aw = (n-j+2).^(a+1) + (n-j).^(a+1) - 2*(n-j+1).^(a+1);  % corrector
        aw(1) = n^(a+1) - (n-a)*(n+1)^a;
        psi(n+2,:) = hst^a/gamma(a+2) * (Ff(pred) + aw.' * Fv(1:n+1,:));
        Fv(n+2,:) = Ff(psi(n+2,:));
    end
    s = k*hst;
    g0 = s0 + bet*max(T-s,0).^a/gamma(a+1);
    wts = hst*ones(N+1,1); wts([1 end]) = hst/2;
    lp = (wts.*g0).' * Fv;
end

function iv = IVcall(F, K, T, C)
    lo = 1e-4; hi = 3;
    for it = 1:100
        mid = 0.5*(lo+hi);
        d1 = (log(F/K) + 0.5*mid^2*T)/(mid*sqrt(T)); d2 = d1 - mid*sqrt(T);
        pr = F*0.5*erfc(-d1/sqrt(2)) - K*0.5*erfc(-d2/sqrt(2));
        if pr > C, hi = mid; else, lo = mid; end
    end
    iv = 0.5*(lo+hi);
end

function out = ifelse(c, a, b)
    if c, out = a; else, out = b; end
end
