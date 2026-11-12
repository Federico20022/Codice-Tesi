clearvars; clc;
% VALIDA_RHKAPPA  Controlli di CharFunRHKappa.
%  1. kappa = 0 contro Pade e Adams.
%  2. H = 1/2 contro la forma chiusa.
%  3. Martingala.
%  4. Differenza di prezzo in bp di IV a kappa = 0.
% Output: Results/valida_rhkappa.mat.
root = fileparts(mfilename('fullpath')); if isempty(root), root = pwd; end
addpath(fullfile(root,'Functions'));
T = 1.0; n = 400; ti = ((0:n-1)+0.5)/n*T;
xf = @(t) 0.03 + 0.02*sqrt(max(t,0)/T);  xi = xf(ti);
u  = [2 5 10 20];
R  = struct();

p0 = [0.087 0.34 -0.676 0];
lmf = log(charfun('rough Heston kappa', p0, xi, n, u, T));
lpd = log(charfun('rough Heston', p0, xi, n, u, T));
lad = LogPhiAdamsRH(p0, u, T, 2000, xf);
R.k0 = [max(abs(lmf-lad)), max(abs(lpd-lad))];
fprintf('kappa=0: max|log phi - Adams|  multi-fattore %.2e   Pade %.2e\n', R.k0);

ph = [0.5 1.4 -0.69 4.9];
R.h12 = max(abs(log(charfun('rough Heston kappa', ph, xi, n, u, T)) - ...
                log(charfun('Heston curve kappa', ph([4 2 3]), xi, n, u, T))));
fprintf('H=1/2:   max|log phi - forma chiusa| = %.2e\n', R.h12);

R.mart = abs(charfun('rough Heston kappa', [0.2 0.8 -0.68 2], xi, n, -1i, T) - 1);
fprintf('martingala |phi(-i)-1| = %.2e\n', R.mart);

K = 100*[0.9 0.95 1 1.05 1.1];
Pm = SINC_fastFT(100, T, K, 0, 0, 'rough Heston kappa', p0, xi, n, 12, 0, 12, 1);
Pp = SINC_fastFT(100, T, K, 0, 0, 'rough Heston',       p0, xi, n, 12, 0, 12, 1);
iv = @(P) arrayfun(@(i) IVcall(100, K(i), T, P(i)), 1:numel(K));
R.dIVbp = 1e4*(iv(Pm) - iv(Pp));
fprintf('kappa=0: IV multi-fattore - IV Pade (bp) = %s\n', mat2str(round(10*R.dIVbp)/10));
save(fullfile(root,'Results','valida_rhkappa.mat'), 'R');

function lp = LogPhiAdamsRH(p, u, T, N, xifun)
    a = p(1)+0.5; nu = p(2); rho = p(3); kap = p(4); h = T/N;
    F = @(v) -0.5*(u.^2+1i*u) + 1i*rho*nu*u.*v + 0.5*nu^2*v.^2;
    G = @(v) F(v) - kap*v;  % lato destro della Volterra
    psi = zeros(N+1,numel(u)); Gv = zeros(N+1,numel(u)); Fv = Gv;
    Gv(1,:) = G(psi(1,:)); Fv(1,:) = F(psi(1,:));
    for m = 0:N-1
        j  = (0:m).';
        bw = h^a/a*((m+1-j).^a - (m-j).^a);  pred = (bw.'*Gv(1:m+1,:))/gamma(a);
        aw = (m-j+2).^(a+1) + (m-j).^(a+1) - 2*(m-j+1).^(a+1);
        aw(1) = m^(a+1) - (m-a)*(m+1)^a;
        psi(m+2,:) = h^a/gamma(a+2)*(G(pred) + aw.'*Gv(1:m+1,:));
        Gv(m+2,:) = G(psi(m+2,:)); Fv(m+2,:) = F(psi(m+2,:));
    end
    s = (0:N)'*h; w = h*ones(N+1,1); w([1 end]) = h/2;
    lp = (w.*xifun(T-s)).' * Fv;
end

function v = IVcall(F, K, T, C)
    lo = 1e-4; hi = 3;
    for it = 1:100
        mid = 0.5*(lo+hi);
        d1 = (log(F/K) + 0.5*mid^2*T)/(mid*sqrt(T)); d2 = d1 - mid*sqrt(T);
        pr = F*0.5*erfc(-d1/sqrt(2)) - K*0.5*erfc(-d2/sqrt(2));
        if pr > C, hi = mid; else, lo = mid; end
    end
    v = 0.5*(lo+hi);
end
