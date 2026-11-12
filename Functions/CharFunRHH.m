function phi = CharFunRHH(params, u, t, nsteps, nojump)
% CHARFUNRHH  CF del rough Hawkes Heston (Bondi, Pulido, Scotti 2022).
%   params = [alpha, b, c, rho, Lambda, sigma0_sq, beta]
%   kernel K(s) = s^(alpha-1)/Gamma(alpha), salti esponenziali di rate 1,
%   curva g0(s) = sigma0^2 + beta*s^alpha/Gamma(alpha+1).
% Riccati-Volterra con schema multi-fattore (Abi Jaber, El Euch 2019, primo
% intervallo da 0) ed Eulero linearizzato-implicito.
% nojump = true spegne i salti (solo per il test di Heston in Valida_RHH).
if nargin < 5, nojump = false; end
alpha = params(1); b   = params(2); c    = params(3); rho = params(4);
Lam   = params(5); sg0 = params(6); bet  = params(7);
if nojump, Lam = 0; end
J = double(~nojump);  % interruttore dei termini di salto

w  = 1i*u(:).';  % E[e^{wX}], Re w in [0,1]
nw = length(w);

% substepping adattivo: la stiffness della Riccati scala come sqrt(c)*|u|*T;
% senza questo, corner con c grande x T lungo x frequenze SINC alte fanno
% oscillare lo schema.
nsteps = max(nsteps, ceil(0.75*sqrt(c)*max(abs(u(:)))*t));

% nodi multi-factor
NN = 40;
if alpha > 0.9999
    m = 1; x = 0;  % alpha=1: kernel costante, un solo fattore
else
    edges = logspace(log10(1e-4/t), log10(20*nsteps/t), NN+1);
    edges(1) = 0;  % il primo fattore copre [0, e1]
    cst = 1/(gamma(alpha)*gamma(1-alpha));
    p1  = (edges(2:end).^(1-alpha) - edges(1:end-1).^(1-alpha))/(1-alpha);
    p2  = (edges(2:end).^(2-alpha) - edges(1:end-1).^(2-alpha))/(2-alpha);
    m   = cst*p1;  % pesi
    x   = p2./p1;  % rate
end
m = m(:); x = x(:);

% Eulero linearizzato-implicito sulla Riccati-Volterra
dt   = t/nsteps;
A    = 0.5*(w.^2 - w) + J*(w*Lam/(1+Lam) - 1);  % costante + salti(v=0)
B    = b + rho*sqrt(c)*w - J;  % coeff. lineare in v
inv_ = 1./(1 + x*dt);  % (NN x 1)
S2   = sum(m.*inv_);
minv = (m.*inv_).';  % (1 x NN)

psij = zeros(length(x), nw);
psi  = zeros(1, nw);
Fk   = zeros(nsteps+1, nw);
Fk(1,:) = FfunRHH(w, psi, b, c, rho, Lam, J);
for k = 1:nsteps
    rat     = J./(1 - psi + Lam*w);  % salti (esplicito, bounded)
    rhs     = A + rat;  % (1 x nw)
    Btot    = B + 0.5*c*psi;  % quadratica linearizzata
    S1      = minv * (psij + dt*rhs);  % (1 x nw)
    psi_new = S1 ./ (1 - dt*Btot*S2);
    psij    = inv_ .* (psij + dt*(rhs + Btot.*psi_new));
    psi     = psi_new;
    Fk(k+1,:) = FfunRHH(w, psi, b, c, rho, Lam, J);
end

% trasformata: int_0^T F(w, psi(s)) g0(T-s) ds (trapezi)
s   = (0:nsteps)*dt;
g0  = sg0 + bet*max(t-s,0).^alpha/gamma(alpha+1);  % clip: t-s puo' fare -1ulp
wts = dt*ones(1, nsteps+1); wts(1) = wts(1)/2; wts(end) = wts(end)/2;
expo = (wts.*g0) * Fk;
% clamp teorico: sulla striscia Re w in [0,1] vale |Phi| <= 1, quindi
% Re(exponent) <= 0. Attivo solo su residui numerici a parametri estremi; ai
% riferimenti non tocca nulla.
expo = min(real(expo), 0) + 1i*imag(expo);
phi  = exp(expo);
end

function out = FfunRHH(w, v, b, c, rho, Lam, J)
% F(u,v) del paper (eq. 7) con nu(dz)=e^{-z}dz in forma chiusa: int
% [e^{(v-Lam*u)z} - u(e^{-Lam z}-1) - 1 - vz] e^{-z} dz = 1/(1-v+Lam*u) +
% u*Lam/(1+Lam) - 1 - v Sanity: F(1,0)=0 (martingala).
out = 0.5*(w.^2 - w) + (b + rho*sqrt(c)*w).*v + 0.5*c*v.^2 ...
    + J*(1./(1 - v + Lam*w) + w*Lam/(1+Lam) - 1 - v);
end
