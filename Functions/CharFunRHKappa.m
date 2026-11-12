function phi = CharFunRHKappa(params, xi, n, u, t)
% CHARFUNRHKAPPA  CF del rough Heston con mean reversion kappa, in forma
% forward variance.
%   params = [H, nu, rho, kappa],  xi = xi0 sui punti medi (riga, n)
% Schema multi-fattore + Eulero implicito + estrapolazione di Richardson.
% Limiti: kappa = 0 -> 'rough Heston', H = 1/2 -> 'Heston curve kappa'
% (verificati in Valida_RHKappa).
H = params(1); nu = params(2); rho = params(3); kap = params(4);
alpha = H + 0.5;
u = u(:).'; nw = numel(u);
xi = xi(:).';

r  = max(1, ceil(0.75*nu*max(abs(u))*t / n));  % sottopassi per nodo
% Eulero e' del primo ordine: si combinano due risoluzioni (estrapolazione
% di Richardson, 2*e(2r) - e(r)) e l'errore diventa del secondo ordine.
% Senza questo passo l'errore a 400 nodi arriva a ~10 bp di IV a T = 1.
expo = 2*Esponente(alpha, nu, rho, kap, xi, n, 2*r, u, t) ...
       - Esponente(alpha, nu, rho, kap, xi, n, r, u, t);
expo = min(real(expo), 0) + 1i*imag(expo);  % |phi| <= 1 sulla striscia
phi  = exp(expo);
end

function expo = Esponente(alpha, nu, rho, kap, xi, n, r, u, t)
% Esponente log phi con n*r passi di Eulero e xi0 costante a tratti.
nw = numel(u);
ns = n*r;
xis = kron(xi, ones(1,r));  % xi0 sui punti medi fini
NN = 80;  % 40 nodi danno ~1e-2 su log phi a u = 20 contro Adams
if alpha > 0.9999
    m = 1; x = 0;
else
    edges = logspace(log10(1e-4/t), log10(20*ns/t), NN+1);
    edges(1) = 0;
    cst = 1/(gamma(alpha)*gamma(1-alpha));
    p1  = (edges(2:end).^(1-alpha) - edges(1:end-1).^(1-alpha))/(1-alpha);
    p2  = (edges(2:end).^(2-alpha) - edges(1:end-1).^(2-alpha))/(2-alpha);
    m   = cst*p1; x = p2./p1;
end
m = m(:); x = x(:);
dt   = t/ns;
A    = -0.5*(u.^2 + 1i*u);
B    = 1i*rho*nu*u - kap;
inv_ = 1./(1 + x*dt);
S2   = sum(m.*inv_);
minv = (m.*inv_).';
psij = zeros(numel(x), nw);
psi  = zeros(1, nw);
Fprev = A;  % F(u, h(0)) con h(0) = 0
expo  = zeros(1, nw);
xflip = flip(xis);  % xi0(T - s) sui punti medi
for k = 1:ns
    Btot    = B + 0.5*nu^2*psi;
    S1      = minv * (psij + dt*A);
    psi_new = S1 ./ (1 - dt*Btot*S2);
    psij    = inv_ .* (psij + dt*(A + Btot.*psi_new));
    psi     = psi_new;
    Fnew    = A + 1i*rho*nu*u.*psi + 0.5*nu^2*psi.^2;
    expo    = expo + xflip(k) * 0.5*(Fprev + Fnew) * dt;
    Fprev   = Fnew;
end
end
