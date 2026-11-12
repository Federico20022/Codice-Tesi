function phi = charfun(model,params,xi,n,u,t)
% CHARFUN  CF del log-forward (martingala: phi(-i) = 1).
%   phi = charfun(model, params, xi, n, u, t)
%   'GBM'                     [sigma]
%   'Heston'                  [lambda, v_bar, eta, v0, rho]
%   'VG'                      [sigma, theta, nu]
%   'CGMY'                    [C, G, M, Y]
%   'rough Heston'            [H, nu, rho, ~]             + xi, n
%   'Heston curve'            [nu, rho]                   + xi, n
%   'Heston curve kappa'      [kappa, nu, rho]            + xi, n
%   'rough Heston kappa'      [H, nu, rho, kappa]         + xi, n
%   'rough Heston CGMY'       [H, nu, rho, ~, C, G, M, Y] + xi, n
%   'rough Heston kappa CGMY' [H, nu, rho, kappa, C, G, M, Y] + xi, n
%   'rough Hawkes Heston'     [alpha, b, c, rho, Lambda, sigma0_sq, beta]
% xi va campionata sui punti medi ((0:n-1)+0.5)/n*t, come vettore riga:
% solo cosi' flip(xi) da' xi0(T-s) senza bias (vedi ControllaGrigliaXi).

if strcmp(model,'GBM')
    sigma = params(1);
    phi = exp(-0.5*1i*u*t*sigma^2 - 0.5*t*(u*sigma).^2);

elseif strcmp(model,'Heston')
    lambda = params(1);
    v_bar  = params(2);
    eta    = params(3);
    v0     = params(4);
    rho    = params(5);
    beta  = lambda - rho*eta*1i*u;
    alpha = -0.5*u.^2 - 0.5*1i*u;
    zeta  = 0.5*eta^2;
    d = sqrt(beta.^2-4*alpha*zeta);
    rm = (beta-d)/eta^2;
    rp = (beta+d)/eta^2;
    g = rm./rp;
    C = rm*t - 2/eta^2*log((1-g.*exp(-d*t))./(1-g));
    D = rm.*(1-exp(-d*t))./(1-g.*exp(-d*t));
    phi = exp(v_bar*lambda*C + v0*D);

elseif strcmp(model,'VG')
    sigma = params(1);
    theta = params(2);
    nu    = params(3);
    fac1 = exp((1i*u*t)/nu * log(1-nu*theta-0.5*(nu*sigma^2)));
    fac2 = (1./(1+2*nu*(0.5*u*sigma).^2-1i*u*nu*theta)).^(t/nu);
    phi = fac1 .* fac2;

elseif strcmp(model,'CGMY')
    C = params(1);
    G = params(2);
    M = params(3);
    Y = params(4);
    % Y = 1 e' un polo di gamma(-Y).
    g = @(u) gamma(-Y)*((M-1i*u).^Y - M^Y + (G+1i*u).^Y - G^Y);
    phi = exp(C*t*(g(u)-1i*u*g(-1i)));

elseif strcmp(model,'rough Heston')
    H   = params(1);
    nu  = params(2);
    rho = params(3);
    alpha = H + 0.5;
    ti  = ((0:(n-1)) + 0.5) / n * t;  % punti medi: vedi nota in testa
    y   = nu * ti.^alpha;
    [rm,rp,p1,p2,p3,q1,q2,q3] = dh_Pade33_coeff(alpha,rho,u);
    m = size(u,2);
    charFun = zeros(1,m);
    for j = 1:m
        dah = dh_Pade33(y,rm(1,j),rp(1,j),p1(1,j),p2(1,j),p3(1,j), ...
                        q1(1,j),q2(1,j),q3(1,j));
        charFun(1,j) = exp(flip(xi) * dah.' * t./n);
    end
    phi = charFun;

elseif strcmp(model,'Heston curve')
    % Limite alpha = 1 (H = 1/2) del rough Heston in formulazione forward
    % variance: Heston classico SENZA mean reversion, con il livello dato
    % dalla curva xi0(t) esattamente come per il caso rough.
    % A alpha = 1 la Riccati frazionaria degenera in una ODE di Riccati a
    % coefficienti costanti,
    %        h'(s) = 1/2 (h - rm)(h - rp),      h(0) = 0,
    %        h(s) = rm (1 - e^{-aa s}) / (1 - g e^{-aa s}),   g = rm/rp,
    nu  = params(1);
    rho = params(2);
    ti  = ((0:(n-1)) + 0.5) / n * t;  % punti medi: vedi nota in testa
    s   = nu * ti;
    a   = u;
    aa  = sqrt(a.*(a+1i) - rho^2.*a.^2);
    rm  = -1i*rho.*a - aa;
    rp  = -1i*rho.*a + aa;
    g   = rm./rp;
    m = size(u,2);
    charFun = zeros(1,m);
    for j = 1:m
        e   = exp(-aa(1,j)*s);
        h   = rm(1,j) * (1 - e) ./ (1 - g(1,j)*e);
        dah = 0.5*(h - rm(1,j)).*(h - rp(1,j));
        charFun(1,j) = exp(flip(xi) * dah.' * t./n);
    end
    phi = charFun;

elseif strcmp(model,'Heston curve kappa')
    % Heston in formulazione forward variance CON mean reversion kappa
    % (modello affine a varianza forward, Gatheral-Keller-Ressel 2019):
    % V_t = xi0(t) + int_0^t nu e^{-kappa(t-s)} sqrt(V_s) dB_s.
    %     D'      = R(u,D) - kappa D,  D(0) = 0,
    kap = params(1);
    nu  = params(2);
    rho = params(3);
    ti  = ((0:(n-1)) + 0.5) / n * t;  % punti medi: vedi nota in testa
    beta  = kap - rho*nu*1i*u;
    alfa  = -0.5*u.^2 - 0.5*1i*u;
    zeta  = 0.5*nu^2;
    d  = sqrt(beta.^2 - 4*alfa*zeta);
    rm = (beta - d)/nu^2;
    rp = (beta + d)/nu^2;
    g  = rm./rp;
    m = size(u,2);
    charFun = zeros(1,m);
    for j = 1:m
        e   = exp(-d(1,j)*ti);
        D   = rm(1,j) * (1 - e) ./ (1 - g(1,j)*e);
        R   = alfa(1,j) + 1i*rho*nu*u(1,j)*D + zeta*D.^2;
        charFun(1,j) = exp(flip(xi) * R.' * t./n);
    end
    phi = charFun;

elseif strcmp(model,'rough Heston CGMY') || strcmp(model,'rHeston CGMY')
    % Componenti indipendenti: entrambe martingale, quindi il prodotto lo
    % e'.
    phi_rH   = charfun('rough Heston', params(1:4), xi, n, u, t);
    phi_CGMY = charfun('CGMY',         params(5:8), xi, n, u, t);
    phi = phi_rH .* phi_CGMY;

elseif strcmp(model,'rough Heston kappa')
    % rough Heston con mean reversion kappa (forward variance).
    phi = CharFunRHKappa(params(1:4), xi, n, u, t);

elseif strcmp(model,'rough Heston kappa CGMY')
    % Come 'rough Heston CGMY', con la parte rough che ha mean reversion.
    phi_rH   = CharFunRHKappa(params(1:4), xi, n, u, t);
    phi_CGMY = charfun('CGMY', params(5:8), xi, n, u, t);
    phi = phi_rH .* phi_CGMY;

elseif strcmp(model,'rough Hawkes Heston')
    % params = [alpha, b, c, rho, Lambda, sigma0_sq, beta] n = numero passi
    % Eulero (>=400 per T lunghe); xi non usata
    phi = CharFunRHH(params, u, t, n);

else
    error('charfun:ModelloIgnoto', 'Modello non riconosciuto: %s', model);
end
end
