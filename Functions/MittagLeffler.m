function E = MittagLeffler(a, b, z)
% MITTAGLEFFLER  E_{a,b}(z) = sum_k z^k / Gamma(a*k + b), per serie.
% Adatta a a in (1/2,1] e |z| <= ~8. Per |z| > 10 perde precisione.
if any(abs(z(:)) > 10)
    warning('MittagLeffler:precisione', ...
        '|z| = %.1f: la serie perde precisione, usare un metodo integrale.', max(abs(z(:))));
end
K = 300;
k = (0:K-1).';
lg = gammaln(a*k + b);  % Gamma(a k + b) > 0 per a,b > 0
E = zeros(size(z));
for i = 1:numel(z)
    zi = z(i);
    if zi == 0, E(i) = 1/gamma(b); continue; end
    % termini in forma logaritmica per evitare overflow di z^k e Gamma
    t = exp(k*log(abs(zi)) - lg) .* sign(zi).^k;
    E(i) = sum(t);
end
end
