function [Xc_sel, pow_sel, INFO] = SelezionaXc(prezzaFn, cfFn, T, K, F, opts)
% SELEZIONAXC  Sceglie il semi-dominio Xc e N = 2^pow del pricer SINC.
%   [Xc, pow, INFO] = SelezionaXc(prezzaFn, cfFn, T, K, F, opts)
% Criterio: primo Xc della scala per cui la IV resta entro TOL_BP dal
% riferimento (il gradino piu' largo). N si ricava da Xc con dx <= DX_MAX.
% La massa di coda e' riportata solo come diagnostica.
%   prezzaFn  @(Xc,pow) -> prezzi forward sugli strike K
%   cfFn      @(u) -> CF del log-forward a T ([] per saltare la diagnostica)
%   opts      .XC_LADDER, .DX_MAX, .POW_MIN, .POW_MAX, .TOL_BP, .CP, .VERBOSE
% Va chiamata una volta a parametri fissati, non dentro lsqnonlin.

    if nargin < 6, opts = struct(); end
    if ~isfield(opts,'XC_LADDER'), opts.XC_LADDER = [6 9 12 18 24 30 40 50]; end
    if ~isfield(opts,'DX_MAX'),    opts.DX_MAX    = 0.008; end
    if ~isfield(opts,'POW_MIN'),   opts.POW_MIN   = 10;    end
    if ~isfield(opts,'POW_MAX'),   opts.POW_MAX   = 16;    end
    if ~isfield(opts,'TOL_BP'),    opts.TOL_BP    = 0.5;   end
    if ~isfield(opts,'CP'),        opts.CP        = -1;    end
    if ~isfield(opts,'VERBOSE'),   opts.VERBOSE   = true;  end

    K   = K(:);
    lad = opts.XC_LADDER(:)';
    nL  = numel(lad);

    % N da Xc: dx = 2*Xc/N <= DX_MAX
    pows = zeros(1,nL);
    for i = 1:nL
        p = ceil(log2(2*lad(i)/opts.DX_MAX));
        pows(i) = min(max(p, opts.POW_MIN), opts.POW_MAX);
    end
    dxs = 2*lad ./ 2.^pows;
    if any(dxs > opts.DX_MAX*(1+1e-9))
        warning(['SelezionaXc: POW_MAX=%d non basta per dx<=%.4g ai gradini ' ...
                 'piu'' larghi. La risoluzione degrada sopra Xc=%g.'], ...
                 opts.POW_MAX, opts.DX_MAX, lad(find(dxs > opts.DX_MAX,1)));
    end

    % prezzi e IV su ogni gradino
    IVs = nan(numel(K), nL);
    for i = 1:nL
        P = prezzaFn(lad(i), pows(i));
        IVs(:,i) = InvertiIV(F, K, T, P(:), opts.CP);
    end

    % Riferimento: il gradino piu' largo.
    ok = all(isfinite(IVs), 2);
    n_persi = sum(~ok);
    dIV_bp = nan(1, nL);
    for i = 1:nL
        if any(~isfinite(IVs(:,i)))
            dIV_bp(i) = Inf;  % prezzo non invertibile: gradino inammissibile
        else
            dIV_bp(i) = 1e4 * max(abs(IVs(ok,i) - IVs(ok,end)));
        end
    end

    % massa di coda (diagnostica)
    pL = nan(1,nL); pR = nan(1,nL);
    if ~isempty(cfFn)
        for i = 1:nL
            [pL(i), pR(i)] = MassaCoda(cfFn, lad(i), pows(i));
        end
    end

    % selezione
    idx = find(dIV_bp <= opts.TOL_BP, 1, 'first');
    if isempty(idx)
        idx = nL;
        warning(['SelezionaXc: nessun gradino sotto %.2f bp. Selezionato il ' ...
                 'piu'' largo (Xc=%g). Estendere XC_LADDER.'], opts.TOL_BP, lad(idx));
    end
    Xc_sel  = lad(idx);
    pow_sel = pows(idx);

    if opts.VERBOSE
        fprintf('  %-6s %-5s %-9s %-11s %-11s %-12s\n', ...
            'Xc','pow','dx','coda sx','coda dx','max dIV (bp)');
        for i = 1:nL
            sel = ''; if i == idx, sel = '  <-- scelto'; end
            if isinf(dIV_bp(i)), sd = '     n.d.';
            else,                sd = sprintf('%9.3f', dIV_bp(i)); end
            fprintf('  %-6g %-5d %-9.5f %-11.3e %-11.3e %s%s\n', ...
                lad(i), pows(i), dxs(i), pL(i), pR(i), sd, sel);
        end
        if n_persi > 0
            fprintf(['  [%d strike con IV non invertibile su almeno un gradino: ' ...
                     'esclusi dal max, ma il gradino e'' marcato n.d.]\n'], n_persi);
        end
    end

    INFO = struct('Xc',lad,'pow',pows,'dx',dxs,'coda_sx',pL,'coda_dx',pR, ...
                  'dIV_bp',dIV_bp,'idx',idx,'IV',IVs,'tol_bp',opts.TOL_BP, ...
                  'n_non_invertibili',n_persi);
end

function [pL, pR] = MassaCoda(cfFn, Xc, pow)
% Q(X < -Xc) e Q(X > +Xc) lette dall'array digitale che il SINC produce.
% Stessa costruzione di SINC_fastFT con Xm = 0 e r = q = 0: d1 e' la
% funzione di ripartizione del log-forward sui nodi della griglia.
    N   = 2^pow;
    wid = 2*Xc;
    n   = 1:2:(N/2);
    wn  = n / wid;
    Il  = 2 ./ n;

    f1l = cfFn(2*pi*wn);
    q1l = zeros(1, N/2);  q1l(2:2:end) = f1l .* Il;
    Iu  = -flip(Il);
    f1u = flip(conj(f1l));
    q1u = zeros(1, N/2);  q1u(1:2:end) = f1u .* Iu;

    q1  = 1i/(2*pi) * [q1l, 0, q1u(1:N/2-1)];
    d1  = fft(q1) + 0.5;
    dd1 = [d1(N/2+1:end), d1(1:N/2)];

    pL = max(real(dd1(1)),   0);  % Q(X < -Xc)
    pR = max(1 - real(dd1(end)), 0);  % Q(X > +Xc - dx)
end

function IV = InvertiIV(F, K, T, P, CP)
% Inversione in spazio forward.
    IV = nan(size(K));
    for i = 1:numel(K)
        try
            if CP == -1
                if K(i) <= F
                    IV(i) = blkimpv(F, K(i), 0, T, P(i), 'Class','put');
                else
                    IV(i) = blkimpv(F, K(i), 0, T, ...
                                    max(P(i) + F - K(i), 1e-12), 'Class','call');
                end
            else
                if K(i) >= F
                    IV(i) = blkimpv(F, K(i), 0, T, P(i), 'Class','call');
                else
                    IV(i) = blkimpv(F, K(i), 0, T, ...
                                    max(P(i) - F + K(i), 1e-12), 'Class','put');
                end
            end
        catch
            IV(i) = NaN;
        end
    end
    % isreal() NON e' element-wise: restituisce un solo logico per l'intero
    % array, quindi ~isreal(IV) propagherebbe un singolo valore complesso su
    % tutto il vettore e azzererebbe l'intera scadenza. imag() invece opera
    % elemento per elemento.
    IV(imag(IV) ~= 0 | real(IV) <= 0) = NaN;
    IV = real(IV);
end
