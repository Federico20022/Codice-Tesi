function out = SINC_fastFT(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,pow,CP)
% SINC_FASTFT  Prezzi europei (o digitali) con il metodo SINC.
%   out = SINC_fastFT(S, t, K, IR, DY, model, params, xi, j, Xc, Xm, pow, CP)
% Dominio [Xm-Xc, Xm+Xc], N = 2^pow nodi.
% CP: 1 call, -1 put, 0 digitale cash-or-nothing call.
    % Trasforma K in colonna se non lo è
    if ~isrow(K)
        K = K';
    end

    N   = 2^pow;
    Xl  = -Xc+Xm;
    Xh  =  Xc+Xm;
    wid =  Xh-Xl;

    n   = 1:2:(N/2);
    wn  = n/wid;
    Il  = 2 ./ n;

    f1l = exp(-1i*2*pi*wn*Xm) .* charfun(model,params,xi,j,2*pi*wn    ,t);
    f2l = exp(-1i*2*pi*wn*Xm) .* charfun(model,params,xi,j,2*pi*wn-1i,t);

    q1l = zeros(1,N/2);
    q1l(2:2:end) = f1l .* Il;
    q2l = zeros(1,N/2);
    q2l(2:2:end) = f2l .* Il;

    Iu  = -flip(Il);
    f1u = flip(conj(f1l));
    f2u = flip(conj(f2l));

    q1u = zeros(1,N/2);
    q1u(1:2:end) = f1u .* Iu;
    q2u = zeros(1,N/2);
    q2u(1:2:end) = f2u .* Iu;

    q1  = 1i / (2*pi) * [q1l 0 q1u(1:(N/2-1))];
    q2  = 1i / (2*pi) * [q2l 0 q2u(1:(N/2-1))];

    d1  = fft(q1) + 0.5;
    d2  = fft(q2) + 0.5;

    d11 = d1(1:N/2);
    d12 = d1(N/2+1:N);
    d21 = d2(1:N/2);
    d22 = d2(N/2+1:N);

    dd1 = [d12 d11];
    dd2 = [d22 d21];

    m   = (-N/2):(N/2-1);
    km  = m * wid / N + (IR-DY)*t + Xm;
    Km  = S * exp(km);

    conP = exp(-IR*t).*real(Km .* dd1);
    aonP = exp(-DY*t).*real(S * dd2);

    conP = interp1(Km,conP,K,'spline');  % Usa spline per maggiore regolarità
    aonP = interp1(Km,aonP,K,'spline');

    if CP == 0
        % Digitale cash-or-nothing call, scontata: e^{-rT} Q(S_T > K). dd1
        % sulla griglia e' Q(S_T < K_m): e' la gamba cash-or-nothing della
        % decomposizione del put, esatta, senza differenze finite.
        Pi0 = interp1(Km, real(dd1), K, 'spline');
        out = exp(-IR*t) * (1 - Pi0(:));
        return
    end

    put = conP - aonP;
    call = put+S*exp(-DY*t)-K*exp(-IR*t);

    if CP == 1
        out = call(:);  % Forza vettore colonna
    elseif CP == -1
        out = put(:);  % Forza vettore colonna
    end
end
