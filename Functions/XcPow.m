function [Xc, pow] = XcPow(modello, stadio)
% XCPOW  Dominio SINC (Xc, pow) per modello e stadio. Fonte unica per
% calibrazioni, diagnostiche e PricingSuite.
%   [Xc, pow] = XcPow(modello)          % stadio 'exact'
%   [Xc, pow] = XcPow(modello, stadio)  % 'fast' | 'exact' | 'guard'
%   modello      fast      exact     guard
%   Heston       6/10      12/12     30/13
%   rHeston      6/10      12/12     30/13
%   CGMY         50/14     50/14     60/15
%   rH+CGMY      18/13     18/13     30/13
%   rHawkesH     12/12     18/13     30/14
% Valori scelti con Diagnostica_Xc (tolleranza 0.5 bp di IV).

    if nargin < 2, stadio = 'exact'; end

    switch modello
        case {'Heston','rHeston'}
            V = struct('fast',[6 10],  'exact',[12 12], 'guard',[30 13]);
        case 'CGMY'
            V = struct('fast',[50 14], 'exact',[50 14], 'guard',[60 15]);
        case {'rH+CGMY','rHCGMY'}
            V = struct('fast',[18 13], 'exact',[18 13], 'guard',[30 13]);
        case 'rHawkesH'
            V = struct('fast',[12 12], 'exact',[18 13], 'guard',[30 14]);
        otherwise
            error('XcPow:modello', ...
                  ['Modello non riconosciuto: %s. Attesi: Heston, rHeston, ' ...
                   'CGMY, rH+CGMY, rHawkesH.'], modello);
    end

    if ~isfield(V, stadio)
        error('XcPow:stadio', ...
              'Stadio non riconosciuto: %s. Attesi: fast, exact, guard.', stadio);
    end

    v   = V.(stadio);
    Xc  = v(1);
    pow = v(2);
end
