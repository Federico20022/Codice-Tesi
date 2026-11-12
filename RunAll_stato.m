function out = RunAll_stato(azione, arg)
% RUNALL_STATO  Stato di Run_All in variabili persistent, immuni ai
% clearvars degli script di calibrazione.
%   RunAll_stato('init', root)  path ed elenco script
%   RunAll_stato('avanza')      passa allo step successivo
%   RunAll_stato('k' | 'n' | 'nome' | 'get')
%   RunAll_stato('fine')        tempo ed esito (controlla il .mat)
%   RunAll_stato('errore', msg)
% Le calibrazioni vanno per prime: gli altri script leggono i loro .mat.

persistent S

switch azione

    case 'init'
        rp  = arg;
        S   = struct();
        S.res = fullfile(rp, 'Results');
        S.scripts = { ...
            'Heston_SPX_20260515_def',       'calib_heston_7x15z.mat'; ...
            'Check_Heston_Curve',            'calib_hestoncurve_7x15z.mat'; ...
            'Heston_Curve_Kappa',            'calib_hestoncurvekappa_7x15z.mat'; ...
            'rough_heston_def',              'calib_rheston_7x15z_curva.mat'; ...
            'CMGY_def',                      'calib_cgmy_7x15z.mat'; ...
            'Rough_Heston_CMGY_def',         'calib_rhcgmy_7x15z.mat'; ...
            'RoughHawkesHeston_SPX_20260515','calib_rhh_7x15z.mat'; ...
            'Rough_Heston_Kappa',            'calib_rhestonkappa_7x15z.mat'; ...
            'Rough_Heston_CMGY_Kappa',       'calib_rhcgmykappa_7x15z.mat'; ...
            'Profile_G',                     'profile_G.mat'; ...
            'Multistart',                    'multistart.mat'; ...
            'Diagnostica_Xc',                'diagnostica_Xc.mat'; ...
            'Valida_Pricer',                 'valida_pricer.mat'; ...
            'Valida_RHH',                    'valida_rhh.mat'; ...
            'Valida_RHKappa',                'valida_rhkappa.mat' };

        % Controllo di coerenza: se uno script non e' sul path, meglio
        % accorgersene adesso che dopo un'ora di calibrazioni.
        mancanti = {};
        for i = 1:size(S.scripts,1)
            if isempty(which(S.scripts{i,1}))
                mancanti{end+1} = S.scripts{i,1};  %#ok<AGROW>
            end
        end
        if ~isempty(mancanti)
            error('RunAll_stato:scriptMancanti', ...
                  'Non trovo sul path: %s', strjoin(mancanti, ', '));
        end

        % +1 slot per la PricingSuite, che gira fuori dal ciclo
        S.tempi = zeros(size(S.scripts,1)+1, 1);
        S.esito = strings(size(S.scripts,1)+1, 1);
        S.k     = 0;
        S.t0    = now;
        S.tstart = now;
        out = [];

    case 'avanza'
        S.k = S.k + 1;
        S.tstart = now;
        out = S.k;

    case 'k'
        out = S.k;

    case 'n'
        out = size(S.scripts,1);

    case 'nome'
        out = S.scripts{min(S.k, size(S.scripts,1)), 1};

    case 'get'
        out = S;

    case 'errore'
        S.tempi(S.k) = (now - S.tstart)*86400;
        S.esito(S.k) = "ERRORE: " + string(arg);
        out = [];

    case 'fine'
        S.tempi(S.k) = (now - S.tstart)*86400;
        if S.k > size(S.scripts,1)
            S.esito(S.k) = "OK";  % pricing suite
        else
            f = fullfile(S.res, S.scripts{S.k,2});  % .mat aggiornato?
            if isfile(f)
                d = dir(f);
                if (now - d.datenum)*86400 < S.tempi(S.k) + 120
                    S.esito(S.k) = "OK";
                else
                    S.esito(S.k) = "MAT NON aggiornato (?)";
                end
            else
                S.esito(S.k) = "MAT MANCANTE";
            end
        end
        out = [];

    otherwise
        error('RunAll_stato:azione', 'Azione non riconosciuta: %s', azione);
end
end
