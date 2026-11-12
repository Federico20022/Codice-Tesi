function Run_Data(csv, data_str, varargin)
% RUN_DATA  Ripete tutta la pipeline della tesi su un'altra data.
%   Run_Data('20261014SPX_merged.csv', '2026-10-14')
% Il CSV va in Data/ con le colonne
%   Instrument,Root,PutCall,CF_LAST,CF_BID,CF_ASK,STRIKE_PRC,EXPIR_DATE
% e una riga con Root = SPOT.
% Ogni script di calibrazione viene copiato in Tmp_AAAAMMGG/ con CSV, data e
% cartella dei risultati sostituiti. Results/ e Figure/ non vengono toccati.
% Output in Results_AAAAMMGG/ e Figure_AAAAMMGG/.
% Opzioni: 'calibra' (true), 'suite' (true), 'skew' (true).

p = inputParser;
addParameter(p, 'calibra', true);
addParameter(p, 'suite',   true);
addParameter(p, 'skew',    true);
parse(p, varargin{:});
O = p.Results;

root = fileparts(mfilename('fullpath'));
tag  = strrep(data_str, '-', '');
assert(numel(tag) == 8, 'La data va scritta come ''AAAA-MM-GG''.');
datenum(data_str);  % controlla il formato
assert(isfile(fullfile(root, 'Data', csv)), 'Non trovo Data/%s', csv);

resdir = fullfile(root, ['Results_' tag]);
figdir = fullfile(root, ['Figure_' tag]);
tmpdir = fullfile(root, ['Tmp_' tag]);
for d = {resdir, figdir, tmpdir}
    if ~isfolder(d{1}), mkdir(d{1}); end
end
addpath(root, fullfile(root, 'Functions'), fullfile(root, 'Models'), fullfile(root, 'Data'), tmpdir);

logf = fullfile(resdir, ['log_' tag '.txt']);
diary off; if isfile(logf), delete(logf); end
diary(logf);
fig_vis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
pulizia = onCleanup(@() FineLavori(tmpdir, fig_vis));

fprintf('==================================================================\n');
fprintf('RUN DATA %s  |  file %s  |  start %s\n', data_str, csv, datestr(now));
fprintf('==================================================================\n');

SCRIPT = {'Heston_SPX_20260515_def', 'Check_Heston_Curve', 'Heston_Curve_Kappa', ...
          'rough_heston_def', 'Rough_Heston_Kappa', 'CMGY_def', ...
          'Rough_Heston_CMGY_def', 'Rough_Heston_CMGY_Kappa', ...
          'RoughHawkesHeston_SPX_20260515'};
ESITO = cell(0, 3);

% 1. calibrazioni
if O.calibra
    for i = 1:numel(SCRIPT)
        nome = SCRIPT{i};
        fprintf('\n\n########## [%d/%d] %s  (%s) ##########\n', i, numel(SCRIPT), nome, data_str);
        t0 = tic;
        try
            fn = PreparaScript(nome, tmpdir, tag, csv, data_str, resdir, figdir, true);
            feval(fn);
            ESITO(end+1,:) = {nome, toc(t0), 'OK'};  %#ok<AGROW>
        catch ME
            ESITO(end+1,:) = {nome, toc(t0), ['ERRORE: ' ME.message]};  %#ok<AGROW>
            fprintf(2, '\n>>> %s FALLITO: %s\n', nome, ME.message);
        end
        close all force
    end
end

% 2. pricing suite con le varianti kappa
if O.suite
    fprintf('\n\n########## PricingSuite_Kappa (%s) ##########\n', data_str);
    t0 = tic;
    try
        fn = PreparaScript('PricingSuite_Kappa', tmpdir, tag, csv, data_str, resdir, figdir, false);
        diary off; feval(fn); diary(logf);  % la suite apre un suo diary
        ESITO(end+1,:) = {'PricingSuite_Kappa', toc(t0), 'OK'};
    catch ME
        diary(logf);
        ESITO(end+1,:) = {'PricingSuite_Kappa', toc(t0), ['ERRORE: ' ME.message]};
        fprintf(2, '\n>>> PricingSuite_Kappa FALLITA: %s\n', ME.message);
    end
    close all force
end

% 3. skew sotto le due settimane e superficie intera
if O.skew
    fprintf('\n\n########## Test_Fuori_Campione (%s) ##########\n', data_str);
    t0 = tic;
    try
        fn = PreparaScript('Test_Fuori_Campione', tmpdir, tag, csv, data_str, resdir, figdir, false);
        feval(fn);
        ESITO(end+1,:) = {'Test_Fuori_Campione', toc(t0), 'OK'};
    catch ME
        ESITO(end+1,:) = {'Test_Fuori_Campione', toc(t0), ['ERRORE: ' ME.message]};
        fprintf(2, '\n>>> Test_Fuori_Campione FALLITO: %s\n', ME.message);
    end
    close all force
end

% 4. riepilogo
fprintf('\n\n==================================================================\n');
fprintf('ESITO %s\n', data_str);
for i = 1:size(ESITO,1)
    fprintf('  %-34s %8.1f s   %s\n', ESITO{i,1}, ESITO{i,2}, ESITO{i,3});
end
diary off
Riepilogo_Data(resdir, data_str, csv, ESITO);
fprintf('\nLog completo: %s\n', logf);
end

function fn = PreparaScript(nome, tmpdir, tag, csv, data_str, resdir, figdir, calib)
% Copia lo script come FUNZIONE in tmpdir con data, CSV e cartelle della
% nuova data.
    src = which(nome);
    assert(~isempty(src), 'Non trovo %s sul path.', nome);
    txt = fileread(src);
    fn  = ['T' tag '_' nome];

    if calib
        txt = Sostituisci(txt, 'CSV_PATH\s*=\s*''[^'']*''\s*;', ...
            sprintf('CSV_PATH = ''%s'';', csv), 1, nome);
        txt = Sostituisci(txt, 'QUOTE_DATE\s*=\s*datenum\(''[^'']*''\)\s*;', ...
            sprintf('QUOTE_DATE = datenum(''%s'');', data_str), 1, nome);
        txt = Sostituisci(txt, '''\.\.''\s*,\s*''Results''\s*\)', ...
            sprintf('''..'', ''Results_%s'')', tag), -1, nome);
    elseif strcmp(nome, 'PricingSuite_Kappa')
        txt = Sostituisci(txt, 'fullfile\(root_path,\s*''Results''\)', ...
            sprintf('''%s''', resdir), -1, nome);
        txt = Sostituisci(txt, 'fullfile\(root_path,\s*''Figure''\)', ...
            sprintf('''%s''', figdir), -1, nome);
    elseif strcmp(nome, 'Test_Fuori_Campione')
        txt = Sostituisci(txt, 'readtable\(''20260515SPX_merged\.csv''', ...
            sprintf('readtable(''%s''', csv), 1, nome);
        txt = Sostituisci(txt, 'datenum\(''2026-05-15''\)', ...
            sprintf('datenum(''%s'')', data_str), -1, nome);
        txt = Sostituisci(txt, 'fullfile\(root,\s*''Results''\)', ...
            sprintf('''%s''', resdir), -1, nome);
        txt = Sostituisci(txt, 'fullfile\(root,\s*''Figure''\)', ...
            sprintf('''%s''', figdir), -1, nome);
    end

    % corpo dello script = tutto prima della prima funzione locale
    k = regexp(txt, '(^|\n)function\s', 'once');
    if isempty(k), corpo = txt; locali = '';
    else,          corpo = txt(1:k); locali = txt(k+1:end); end

    extra = '';
    if calib
        extra = sprintf([ ...
            '\n%% --- aggiunto da Run_Data: sintesi per il riepilogo ---\n' ...
            'RIASS = struct(''script'', ''%s'', ''theta'', theta_ottimi(:).'');\n' ...
            'vv = {''RMSE_paper'',''RMSE_Vol_abs'',''T_uniq'',''V_fair'',''Kvar_mkt'',''PROF'',''H_GRID'',''KAPPA_GRID'',''S0_feed'',''valid''};\n' ...
            'for iv_ = 1:numel(vv)\n' ...
            '    if exist(vv{iv_}, ''var''), RIASS.(vv{iv_}) = eval(vv{iv_}); end\n' ...
            'end\n' ...
            'save(fullfile(''%s'', ''riass_%s.mat''), ''RIASS'');\n'], nome, resdir, nome);
    end
    out = sprintf('function %s()\n%% Copia temporanea di %s generata da Run_Data per il %s.\n', fn, nome, data_str);
    out = [out corpo extra sprintf('\nend\n\n') locali];
    f = fopen(fullfile(tmpdir, [fn '.m']), 'w');
    fwrite(f, out); fclose(f);
    rehash;
end

function txt = Sostituisci(txt, pat, rep, attese, nome)
% attese = numero esatto di occorrenze, -1 = almeno una.
    n = numel(regexp(txt, pat));
    if attese > 0
        assert(n == attese, 'In %s il pattern %s compare %d volte invece di %d.', nome, pat, n, attese);
    else
        assert(n >= 1, 'In %s non trovo il pattern %s.', nome, pat);
    end
    txt = regexprep(txt, pat, regexprep(rep, '\\', '\\\\'));
end

function FineLavori(tmpdir, fig_vis)
    set(0, 'DefaultFigureVisible', fig_vis);
    diary off
    if isfolder(tmpdir)
        rmpath(tmpdir);
        rmdir(tmpdir, 's');
    end
end
