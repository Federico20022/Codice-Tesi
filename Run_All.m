% RUN_ALL  Rigenera da zero tutti i risultati del 15 maggio 2026.
% Nove calibrazioni, diagnostiche, validazioni e PricingSuite, poi un
% riepilogo con tempi ed esiti. Imposta da solo i path. Log in runall.txt.
% Lo stato vive in RunAll_stato, perche' gli script fanno clearvars.
clc;
diary off
if isfile('runall.txt'), delete('runall.txt'); end
diary('runall.txt')
% 1. Path e stato
RunAll_stato('init', RunAll_path());

fprintf('==================================================================\n');
fprintf('RUN ALL  |  start %s\n', datestr(now));
fprintf('==================================================================\n');

% 2. Calibrazioni, diagnostiche e validazioni
while RunAll_stato('k') < RunAll_stato('n')
    RunAll_stato('avanza');
    fprintf('\n\n########## [%d/%d] %s ##########\n', ...
        RunAll_stato('k'), RunAll_stato('n'), RunAll_stato('nome'));
    try
        run(RunAll_stato('nome'));  % lo script salva da se' il suo .mat
        RunAll_stato('fine');
    catch ME
        RunAll_stato('errore', ME.message);
        fprintf(2, '\n>>> FALLITO: %s\n', ME.message);
    end
end

% 3. PricingSuite
RunAll_stato('avanza');
fprintf('\n\n########## [suite] PricingSuite_SPX_20260515 ##########\n');
try
    run('PricingSuite_SPX_20260515');
    RunAll_stato('fine');
catch ME
    RunAll_stato('errore', ME.message);
    fprintf(2, '\n>>> PricingSuite FALLITA: %s\n', ME.message);
end

% 4. Riepilogo
RunAll_riepilogo();
diary off
fprintf('\nLog scritto in: %s\n', fullfile(pwd, 'runall.txt'));
