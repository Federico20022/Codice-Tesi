% RUN_CALENDARIO  Ripete l'estensione (Sez. 7.8) sulle quattro date con il
% tempo in giorni di borsa (NYSE / 252) e confronta con il calendario.
% Da lanciare dopo Run_Estensione_Date. Usa Estensione_Borsa.
% Output: Results_AAAAMMGG/*_borsa.*, confronto_calendario.txt.
% Tempo: 4-6 ore in tutto.
DATE = { '20260327SPX_merged.csv', '2026-03-27'
         '20260515SPX_merged.csv', '2026-05-15'
         '20260727SPX_merged.csv', '2026-07-27'
         '20261014SPX_merged.csv', '2026-10-14' };
FAI_PROFILO    = false;  % true aggiunge circa 25 minuti per data
SOLO_CONFRONTO = false;  % true: non ricalibra, rifa' solo il confronto

root = fileparts(mfilename('fullpath')); if isempty(root), root = pwd; end
ESITO = cell(0,3);
if ~SOLO_CONFRONTO
    for i = 1:size(DATE,1)
        t0 = tic;
        try
            Estensione_Borsa(DATE{i,1}, DATE{i,2}, 'tempo', 'borsa', 'profilo', FAI_PROFILO);
            ESITO(end+1,:) = {DATE{i,2}, toc(t0), 'OK'};  %#ok<SAGROW>
        catch ME
            diary off
            ESITO(end+1,:) = {DATE{i,2}, toc(t0), ['ERRORE: ' ME.message]};  %#ok<SAGROW>
            fprintf(2, '\n>>> %s FALLITA: %s\n', DATE{i,2}, ME.message);
        end
        close all force
    end
end

% confronto calendario / borsa
nd = size(DATE,1); V = nan(14, nd);
for j = 1:nd
    tag = strrep(DATE{j,2},'-','');
    fc = fullfile(root, ['Results_' tag], ['estensione_' tag '.mat']);
    if ~isfile(fc) && strcmp(DATE{j,2}, '2026-05-15'), fc = fullfile(root, 'Results', 'estensione.mat'); end  % 15 maggio della tesi
    fb = fullfile(root, ['Results_' tag], ['estensione_' tag '_borsa.mat']);
    if ~isfile(fc) || ~isfile(fb), fprintf('  [salto %s: manca un file]\n', DATE{j,2}); continue; end
    C = load(fc); B = load(fb);
    V(1,j)  = C.R.rhk.theta(1);  V(2,j)  = B.R.rhk.theta(1);
    V(3,j)  = C.R.rhk.theta(4);  V(4,j)  = B.R.rhk.theta(4);
    V(5,j)  = C.R.rhk.M1.tot;    V(6,j)  = B.R.rhk.M1.tot;
    V(7,j)  = C.R.rh.theta(1);   V(8,j)  = B.R.rh.theta(1);
    if isfield(C.R,'rhh') && isfield(B.R,'rhh')
        V(9,j) = C.R.rhh.theta(1) - 0.5; V(10,j) = B.R.rhh.theta(1) - 0.5;
    end
    V(11,j) = C.cf(1);           V(12,j) = B.cf(1);
    okc = ~isnan(C.psiM) & C.psiM < 0; okb = ~isnan(B.psiM) & B.psiM < 0;
    V(13,j) = C.psiM(find(okc,1)) / C.psiM(find(okc,1,'last'));
    V(14,j) = B.psiM(find(okb,1)) / B.psiM(find(okb,1,'last'));
end
RIGHE = { 'rHk: H, calendario',             '%.3f'
          'rHk: H, borsa',                  '%.3f'
          'rHk: kappa, calendario',         '%.2f'
          'rHk: kappa, borsa',              '%.2f'
          'rHk: RMSE %, calendario',        '%.2f'
          'rHk: RMSE %, borsa',             '%.2f'
          'rH: H, calendario',              '%.3f'
          'rH: H, borsa',                   '%.3f'
          'rHH: H_eff, calendario',         '%.3f'
          'rHH: H_eff, borsa',              '%.3f'
          'mercato: pendenza < 2w, cal.',   '%+.2f'
          'mercato: pendenza < 2w, borsa',  '%+.2f'
          'mercato: psi corta/2w, cal.',    '%.2f'
          'mercato: psi corta/2w, borsa',   '%.2f' };
Fmt = @(f,x) subsref({sprintf(f,x), '---'}, struct('type','{}','subs',{{1+isnan(x)}}));

fid = fopen(fullfile(root,'confronto_calendario.txt'),'w');
fprintf(fid, 'ESTENSIONE IN TEMPO DI CALENDARIO E DI BORSA  (%s)\n', datestr(now));
fprintf(fid, '%s\n%-34s', repmat('=',1,34+14*nd), '');
for j = 1:nd, fprintf(fid, '%14s', DATE{j,2}); end
fprintf(fid, '\n');
for r = 1:size(RIGHE,1)
    fprintf(fid, '%-34s', RIGHE{r,1});
    for j = 1:nd, fprintf(fid, '%14s', Fmt(RIGHE{r,2}, V(r,j))); end
    fprintf(fid, '\n');
    if mod(r,2) == 0, fprintf(fid, '\n'); end
end
fprintf(fid, 'Escursione di H fra le date, rough Heston-kappa: calendario %.3f, borsa %.3f\n', ...
    max(V(1,:)) - min(V(1,:)), max(V(2,:)) - min(V(2,:)));
fprintf(fid, 'Escursione di H_eff fra le date, rough Hawkes:   calendario %.3f, borsa %.3f\n', ...
    max(V(9,:)) - min(V(9,:)), max(V(10,:)) - min(V(10,:)));
if ~isempty(ESITO)
    fprintf(fid, '\nESITO\n');
    for i = 1:size(ESITO,1), fprintf(fid, '  %s  %8.0f s  %s\n', ESITO{i,:}); end
end
fclose(fid);
type(fullfile(root,'confronto_calendario.txt'));
