% RUN_ESTENSIONE_DATE  Estensione alle scadenze corte (Sez. 7.8) sulle altre
% tre date e confronto con il 15 maggio. Da lanciare dopo Run_Data.
% Output: Results_AAAAMMGG/estensione_*, confronto_estensione.txt,
% Figure/tab_estensione_date.tex. Tempo: circa 4.5 ore.
DATE = { '20260327SPX_merged.csv', '2026-03-27'
         '20260727SPX_merged.csv', '2026-07-27'
         '20261014SPX_merged.csv', '2026-10-14' };
SOLO_CONFRONTO = false;  % true: non ricalibra, rifa' solo il confronto

root = fileparts(mfilename('fullpath')); if isempty(root), root = pwd; end
Fmt    = @(f,x) subsref({sprintf(f,x), '---'}, struct('type','{}','subs',{{1+isnan(x)}}));
FmtTex = @(f,x) subsref({strrep(sprintf(f,x),'-','$-$'), '---'}, struct('type','{}','subs',{{1+isnan(x)}}));
ESITO = cell(0,3);
if ~SOLO_CONFRONTO
    for i = 1:size(DATE,1)
        t0 = tic;
        try
            Estensione_Data(DATE{i,1}, DATE{i,2});
            ESITO(end+1,:) = {DATE{i,2}, toc(t0), 'OK'};  %#ok<SAGROW>
        catch ME
            diary off
            ESITO(end+1,:) = {DATE{i,2}, toc(t0), ['ERRORE: ' ME.message]};  %#ok<SAGROW>
            fprintf(2, '\n>>> %s FALLITA: %s\n', DATE{i,2}, ME.message);
        end
        close all force
    end
end

% confronto
F = {fullfile(root,'Results','estensione.mat'), '2026-05-15'};
for i = 1:size(DATE,1)
    tag = strrep(DATE{i,2},'-','');
    F(end+1,:) = {fullfile(root,['Results_' tag],['estensione_' tag '.mat']), DATE{i,2}};  %#ok<SAGROW>
end
[~, o] = sort(datenum(F(:,2))); F = F(o,:);
C = struct('data',{},'S',{});
for i = 1:size(F,1)
    if isfile(F{i,1}), C(end+1) = struct('data',F{i,2},'S',load(F{i,1}));  %#ok<SAGROW>
    else, fprintf('  [salto %s: manca %s]\n', F{i,2}, F{i,1}); end
end
assert(~isempty(C), 'Nessun file di estensione trovato.');

nd = numel(C); V = nan(18, nd);
for j = 1:nd
    S = C(j).S; R = S.R; ok = ~isnan(S.psiM) & S.psiM < 0;
    V(1,j)  = numel(S.T_cor);
    V(2,j)  = S.cf(1);
    V(3,j)  = S.psiM(find(ok,1)) / S.psiM(find(ok,1,'last'));
    V(4,j)  = R.rhk.theta0(1);  V(5,j) = R.rhk.theta0(4);
    V(6,j)  = R.rhk.theta(1);   V(7,j) = R.rhk.theta(4);
    V(8,j)  = R.rhk.M1.tot;     V(9,j) = R.rhk.M1.cor;
    V(10,j) = R.rh.theta(1);    V(11,j) = R.rh.M1.tot;
    V(12,j) = R.hck.theta(1);   V(13,j) = R.hck.M1.tot;
    if isfield(R,'rhh'), V(14,j) = R.rhh.theta(1) - 0.5; V(15,j) = R.rhh.M1.tot; end
    if isfield(S,'PRO') && any(~isnan(S.PRO(:,5)))
        [mn, im] = min(S.PRO(:,5)); V(16,j) = S.PRO(im,1);
        V(17,j) = 100*(S.PRO(S.PRO(:,1)==0.50,5)/mn - 1);
        V(18,j) = 100*(S.PRO(S.PRO(:,1)==0.10,5)/mn - 1);
    end
end
RIGHE = { 'scadenze corte aggiunte',               '%.0f',  'Short maturities added'
          'mercato: pendenza skew < 2 sett.',      '%+.2f', 'Market skew slope below two weeks'
          'mercato: rapporto psi corta / psi 2w',  '%.2f',  'Market ratio $\psi$ shortest over two weeks'
          'rHk standard: H',                       '%.3f',  '$H$, rough Heston-$\kappa$, seven maturities'
          'rHk standard: kappa',                   '%.2f',  '$\kappa$, rough Heston-$\kappa$, seven maturities'
          'rHk esteso: H',                         '%.3f',  '$H$, rough Heston-$\kappa$, extended'
          'rHk esteso: kappa',                     '%.2f',  '$\kappa$, rough Heston-$\kappa$, extended'
          'rHk esteso: RMSE tutte %',              '%.2f',  'RMSE, rough Heston-$\kappa$, all (\%)'
          'rHk esteso: RMSE corte %',              '%.2f',  'RMSE, rough Heston-$\kappa$, short (\%)'
          'rH esteso: H',                          '%.3f',  '$H$, rough Heston, extended'
          'rH esteso: RMSE tutte %',               '%.2f',  'RMSE, rough Heston, all (\%)'
          'Hck esteso: kappa',                     '%.2f',  '$\kappa$, Heston-curve-$\kappa$, extended'
          'Hck esteso: RMSE tutte %',              '%.2f',  'RMSE, Heston-curve-$\kappa$, all (\%)'
          'rHH esteso: H_eff',                     '%.3f',  '$H_{\mathrm{eff}}$, rough Hawkes Heston, extended'
          'rHH esteso: RMSE tutte %',              '%.2f',  'RMSE, rough Hawkes Heston, all (\%)'
          'profilo esteso: H al minimo',           '%.2f',  'Profile: $H$ at the minimum'
          'profilo esteso: H=0.5 sopra il min %',  '%.0f',  'Profile: excess at $H = 1/2$ (\%)'
          'profilo esteso: H=0.1 sopra il min %',  '%.0f',  'Profile: excess at $H = 0.10$ (\%)' };

fid = fopen(fullfile(root,'confronto_estensione.txt'),'w');
fprintf(fid, 'ESTENSIONE CON LE SCADENZE CORTE: CONFRONTO FRA LE DATE  (%s)\n', datestr(now));
fprintf(fid, '%s\n%-40s', repmat('=',1,40+14*nd), '');
for j = 1:nd, fprintf(fid, '%14s', C(j).data); end
fprintf(fid, '\n');
for r = 1:size(RIGHE,1)
    fprintf(fid, '%-40s', RIGHE{r,1});
    for j = 1:nd, fprintf(fid, '%14s', Fmt(RIGHE{r,2}, V(r,j))); end
    fprintf(fid, '\n');
end
for j = 1:nd
    fprintf(fid, '\n%s, profilo esteso (H, kappa, RMSE):', C(j).data);
    if isfield(C(j).S,'PRO'), fprintf(fid, '  %.2f/%.2f/%.2f', C(j).S.PRO(:,[1 4 5])'); end
end
if ~isempty(ESITO)
    fprintf(fid, '\n\nESITO\n');
    for i = 1:size(ESITO,1), fprintf(fid, '  %s  %8.0f s  %s\n', ESITO{i,:}); end
end
fclose(fid);
type(fullfile(root,'confronto_estensione.txt'));

if ~isfolder(fullfile(root,'Figure')), mkdir(fullfile(root,'Figure')); end
fid = fopen(fullfile(root,'Figure','tab_estensione_date.tex'),'w');
fprintf(fid, '%% Generato da Run_Estensione_Date.m il %s\n', datestr(now));
fprintf(fid, '\\begin{tabular}{l%s}\n\\toprule\n', repmat('r',1,nd));
for j = 1:nd, fprintf(fid, ' & %s', datestr(datenum(C(j).data),'dd mmm')); end
fprintf(fid, '\\\\\n\\midrule\n');
for r = 1:size(RIGHE,1)
    if any(r == [4 6 10 16]), fprintf(fid, '\\addlinespace\n'); end
    fprintf(fid, '%s', RIGHE{r,3});
    for j = 1:nd, fprintf(fid, ' & %s', FmtTex(RIGHE{r,2}, V(r,j))); end
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n'); fclose(fid);
fprintf('\nScritti confronto_estensione.txt e Figure/tab_estensione_date.tex\n');
