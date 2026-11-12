function Confronto_Date(cartelle)
% CONFRONTO_DATE  Tabella di confronto fra le date (Sez. 6.7).
%   Confronto_Date                      % tutte le cartelle Results_*
%   Confronto_Date({'Results_20260515','Results_20261014'})
% Legge riepilogo_AAAAMMGG.mat di ogni cartella.
% Scrive confronto_date.txt e Figure/tab_date.tex.
root = fileparts(mfilename('fullpath')); if isempty(root), root = pwd; end
if nargin < 1
    d = dir(fullfile(root, 'Results_*'));
    cartelle = sort({d([d.isdir]).name});
end
RR = {};
for i = 1:numel(cartelle)
    tag = strrep(cartelle{i}, 'Results_', '');
    f = fullfile(root, cartelle{i}, ['riepilogo_' tag '.mat']);
    if isfile(f), S = load(f); RR{end+1} = S.R;  %#ok<AGROW>
    else, fprintf('  [salto %s: manca riepilogo_%s.mat]\n', cartelle{i}, tag); end
end
assert(~isempty(RR), 'Nessun riepilogo trovato: lancia prima Run_Data.');
nd = numel(RR);

% righe della tabella: {etichetta testo, etichetta LaTeX, funzione, formato}
rm = @(R, k) Campo(R.(k), 'RMSE_paper');
righe = {
 'spot',                         'Spot',                                         @(R) R.mk.S0,               '%.0f'
 'sqrt(Kvar) prima scadenza',    '$\sqrt{K_{\mathrm{var}}}$, first maturity (\%)', @(R) Primo(R.mk.vs),      '%.1f'
 'sqrt(Kvar) a circa 1 anno',    '$\sqrt{K_{\mathrm{var}}}$, one year (\%)',     @(R) UnAnno(R.mk),          '%.1f'
 'RMSE Heston',                  'Heston',                                       @(R) rm(R,'hes'),           '%.2f'
 'RMSE Heston-curve',            'Heston-curve',                                 @(R) rm(R,'hcur'),          '%.2f'
 'RMSE Heston-curve-kappa',      'Heston-curve-$\kappa$',                        @(R) rm(R,'hck'),           '%.2f'
 'RMSE rough Heston',            'rough Heston',                                 @(R) rm(R,'rh'),            '%.2f'
 'RMSE rough Heston-kappa',      'rough Heston-$\kappa$',                        @(R) rm(R,'rhk'),           '%.2f'
 'RMSE CGMY',                    'CGMY',                                         @(R) rm(R,'cgmy'),          '%.2f'
 'RMSE composto',                'rough Heston + CGMY',                          @(R) rm(R,'rhc'),           '%.2f'
 'RMSE composto-kappa',          'rough Heston-$\kappa$ + CGMY',                 @(R) rm(R,'rhck'),          '%.2f'
 'RMSE rough Hawkes Heston',     'rough Hawkes Heston',                          @(R) rm(R,'rhh'),           '%.2f'
 'H rough Heston',               '$H$, rough Heston',                            @(R) R.D.H_rh,              '%.3f'
 'H composto',                   '$H$, rough Heston + CGMY',                     @(R) R.D.H_rhc,             '%.3f'
 'H rough Heston-kappa',         '$H$, rough Heston-$\kappa$',                   @(R) R.D.H_rhk,             '%.3f'
 'kappa rough Heston-kappa',     '$\kappa$, rough Heston-$\kappa$',              @(R) R.D.k_rhk,             '%.2f'
 'kappa Heston-curve-kappa',     '$\kappa$, Heston-curve-$\kappa$',              @(R) R.D.k_hck,             '%.2f'
 'H composto-kappa',             '$H$, rough Heston-$\kappa$ + CGMY',            @(R) R.D.H_rhck,            '%.3f'
 'H_eff rough Hawkes Heston',    '$H_{\mathrm{eff}}$, rough Hawkes Heston',      @(R) R.D.Heff,              '%.3f'
 'riduzione da ruvidezza k=0 %', 'Reduction from roughness, $\kappa = 0$ (\%)',  @(R) R.D.rid_rough,         '%.0f'
 'riduzione da mean reversion %','Reduction from mean reversion, $H = 1/2$ (\%)',@(R) R.D.rid_mr,            '%.0f'
 'riduzione da ruvidezza k lib %','Reduction from roughness, $\kappa$ free (\%)',@(R) R.D.rid_rough_k,       '%.0f'
 'riduzione dai salti %',        'Reduction from independent jumps (\%)',        @(R) R.D.rid_jump,          '%.0f'
 'spostamento di H dai salti',   'Change in $H$ from independent jumps',         @(R) R.D.dH_jump,           '%+.3f'
 'quota QV salti composto %',    'Jump share of quadratic variation (\%)',       @(R) R.D.qv_rhc,            '%.1f'
 };
nr = size(righe, 1);
V = nan(nr, nd);
for j = 1:nd
    for i = 1:nr
        try, V(i,j) = righe{i,3}(RR{j}); catch, V(i,j) = NaN; end
    end
end

% testo
out = fullfile(root, 'confronto_date.txt');
fid = fopen(out, 'w');
fprintf(fid, 'CONFRONTO FRA LE DATE  (%s)\n%s\n', datestr(now), repmat('=', 1, 80));
fprintf(fid, '%-32s', ''); for j = 1:nd, fprintf(fid, '%14s', RR{j}.data); end; fprintf(fid, '\n');
for i = 1:nr
    fprintf(fid, '%-32s', righe{i,1});
    for j = 1:nd, fprintf(fid, '%14s', sprintf(righe{i,4}, V(i,j))); end
    fprintf(fid, '\n');
end
fprintf(fid, '\nCLASSIFICA PER ERRORE RELATIVO\n');
for j = 1:nd
    fprintf(fid, '  %s: %s\n', RR{j}.data, strjoin(RR{j}.classifica, ' < '));
end
fclose(fid);
type(out);

% LaTeX
figdir = fullfile(root, 'Figure'); if ~isfolder(figdir), mkdir(figdir); end
ftex = fullfile(figdir, 'tab_date.tex');
fid = fopen(ftex, 'w');
fprintf(fid, '%% Generato da Confronto_Date.m il %s\n', datestr(now));
fprintf(fid, '\\begin{tabular}{l%s}\n\\toprule\n', repmat('r', 1, nd));
fprintf(fid, ' '); for j = 1:nd, fprintf(fid, ' & %s', datestr(datenum(RR{j}.data), 'dd mmm yyyy')); end
fprintf(fid, '\\\\\n\\midrule\n');
blocchi = [3 12 19 25];  % \midrule dopo queste righe
for i = 1:nr
    fprintf(fid, '%s', righe{i,2});
    for j = 1:nd
        if isnan(V(i,j)), s = '---';
        else, s = strrep(sprintf(righe{i,4}, V(i,j)), '-', '$-$'); end
        fprintf(fid, ' & %s', s);
    end
    fprintf(fid, '\\\\\n');
    if any(i == blocchi(1:end-1)), fprintf(fid, '\\midrule\n'); end
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
fprintf('\nTabella LaTeX: %s\n', ftex);
end

function v = Campo(S, nome)
    if isempty(S) || ~isfield(S, nome), v = NaN; else, v = S.(nome); end
end
function v = Primo(x)
    if isempty(x), v = NaN; else, v = x(1); end
end
function v = UnAnno(mk)
    if isempty(mk.T), v = NaN; return; end
    [~, i] = min(abs(mk.T - 1)); v = mk.vs(i);
end
