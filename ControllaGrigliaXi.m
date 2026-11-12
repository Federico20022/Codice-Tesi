function ok = ControllaGrigliaXi(root)
% CONTROLLAGRIGLIAXI  Verifica che ogni file costruisca la griglia di xi0
% sui punti medi ti = ((0:n-1)+0.5)/n*T.
%   ok = ControllaGrigliaXi          % cartella di questo file
%   ok = ControllaGrigliaXi(root)
% Con gli estremi sinistri flip(xi) sbaglia di mezzo passo (fino a 10 bp
% di IV). Il controllo legge i sorgenti, quindi e' deterministico.

if nargin < 1 || isempty(root)
    root = fileparts(mfilename('fullpath'));
end

files = { 'Functions/charfun.m'
          'Models/rough_heston_def.m'
          'Models/Rough_Heston_CMGY_def.m'
          'Models/Check_Heston_Curve.m'
          'PricingSuite_SPX_20260515.m' };

verdetto = {'  <-- MANCA IL +0.5', '  ok'};
ok       = true;
trovate  = 0;

fprintf('\n=== Controllo convenzione griglia xi (punti medi) ===\n');
for k = 1:numel(files)
    f = fullfile(root, files{k});
    if ~isfile(f)
        fprintf('  [assente] %s\n', files{k});
        continue
    end
    testo = strsplit(fileread(f), newline);
    n_qui = 0;
    for i = 1:numel(testo)
        L = testo{i};
        if isempty(regexp(L, '^\s*ti\s*=', 'once')), continue; end
        n_qui   = n_qui + 1;
        trovate = trovate + 1;
        buona   = ~isempty(strfind(L, '0.5'));  %#ok<STREMP>
        ok      = ok && buona;
        fprintf('  %-34s r.%-5d %s%s\n', files{k}, i, strtrim(L), verdetto{buona+1});
    end
    if n_qui == 0
        fprintf('  %-34s nessuna assegnazione di ti trovata\n', files{k});
    end
end

fprintf('\n');
if trovate == 0
    fprintf(2, 'Nessuna griglia trovata: controllare il percorso root.\n');
    ok = false;
elseif ok
    fprintf('Tutte le %d griglie usano i punti medi: convenzione COERENTE.\n', trovate);
else
    fprintf(2, ['Almeno una griglia usa ancora gli estremi sinistri.\n' ...
                'Sostituire  ti = (0:N-1)/N*T  con  ti = ((0:N-1)+0.5)/N*T.\n']);
end
fprintf('\n');
end
