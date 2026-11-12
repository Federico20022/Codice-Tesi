function R = Riepilogo_Data(resdir, data_str, csv, ESITO)
% RIEPILOGO_DATA  Sintesi di una data, nell'ordine del Cap. 6.
%   R = Riepilogo_Data('Results_20261014', '2026-10-14')
% Legge i riass_*.mat di Run_Data e scrive riepilogo_AAAAMMGG.txt.
if nargin < 3, csv = ''; end
if nargin < 4, ESITO = {}; end
tag = strrep(data_str, '-', '');

M = { 'Heston_SPX_20260515_def',        'hes',  'Heston',                     5; ...
      'Check_Heston_Curve',             'hcur', 'Heston-curve',               2; ...
      'Heston_Curve_Kappa',             'hck',  'Heston-curve-kappa',         3; ...
      'rough_heston_def',               'rh',   'rough Heston',               3; ...
      'Rough_Heston_Kappa',             'rhk',  'rough Heston-kappa',         4; ...
      'CMGY_def',                       'cgmy', 'CGMY',                       3; ...
      'Rough_Heston_CMGY_def',          'rhc',  'rough Heston + CGMY',        6; ...
      'Rough_Heston_CMGY_Kappa',        'rhck', 'rough Heston-kappa + CGMY',  7; ...
      'RoughHawkesHeston_SPX_20260515', 'rhh',  'rough Hawkes Heston',        7 };

R = struct('data', data_str, 'tag', tag, 'csv', csv);
for i = 1:size(M,1)
    f = fullfile(resdir, ['riass_' M{i,1} '.mat']);
    if isfile(f)
        S = load(f); R.(M{i,2}) = S.RIASS;
    else
        R.(M{i,2}) = [];
    end
end

rm  = @(k) Campo(R.(k), 'RMSE_paper');
ra  = @(k) Campo(R.(k), 'RMSE_Vol_abs');
th  = @(k, j) Elem(R.(k), j);
rid = @(a, b) 100*(1 - b/a);

% quantita' derivate (Tab. 6.3 e 6.7)
D = struct();
D.H_rh   = th('rh', 1);    D.H_rhk = th('rhk', 1);  D.k_rhk = th('rhk', 4);
D.k_hck  = th('hck', 1);   D.H_rhc = th('rhc', 1);  D.H_rhck = th('rhck', 1);
D.k_rhck = th('rhck', 8);  D.Heff  = th('rhh', 1) - 0.5;  D.b_rhh = th('rhh', 2);
D.rid_rough   = rid(rm('hcur'), rm('rh'));  % ruvidezza, kappa = 0
D.rid_mr      = rid(rm('hcur'), rm('hck'));  % mean reversion, H = 1/2
D.rid_rough_k = rid(rm('hck'),  rm('rhk'));  % ruvidezza, kappa libero
D.rid_mr_r    = rid(rm('rh'),   rm('rhk'));  % mean reversion sopra la ruvidezza
D.rid_jump    = rid(rm('rh'),   rm('rhc'));  % salti, kappa = 0
D.rid_jump_k  = rid(rm('rhk'),  rm('rhck'));  % salti, kappa libero
D.dH_jump     = D.H_rhc  - D.H_rh;
D.dH_jump_k   = D.H_rhck - D.H_rhk;
[D.qv_rhc,  D.vsl_rhc ] = QuotaSalti(R.rhc,  1.0);
[D.qv_rhck, D.vsl_rhck] = QuotaSalti(R.rhck, 1.0);
if ~isempty(R.rhh), L = R.rhh.theta(5); D.qv_rhh = 100*2*L^2/(1+2*L^2); else, D.qv_rhh = NaN; end
R.D = D;

% classifica per errore relativo
err = cellfun(@(k) rm(k), M(:,2));
[~, ord] = sort(err); ord = ord(~isnan(err(ord)));
R.classifica = M(ord, 3);

% descrizione del mercato
mk = struct('T', [], 'vs', [], 'S0', NaN);
if ~isempty(R.rh) && isfield(R.rh, 'T_uniq') && isfield(R.rh, 'V_fair')
    mk.T = R.rh.T_uniq(:).'; mk.vs = 100*sqrt(R.rh.V_fair(:).');
end
if ~isempty(R.rh) && isfield(R.rh, 'S0_feed'), mk.S0 = R.rh.S0_feed; end
R.mk = mk;

% fuori campione
R.ps = [];
f = fullfile(resdir, 'PricingSuite_kappa.mat');
if isfile(f), R.ps = load(f); end
R.tfc = [];
f = fullfile(resdir, 'test_fuori_campione.mat');
if isfile(f), R.tfc = load(f); end

% SCRITTURA
out = fullfile(resdir, ['riepilogo_' tag '.txt']);
fid = fopen(out, 'w');
P = @(varargin) fprintf(fid, varargin{:});
P('RIEPILOGO %s   (file %s, scritto %s)\n', data_str, csv, datestr(now));
P('%s\n\n', repmat('=', 1, 72));

P('MERCATO\n');
P('  spot (feed)         %.2f\n', mk.S0);
if ~isempty(mk.T)
    P('  scadenze (giorni)  %s\n', sprintf('%7.0f', 365*mk.T));
    P('  sqrt(Kvar) (%%)     %s\n', sprintf('%7.2f', mk.vs));
    P('  pendenza sqrt(Kvar) ultima - prima: %+.2f punti di vol\n', mk.vs(end) - mk.vs(1));
end

P('\nERRORI DI CALIBRAZIONE (Tab. 6.2 e 6.3)\n');
P('  %-28s %6s %12s %12s   %s\n', 'modello', 'param', 'RMSE rel %', 'RMSE abs %', 'parametri');
for i = 1:size(M,1)
    k = M{i,2};
    if isempty(R.(k)), P('  %-28s %6d %12s\n', M{i,3}, M{i,4}, 'MANCANTE'); continue; end
    P('  %-28s %6d %12.3f %12.3f   %s\n', M{i,3}, M{i,4}, rm(k), ra(k), mat2str(R.(k).theta, 4));
end
P('  classifica: %s\n', strjoin(R.classifica, ' < '));

P('\nRUVIDEZZA, MEAN REVERSION E SALTI\n');
P('  H rough Heston (kappa = 0)              %.4f\n', D.H_rh);
P('  H, kappa rough Heston-kappa             %.4f, %.3f\n', D.H_rhk, D.k_rhk);
P('  kappa Heston-curve-kappa                %.3f\n', D.k_hck);
P('  H composto, H composto-kappa            %.4f, %.4f  (kappa %.3f)\n', D.H_rhc, D.H_rhck, D.k_rhck);
P('  H_eff rough Hawkes Heston, b            %.4f, %.3f\n', D.Heff, D.b_rhh);
P('  riduzione da ruvidezza, kappa = 0       %5.1f%%   (Heston-curve -> rough Heston)\n', D.rid_rough);
P('  riduzione da mean reversion, H = 1/2    %5.1f%%   (Heston-curve -> Heston-curve-kappa)\n', D.rid_mr);
P('  riduzione da ruvidezza, kappa libero    %5.1f%%   (Heston-curve-kappa -> rough Heston-kappa)\n', D.rid_rough_k);
P('  riduzione da kappa sopra la ruvidezza   %5.1f%%   (rough Heston -> rough Heston-kappa)\n', D.rid_mr_r);
P('  riduzione dai salti, kappa = 0          %5.1f%%   (H si sposta di %+.4f)\n', D.rid_jump, D.dH_jump);
P('  riduzione dai salti, kappa libero       %5.1f%%   (H si sposta di %+.4f)\n', D.rid_jump_k, D.dH_jump_k);
P('  quota QV dei salti a 1 anno: composto %.1f%%, composto-kappa %.1f%%, rHH %.1f%%\n', D.qv_rhc, D.qv_rhck, D.qv_rhh);
P('  VS/L della sola componente di salto: composto %.3f, composto-kappa %.3f\n', D.vsl_rhc, D.vsl_rhck);

if ~isempty(R.rhk) && isfield(R.rhk, 'PROF') && isfield(R.rhk, 'H_GRID')
    P('\nPROFILO IN H CON KAPPA LIBERO (Tab. 6.4, destra)\n');
    P('  %6s %8s %8s %8s %10s\n', 'H', 'nu', 'rho', 'kappa', 'RMSE %');
    for i = 1:numel(R.rhk.H_GRID)
        P('  %6.2f %8.4f %8.4f %8.4f %10.3f\n', R.rhk.H_GRID(i), R.rhk.PROF(i,:));
    end
end
if ~isempty(R.hck) && isfield(R.hck, 'PROF') && isfield(R.hck, 'KAPPA_GRID')
    P('\nPROFILO IN KAPPA A H = 1/2 (Tab. 6.4, sinistra)\n');
    P('  %6s %8s %8s %10s\n', 'kappa', 'nu', 'rho', 'RMSE %');
    for i = 1:numel(R.hck.KAPPA_GRID)
        P('  %6.1f %8.4f %8.4f %10.3f\n', R.hck.KAPPA_GRID(i), R.hck.PROF(i,:));
    end
end

if ~isempty(R.ps)
    S = R.ps; nm = S.MODELLI;
    P('\nFUORI CAMPIONE (PricingSuite_Kappa)\n');
    [it, ik] = deal(find(abs(S.T_call - 0.5) < 1e-9), find(abs(S.Kf - 1.10) < 1e-9));
    c = squeeze(S.CallPrices(it, ik, :));
    P('  call T = 0.5, K/F = 1.10: %s\n', Coppie(nm, c, '%.4f'));
    [it, ik] = deal(find(abs(S.T_dig - 0.25) < 1e-9), find(abs(S.Kf_d - 1.10) < 1e-9));
    d = squeeze(S.Digitals(it, ik, :));
    P('  digitale T = 0.25, K/F = 1.10: %s\n', Coppie(nm, d, '%.4f'));
    P('  pendenza skew su T <= 0.25: %s\n', Coppie(nm, S.Slope, '%+.3f'));
    g5 = 1:5; g0 = [find(strcmp(nm,'Heston')) find(strcmp(nm,'rHeston')) find(strcmp(nm,'CGMY')) ...
                    find(strcmp(nm,'rH+CGMY')) find(strcmp(nm,'rHawkesH'))];
    sp = @(v, g) 100*(max(v(g)) - min(v(g)))/min(v(g));
    P('  spread call (gruppo con kappa / cinque del Cap. 7): %.1f%% / %.1f%%\n', sp(c, g5), sp(c, g0));
    P('  spread digitale (gruppo con kappa / cinque del Cap. 7): %.1f%% / %.1f%%\n', sp(d, g5), sp(d, g0));
end

if ~isempty(R.tfc)
    P('\nTEST FUORI CAMPIONE (Test_Fuori_Campione)\n');
    T = R.tfc;
    if isfield(T, 'R_all') && isfield(T, 'MOD')
        P('  RMSE rel su tutta la superficie filtrata:\n');
        for j = 1:numel(T.MOD)
            P('    %-28s %7.2f%%\n', T.MOD(j).nome, T.R_all(j));
        end
    end
    if isfield(T, 'Tsk') && isfield(T, 'psiM')
        ok = ~isnan(T.psiM) & T.psiM < 0 & T.Tsk < 14/365;
        if sum(ok) >= 3
            cf = polyfit(log(T.Tsk(ok)), log(abs(T.psiM(ok))), 1);
            P('  pendenza log-log dello skew di mercato sotto le due settimane: %+.3f (%d scadenze)\n', cf(1), sum(ok));
        else
            P('  meno di tre scadenze sotto le due settimane: pendenza non stimata\n');
        end
    end
end

if ~isempty(ESITO)
    P('\nESITO DEGLI SCRIPT\n');
    for i = 1:size(ESITO,1), P('  %-34s %8.1f s   %s\n', ESITO{i,1}, ESITO{i,2}, ESITO{i,3}); end
end
fclose(fid);
save(fullfile(resdir, ['riepilogo_' tag '.mat']), 'R');
fprintf('Riepilogo scritto in %s\n', out);
type(out);
end

function v = Campo(S, nome)
    if isempty(S) || ~isfield(S, nome), v = NaN; else, v = S.(nome); end
end

function v = Elem(S, j)
    if isempty(S) || numel(S.theta) < j, v = NaN; else, v = S.theta(j); end
end

function [qv, vsl] = QuotaSalti(S, T)
% Quota della variazione quadratica dovuta ai salti a T e VS/L della sola
% componente di salto, per i modelli composti (Sez. 6.4).
    qv = NaN; vsl = NaN;
    if isempty(S) || ~isfield(S, 'Kvar_mkt') || ~isfield(S, 'T_uniq'), return; end
    t = S.theta; C = t(4); G = t(5); Mm = t(6); Y = t(7);
    Lcg = 2*C*( gamma(-Y)*((Mm-1)^Y - Mm^Y + (G+1)^Y - G^Y) - gamma(1-Y)*(Mm^(Y-1) - G^(Y-1)) );
    Vcg = C * gamma(2-Y) * (Mm^(Y-2) + G^(Y-2));
    W  = pchip([0; S.T_uniq(:)], [0; S.T_uniq(:).*S.Kvar_mkt(:)]);
    N  = 400; ti = ((0:N-1) + 0.5)/N * T; h = 1e-5;
    xi = max((ppval(W, ti + h) - ppval(W, max(ti - h, 0))) ./ (ti + h - max(ti - h, 0)), 1e-8);
    qv  = 100 * Vcg / (mean(max(xi - Lcg, 1e-6)) + Vcg);
    vsl = Vcg / Lcg;
end

function s = Coppie(nm, v, fmt)
    s = '';
    for j = 1:numel(nm), s = [s sprintf(['%s ' fmt '  '], nm{j}, v(j))]; end  %#ok<AGROW>
end
