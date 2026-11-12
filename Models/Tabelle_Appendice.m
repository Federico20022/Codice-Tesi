% TABELLE_APPENDICE  Scrive in Figure/ le tabelle dell'Appendice C da
% Results/valida_pricer.mat e Results/multistart.mat. Non calcola nulla.
% Si lancia da solo dopo Figure_Tesi, oppure da Figure_Tesi prima delle
% funzioni locali. Niente clearvars, per non svuotare il chiamante.

if ~exist('figdir','var') || ~exist('resdir','var')
    root_app = fileparts(mfilename('fullpath'));
    if isempty(root_app), root_app = pwd; end
    resdir = fullfile(root_app, '..', 'Results');
    figdir = fullfile(root_app, '..', 'Figure');
end
if ~exist(figdir,'dir'), mkdir(figdir); end

fprintf('\nTabelle dell''appendice C...\n');
Wv = load(fullfile(resdir,'valida_pricer.mat'));
VA = Wv.VALIDA;

% C.1.1  SINC contro l'inversione di Lewis 24 combinazioni strike/scadenza
% sotto Heston.
fid = ApriTabApp(figdir,'tab_app_lewis');
fprintf(fid,'\\begin{tabular}{rrrrr}\n\\toprule\n');
fprintf(fid,'$T$ & $K/F$ & Lewis & SINC & Rel. error\\\\\n\\midrule\n');
for i = 1:numel(VA.A.T)
    fprintf(fid,'%.4f & %.2f & %.8f & %.8f & %.2e\\\\\n', ...
        VA.A.T(i), VA.A.KF(i), VA.A.lewis(i), VA.A.sinc(i), VA.A.err_rel(i));
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{4}{l}{Maximum relative error} & %.2e\\\\\n', VA.A.err_rel_max);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.1.2  convergenza interna in N
% Non dipende da alcun riferimento esterno: se pow=14 e pow=16 coincidono,
% la configurazione accurata e' convergente.
fid = ApriTabApp(figdir,'tab_app_convergenza');
fprintf(fid,'\\begin{tabular}{rrrrrr}\n\\toprule\n');
fprintf(fid,'$T$ & $K/F$ & $N=2^{12}$ & $N=2^{14}$ & $N=2^{16}$ & $|2^{12}-2^{16}|$\\\\\n');
fprintf(fid,'\\midrule\n');
for i = 1:numel(VA.Abis.T)
    fprintf(fid,'%.4f & %.2f & %.8f & %.8f & %.8f & %.2e\\\\\n', ...
        VA.Abis.T(i), VA.Abis.KF(i), VA.Abis.pow12(i), VA.Abis.pow14(i), ...
        VA.Abis.pow16(i), abs(VA.Abis.pow12(i)-VA.Abis.pow16(i)));
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{5}{l}{Maximum difference} & %.2e\\\\\n', VA.Abis.max_diff_12_16);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.1.3  errore della configurazione coarse
% Una colonna con "n.a." segnala un prezzo non positivo: quella
% configurazione non e' utilizzabile in una loss pesata.
fid = ApriTabApp(figdir,'tab_app_coarse');
nc = size(VA.B.cfg,1);
fprintf(fid,'\\begin{tabular}{rr%s}\n\\toprule\n',repmat('r',1,nc));
fprintf(fid,'$T$ & $K/F$');
for j = 1:nc
    fprintf(fid,' & $X_c=%g$, $2^{%g}$', VA.B.cfg(j,1), VA.B.cfg(j,2));
end
fprintf(fid,'\\\\\n\\midrule\n');
for i = 1:numel(VA.B.T)
    fprintf(fid,'%.4f & %.2f', VA.B.T(i), VA.B.KF(i));
    for j = 1:nc
        if isnan(VA.B.dIV_bp(i,j)), fprintf(fid,' & n.a.');
        else,                       fprintf(fid,' & %.1f', VA.B.dIV_bp(i,j)); end
    end
    fprintf(fid,'\\\\\n');
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.1.4  Pade al bordo H = 1/2
% Questo NON e' un test di convergenza di griglia: il Pade converge al
% proprio valore, che ad alpha = 1 non e' la soluzione esatta.
fid = ApriTabApp(figdir,'tab_app_pade');
fprintf(fid,'\\begin{tabular}{rrrrr}\n\\toprule\n');
fprintf(fid,'$T$ & $K/F$ & Closed-form Heston & Pad\\''e $(3,3)$ & $\\Delta$IV (bp)\\\\\n');
fprintf(fid,'\\midrule\n');
for i = 1:numel(VA.C.T)
    fprintf(fid,'%.4f & %.2f & %.6f & %.6f & %+.1f\\\\\n', ...
        VA.C.T(i), VA.C.KF(i), VA.C.heston_chiuso(i), VA.C.pade(i), VA.C.dIV_bp(i));
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{4}{l}{Maximum absolute discrepancy} & %.1f\\\\\n', VA.C.max_abs_bp);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.2.1  martingala e put-call parity
% Per Heston la martingala e' esatta; per il ramo Pade il residuo e' un
% limite e scala linearmente con l'offset usato per valutare phi(-i).
fid = ApriTabApp(figdir,'tab_app_martingala');
fprintf(fid,'\\begin{tabular}{lrrr}\n\\toprule\n');
fprintf(fid,'Model & $T$ & $|\\varphi_T(-i)-1|$ & $\\max|C-P-(F-K)|$\\\\\n\\midrule\n');
nomiE = {'Heston','Heston','rough Heston','rough Heston'};
for i = 1:numel(VA.E.T)
    fprintf(fid,'%s & %.4f & %.2e & %.2e\\\\\n', ...
        nomiE{i}, VA.E.T(i), VA.E.martingala(i), VA.E.parity(i));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.2.2  convenzione della griglia di Riccati
% Con i punti medi l'errore deve scalare come 2^(1+alpha): un rapporto
% vicino a 2 segnalerebbe griglie disallineate.
fid = ApriTabApp(figdir,'tab_app_griglia');
fprintf(fid,'\\begin{tabular}{rrrr}\n\\toprule\n');
fprintf(fid,'$T$ & $N$ & ATM price & $|p(N)-p(2N)|$\\\\\n\\midrule\n');
for r = 1:size(VA.G2.prezzi,1)
    for j = 1:numel(VA.G2.N)
        if j==1, et=sprintf('%.4f',VA.G2.T(r)); else, et=''; end
        if j < numel(VA.G2.N)
            fprintf(fid,'%s & %d & %.8f & %.3e\\\\\n', et, VA.G2.N(j), ...
                VA.G2.prezzi(r,j), abs(VA.G2.prezzi(r,j)-VA.G2.prezzi(r,j+1)));
        else
            fprintf(fid,'%s & %d & %.8f & ---\\\\\n', et, VA.G2.N(j), VA.G2.prezzi(r,j));
        end
    end
    fprintf(fid,'\\multicolumn{3}{l}{measured ratio} & %.2f\\\\\n', VA.G2.rapporto_misurato(r));
    if r < size(VA.G2.prezzi,1), fprintf(fid,'\\addlinespace\n'); end
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{3}{l}{expected ratio $2^{1+\\alpha}$, $\\alpha=%.4f$} & %.2f\\\\\n', ...
    VA.G2.alpha, VA.G2.rapporto_atteso);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.2.3  sensibilita' del caso nested a H = 1/2 - eps
fid = ApriTabApp(figdir,'tab_app_eps');
fprintf(fid,'\\begin{tabular}{r%s}\n\\toprule\n$\\varepsilon$',repmat('r',1,numel(VA.F.T)));
for j = 1:numel(VA.F.T), fprintf(fid,' & ATM $T=%.3f$', VA.F.T(j)); end
fprintf(fid,'\\\\\n\\midrule\n');
for i = 1:numel(VA.F.eps)
    fprintf(fid,'%.0e', VA.F.eps(i));
    for j = 1:numel(VA.F.T), fprintf(fid,' & %.6f', VA.F.prezzi(i,j)); end
    fprintf(fid,'\\\\\n');
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.3.1  variance swap estratte
fid = ApriTabApp(figdir,'tab_app_kvar');
fprintf(fid,'\\begin{tabular}{rrr}\n\\toprule\n');
fprintf(fid,'$T$ (yr) & $\\mathcal{K}_{\\mathrm{var}}(T)$ & $\\sqrt{\\mathcal{K}_{\\mathrm{var}}}$ (\\%%)\\\\\n');
fprintf(fid,'\\midrule\n');
for i = 1:numel(VA.T_uniq)
    fprintf(fid,'%.4f & %.5f & %.3f\\\\\n', VA.T_uniq(i), VA.V_fair(i), ...
        100*sqrt(VA.V_fair(i)));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.3.2  curva PCHIP contro livelli flat
% Stessa varianza totale a ogni scadenza, forma diversa: la tabella misura
% quanto vale la differenza in punti base di IV.
fid = ApriTabApp(figdir,'tab_app_pchipflat');
fprintf(fid,'\\begin{tabular}{rrrrr}\n\\toprule\n');
fprintf(fid,'$T$ & $K/F$ & IV pchip (\\%%) & IV flat (\\%%) & $\\Delta$IV (bp)\\\\\n');
fprintf(fid,'\\midrule\n');
for i = 1:numel(VA.D.T)
    fprintf(fid,'%.4f & %.2f & %.3f & %.3f & %+.1f\\\\\n', ...
        VA.D.T(i), VA.D.KF(i), 100*VA.D.iv_pchip(i), 100*VA.D.iv_flat(i), VA.D.dIV_bp(i));
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{4}{l}{Maximum absolute discrepancy} & %.1f\\\\\n', VA.D.max_abs_bp);
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

% C.4  multistart
% Una riga per modello: quante partenze, quanti bacini, accordo
% sull'obiettivo e escursione dei parametri entro il bacino.
if isfile(fullfile(resdir,'multistart.mat'))
    Ms = load(fullfile(resdir,'multistart.mat'));
    mm  = {'rheston','rhcgmy','cgmy'};
    eti = {'rough Heston','rough Heston + CGMY','CGMY'};
    fid = ApriTabApp(figdir,'tab_app_multistart');
    fprintf(fid,'\\begin{tabular}{llrrl}\n\\toprule\n');
    fprintf(fid,'Model & Parameter & Restarts & Basins & Range across restarts\\\\\n');
    fprintf(fid,'\\midrule\n');
    for k = 1:numel(mm)
        S = Ms.OUT.(mm{k});
        for j = 1:numel(S.nomi)
            if j==1
                et1 = eti{k};
                et2 = sprintf('%d', numel(S.loss));
                et3 = sprintf('%d', S.n_bacini);
            else
                et1=''; et2=''; et3='';
            end
            fprintf(fid,'%s & $%s$ & %s & %s & $[%.6g,\\; %.6g]$\\\\\n', ...
                et1, S.nomi{j}, et2, et3, min(S.theta(:,j)), max(S.theta(:,j)));
        end
        fprintf(fid,'\\addlinespace\n');
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);

    fprintf('  tab_app_multistart\n');
end

fprintf('Tabelle dell''appendice scritte in %s\n', figdir);

function fid = ApriTabApp(cartella, nome)
    fid = fopen(fullfile(cartella,[nome '.tex']),'w');
    fprintf(fid,'%% Generato da Tabelle_Appendice.m il %s -- non modificare a mano\n', ...
        datestr(now));
    fprintf('  %s\n', nome);
end
