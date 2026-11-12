clearvars; clc;
% FIGURA_ESTENSIONE  Ridisegna Fig. 7.6 e Tab. 7.7 da Results/estensione.mat.
% Nessun ricalcolo, pochi secondi.
root = fileparts(mfilename('fullpath')); if isempty(root), root = pwd; end
resdir = fullfile(root,'Results'); figdir = fullfile(root,'Figure');
S = load(fullfile(resdir,'estensione.mat'));
R = S.R; PRO = S.PRO; MOD = S.MOD;

% figura: profilo in H con kappa libero
H6 = [0.05 0.10 0.20 0.30 0.40 0.50]; R6 = [4.034 3.859 3.553 3.334 3.240 3.294];  % Tab. 6.4
ok = ~isnan(PRO(:,1));
f = figure('Color','w','Units','centimeters','Position',[2 2 11 7]);
hold on; grid on; box on
plot(H6, R6, 'o--', 'Color',[0.45 0.45 0.45], 'LineWidth',1.2, 'MarkerSize',5, ...
    'DisplayName','Seven calibration maturities');
plot(PRO(ok,1), PRO(ok,5), 's-', 'Color',[0.64 0.08 0.18], 'LineWidth',1.5, 'MarkerSize',5, ...
    'MarkerFaceColor',[0.64 0.08 0.18], 'DisplayName','With the three short maturities');
Hst = R.rhk.theta(1); Rst = R.rhk.M1.tot;  % ottimo libero di rough Heston-kappa (Tab. 7.7)
plot(Hst, Rst, 'p', 'Color','k', 'MarkerFaceColor','k', 'MarkerSize',9, 'HandleVisibility','off');
xlim([0 0.52]); ylim([3.1 4.8]);
set(gca,'FontSize',8,'XTick',0:0.1:0.5);
xlabel('$H$','Interpreter','latex','FontSize',9);
ylabel('Relative RMSE (\%)','Interpreter','latex','FontSize',9);
legend('Location','north','Box','off','FontSize',8);
try, exportgraphics(f, fullfile(figdir,'fig_estensione.pdf'),'ContentType','vector');
catch, set(f,'PaperPositionMode','auto'); print(f, fullfile(figdir,'fig_estensione.pdf'),'-dpdf','-painters'); end
close(f);

% tabella
NOME = struct('hck','Heston-curve-$\kappa$','rh','rough Heston','rhk','rough Heston-$\kappa$', ...
              'rhh','rough Hawkes Heston');
fid = fopen(fullfile(figdir,'tab_estensione.tex'),'w');
fprintf(fid,'%% Generato da Figura_Estensione.m il %s\n', datestr(now));
fprintf(fid,'\\begin{tabular}{lrrrrrrr}\n\\toprule\n');
fprintf(fid,' & & & & & \\multicolumn{3}{c}{RMSE (\\%%)}\\\\\n\\cmidrule(lr){6-8}\n');
fprintf(fid,'Model & $H$ & $\\kappa$ & $\\nu$ & $\\rho$ & All & Short & Standard\\\\\n\\midrule\n');
for im = 1:numel(MOD)
    nm = MOD{im}; th = R.(nm).theta; M = R.(nm).M1;
    switch nm
        case 'hck', v = [0.5,       th(1),  th(2),       th(3)];
        case 'rh',  v = [th(1),     0,      th(2),       th(3)];
        case 'rhk', v = [th(1),     th(4),  th(2),       th(3)];
        case 'rhh', v = [th(1)-0.5, -th(2), sqrt(th(3)), th(4)];
    end
    riga = sprintf('%s & %.3f & %.2f & %.3f & %.3f & %.2f & %.2f & %.2f', NOME.(nm), v, M.tot, M.cor, M.std);
    riga = strrep(riga, '& -', '& $-$');  % segno meno tipografico
    fprintf(fid, '%s\\\\\n', riga);
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);
fprintf('Scritti Figure/fig_estensione.pdf e Figure/tab_estensione.tex\n');
