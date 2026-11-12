clearvars; clc;
% RIDISEGNA_SKEW  Ridisegna Figure/fig_skew_mkt.pdf (Fig. 7.5) da
% Results/test_fuori_campione.mat. Pochi secondi.
root = fileparts(mfilename('fullpath')); if isempty(root), root = pwd; end
addpath(fullfile(root,'Data')); resdir = fullfile(root,'Results'); figdir = fullfile(root,'Figure');
R = load(fullfile(resdir,'test_fuori_campione.mat'));
Tsk = R.Tsk(:); psiM = R.psiM(:); psi = R.psi; MOD = R.MOD;

% dati (stesso blocco di Test_Fuori_Campione)
raw = readtable('20260515SPX_merged.csv','Delimiter',',','VariableNamingRule','preserve');
is_spot = strcmp(raw.Root,'SPOT'); S0_feed = raw.CF_LAST(is_spot); raw = raw(~is_spot,:);
is_call = strcmp(raw.PutCall,'C');
kp = ~isnan(raw.STRIKE_PRC) & raw.CF_BID>0 & raw.CF_ASK>raw.CF_BID & raw.CF_BID>=0.375;
raw = raw(kp,:); is_call = is_call(kp);
K_all = raw.STRIKE_PRC; mid_all = 0.5*(raw.CF_BID+raw.CF_ASK);
T_all = (datenum(raw.EXPIR_DATE)-datenum('2026-05-15'))/365;
[~,ia] = unique([T_all,K_all,double(is_call)],'rows','stable');
T_all = T_all(ia); K_all = K_all(ia); mid_all = mid_all(ia); is_call = is_call(ia);
psL = nan(size(Tsk)); psH = nan(size(Tsk));
for i = 1:numel(Tsk)
    Ti = Tsk(i); mC = (T_all==Ti)&is_call; mP = (T_all==Ti)&~is_call;
    [Kc,iC,iP] = intersect(K_all(mC),K_all(mP)); cI = find(mC); pI = find(mP);
    d = mid_all(cI(iC)) - mid_all(pI(iP)); nr = abs(Kc/S0_feed-1) <= 0.25;
    if sum(nr) >= 5, Kr = Kc(nr); dr = d(nr); else, Kr = Kc; dr = d; end
    co = [ones(size(Kr)),Kr]\dr; Di = -co(2); Fi = co(1)/Di;
    m = (T_all==Ti); Ki = K_all(m); mi = mid_all(m); ci = is_call(m); Pf = nan(size(Ki));
    for j = 1:numel(Ki)
        if Ki(j) < Fi && ~ci(j),     Pf(j) = mi(j)/Di;
        elseif Ki(j) >= Fi && ci(j), Pf(j) = mi(j)/Di - Fi + Ki(j); end
    end
    o = ~isnan(Pf) & Pf > 1e-4; Ki = Ki(o); Pf = Pf(o);
    iv = nan(size(Ki));
    for j = 1:numel(Ki)
        try
            if Ki(j) <= Fi, iv(j) = blkimpv(Fi,Ki(j),0,Ti,Pf(j),'Class','put');
            else,           iv(j) = blkimpv(Fi,Ki(j),0,Ti,Pf(j)+Fi-Ki(j),'Class','call'); end
        catch, end
    end
    k = log(Ki/Fi); g = ~isnan(iv) & iv>0.01 & iv<2 & abs(k)<0.5; k = k(g); iv = iv(g);
    [~,a] = min(abs(k)); kw = max(0.75*iv(a)*sqrt(Ti), 0.006);
    psL(i) = SkewMercato(k, iv, 0.6*kw); psH(i) = SkewMercato(k, iv, 1.4*kw);
end

% figura
CC = { [0.00 0.45 0.74], [0.30 0.60 0.85], [0.10 0.60 0.55], [0.85 0.33 0.10], ...
       [0.64 0.08 0.18], [0.47 0.67 0.19], [0.49 0.18 0.56], [0.75 0.40 0.75], ...
       [0.93 0.69 0.13] };
STY = {'-','--','-','-','--','-.','-','--','-'};
okM = ~isnan(psiM) & psiM < 0;
f = figure('Color','w','Units','centimeters','Position',[2 2 16 7.2]);
tit = {'Without mean reversion','With mean reversion'};
yl = [0.13 3.6]; xl = [0.007 1.75];
for g = 0:1
    ax = subplot(1,2,g+1); hold on; box on; grid on
    set(ax,'XScale','log','YScale','log','FontSize',8, ...
        'XTick',[0.01 0.02 0.05 0.1 0.2 0.5 1],'XTickLabel',{'0.01','0.02','0.05','0.1','0.2','0.5','1'}, ...
        'YTick',[0.15 0.2 0.3 0.5 1 2 3],'YTickLabel',{'0.15','0.2','0.3','0.5','1','2','3'});
    patch([xl(1) 0.0384 0.0384 xl(1)], [yl(1) yl(1) yl(2) yl(2)], [0.92 0.92 0.92], ...
        'EdgeColor','none','HandleVisibility','off');
    a = abs(psiM); lo = min(abs(psL),abs(psH)); hi = max(abs(psL),abs(psH));
    e1 = max(a - lo, 0); e2 = max(hi - a, 0); e1(isnan(e1)) = 0; e2(isnan(e2)) = 0;
    errorbar(Tsk(okM), a(okM), e1(okM), e2(okM), 'ko', 'MarkerSize',3, ...
        'MarkerFaceColor','k', 'LineStyle','none', 'CapSize',2, 'DisplayName','Market');
    for j = find([MOD.gruppo] == g)
        plot(Tsk, abs(psi(:,j)), STY{j}, 'Color', CC{j}, 'LineWidth', 1.3, 'DisplayName', MOD(j).nome);
    end
    xlim(xl); ylim(yl);
    xlabel('$T$ (years)','Interpreter','latex','FontSize',9);
    if g == 0
        ylabel('$|\partial\sigma_{\mathrm{BS}}/\partial k|$ at the money','Interpreter','latex','FontSize',9);
    end
    title(tit{g+1},'FontSize',9,'FontWeight','normal');
    legend('Location','southwest','Interpreter','tex','Box','off','FontSize',7);
end
try
    exportgraphics(f, fullfile(figdir,'fig_skew_mkt.pdf'),'ContentType','vector');
catch
    set(f,'PaperPositionMode','auto');
    print(f, fullfile(figdir,'fig_skew_mkt.pdf'), '-dpdf', '-painters');
end
close(f); fprintf('Scritto Figure/fig_skew_mkt.pdf\n');

function s = SkewMercato(k, iv, kw)
    m = abs(k) <= kw;
    if sum(m) < 5 || sum(k(m) < 0) < 2 || sum(k(m) > 0) < 2, s = NaN; return; end
    p = polyfit(k(m), iv(m), 2); s = p(2);
end
