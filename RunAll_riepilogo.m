function RunAll_riepilogo()
% Stampa il riepilogo finale leggendo lo stato persistent.
S = RunAll_stato('get');
fprintf('\n\n==================================================================\n');
fprintf('RIEPILOGO RUN ALL  (totale %.1f min)\n', (now - S.t0)*1440);
fprintf('==================================================================\n');
for k = 1:size(S.scripts,1)
    fprintf('  %-34s %7.1f s   %s\n', S.scripts{k,1}, S.tempi(k), S.esito(k));
end
fprintf('  %-34s %7.1f s   %s\n', 'PricingSuite', S.tempi(end), S.esito(end));
fprintf('==================================================================\n');
if all(S.esito == "OK")
    fprintf('Tutto rigenerato correttamente. I numeri della tesi sono coerenti.\n');
else
    fprintf(2, 'Attenzione: alcuni step non OK, controlla il riepilogo sopra.\n');
end
end
