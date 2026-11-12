# SINC Fourier Pricing for Rough Volatility and Jump Models

MATLAB code for the master's thesis *SINC Fourier Pricing for Rough Volatility and Jump Models: a Comparative Calibration on the SPX Surface*, by Federico Fabbri and Federico Savini (Politecnico di Milano, advisor Prof. Daniele Marazzina).

The code calibrates nine models to the SPX implied volatility surface with the SINC Fourier pricer. The same sample, weights and pricer are used for every model. The calibrated models are then compared on instruments outside the calibration sample.

| Family | Models |
|---|---|
| Diffusive | Heston, Heston-curve, Heston-curve-κ |
| Rough | rough Heston (Padé 3,3), rough Heston-κ (multi-factor) |
| Jumps | CGMY, rough Heston + CGMY, rough Heston-κ + CGMY |
| Rough with jumps | rough Hawkes Heston |

## Data

The option quotes come from LSEG Workspace (Refinitiv). **They are under licence and are not included in this repository.** The same applies to the `.mat` files, which contain market quotes and calibrated surfaces.

To run the code, put one CSV per date in `Data/`, named `AAAAMMGGSPX_merged.csv`, with these columns:

```
Instrument,Root,PutCall,CF_LAST,CF_BID,CF_ASK,STRIKE_PRC,EXPIR_DATE
```

One row must have `Root = SPOT` and the index level in `CF_LAST`. The thesis uses 27 March, 15 May, 27 July and 14 October 2026.

## Requirements

- MATLAB R2024b
- Optimization Toolbox (`lsqnonlin`)
- Global Optimization Toolbox (`ga`)
- Financial Toolbox (`blkprice`, `blkimpv`, `busdays`)

The validation scripts `Valida_RHH` and `Valida_RHKappa` need no toolbox.

## Structure

| Folder / file | Content |
|---|---|
| `Functions/` | Pricer and characteristic functions: `SINC_fastFT`, `charfun`, `CharFunRHH`, `CharFunRHKappa`, Padé coefficients, `XcPow`, `SelezionaXc` |
| `Models/` | One calibration script per model, plus diagnostics (`Multistart`, `Profile_G`, `Diagnostica_Xc`, `Robustezza_H`, `Sensibilita_SVI`) and `Figure_Tesi` |
| Root scripts | Pipelines (`Run_All`, `Run_Data`), pricing suites, out-of-sample tests, short-maturity extension, validations |
| `Figure/` | Figures (`fig_*.pdf`) and LaTeX tables (`tab_*.tex`) of the thesis, 15 May 2026 |
| `Figure_AAAAMMGG/` | Figures and tables for each of the four dates |
| `Results_AAAAMMGG/` | Text summaries for each date (`riepilogo_*.txt`, pricing suite and extension logs) |
| `confronto_*.txt` | Comparisons across dates, short-maturity extension and trading-day clock |

Comments in the code are in Italian. Text inside figures and tables is in English.

## Reproducing the results

Run everything from the repository root.

**15 May 2026 (Chapters 6 and 7)**

1. `Run_All`: nine calibrations, diagnostics, validations and pricing suite (about 4 hours).
2. `Figure_Tesi`: figures and tables of Chapters 6, 7 and Appendix C.
3. `PricingSuite_Kappa`: pricing of the κ variants.
4. `Test_Fuori_Campione`: market skew below two weeks and error on the full surface.
5. `Estensione`, then `Figura_Estensione`: calibration with the short maturities.

**Other dates**

1. `Run_Data('AAAAMMGGSPX_merged.csv', 'AAAA-MM-GG')` for each date (4–5 hours per date).
2. `Confronto_Date`: table across the dates.
3. `Run_Estensione_Date`: short-maturity extension on the other three dates.
4. `Run_Calendario`: the same extension with time measured in trading days.

**Numerical checks**

`Valida_Pricer`, `Valida_RHH`, `Valida_RHKappa` and `ControllaGrigliaXi` reproduce the tests of Appendix C.

## Main result (15 May 2026)

Relative IV RMSE on 105 quotes (7 maturities × 15 strikes):

| Model | Parameters | RMSE (%) |
|---|---|---|
| rough Hawkes Heston | 7 | 2.22 |
| Heston | 5 | 2.90 |
| rough Heston-κ + CGMY | 7 | 3.05 |
| rough Heston-κ | 4 | 3.24 |
| Heston-curve-κ | 3 | 3.29 |
| rough Heston + CGMY | 6 | 3.76 |
| rough Heston | 3 | 4.01 |
| CGMY | 3 | 5.59 |
| Heston-curve | 2 | 10.35 |

The full numbers for every date are in `Results_AAAAMMGG/riepilogo_AAAAMMGG.txt` and `confronto_date.txt`.
