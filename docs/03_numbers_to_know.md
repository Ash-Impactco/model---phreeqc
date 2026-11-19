# 15 numbers to know by heart

All are model outputs from **placeholder** soil and rock values (literature soil, Marques Britas basalt). Say "placeholder" every time you quote one.

| # | Number | Where it comes from | How sure |
|---|---|---|---|
| 1 | **16 plots, 4 blocks, 10.24 ha** (0.64 ha each, 80 x 80 m) | site.py fitted to the open ground of U3 to U5, avoiding trees + 5 m, drainage + 10 m, slopes over 20% | Solid: built from the QGIS layers. Check on the ground before marking out |
| 2 | **0, 20, 40, 80 t/ha** | Design choice; spread wide enough to fit a dose response and test curvature | Design, open for discussion |
| 3 | **358 t dry; about 430 t delivered (18 loads)** | 16 plots x rate x 0.64 ha, +10% contingency, 8% moisture | Solid arithmetic; moisture is a placeholder |
| 4 | **320 kg CO2 per t** | XRF oxides to charge (8.5 eq/kg), minus rock carbon and sulfur, x 0.86 downstream retention. Ca + Mg only: 263 | Ceiling, not a forecast. Changes when Tecnovia XRF arrives |
| 5 | **About 7% dissolved at 30 months** (honest range 3 to 13%) | Shrinking-core kinetics, central field factor 1 | Weak: the field factor alone moves it from 2.6% to 11.5% |
| 6 | **70% of early dissolution comes from grains under 63 micrometres** (half from under 10) | Grain-size classes in geochem.py | The share is model-solid; the fines shape is assumed. Laser diffraction needed |
| 7 | **About 10% of released cations leave as bicarbonate at pH 5.4** | 1/(1+10^(6.35 - pH)); 22% at 5.8, 41% at 6.2 | Rule of thumb; cross-check against soil CO2 gives the same order |
| 8 | **4.3 kg CO2e per t** for haul, loading and spreading (+1.8 rock carbon, +0.5 sulfur) | 26 km x 0.12 kg/tkm + 0.4 + 0.8 | Placeholder haul data; about 6% low because we haul wet tonnes |
| 9 | **Whole trial net about -1.3 t CO2 at 30 months** (range -2.1 to +0.2) | Gross 1.1 t minus about 2.4 t of deductions | Direction is robust: in every central-ish case the grant period is near zero or negative |
| 10 | **Positive only after years, at 40 and 80 t/ha** (10 years: +0.02 and +0.25 t/ha) | Same model run to 120 months | Weak in size, robust in shape |
| 11 | **pH +0.44 vs control at 40 t/ha; exchangeable Ca + Mg +1.5 cmolc/kg** (30 months) | Soil block of geochem.py (+0.34 vs own baseline, controls drift -0.10) | Depends on buffer capacity and field factor (+0.23 at field factor 0.3) |
| 12 | **Lysimeter alkalinity +0.07 (40 t/ha) and +0.20 meq/L (80 t/ha)** vs lab limit 0.05 | Bicarbonate flux / 1,200 mm drainage | Upper bound: subsoil below 10 cm may consume some |
| 13 | **Soil mass balance: smallest detectable 0.17 to 0.35, expected 0.07** | Simulated cores, mrv.py | Conclusion robust; actually worse than shown (exchangeable cations still counted) |
| 14 | **Power: 16 plots about 0.95, 12 plots about 0.8** (medium variability) | 300 simulated trials per design, stats.py | Conditional on the central field factor and placeholder variability |
| 15 | **Field factor, starting pH, fines: rank correlations 0.83, 0.35, 0.31** | Monte Carlo, 600 runs, 9 inputs | Ranking robust; exact values not. These are the three things to measure first |

Bonus for Luis: metals at 80 t/ha raise soil Ni by about 10% and Cr by about 20% over (placeholder) background; the background itself already exceeds the sludge-law screening limits, which is common in volcanic soils and has to be confirmed in the December screening.
