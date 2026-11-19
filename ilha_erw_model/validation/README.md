# Benchmark against published basalt trials

19 November 2026. Which of the two Lomba da Maia models agrees with what has actually been measured in the real world?

- **Trial model**: `ilha_erw_model/` (this folder's parent), run here under each study's conditions with `validation/benchmark_literature.py`. The model code was not changed.
- **Earlier PHREEQC model**: repository root (`app.py`, `erw/model.py`). Not re-run; its own reported results (`data/lab/15_model_results_10yr.csv`) are used.

There are no measurements from Lomba da Maia yet. The benchmarks are five published studies: two of them measure how fast basalt dissolves, three measure how much CO2 is actually removed.

## Results

| Study | What was measured | Measured | Trial model (field factor 1) | Trial model (field factor 0.1) | Earlier PHREEQC model |
|---|---|---|---|---|---|
| Reershemius et al. 2023, UK mesocosm, 50 t/ha, p80 35 um, soil pH 6.45, 235 days | Share of rock cations dissolved | **15.7 +/- 3.1 %** | 10.7 % | 1.4 % | No matching run; even milled 100 um rock only 2.7 % in 10 years |
| Beerling et al. 2024, Illinois field, 4 x 50 t/ha, p80 267 um | Ca + Mg lost from the basalt | **16 +/- 6 %** | 6.1 % | 1.0 % | 0.67 % in 10 years (0/4 mm); 2.7 % (100 um) |
| Holden et al. 2024, Australian field, acidic soil pH 5.8, <5 mm basalt, 50 t/ha/yr for 5 years | Extra bicarbonate leaving the soil, t CO2/ha/yr | **0.026, not significant** | 0.69 | 0.13 | About 0.0003 (per 40 t/ha one-off) |
| Dupla et al. 2025, Swiss field trials, 3 years | Realised CO2 removal, t CO2/ha/yr | **0.10 +/- 0.03** | 0.04 (40 t/ha), 0.10 (80 t/ha) | 0.007, 0.015 | 0.0003; 0.026 in its optimistic scenario |
| Vienne et al. 2026, mesocosm, up to 200 t/ha, 101 days | Realised CO2 removal | **None detected**; cations held in the soil | 0.04 t/ha at 50 t/ha | 0.004 | Near zero |

The full table with every assumption is in `benchmark_literature.csv`.

## What it says

**1. How fast the rock dissolves: the trial model is close, the earlier model is far too slow.**
Both direct measurements (Reershemius, Beerling) show 16% of the cations gone within months to a few years. The trial model at field factor 1 gives 6 to 11%, so it is if anything a little slow (by a factor of about 1.5 to 2.5), even though the Illinois rock has fewer fast-weathering minerals than ours. The earlier PHREEQC model gives under 3% in ten years, 10 to 100 times slower than measured. Its field factor of 0.1 is contradicted by both studies.

**2. How much CO2 is actually removed in acidic soil: the trial model is too optimistic.**
The closest real-world analogue to Lomba da Maia is the Australian trial (acidic soil, coarse unmilled basalt). After five years and 250 t/ha it measured an extra 0.026 t CO2/ha/yr, not statistically different from zero, because strong acids (soil acidity, fertiliser nitrate, organic acids) did most of the weathering. The trial model predicts 0.69 under those conditions: about 25 times too high. In the Swiss trials it is about right (0.04 to 0.10 vs 0.10). In the 101-day mesocosm it predicts some removal where none was seen, because it lets bicarbonate leave the top 10 cm at once with no lag. The earlier model predicts almost no removal, which matches acidic soils, but only because it barely dissolves the rock: it gets the right answer for the wrong reason, and is 300 times below the Swiss result.

**3. Where the cations go: both models point the right way.**
Real soils keep most released cations (exchange sites and iron oxide coatings; Vienne et al. found the oxide-bound pool about 10 times the exchangeable one). The trial model keeps about half on exchange sites and 10% in secondary minerals; the secondary share is probably too small.

## Verdict

- **Use the trial model as the reference.** It reproduces measured dissolution within a factor of about 2, and dissolution drives everything the trial will actually measure in 30 months (soil pH, exchangeable Ca + Mg, total Ca + Mg, lysimeter cations).
- **Treat its carbon removal as an upper bound in acidic soil.** The real-world evidence for acidic soils (Holden) says realised removal in the first years is close to zero. That strengthens the existing conclusion: no carbon claims for the grant period.
- **Do not use the earlier PHREEQC model's dissolution rates** (field factor 0.1) without recalibration. Its soil exchange chemistry is still useful as a reference for where the cations go.

## Proposed model changes (not applied; need approval)

1. Central field factor: the benchmarks support 1 or slightly above, not lower. Keep 1 as central, narrow the low case from 0.2 to about 0.5.
2. Bicarbonate export: add a lag or a subsoil layer so alkalinity is not counted until it passes the rooting zone, and test the bicarbonate share against the Holden result.
3. Secondary / oxide-bound retention: raise the 10% placeholder range to include much higher values (Vienne et al.).

## Limits of this check
Each study has its own rock, soil, climate and crop. Particle sizes were rebuilt from a single p80 value; temperatures for the mesocosms and the Illinois soil pH were assumed; repeated annual applications were approximated by adding up independent yearly batches. The Dupla comparison uses Lomba da Maia conditions because their rate and soil were not retrieved. Two papers (Holden, Dupla) were read from their abstracts only. This is an order-of-magnitude check, not a calibration.

## Sources
- Reershemius, T. et al. (2023) Initial validation of a soil-based mass-balance approach for empirical monitoring of enhanced rock weathering rates. Environ. Sci. Technol. 57, 19497-19507. https://arxiv.org/abs/2302.05004
- Beerling, D.J. et al. (2024) Enhanced weathering in the US Corn Belt delivers carbon removal with agronomic benefits. PNAS 121, e2319436121. https://doi.org/10.1073/pnas.2319436121
- Holden, F. et al. (2024) In-field carbon dioxide removal via weathering of crushed basalt applied to acidic tropical agricultural soil. Sci. Total Environ. https://www.sciencedirect.com/science/article/pii/S004896972406724X
- Dupla, X., Bertagni, M.B., Grand, S. (2025) Three years of field trials indicate a sustained enhanced rock weathering signal with limited CO2 removal. Environ. Sci. Technol. 59, 25751-25764. https://doi.org/10.1021/acs.est.5c09820
- Vienne, A. et al. (2026) Weathering without realizing inorganic CO2 removal revealed through base cation monitoring. SOIL.
