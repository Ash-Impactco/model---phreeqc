# Model map: what each part does

Model: `ilha_erw_model/`, version 1, first run on placeholder values (19 November 2026).
Every soil and rock number below is a placeholder until the December screening and the Tecnovia fines analyses are loaded.

## run_all.py (the conductor)
Reads `config.json`, then checks `inputs/` for filled-in lab CSVs and lets any it finds overwrite the placeholders (and marks them "measured"). Then calls the modules below in order and builds the report. **In:** config.json, inputs/*.csv. **Out:** everything in outputs/ and report/. **Assumption that matters most:** when a lab CSV exists, the model uses the simple average of all rows (for example, mean soil pH across all screening composites) for the whole trial. Plot-to-plot differences in starting soil are not carried into the chemistry, only into the noise used for power.

## erw/common.py (plumbing)
Loads config.json, keeps track of which values were overwritten, writes the parameter status table, and holds the black-and-grey plotting style. No science in here.

## erw/site.py (where the plots go)
Takes the QGIS layers (units U1 to U16, trees + 5 m no-spread zones, drainage lines, DEM, slope). Keeps only units U3 to U5, cuts out no-spread areas and a 10 m buffer around drainage lines, then tries a grid of 80 x 80 m plots with 8 m alleys at several rotations and offsets along the contour and keeps the grid that fits the most plots. Drops plots steeper than 20% on average, picks the 16 that form the most compact groups of 4 (the blocks), orders blocks by elevation, then randomises rates within each block, redrawing until rate is not lined up with slope or position on the slope. Works out tonnes per plot. **In:** site_data/, design settings. **Out:** trial_plots.csv / .geojson / .gpkg, layout map, tonnage summary. **Assumption that matters most:** that grouping plots into blocks by distance and elevation captures the main soil differences. If the soil changes in a different direction (say, with old field boundaries or past manure spreading), the blocks will not soak up that variation and power will be lower than modelled.

## erw/geochem.py (the chemistry engine, the heart of the model)
Five jobs:
1. **Rock potential.** Turns the rock's oxide analysis (XRF) into how many "charge units" of calcium, magnesium, sodium and potassium a tonne holds (equivalents per kg), and so how much CO2 it could capture at most (1 mol CO2 per equivalent), minus rock carbonate and sulfur, times the share that stays stored downstream (0.86). Result: about 320 kg CO2 per tonne.
2. **Minerals.** Splits the rock into olivine, augite, plagioclase (labradorite), glass and inert phases, and gives glass whatever cations are left over after the crystals.
3. **Dissolution over time.** Each mineral in each grain-size class shrinks from the outside in (shrinking core), at a lab-measured speed that depends on soil pH and temperature, multiplied by surface roughness (7) and by the field factor. Monthly time step.
4. **Soil and carbon.** Each month, the cations released are split: 10% lost to new clay minerals, 5% taken up by grass, then the remainder divided into bicarbonate (carbon removed) by the share `1/(1+10^(6.35 - pH))`, and the rest neutralising soil acidity (70% held on exchange sites, 30% leaving with nitrate). The acid-neutralising part raises pH through the buffer capacity; the fertiliser acid load pushes pH down on every plot. The carbon account then subtracts rock carbon, sulfur and trucking/spreading emissions.
5. **Monte Carlo and metals.** Draws 600 combinations of the uncertain inputs to give P10/P50/P90 and a ranking of what matters; and a simple mixing calculation of trace metals in the top 10 cm against Portuguese sludge-law limits.

**In:** feedstock, kinetics, soil, emissions, limits sections of config.json. **Out:** feedstock_potential.csv, feedstock_mineral_balance.csv, carbon_account_by_scenario.csv, trial_total_CDR_30m.csv, mc_*.csv, trace_metal_mass_balance.csv, figures 2 to 4. **Assumption that matters most:** the field factor (how fast rock weathers in the field compared with a lab beaker). It alone explains most of the spread in the carbon result (rank correlation 0.83). Close second, hidden inside the grain sizes: about 70% of everything released in 30 months comes from the grains finer than 63 micrometres, and half from grains finer than 10 micrometres, so the shape of the finest sieve fraction matters enormously.

## erw/mrv.py (what we can actually measure)
Asks "would we see it?" for each measurement:
- **Tracer contrast:** is the rock different enough from the soil in Ti and Zr to use them as markers of how much rock is in a core? (No: Ti signal-to-noise 0.25 to 0.93; Zr goes the wrong way because soil Zr is higher than rock.)
- **Soil mass balance:** simulates cores before and after spreading and estimates the dissolved fraction two ways (from the Ti tracer, or from the dose measured in catch trays). Reports the smallest dissolved fraction each could detect.
- **Water:** turns the bicarbonate leaving each plot into a concentration in drainage water at plot scale (lysimeters) and at the six stream points W1 to W6.
- **Mesocosms:** soil columns in the lab, 0 to 160 t/ha, 12 months, at 20 C without fertiliser acid; expected leachate alkalinity against the lab detection limit.
- **Sampling plan:** counts of samples per round.

**In:** central run from geochem, soil heterogeneity and lab precision, water points. **Out:** tracer_contrast.csv, tau_detectability.csv, water_signal.csv, mesocosm_design.csv, sampling_schedule_and_counts.csv, figure 5. **Assumption that matters most:** that core-to-core variability is 25% CV for CaO and 20% for MgO, and that spreading is 25% uneven at core scale. These are placeholders, and the whole "can't see dissolution in soil" conclusion rests on them.

## erw/stats.py (how many plots, and the analysis we will run)
Takes the expected change at 30 months for three responses (soil pH, exchangeable Ca + Mg, total Ca + Mg) from the central chemistry run, adds realistic noise (plot-to-plot differences plus sampling and lab error), and simulates the trial 300 times for each of four designs (10, 12, 15, 16 plots) at low, medium and high variability. Power = share of simulated trials where the dose-response slope is significant at p < 0.05. Also contains the analysis that will be run on the real data: regression on rate with blocks, curvature test, MANOVA, each rate vs control with Holm correction, and a mixed model over sampling rounds. Runs it on a simulated example so the code is tested. **In:** central chemistry, soil noise placeholders, design list. **Out:** power_by_design.csv, expected_effects_30m.csv, plot_level_noise.csv, example analysis and figures 6 and 7. **Assumption that matters most:** that the true effect is the central chemistry effect. Power is calculated as if we already knew the field factor is 1. If it is 0.3, the pH effect at 40 t/ha drops from 0.34 to 0.13 and power falls a long way.

## erw/report.py (the PDF)
Writes the 14-page HTML report from the tables and figures and converts it to PDF with Playwright's Chromium. No science. One plumbing bug: it calls `python3` by name for the PDF step, so it fails inside a virtual environment unless that environment is first on PATH.

## config.json (the single source of truth)
Every number in the model, each with a status: placeholder (replace with site or lab data), literature (published value), design (our choice), measured (our own data). 68 entries are tracked in outputs/parameter_status.csv: 32 placeholders, 31 design choices, 5 literature values, 0 measured.
