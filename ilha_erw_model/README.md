# Ilha Carbon: Lomba da Maia basalt rate trial, soil and geological model

Fresh base, 19 November 2026. The model code is ready; the soil and rock data are not in yet.

## Status
- **Site and design: ready.** Plot layout, blocks, randomisation and dry tonnage run from the QGIS layers in `site_data/`.
- **Soil and rock: waiting for data.** Every soil and rock value in `config.json` is empty (status `missing`). Until the December soil screening and the Tecnovia fines analyses are in `inputs/`, a run lays out the plots and then stops, listing what is missing in `outputs/missing_inputs.csv`.
- **Model parameters still assumed** (status `placeholder` in `config.json`): field factor, exchange retention, secondary-mineral and plant-uptake losses, climate, haulage and emission factors. These are calibrated later from the mesocosms and first field rounds, or replaced when known.
- **Not yet reviewed** by the MRV lead. `validation/` compares the model with published basalt trials.

## Run
    pip install numpy pandas scipy statsmodels matplotlib geopandas rasterio shapely pyproj playwright
    python3 run_all.py            # full run once inputs are complete, about 3 minutes; writes outputs/ and report/
    python3 run_all.py --quick    # fewer simulations, for checking changes
    python3 run_all.py analyse inputs/plot_results.csv   # analyse trial results only

The PDF report step uses Playwright's Chromium (`python3 -m playwright install chromium` once).

## Adding data
1. Copy the matching file from `templates/` into `inputs/` (same name).
2. Fill it with lab or field results. `templates/README_templates.txt` explains each file.
3. Run `python3 run_all.py`. Values read from `inputs/` are marked `measured` in `outputs/parameter_status.csv`.

| File | Source | Fills |
|---|---|---|
| `feedstock_xrf.csv` | Tecnovia composites, fusion XRF and ICP-MS | rock oxides, trace metals |
| `feedstock_xrd.csv` | quantitative XRD, glass by difference | mineral mode |
| `feedstock_psd.csv` | sieve plus laser diffraction | grain-size distribution |
| `feedstock_other.csv` | moisture, inorganic carbon, total sulfur | deductions, wet tonnage |
| `soil_screening.csv` | December screening composites | soil pH, buffer capacity, bulk density, CEC, base saturation, total elements, plot-to-plot variability |
| `soil_cores.csv` | individually analysed cores | core-to-core variability |
| `site_other.csv` | farm N and lime records, lab QC | background acid load, lab precision, optional climate |

## Folders
    config.json   every parameter, with value, status (missing / placeholder / literature / design / measured) and a note
    site_data/    QGIS layers (units, trees + 5 m, drainage, water points, DEM, slope)
    erw/          site.py (layout), geochem.py (feedstock, kinetics, soil, carbon, metals, Monte Carlo),
                  mrv.py (measurement design, mesocosms, sampling), stats.py (power and analysis), report.py
    templates/    empty input forms
    inputs/       filled forms go here
    outputs/      layout files now (trial_plots.geojson, .gpkg, .csv, map); all tables and figures once inputs are complete
    validation/   comparison with published basalt trials
    report/       PDF report, built once inputs are complete

## Change the design
Edit `config.json`: candidate units, plot size, alley, rates, number of blocks, seed. Layout, tonnage and every downstream number are recomputed. Open `outputs/trial_plots_EPSG5015.gpkg` in QGIS to check plots on the ground before marking them out.
