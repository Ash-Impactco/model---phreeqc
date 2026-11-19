# Ilha Carbon: Lomba da Maia basalt rate trial, soil and geological model

Version 1, 19 November 2026. First run on placeholder values.

## Run
    pip install numpy pandas scipy statsmodels matplotlib geopandas rasterio shapely pyproj playwright
    python3 run_all.py            # full run, about 3 minutes; writes outputs/ and report/
    python3 run_all.py --quick    # fewer simulations, for checking changes
    python3 run_all.py analyse inputs/plot_results.csv   # analyse real trial data only

The PDF step uses Playwright's Chromium (`python3 -m playwright install chromium` once). If it is missing, the HTML report is still written.

## Folders
    config.json   every parameter, with value, status (placeholder / literature / design / measured) and a note
    site_data/    QGIS layers from Lomba_da_Maia_Plots_v2 (units, trees + 5 m, drainage, water points, DEM, slope)
    erw/          site.py (layout), geochem.py (feedstock, kinetics, soil, carbon, metals, Monte Carlo),
                  mrv.py (measurement design, mesocosms, sampling), stats.py (power and analysis), report.py
    templates/    one CSV per data source; fill and copy to inputs/
    inputs/       filled CSVs here override the placeholders (empty in this version)
    outputs/      tables (.csv), figures (.png), trial_plots.geojson (WGS84) and trial_plots_EPSG5015.gpkg
    report/       Ilha_ERW_Model_Report.pdf and .html

## Change the design
Edit config.json: candidate units, plot size, alley, rates, number of blocks, seed. The layout, tonnage,
power and every downstream number are recomputed. Open outputs/trial_plots_EPSG5015.gpkg in the QGIS project
to check plots on the ground before marking them out.

## Main caveat
Until the December soil screening and the Tecnovia fines results are in inputs/, every number is illustrative.
outputs/parameter_status.csv lists what is still a placeholder.
