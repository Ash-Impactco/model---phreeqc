How to use the templates
1. Copy a template into ../inputs/ with the same file name.
2. Replace the example rows with real data (keep the column names).
3. Run python3 run_all.py. Any file in inputs/ overrides the placeholder in config.json and is listed as measured.

feedstock_xrf.csv     lab XRF and ICP results on the Tecnovia composites (averaged across samples)
feedstock_xrd.csv     quantitative mineralogy; glass by difference
feedstock_psd.csv     sieve plus laser diffraction, mass fraction per size bin
feedstock_other.csv   moisture, inorganic carbon, sulfur
soil_screening.csv    one row per screening composite (December), duplicate flagged with is_duplicate = 1
soil_cores.csv        optional: individually analysed cores, grouped, to measure core-to-core variability
plot_results.csv      trial results per plot and round (raw values; the model computes changes from month 0)
                      run on its own with: python3 run_all.py analyse inputs/plot_results.csv
