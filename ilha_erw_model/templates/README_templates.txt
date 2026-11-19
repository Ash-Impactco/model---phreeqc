How to use the templates
1. Copy a template into ../inputs/ with the same file name.
2. Fill it with real results (keep the column names; leave a cell empty if not measured).
3. Run python3 run_all.py. Until every required value is filled, the run lays out the plots and
   then stops, listing what is still missing (outputs/missing_inputs.csv).

feedstock_xrf.csv     Tecnovia composites: one row per analyte per sample, averaged by the model.
                      Oxides in wt%: SiO2 Al2O3 Fe2O3 CaO MgO TiO2 Na2O K2O
                      Trace in ppm:  Zr Ni Cr Co Cu Zn Pb As Cd Hg V
feedstock_xrd.csv     quantitative mineralogy (Rietveld), glass by difference; mineral names as listed
feedstock_psd.csv     sieve plus laser diffraction: mass fraction per size bin (add finer bins if measured)
feedstock_other.csv   moisture (wet basis), inorganic carbon, total sulfur
soil_screening.csv    one row per December screening composite; field duplicate flagged is_duplicate = 1.
                      pH_buffer_cmolc_kg_pH = titration on the composites. 4 or more composites are
                      needed to estimate plot-to-plot variability.
soil_cores.csv        individually analysed cores (5 to 6 per composite, 2 or more composites), with
                      'group' = composite id: gives core-to-core variability
site_other.csv        background acid load (from farm N and lime records) and lab precision;
                      optional site climate
plot_results.csv      trial results per plot and sampling round (raw values; month 0 = baseline).
                      Analyse on its own with: python3 run_all.py analyse inputs/plot_results.csv
