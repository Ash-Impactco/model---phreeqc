# Prompt for extending the model (use in Claude / Fable, or Claude Code inside VS Code)

Open the `erw-phreeqc` folder, then paste everything below the line.

---

You are working in my repository `erw-phreeqc`. It holds a live PHREEQC model (phreeqpython + llnl.dat) of enhanced rock weathering on Azorean andosol pasture, with a Streamlit app (`app.py`), the model (`erw/model.py`), all inputs (`data/inputs.json`), a batch runner (`scripts/run_batch.py`) and a GitHub Actions workflow.

Start by reading `README.md`, `erw/model.py` and `data/inputs.json`. Then run `python scripts/run_batch.py` and confirm you reproduce `outputs/results.csv`: base case about 0.10 tCO2e/ha released and about 0.003 t/ha removed over 10 years. Do not change any results until that matches.

The audience is a university geology and lab team (INOVA, Azores), so every change must be cited, documented in the README and adjustable from `data/inputs.json` or the app sidebar. This remains a scenario model, not a prediction. British English throughout.

Make these improvements one at a time. After each, run the batch, show me the before/after numbers, and commit:

1. **Andosol acidity.** Replace the single calibrated HX pool with exchangeable Al (AlX3, with gibbsite or amorphous Al(OH)3 control) plus pH-dependent charge (a SURFACE block, or a second exchanger with a pKa). Inputs: exchangeable acidity (KCl), effective CEC, CEC at pH 7. Keep the old behaviour as an option.
2. **Mineralogy input.** Let `inputs.json` take measured XRD/Rietveld fractions (plagioclase composition, pyroxene, olivine, glass, oxides) instead of the norm, with a switch in the app.
3. **Glass kinetics.** Implement the full Gislason & Oelkers (2003) rate law, with the Al3+ dependence, for basaltic glass.
4. **Seasonal hydrology.** Monthly drainage for north-east São Miguel (cite the rainfall and evapotranspiration source). Run TRANSPORT in monthly chunks with varying time steps.
5. **Plant uptake and grazing.** An annual removal of Ca, Mg and K in harvested grass, using literature values for Azorean pasture, as an adjustable sink.
6. **Trace metals.** Add Ni, Cr and Co to the rock as trace components, released in proportion to olivine and pyroxene dissolution. Report topsoil accumulation against EU sewage sludge directive limits. Leave the ppm inputs blank until measured.
7. **Detectability panel.** A new app tab estimating the minimum detectable change in exchangeable Ca+Mg and in a (Ca+Mg)/Ti ratio for a given number of composite samples and a field coefficient of variation (inputs), so we can design the sampling plan.
8. **Sensitivity.** A batch mode that varies each key parameter one at a time and saves a tornado chart of 10-year CO2 removed.
9. **Documentation.** Update the README "How the model works" and "Main limitations" sections, and write `docs/methods.md` (3–4 pages) with the equations, parameter table (measured / stand-in / assumed), results and questions for INOVA.

Do not hard-code results, invent data, or present stand-in values as measurements. If something cannot be done reliably, say so and stop.
