# Ilha Carbon: live PHREEQC model of basalt weathering on Azorean pasture

An interactive enhanced rock weathering (ERW) model. You change the rock, soil and field settings, press **Run PHREEQC**, and the charts update year by year while PHREEQC runs. A 10-year run takes about 15–20 seconds.

**This is a scenario model, not a validated prediction.** The rock data comes from Marques Britas fines and the soil data from two Sanguinho pastures. Both are stand-ins until TECHNOVIA's rock analysis and the Lomba da Maia screening are in. The mineralogy is estimated from bulk chemistry, because there is no XRD yet.

## What is in here

| Path | What it is |
|---|---|
| `app.py` | The live app (Streamlit + Plotly) |
| `erw/model.py` | The PHREEQC model: mineral norm, kinetics, exchange, transport, results |
| `data/inputs.json` | **All input data**: rock oxides, grain size, soils, rate constants, defaults, references |
| `data/scenarios.json` | Scenarios for the batch run |
| `data/lab/` | **All project data**: 18 CSVs (rock, soil, site, climate, model, Isometric terms, emission factors) plus the same in one Excel file |
| `data/llnl.dat` | PHREEQC thermodynamic database (USGS, public domain) |
| `scripts/run_batch.py` | Runs every scenario and writes `outputs/results.csv`, `outputs/scenarios.png` and `.pqi` files |
| `.github/workflows/run-scenarios.yml` | Runs the batch on GitHub and lets you download the results |
| `.devcontainer/` | Lets the whole thing run in GitHub Codespaces with no installation |
| `erw/groundwater.py` | Groundwater step: MODFLOW 6 flow and transport coupled to PHREEQC (via mf6rtm / PhreeqcRM) |
| `scripts/get_modflow.py` | Downloads the MODFLOW 6 executable and library into `bin/` (run once) |
| `scripts/run_groundwater.py` | Runs the soil-to-stream chain for two scenarios, writing `outputs/groundwater.csv` |
| `FABLE_PROMPT.md` | Prompt for extending the model with Claude |

## Run it on your Mac (VS Code)

1. Open this folder in VS Code (File › Open Folder).
2. Open a terminal in VS Code (Terminal › New Terminal) and run:
   ```
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   python scripts/get_modflow.py
   streamlit run app.py
   ```
3. Your browser opens the app at http://localhost:8501.

From then on, `source .venv/bin/activate` then `streamlit run app.py` is all you need. You can also use the Run and Debug panel: choose **Run live app (Streamlit)**.

If `pip install phreeqpython` fails on your Mac, use Codespaces (next section). It runs on Linux, where phreeqpython always works.

## Run it on GitHub (no installation)

1. Create a new private repository on github.com, for example `erw-phreeqc`, and upload this folder. You can drag and drop it, or use `git init`, `git add .`, `git commit`, `git push`.
2. **Codespaces (live app in the browser):** on the repo page, click Code › Codespaces › Create codespace. It installs everything and opens the app on port 8501. The free monthly allowance is plenty for this.
3. **Actions (batch runs):** open the Actions tab › Run ERW scenarios › Run workflow. When it finishes, download `erw-results` (CSV, chart and `.pqi` files). It also runs automatically whenever you change anything in `data/`, `erw/` or `scripts/`.
4. **Share a link (optional):** at share.streamlit.io, sign in with GitHub, choose the repo and `app.py`, and deploy. INOVA can then use the app in a browser. Keep the repo private if the data is sensitive.

## Run it in the PHREEQC program itself

In the app's **Download** tab, save `erw_basalt.pqi` and `erw_control.pqi`, or take them from `outputs/pqi/` after a batch run. Then run:
```
~/Documents/phreeqc-3.5.0-14000/bin/phreeqc erw_basalt.pqi erw_basalt.out ~/Documents/phreeqc-3.5.0-14000/database/llnl.dat
```
Always use **llnl.dat**. phreeqc.dat does not contain the basalt silicates.

## Updating the data

Edit `data/inputs.json` when the real results arrive:
- **TECHNOVIA rock:** `feedstock.oxides_wt_pct`, plus `grain_size_bins` from the laser PSD
- **Lomba da Maia soil:** add a new entry under `soils` (pH, CEC, exchangeable Ca, Mg, K)
- **XRD / glass:** once measured, set the glass fraction in the app, or extend the model to take mineral fractions directly (see `FABLE_PROMPT.md`)

## How the model works

- **Column:** 0–50 cm in 10 cm cells, with basalt mixed into 0–20 cm. Steady drainage, fixed soil CO2 and O2. Gibbsite, kaolinite, amorphous silica, calcite and Fe(OH)3 can precipitate.
- **Rock:** a cation-conserving norm from the oxides (plagioclase, K-feldspar, diopside, forsterite, fayalite, inert Fe-Ti oxides), plus optional basaltic glass.
- **Dissolution rates:** Palandri & Kharaka (2004), acid and neutral terms, Arrhenius-corrected. Area is geometric surface area × roughness × field factor, with shrinking-core loss. Glass follows an approximate Gislason & Oelkers (2003) rate.
- **Exchange:** Gaines-Thomas. The soil's acidity is a proton pool calibrated to the measured base saturation, either the "CEC at pH 7" version or a small "effective CEC" version.
- **CO2 removed:** extra alkalinity leaving the column base, compared with an identical no-rock column. Cations still held on exchange sites are reported separately and not counted.

## Isometric alignment (EWA v1.3)

Every run reports an **Isometric accounting** tab, laid out like Equations 1, 2 and 4 of Enhanced Weathering in Agriculture v1.3, per hectare and per reporting period:

- **Terms:** Weathering, BiomassLoss, NetNewCarbonate, NetNewSilicate, NetSorption and NonCarbonicNeut, then export from the **Near-Field Zone (0–20 cm)**, × η RiverOcean, then Stored − Emissions = **Removal**. The counterfactual is subtracted throughout, because every term is a rock column minus a no-rock control column.
- **Removal is counted at 20 cm** (the NFZ boundary), as the protocol requires. It is not counted at 50 cm.
- **Placeholders:** project emissions use placeholder factors (`data/lab/18_emission_factors.csv`), and η RiverOcean is set to 0.85 until replaced.
- **Forecast only:**
  - Under v1.3, credits must be quantified from **soil or porewater measurements** (List 1).
  - Reactive transport models are not yet accepted as a quantification approach. They can support the sorption time-lag term (guidance G-V660-0, with the model and its data provided), sampling design and the PDD.
  - Credits are issued at about the 30th percentile (40th with validation) of the measured distribution, so credited figures will be lower than this central forecast.
- **What the forecast shows so far:** with stand-in data and placeholder emissions, the as-supplied 0/4 mm fines give a **net negative** result over 10 years (about −0.16 t CO2e/ha), because transport and spreading emissions exceed the small removal. Only the optimistic milled-and-glassy scenario is net positive (about +0.07 t/ha over 10 years).

## Groundwater step (MODFLOW 6 + PHREEQC)

The soil model gives the chemistry of water draining below the root zone. The groundwater step routes that drainage through a shallow basalt aquifer to a stream, using MODFLOW 6 for flow and transport and PHREEQC for the chemistry. The two are coupled through mf6rtm, which uses modflowapi and PhreeqcRM. It shows how long the extra alkalinity takes to reach the stream, how much is diluted, and how much is lost on the way, for example to calcite precipitation.

- **Set-up:** a one-layer flow line from the upslope edge of the field to a stream. Recharge equals the soil drainage. Field cells receive basalt-amended drainage; the rest receive control drainage. It runs as a pair: rock and control.
- **All hydrogeology is ASSUMED:** distance to stream 1,500 m, flow path under the field 600 m, thickness 20 m, hydraulic conductivity 10 m/day, porosity 0.15, background groundwater chemistry. Replace these with Azorean data, for example from the regional water authority (DROTRH) or published hydrogeology of São Miguel.
- **Simplification:** the final soil year's drainage is applied as constant recharge.
- **What it tells you so far:** in the optimistic soil scenario, about 96% of the extra alkalinity reaches the stream, after accounting for dilution by recharge outside the field, with no calcite loss. In the base case the signal (about 1 µeq/L) is too small to resolve.
- **In the app:** after a soil run, use the "Groundwater step" section. **In batch:** `python scripts/run_groundwater.py`.

## Main limitations

Mineralogy is normative. Rates are lab rates scaled by an assumed field factor. Andosol acid buffering is simplified to one pool. Flow is steady, with no seasons. There is no plant uptake or grazing removal, and no trace metals.
