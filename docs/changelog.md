# Changelog

Every change to the model, the data pack or the tooling, with its effect on the base case
(**S1**: Marques Britas 0/4 mm fines on the Moreira soil, 40 t/ha, 10 years). All figures are
model forecasts on stand-in data, in t CO2e per hectare. They are not creditable tonnes.

| Step | Change | Cations released | Leaving at 50 cm | Leaving the 20 cm NFZ | Isometric net removal |
|---|---|---|---|---|---|
| 0 | Model as supplied (`data/inputs.json`) | 0.0960 | 0.0027 | 0.0081 | −0.1635 |
| 1 | Inputs read from `data/lab`; grain-size bins interpolated from table 03 | 0.1037 | 0.0029 | 0.0088 | −0.1630 |
| 2.1 | Andosol acidity added as an option (base case unchanged) | 0.1037 | 0.0029 | 0.0088 | −0.1630 |
| | *S1 with the Andosol option switched on (scenario S9)* | *0.1055* | *0.0029* | *0.0095* | *−0.1623* |

The Step 0 numbers stay reproducible with `python scripts/run_batch.py --source json`, and a test
checks them on every run (`tests/reference/step0_results.csv`).

## Step 2.1: Andosol acidity (exchangeable Al plus pH-dependent charge)

**What changed**

A third acidity option, `acidity_mode = "andosol"`, next to the two existing pools (which are
unchanged and still the default):

- **Permanent-charge sites (X)** hold Ca, Mg, K, Na and **exchangeable Al** (AlX3). The Al in the
  soil water is held by an Al hydroxide phase, either amorphous Al(OH)3 (log K 10.8 at 25 °C, as in
  phreeqc.dat and wateq4f.dat; the default) or gibbsite (llnl.dat). When the pH rises, the
  exchangeable Al comes off and precipitates, using up three OH⁻ per Al.
- **Variable-charge sites (Xv)**, a second exchanger with one proton-binding constant: protonated
  and uncharged at low pH, taking up cations as the pH rises (allophane, imogolite, organic matter).
- **Calibration, from the soil report only:** charge at the soil pH = effective CEC (bases plus
  KCl acidity); charge at pH 7 = CEC at pH 7; exchangeable Al = KCl acidity. The exchange constant
  of AlX3 is fitted so the exchangeable Al is in equilibrium with the Al phase at the soil pH.
- **Still assumed:** the KCl acidity (1.0 cmol(+)/kg, not measured at Sanguinho), the half-point pH
  of the variable charge (midway between the soil pH and 7), and the choice of Al phase. All three
  are rows in `14_model_settings.csv` and controls in the app.
- A SURFACE block (surface complexation) was not used: the second exchanger does the same job,
  keeps the charge balance simple and adds no run time (16 s for 10 years).
- llnl.dat contains the Al13 polymer, which would reach molar concentrations at amorphous Al(OH)3
  saturation below pH 6. That is not what soil solutions show, and phreeqc.dat and wateq4f.dat
  leave the species out. It is switched off for the amorphous option only (written out in the
  `.pqi` file). Nothing else in llnl.dat is changed.

**Effect on results**

The base case keeps the large pool, so S1 does not move. Two scenarios were added to show the
effect (10 years, t CO2e/ha):

| Scenario | Released | Leaving the NFZ | Held on exchange (NFZ) | Net removal |
|---|---|---|---|---|
| S1 Base, large pool | 0.1037 | 0.0088 | 0.0967 | −0.1630 |
| S6 Base, small pool | 0.1144 | 0.0226 | 0.0940 | −0.1512 |
| **S9 Base, Andosol acidity** | 0.1055 | 0.0095 | 0.0972 | −0.1623 |
| S7 Optimistic, small pool | 2.1391 | 0.9488 | 1.1841 | +0.0961 |
| **S10 As S7, Andosol acidity** | 2.2176 | 0.2454 | 1.9646 | −0.5018 |

**What this shows**

- With the measured CEC at pH 7 honoured, the Andosol model behaves like the large pool, not the
  small one. The small pool ignores the charge that develops between the soil pH and pH 7.
- The optimistic scenario (S7) is net positive **only** under the small pool. With the same rock,
  grain size, glass and rate factor but Andosol acidity (S10), most of the released cations stay
  on the exchanger for the 10 years and the result is net negative. The acidity pool is therefore
  as important as the field rate factor, and KCl acidity and ECEC at Lomba da Maia are needed to
  choose between the pools.
- Fitted Al exchange constants (log K) for the four Sanguinho soils: 0.90, −0.40, 1.16 and 1.33
  with amorphous Al(OH)3, against 4.00, 2.72, 4.26 and 4.43 with gibbsite. The llnl.dat value is
  0.41. Amorphous control is therefore the more consistent choice for an assumed 1 cmol(+)/kg of
  exchangeable Al, which is why it is the default.
- In every acidity option the modelled soil-water pH drifts from the measured 5.8 to about 6.1-6.2
  in the first year, in the control as well. This comes from the starting soil water (0.1 mmol/L
  Ca, assumed) being more concentrated than the rain that replaces it. It was already in the
  model and is not changed here; it cancels in rock minus control only to first order.

## Step 1: data integration (lab files feed the model)

**What changed**

- `erw/data_io.py` loads every CSV in `data/lab/`. It reads values such as `<0.05`, `2850-3000` and
  `~11`, keeps blanks blank, and keeps the `status` of every row. A measurement on a rock or soil
  that is not the project's own (Marques Britas, Sanguinho) is reported as a **stand-in**, never
  as measured.
- The model now takes its inputs from the lab files: rock from table 01 (density from 02, grain
  size from 03, mineralogy from 05 when measured), soil from 08 or 09, site and climate from 11 and
  12, rate constants from 13, settings from 14, emission factors from 18. Each input carries its
  status and its source file into the run (`cfg["provenance"]`, shown in the app and written to
  `outputs/inputs_S1.csv`).
- `data/inputs.json` is now only a fallback for values the lab files do not hold, and the frozen
  Step 0 reference (`default_config(source="json")`).
- `erw/validate.py` writes `outputs/validation.md` (rules listed in the report).
- `erw/samples.py` and `scripts/run_samples.py` run the model once per soil composite in parallel,
  plus a site-mean run. Table 09 is still empty, so nothing runs for Lomba da Maia yet.
- `tests/` (pytest, 89 tests): loader, each validation rule, CO2 potential, regressions on S1.

**Effect on results, and why**

Only one thing moved the numbers: the grain-size bins. The original bins were typed in by hand
with round diameters (30, 150, 600, 1500, 3000 µm). They are now built from the PSD points in
table 03 by interpolating the grading curve as a straight line in log(size), and each bin gets
the diameter that keeps its surface area, ln(d2/d1) / (1/d1 − 1/d2). That gives 30, 116, 462,
1386 and 2773 µm for the same mass split: about 8 % more surface area, so about 8 % more
dissolution in S1. Everything else read from the lab files is identical to the original inputs
(a test checks this).

| Scenario | Released, before → after | Leaving the NFZ, before → after | Net removal, before → after |
|---|---|---|---|
| S1 Base | 0.0960 → 0.1037 | 0.0081 → 0.0088 | −0.1635 → −0.1630 |
| S2 Soil B | 0.1261 → 0.1361 | 0.0111 → 0.0120 | −0.1610 → −0.1602 |
| S3 Milled P80 1 mm | 0.1965 → 0.2044 | 0.0166 → 0.0173 | −0.2463 → −0.2457 |
| S4 Milled P80 100 µm | 0.3901 → 0.3997 | 0.0330 → 0.0339 | −0.6823 → −0.6816 |
| S5 Base + 20 % glass | 0.2054 → 0.2222 | 0.0265 → 0.0280 | −0.1478 → −0.1466 |
| S6 Small acidity pool | 0.1060 → 0.1144 | 0.0210 → 0.0226 | −0.1526 → −0.1512 |
| S7 Optimistic | 2.0917 → 2.1391 | 0.9125 → 0.9488 | +0.0652 → +0.0961 |
| S8 Upper bound | 4.3135 → 4.4036 | 3.0943 → 3.1879 | +1.9197 → +1.9993 |

Groundwater step (S7): share of the extra alkalinity reaching the stream 95.8 % → 95.3 %.

**Changes to the data pack** (all labelled, none are new measurements)

- `09_soil_lomba_da_maia_screening.csv`: every data row had one cell more than the header, so in
  a spreadsheet `TO MEASURE` sat under `lab`, and `15` (subsamples) under `bulk_density_g_cm3`.
  One empty cell was removed from each row. No value was lost. The loader now reports any row
  whose length does not match the header.
- `02_rock_physical_properties.csv`: one row added, the particle density the model already used
  for the fines (2900 kg/m³, **assumed**, within the 2850-3000 measured on the olivine basalt block).
- `03_rock_grain_size.csv`: two **assumed** grading curves added for the milling scenarios
  (`MILL-P80-1MM`, `MILL-P80-100UM`), so they are no longer hidden in `inputs.json`.
- `13_model_rate_parameters.csv`: a `status` column added (all **assumed**: literature lab rates).
- `14_model_settings.csv`: seven settings that were hidden in the code or in `inputs.json` are
  now rows with a status (glass fraction, assumed KCl acidity, assumed Na share, run length, the
  30 µm diameter assumed for the fraction finer than 63 µm, NFZ depth, η RiverOcean).
- `15_model_results_10yr.csv`: rewritten by the batch (`python scripts/run_batch.py --update-lab`).
- `Ilha_Carbon_ERW_Data.xlsx` was **not** regenerated and is now behind the CSVs for these rows.

**Findings worth a look**

- The diameter assumed for the fraction finer than 63 µm (30 µm) carries about two thirds of the
  rock's surface area in the model. It is the single most useful number a laser PSD will replace.
- The validation report fails `MB-POLISHED` (oxide total 93.2 %) and keeps the potassium conflict
  and the Mn2O4 question open. It cannot settle the potassium conflict from the data in hand.

## Step 0: set-up and reproduction

- Project unpacked into the repository root. All eight scenarios reproduce
  `data/lab/15_model_results_10yr.csv` exactly; the S7 coupled groundwater run delivers 95.8 %.
- `scripts/get_modflow.py`: falls back to the direct release download when the GitHub API is
  refused. No effect on results.
