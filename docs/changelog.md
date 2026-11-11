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
| 2.2 | Glass: full Gislason & Oelkers (2003) rate law (S1 has no glass, so unchanged) | 0.1037 | 0.0029 | 0.0088 | −0.1630 |
| 2.3 | Correction: the year-end state was counted twice in the export sum | 0.1037 | 0.0028 | 0.0083 | −0.1633 |
| | *S1 with seasonal drainage switched on (scenario S11, assumed pattern)* | *0.1039* | *0.0027* | *0.0081* | *−0.1635* |
| 2.4 | Plant uptake sink added, zero until the share leaving the field is known | 0.1037 | 0.0028 | 0.0083 | −0.1633 |
| | *S1 with all grass carried off (scenario S12, upper bound)* | *0.1185* | *0.0029* | *0.0087* | *−0.1630* |
| 2.5 | Trace metals (post-processing only; no effect on the chemistry) | 0.1037 | 0.0028 | 0.0083 | −0.1633 |

The Step 0 numbers stay reproducible with `python scripts/run_batch.py --source json`, and a test
checks them on every run (`tests/reference/step0_results.csv`).

## Step 2.5: trace metals

**What changed**

- `erw/metals.py` and a "Trace metals" tab in the app. For each element with a value in
  `06_rock_trace_metals.csv` it reports the rise in the topsoil's total content (all the metal in
  the rock, mixed into 0-20 cm) and, for Ni, Cr and Co, the part released as olivine and pyroxene
  dissolve. Both are in mg per kg of soil, set against the baseline and the limit in the new
  `20_trace_metal_limits.csv`.
- The model now records the rock mass left in each mineral and each grain-size bin every year
  (also needed for the dissolution charts in Step 3).
- **All inputs are blank**: table 06 has no values yet and table 20 waits for your limits. With
  blank inputs the tab says so and shows no numbers. Nothing was assumed.

**Assumptions built in** (both on the high side for accumulation): all Ni, Cr and Co is taken to
sit in olivine and pyroxene, and released metals stay in the topsoil.

**Effect on results:** none. For scale only: at 40 t/ha into 1,600 t/ha of topsoil, every
100 mg/kg in the rock adds about 2.4 mg/kg to the soil's total content. In S1 only 0.66 % of the
olivine and pyroxene dissolves in 10 years, so the released part would be under 1 % of that.

## Step 2.4: plant uptake and grazing

**What changed**

- A sink that removes Ca, Mg and K from the root zone (0-20 cm) every year, in the rock column and
  in the control column. Each cation is taken up in exchange for H⁺, as roots do, so removal
  acidifies the soil. Uptake is demand-driven and supply-limited (Michaelis-Menten), so it takes
  close to the stated amount unless the soil water runs short.
- Inputs are in the new table `19_plant_uptake_and_grazing.csv`. Net removal = pasture yield ×
  herbage content × share leaving the field, plus milk sold × milk content.
- The Isometric term **BiomassLoss** is now calculated: extra cations removed in biomass, rock
  minus control.

**Status of the inputs: this could only be done in part**

| Input | Value | Status |
|---|---|---|
| Pasture yield | 6 (5-7) t DM/ha/yr | **stand-in**: São Miguel dairy farms, Morais et al. (2018). Read from a search summary, not from the paper: check it |
| Ca, Mg, K in herbage | 4.0, 2.0, 30.5 g/kg DM | **assumed**: general extension table, not Azorean. I found no Azorean herbage analyses |
| Ca, Mg, K in milk | 1.2, 0.1, 1.3 g/L | **assumed**: general dairy text |
| Share of herbage cations leaving the field | blank | **to find**. 1 for silage or hay taken away, far less under grazing |
| Milk sold per hectare | blank | **to find** (farm records) |

Because the last two are unknown, the sink is **zero in the base case** and BiomassLoss stays 0.
I did not guess a grazing share.

**Effect on results** (10 years, t CO2e/ha)

| Scenario | Released | Leaving the NFZ | Removed in grass, rock / control | BiomassLoss | Net removal |
|---|---|---|---|---|---|
| S1 Base, no removal | 0.1037 | 0.0083 | 0 / 0 | 0 | −0.1633 |
| **S12 All grass carried off** | 0.1185 | 0.0087 | 1.660 / 1.659 | 0.0005 | −0.1630 |

**What this shows**

- Even the upper bound barely changes the removal figures. The grass takes the same amount with
  and without rock, so the difference (BiomassLoss) is close to zero, and the extra acidity from
  harvest raises rock dissolution by 14 %.
- The gross uptake (24 kg Ca, 12 kg Mg, 183 kg K per hectare per year) equals about 0.30 t CO2e/ha
  each year, three times what the rock releases in 10 years in S1. So for detection, the yearly
  cycling of K through grass, dung and urine is a far larger signal than the weathering. This
  matters for the sampling design more than for the accounting.
- In S12 the grass can take only 55 % of the stated amount over 10 years, because the model has no
  fertiliser, slurry or soil-mineral supply and the exchangeable K runs down. Treat S12 as a
  bound, not a forecast.

## Step 2.3: seasonal drainage, and a correction to the export sum

**Correction (changes every scenario's export by about −5 %)**

While adding seasons I found that the model as supplied over-counted the alkalinity leaving the
column. PHREEQC writes the state *before* the first shift of every TRANSPORT block as well as the
state after each shift. The model ran one block per year and summed every row, so the state at
each year end was counted twice: 21 states per year for 20 shifts. The water that leaves at a
shift is the water that was in the cell before it, so the sum now runs over the states at the
start of each shift. With many blocks per year (seasons) the old sum would have been about 30 %
too high. Released cations and cations held on exchange are not affected.

The frozen Step 0 reference keeps the original sum (`export_accounting = "step0"` in
`legacy_config`), so `--source json` still reproduces the supplied numbers. A new test checks that
released = held on exchange + leaving the NFZ + calcite, within 3 %.

**Seasonal drainage (new option, off by default)**

- The yearly drainage is spread over the months by the shares in `12_climate_and_hydrology.csv`.
  A pore-water shift happens each time the cumulative drainage fills one cell, so wet months get
  short time steps (about 12 days) and dry months long ones (about 36 days). Yearly totals and
  the number of shifts are unchanged. The run starts in October (spreading window, table 14).
- **The monthly shares are an assumption**, not station data: a two-level pattern with 75 % of
  the drainage from October to March. I could not read the IPMA climate normals from this session
  (the site is blocked here), so table 12 still says `TO FIND` for rainfall and
  evapotranspiration, and seasonal drainage stays **off** in the base case.
- Temperature is still constant through the year.

**Effect on results** (10 years, t CO2e/ha, before → after the correction)

| Scenario | Leaving at 50 cm | Leaving the NFZ | Net removal |
|---|---|---|---|
| S1 Base | 0.0029 → 0.0028 | 0.0088 → 0.0083 | −0.1630 → −0.1633 |
| S2 Soil B | 0.0039 → 0.0037 | 0.0120 → 0.0113 | −0.1602 → −0.1608 |
| S3 Milled P80 1 mm | 0.0058 → 0.0055 | 0.0173 → 0.0164 | −0.2457 → −0.2465 |
| S4 Milled P80 100 µm | 0.0113 → 0.0107 | 0.0339 → 0.0321 | −0.6816 → −0.6831 |
| S5 Base + 20 % glass | 0.0118 → 0.0112 | 0.0280 → 0.0265 | −0.1466 → −0.1478 |
| S6 Small acidity pool | 0.0075 → 0.0071 | 0.0226 → 0.0215 | −0.1512 → −0.1521 |
| S7 Optimistic | 0.5549 → 0.5268 | 1.4619 → 1.3893 | +0.5322 → +0.4705 |
| S8 Upper bound | 2.3169 → 2.2069 | 3.6242 → 3.4517 | +2.3701 → +2.2236 |
| S9 Andosol acidity | 0.0029 → 0.0028 | 0.0095 → 0.0090 | −0.1623 → −0.1628 |
| S10 As S7, Andosol | 0.0258 → 0.0244 | 0.0989 → 0.0939 | −0.6263 → −0.6306 |
| **S11 Base, seasonal drainage (new)** | 0.0027 | 0.0081 | −0.1635 |

**What this shows**

Seasonal drainage on its own changes little (S11 against S1: 0.0081 against 0.0083 leaving the
NFZ), because the cations are buffered by the exchanger and the temperature is constant. It will
matter more once monthly temperature and measured rainfall are in.

## Step 2.2: basaltic glass, full Gislason & Oelkers (2003) rate law

**What changed**

Glass now dissolves at r = A · exp(−Ea/RT) · (a_H⁺³ / a_Al³⁺)^(1/3), with A = 10^−5.6 mol Si/cm²/s
(10^−1.6 per m²) and Ea = 25.5 kJ/mol, per **geometric** surface area (Gislason & Oelkers 2003,
Geochim. Cosmochim. Acta 67, 3817-3832; fitted by the authors for 6-300 °C and pH 1-11). The
constants are in `13_model_rate_parameters.csv`; the choice of law is a row in
`14_model_settings.csv` and a control in the app. The older constant rate is kept as the option
`approximate`.

- The roughness factor is **not** applied to glass under the full law, because the published
  constant is already per geometric area. The field rate factor still applies.
- The law is the far-from-equilibrium form. No chemical affinity term is included.
- At the first moment of a run there is no dissolved Al, so the Al activity at saturation with the
  soil's Al phase is used until dissolution supplies Al.

**Effect on results** (10 years, t CO2e/ha, before → after)

| Scenario | Released | Leaving the NFZ | Net removal |
|---|---|---|---|
| S1 to S4, S6, S9 (no glass) | unchanged | unchanged | unchanged |
| S5 Base + 20 % glass | 0.2222 → 0.2205 | 0.0280 → 0.0280 | −0.1466 → −0.1466 |
| S7 Optimistic, small pool | 2.1391 → 2.6602 | 0.9488 → 1.4619 | +0.0961 → +0.5322 |
| S8 Upper bound | 4.4036 → 4.7725 | 3.1879 → 3.6242 | +1.9993 → +2.3701 |
| S10 As S7, Andosol acidity | 2.2176 → 1.0603 | 0.2454 → 0.0989 | −0.5018 → −0.6263 |

**What this shows**

- Where gibbsite holds the dissolved Al (S5), the full law gives almost the same rate as the old
  constant: with Al at gibbsite saturation the term (a_H⁺³/a_Al³⁺)^(1/3) does not depend on pH.
- Glass dissolution now depends on **which secondary phase holds the Al**:
  - In S7 and S8 the model lets kaolinite form, which pulls Al³⁺ to 6-15 times below gibbsite
    saturation in the rock layer, so the glass dissolves about twice as fast.
  - In S10, amorphous Al(OH)3 holds about 1000 times more Al³⁺ than gibbsite, so the glass
    dissolves about 10 times more slowly, and 20 % glass becomes a penalty rather than a bonus.
- Equilibrium precipitation of kaolinite within 10 years is an assumption of the model (allophane,
  imogolite or halloysite are more likely in an Andosol). With the full glass law, that assumption
  now matters for the rate as well as for the Si budget. A question for INOVA: which Al phase to
  expect in these soils, and whether soil-water Al and Si can be measured.

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
