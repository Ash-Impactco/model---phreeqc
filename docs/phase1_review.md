# Phase 1: run check and sceptical review

All numbers below come from placeholder soil and rock values. None is a site result.

## 1. Does the model reproduce?

Re-ran on 19 November 2026 (Python 3.9.6, numpy 2.0.2, pandas 2.3.3, scipy 1.13.1, statsmodels 0.14.6, in `.venv`). `--quick` ran in about 1 minute, the full run in 1 min 43 s.

Every output table is byte-identical to the shipped version **except** `power_by_design.csv`, which differs by up to 0.06 in individual cells. That is simulation noise from different library versions (300 simulated trials per cell gives a standard error of 0.02 to 0.03 on a power value), not a change in the model.

| Headline | Shipped claim | Re-run | OK? |
|---|---|---|---|
| Layout | 16 x 0.64 ha = 10.24 ha, 4 blocks, U3/U4/U5, 390 to 445 m, slopes 7 to 17% | same; slopes 6.8 to 16.8% | yes |
| Basalt | 358 t dry, 394 t +10%, 428 t wet, 18 loads; 12 plots 296 t | 358.4, 394.2, 428.6, 17.1 so 18 loads, 295.7 | yes |
| Rock potential | about 320 kg CO2/t | 320.1 | yes |
| Dissolved in 30 months | about 7%, range 3 to 15% | 6.7% at 40 t/ha; field factor 0.2 to 3 gives 2.6 to 11.5%; Monte Carlo P10 to P90 3.4 to 12.9% | yes, but "3 to 13%" is the honest range |
| Whole trial net at 30 months | -1.3 t (-2.1 to +0.2) | -1.31 (-2.10 to +0.20) | yes |
| Positive later at 40 and 80 t/ha | | 10 years: 20 t/ha -0.02, 40 t/ha +0.02, 80 t/ha +0.25 t/ha | yes |
| Biggest drivers | field factor, starting pH, fines | Spearman 0.83, 0.35, 0.31 | yes |
| Lysimeter alkalinity | +0.07 and +0.20 meq/L at 40, 80 t/ha | 0.075, 0.202 | yes |
| Soil mass balance smallest detectable | 0.17 to 0.35 vs 0.07 | 0.174 to 0.355 (dose method, 80 t/ha, 4 plots) | yes |
| Streams | at most 0.01 meq/L | 0.0014 (W1) to 0.0104 (W6) | yes |
| Power | 16 plots 0.90 to 0.97; 12 plots 0.75 to 0.85 | 16 plots 0.94 to 0.95; 12 plots 0.78 to 0.88 | yes (within noise) |
| Metals | Ni, Cr, Cu, Zn background already above limits; +10% Ni, +20% Cr at 80 t/ha | +9.6% Ni, +20.5% Cr | yes |

One plumbing bug: `erw/report.py` calls `python3` by name for the PDF step, so the PDF fails inside a virtual environment. Workaround used: put `.venv/bin` first on PATH. Proposed fix (one line, no science): use `sys.executable`.

## 2. Review as a sceptical geochemist

Graded: **A** = would change a conclusion or a design decision; **B** = changes numbers noticeably but not the conclusion; **C** = tidy-up.

### Checks that passed
- **Units in the CO2 potential.** wt% x 10 / molar mass x charge = eq/kg; x 44.01 = kg CO2/t at 1 mol CO2 per equivalent. CaO 10% gives 3.567 eq/kg, which is right. Rock carbon (0.05% C gives 1.83 kg CO2/t) and sulfur (0.02% S, 2 eq per mol S, gives 0.55 kg CO2/t) are right.
- **Shrinking core.** rate (mol/m2/s) x molar volume (m3/mol) = surface retreat (m/s); remaining mass = ((r0 - retreat)/r0)^3. Correct for spheres. The Palandri and Kharaka constants for forsterite, augite and labradorite match the published table.
- **pH buffer step.** dpH = (acid neutralised - acid added) / (soil mass x buffer capacity / 100). The /100 turns cmolc into eq. Check: 1.5 keq/ha/yr over 30 months into 750 t of soil at 5 cmolc/kg/pH = 0.10 pH drop on the controls, which is what the model gives.
- **Emissions.** 26 km x 0.12 + 0.4 + 0.8 = 4.32 kg CO2e/t. Arithmetic right (but see B2).
- **Power simulation noise.** Noise = plot-to-plot SD combined with sampling and lab error, doubled for before-minus-after; block effects drawn and then removed by the block term. Holm step-down is implemented correctly.

### A1. Power is calculated as if we already knew the field factor
`stats.power` uses the central chemistry effect only. With field factor 0.3 the pH change at 40 t/ha is 0.13 instead of 0.34, exchangeable Ca + Mg 0.81 instead of 1.54 cmolc/kg. Power at 16 plots would then fall well below 80% for pH. The "16 plots gives 0.90 to 0.97" statement is conditional on central kinetics and should say so. **Fix:** add a power run at field factor 0.3 (and 3) so the team sees the downside. I can do this as a precompute for the simulator without changing the model.

### A2. The finest sieve fraction drives the whole early result, and its shape is assumed
71% of all cations released in 30 months come from grains below 63 micrometres, and 50% from grains below 10 micrometres. The model spreads the 2 to 63 micrometre mass evenly on a log scale, so a sixth of it sits at 2 to 3.5 micrometres and dissolves within months. If the real fines were mostly 20 to 63 micrometres, dissolution at 40 t/ha falls from 6.7% to 4.6% and gross removal at 30 months falls by 42%. Grains below 2 micrometres are ignored entirely. **Fix:** laser diffraction on the below-63 fraction of the Tecnovia fines (not just sieving), plus BET. This is the cheapest uncertainty reduction available.

### A3. The soil mass-balance test assumes all dissolved cations leave, but the soil block says half stay
In the chemistry, 51% of released cations at 40 t/ha are still sitting on exchange sites in the top 10 cm at 30 months. Fusion XRF measures total Ca + Mg, which includes those exchangeable cations. So a mass balance on total Ca + Mg would only "see" about half of the dissolution (an apparent 3.3% instead of 6.7%). `mrv.tau_detectability` assumes all dissolved cations are gone, so it is optimistic. The conclusion (cannot detect dissolution in 30 months) gets stronger, not weaker. **Fix:** when we do the mass balance, subtract exchangeable Ca + Mg (measured anyway in the composite) from total, and say so in the protocol. Worth discussing with Juan.

### B1. The bicarbonate share rule is a heuristic, but a reasonable one
`1/(1 + 10^(pKa1 - pH))` is the share of dissolved inorganic carbon present as bicarbonate rather than dissolved CO2. Using it as "the share of released cations that leave as bicarbonate" is not a law of chemistry; the true limit is set by soil-air CO2 and the pH of drainage water. I cross-checked: at 40 t/ha the model exports about 900 eq/ha/yr of bicarbonate at pH 5.74; water in equilibrium with 1 to 3% soil CO2 at that pH and 1.2 m drainage could carry 1,000 to 3,000 eq/ha/yr. Same order, so the rule is not misleading. But it makes the carbon number very pH-sensitive (share is 4% at pH 5.0, 10% at 5.4, 22% at 5.8, 41% at 6.2). **Fix:** none needed now; measure soil-air CO2 and lysimeter DIC so we can replace the rule with data.

### B2. Haul emissions are charged per dry tonne, but we truck wet tonnes
At 8% moisture we haul 1.087 t per dry tonne. Haul emissions rise from 3.12 to 3.39 kg/t; total from 4.32 to 4.59 kg CO2e/t (+6%). Net at 40 t/ha, 30 months, moves from -0.160 to -0.171 t/ha. Also the 10% contingency is hauled but not spread; its emissions belong in the trial footprint. **Fix:** small change in `accounting`; needs your approval.

### B3. Alkalinity is counted as removed when it leaves the top 10 cm
The model has one soil layer (0 to 10 cm). Bicarbonate leaving it is counted as gross removal and is assumed to arrive unchanged at the 30 cm lysimeters. In andic subsoils some of that alkalinity can be consumed by acidity below 10 cm. So gross removal and the lysimeter signal are both upper bounds. A registry would want export below the rooting zone. **Fix:** say "upper bound" in the application; later add a second soil layer or calibrate on lysimeter data.

### B4. The Monte Carlo is narrower than the real uncertainty
It varies 9 inputs, but not rock chemistry, the fines shape (A2), roughness, glass kinetics, mineral mode, or the noise placeholders. P10 to P90 is therefore an under-statement. **Fix:** say so on every uncertainty chart; add fines shape and roughness to the Monte Carlo later.

### B5. Fertiliser acid load lands entirely in the top 10 cm
All 1.5 keq/ha/yr is applied to 750 t of soil, so control pH falls 0.4 units over 10 years. If the acid is spread over the root zone the decline is smaller and the starting pH effect on carbon is milder. **Fix:** split the acid load by depth once Luis's N records are in.

### C1. Mineral mode does not match the XRF element by element
The crystal phases alone imply 7.7% MgO (XRF 6.9%) and 10.4% CaO (XRF 10.0%), and only 1.8% Na2O (XRF 3.4%). The model only checks the total charge, so glass ends up with whatever is left over (3.1 eq/kg, whereas real basaltic glass is closer to the bulk rock). Placeholder issue; the real XRD and microprobe will fix it. Effect on dissolution is small.

### C2. Glass rate may double-count roughness
Gislason and Oelkers (2003) normalised basaltic glass rates to geometric surface. If the config's glass constants come from there, multiplying by roughness 7 counts surface roughness twice. Removing it changes 30-month dissolution from 6.7% to 6.2%. Small, but ask Juan where the glass constants came from.

### C3. Small inconsistencies
- Controls lose pH under the acid load but their base saturation does not fall. Cosmetic.
- Between-plot SD for total Ca + Mg (0.05 eq/kg) is hard-coded in `stats.noise_sd`, not in config.json, so it is not flagged as a placeholder.
- Mesocosm control column shows "5.7% dissolved" at 0 t/ha (it is the per-kg rate with no rock). Harmless.

## 3. What I have not changed
Nothing in the model. All of A1 to C3 are proposals. A1 needs no model change (extra precompute). B2 and the report.py path fix are one-line code changes that need your OK.
