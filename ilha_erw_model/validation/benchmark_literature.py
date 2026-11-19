"""Benchmark the trial model against published basalt ERW results.

Runs the trial model (erw/geochem.py, unchanged) under conditions close to each published
study and compares the predicted dissolved fraction and bicarbonate export with what was
measured. Also tabulates the earlier PHREEQC scenario model (repository root) from its own
results file data/lab/15_model_results_10yr.csv, since that model was not re-run here.

Run from ilha_erw_model/:   python3 validation/benchmark_literature.py
(needs the soil and rock inputs in inputs/; the saved CSV comes from the 19 Nov 2026 stand-in run)
Writes validation/benchmark_literature.csv. Approximations are listed in the 'notes' column.
"""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from erw.common import Config
from erw import geochem as G
from run_all import load_inputs, missing_inputs

HERE = os.path.dirname(os.path.abspath(__file__))

_c = Config()
load_inputs(_c, log=lambda *a: None)
_miss = missing_inputs(_c)
_need = {"feedstock.oxides_pct", "feedstock.mineral_mode_wt_frac", "feedstock.psd_mass_frac", "feedstock.inorganic_C_pct",
         "feedstock.total_S_pct", "soil.pH_H2O", "soil.bulk_density_g_cm3", "soil.CEC_cmolc_kg",
         "soil.base_saturation_frac", "soil.pH_buffer_cmolc_kg_pH", "soil.acid_load_keq_ha_yr"}
if set(_miss.parameter) & _need:
    sys.exit("Benchmark needs the rock and soil inputs first (see outputs/missing_inputs.csv). "
             "The results in validation/benchmark_literature.csv were produced on 19 Nov 2026 with earlier stand-in values.")


def cfg_with(psd=None, pH=None, depth=None):
    c = Config()
    load_inputs(c, log=lambda *a: None)
    if psd:
        c.set("feedstock.psd_bins_mm", [b for b, _ in psd], "benchmark")
        c.set("feedstock.psd_mass_frac", [w for _, w in psd], "benchmark")
    if pH is not None:
        c.set("soil.pH_H2O", pH, "benchmark")
    if depth is not None:
        c.set("soil.mixing_depth_cm", depth, "benchmark")
    return c


def single(c, rate, months, T, ff):
    _, fs = G.feedstock(c)
    _, eq, _ = G.mineral_table(c, fs)
    return G.simulate(c, [0, rate], months, eq, T_C=T, P={"field_factor": ff})


def cohorts(c, rate, n_apps, months, T, ff):
    """Annual applications superposed (each cohort weathers independently from the same starting pH)."""
    s = single(c, rate, months, T, ff)
    rel = s["release"][0, 1]
    bic = s["bicarb"][0, 1]
    fd = s["frac_dissolved"][0, 1]
    tot_rel = np.zeros(months + 1); tot_bic = np.zeros(months + 1)
    for k in range(n_apps):
        start = 12 * k
        if start > months:
            break
        tot_rel[start:] += rel[: months + 1 - start]
        tot_bic[start:] += bic[: months + 1 - start]
    ages = [months - 12 * k for k in range(n_apps) if 12 * k <= months]
    frac_time_integrated = float(np.mean([fd[a] for a in ages]))
    return tot_rel, tot_bic, frac_time_integrated


OUR_PSD = None  # placeholder 0-4 mm quarry fines from config.json
FINE_P80_35 = [([0.002, 0.010], 0.25), ([0.010, 0.035], 0.55), ([0.035, 0.063], 0.15), ([0.063, 0.125], 0.05)]
P80_267 = [([0.002, 0.063], 0.20), ([0.063, 0.125], 0.20), ([0.125, 0.25], 0.38), ([0.25, 0.5], 0.17), ([0.5, 1.0], 0.05)]

rows = []
for ff in (1.0, 0.1):
    lab = f"trial model, field factor {ff}"
    # 1 Reershemius et al. 2023: mesocosm, 50 t/ha into 12 cm, p80 35 um, soil pH 6.45, 235 days
    s = single(cfg_with(FINE_P80_35, 6.45, 12), 50, 8, 20.0, ff)
    rows.append({"benchmark": "Reershemius 2023 mesocosm, 235 d", "metric": "share of rock cations dissolved",
                 "observed": "15.7 +/- 3.1 %", "model": lab, "predicted": f"{100*s['frac_dissolved'][0,1,8]:.1f} %",
                 "notes": "8 months at an assumed 20 C (growth-room temperature not extracted); PSD built to p80 35 um"})
    # 2 Beerling et al. 2024: field, Illinois, 4 annual x 50 t/ha, p80 267 um, near-neutral soil, ~11 C
    _, _, fti = cohorts(cfg_with(P80_267, 6.2, 10), 50, 4, 48, 11.0, ff)
    rows.append({"benchmark": "Beerling 2024 field, 4 yr", "metric": "Ca+Mg lost from basalt (time-integrated over 4 applications)",
                 "observed": "16 +/- 6 %", "model": lab, "predicted": f"{100*fti:.1f} %",
                 "notes": "assumed start pH 6.2 and 11 C; their metabasalt has only ~11% fast minerals, ours ~22% (olivine+glass), so the model should run high"})
    # 3 Holden et al. 2024: field, tropical Australia, acidic Ultisol pH 5.8, <5 mm basalt, 50 t/ha/yr 2018-2022, HCO3 measured 2022-23
    _, bic, _ = cohorts(cfg_with(None, 5.8, 10), 50, 5, 60, 24.0, ff)
    flux = (bic[60] - bic[48]) * 44.01 / 1e6
    rows.append({"benchmark": "Holden 2024 field, acidic soil, year 5", "metric": "extra bicarbonate export, t CO2/ha/yr",
                 "observed": "0.026 (not significant, p = 0.08)", "model": lab, "predicted": f"{flux:.3f}",
                 "notes": "our 0-4 mm PSD for their <5 mm rock; assumed 24 C; model export is at 10 cm, theirs at 1.25 m (model is an upper bound)"})
    # 4 Dupla et al. 2025: Swiss field trials, realised removal ~0.1 t/ha/yr; compare our 40 and 80 t/ha in year 3
    s = single(cfg_with(), 80, 36, 15.0, ff)
    s40 = single(cfg_with(), 40, 36, 15.0, ff)
    y3 = lambda x: (x["bicarb"][0, 1, 36] - x["bicarb"][0, 1, 24]) * 44.01 / 1e6
    rows.append({"benchmark": "Dupla 2025 field (Switzerland), 3 yr", "metric": "realised CO2 removal, t CO2/ha/yr",
                 "observed": "0.10 +/- 0.03", "model": lab, "predicted": f"{y3(s40):.3f} (40 t/ha), {y3(s):.3f} (80 t/ha)",
                 "notes": "site conditions of Lomba da Maia, not theirs (their rate and soil not retrieved): order-of-magnitude check only"})
    # 5 Vienne et al. 2026: mesocosm 101 d, 0-200 t/ha: no realised inorganic CO2 removal
    s = single(cfg_with(), 50, 3, 20.0, ff)
    rows.append({"benchmark": "Vienne 2026 mesocosm, 101 d", "metric": "realised CO2 removal in first ~3 months",
                 "observed": "none detected; cations held in soil (reducible pool ~10x exchangeable)", "model": lab,
                 "predicted": f"{s['bicarb'][0,1,3]*44.01/1e6:.3f} t CO2/ha at 50 t/ha (any value > 0 is an over-prediction)",
                 "notes": "Lomba da Maia placeholders at 20 C"})

old = [
    ("Reershemius 2023 mesocosm, 235 d", "share of rock cations dissolved", "15.7 +/- 3.1 %",
     "no matching run; finest scenario S4 (p80 100 um) dissolves 2.7 % in 10 years"),
    ("Beerling 2024 field, 4 yr", "Ca+Mg lost from basalt", "16 +/- 6 %",
     "S1 base (0/4 mm) 0.67 % in 10 years; S4 (p80 100 um) 2.7 % in 10 years"),
    ("Holden 2024 field, acidic soil, year 5", "extra bicarbonate export, t CO2/ha/yr", "0.026 (not significant)",
     "S1 0.0003 per yr at 40 t/ha one-off (0.0027 over 10 years)"),
    ("Dupla 2025 field (Switzerland), 3 yr", "realised CO2 removal, t CO2/ha/yr", "0.10 +/- 0.03",
     "S1 0.0003 per yr; S7 optimistic 0.026 per yr"),
    ("Vienne 2026 mesocosm, 101 d", "realised CO2 removal in first ~3 months", "none detected",
     "near zero (consistent)"),
]
for b, m, o, p in old:
    rows.append({"benchmark": b, "metric": m, "observed": o, "model": "earlier PHREEQC model (reported results, not re-run)",
                 "predicted": p, "notes": "from data/lab/15_model_results_10yr.csv on main; field factor 0.1, roughness 5"})

tab = pd.DataFrame(rows).sort_values(["benchmark", "model"])
tab.to_csv(os.path.join(HERE, "benchmark_literature.csv"), index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 70)
print(tab[["benchmark", "model", "observed", "predicted"]].to_string(index=False))
