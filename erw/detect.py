"""
Measurement-plan tools for soil-based quantification (Isometric EWA v1.3, section 10.4).

1. Minimum detectable change: the smallest change in a soil property that the planned number of
   composite samples can tell from zero, given how much the field varies.
2. Tracer mass balance: how much of the feedstock's Ca and Mg has left the sampled layer, from
   baseline and post-application soil analyses and an immobile tracer (Ti or Zr). This is the
   calculation behind soil-based quantification; here it is ready for real data and tested on
   made-up data. It never uses model output.

Both need field data that do not exist yet (the coefficient of variation at Lomba da Maia, total Ca,
Mg and Ti in the soil). Until then the inputs are choices made on the page, labelled as such.
"""
from __future__ import annotations

from statistics import NormalDist

import numpy as np

MOLAR = dict(Ca=40.078, Mg=24.305, Ti=47.867, Zr=91.224, CaO=56.077, MgO=40.304, TiO2=79.866)
OXIDE_TO_ELEMENT = dict(CaO=MOLAR["Ca"] / MOLAR["CaO"], MgO=MOLAR["Mg"] / MOLAR["MgO"],
                        TiO2=MOLAR["Ti"] / MOLAR["TiO2"])


# ---------------------------------------------------------------- 1. minimum detectable change
def _t(prob: float, df: float) -> float:
    """Student t quantile from the normal quantile (Cornish-Fisher, two terms). Good to about 1 %
    for df >= 4, which is enough for sampling design."""
    z = NormalDist().inv_cdf(prob)
    return z + (z ** 3 + z) / (4 * df) + (5 * z ** 5 + 16 * z ** 3 + 3 * z) / (96 * df ** 2)


def minimum_detectable_change(cv_pct: float, n: int, alpha: float = 0.05, power: float = 0.8,
                              paired: bool = False, correlation: float = 0.0, one_sided: bool = True) -> float:
    """Smallest relative change (in % of the mean) detectable with n composites per group.

    Two groups are compared: baseline against post-application, or treated against control.
        unpaired:  MDC = (t_alpha + t_beta) * CV * sqrt(2 / n),              df = 2(n - 1)
        paired:    the same points sampled again; MDC shrinks by sqrt(1 - correlation), df = n - 1
    cv_pct is the coefficient of variation between composites, in per cent.
    """
    if n < 2:
        return float("inf")
    df = (n - 1) if paired else 2 * (n - 1)
    ta = _t(1 - alpha if one_sided else 1 - alpha / 2, df)
    tb = _t(power, df)
    spread = (2 * (1 - correlation if paired else 1.0) / n) ** 0.5
    return (ta + tb) * cv_pct * spread


def composites_needed(cv_pct: float, change_pct: float, alpha: float = 0.05, power: float = 0.8,
                      paired: bool = False, correlation: float = 0.0, n_max: int = 5000):
    """Smallest number of composites per group that can detect a change of change_pct. None if
    more than n_max would be needed."""
    if change_pct <= 0:
        return None
    for n in range(2, n_max + 1):
        if minimum_detectable_change(cv_pct, n, alpha, power, paired, correlation) <= change_pct:
            return n
    return None


def exchangeable_signal(hist, depth_cells: int = 2) -> dict:
    """Modelled change in exchangeable Ca+Mg in the topsoil, rock against control, per year (%)."""
    rel, absolute = [], []
    for p in hist["profiles"]:
        r = float(np.mean(p["excamg_rock"][:depth_cells]))
        c = float(np.mean(p["excamg_ctrl"][:depth_cells]))
        absolute.append(r - c)
        rel.append((r - c) / c * 100 if c else 0.0)
    return dict(year=list(hist["year"]), relative_pct=rel, absolute_cmol_kg=absolute)


def feedstock_camg_ti(oxides: dict) -> dict:
    """Ca, Mg and Ti of the rock in mg/kg from the oxide analysis (wt %)."""
    return {el: oxides[ox] * OXIDE_TO_ELEMENT[ox] * 1e4 for el, ox in (("Ca", "CaO"), ("Mg", "MgO"), ("Ti", "TiO2"))}


def ratio_signal(hist, cfg, soil_camg_mg_kg=None, soil_ti_mg_kg=None) -> dict:
    """Modelled change in total (Ca+Mg)/Ti of the topsoil as the rock loses cations, per year (%).

    The rock keeps its Ti and loses Ca and Mg, so the ratio falls. The size of the fall depends on
    how much Ca, Mg and Ti the soil itself holds, which is not measured yet. Without those two
    numbers the function returns the change in the rock alone, which is the largest it can be.
    The loss used is what has left the Near-Field Zone, because a total analysis still counts the
    cations sitting on exchange sites.
    """
    f = feedstock_camg_ti(cfg["oxides"])
    soil_t = cfg["mix_depth_cm"] / 100.0 * cfg["bulk_density_kg_m3"] * 10.0       # t/ha
    rock_t = float(cfg["application_t_ha"])
    ox = cfg["oxides"]
    # charge of all the base cations in the rock spread on one hectare (eq/ha). The exported charge is
    # taken to be Ca and Mg in the same proportion as in the rock (congruent loss).
    eq_rock = rock_t * 1e6 * (2 * ox["CaO"] / 56.077 + 2 * ox["MgO"] / 40.304 + 2 * ox["Na2O"] / 61.979
                              + 2 * ox["K2O"] / 94.196) / 100.0
    lost = []
    for i in hist["iso"]:
        exported_eq_ha = i["exported_from_nfz"] * 1e6 / 44.01                    # t CO2e/ha -> eq/ha
        lost.append(min(1.0, max(0.0, exported_eq_ha / eq_rock)) if eq_rock else 0.0)
    have_soil = soil_camg_mg_kg is not None and soil_ti_mg_kg is not None and soil_ti_mg_kg > 0
    out = []
    for share in lost:
        if have_soil:
            camg0 = soil_t * soil_camg_mg_kg + rock_t * (f["Ca"] + f["Mg"])
            ti = soil_t * soil_ti_mg_kg + rock_t * f["Ti"]
            camg1 = camg0 - share * rock_t * (f["Ca"] + f["Mg"])
            out.append((camg1 / ti) / (camg0 / ti) * 100 - 100)
        else:
            out.append(-share * 100)
    return dict(year=list(hist["year"]), relative_pct=out, share_of_rock_camg_lost=lost,
                basis="soil and rock mixture" if have_soil else "rock alone (upper bound: soil Ca, Mg and Ti not given)")


# ---------------------------------------------------------------- 2. tracer mass balance
def tracer_mass_balance(baseline: dict, post: dict, feedstock: dict, tracer: str = "Ti", cations=("Ca", "Mg"),
                        soil_mass_t_ha: float | None = None) -> dict:
    """Feedstock cations lost from the sampled layer, from soil analyses and an immobile tracer.

    All inputs are total concentrations in the same unit (for example mg/kg): `baseline` is the soil
    before spreading, `post` the soil after weathering, `feedstock` the rock.

      x          mass fraction of feedstock in the post sample, from the tracer:
                 x = (T_post - T_base) / (T_rock - T_base)
      C_mix      cation concentration expected if nothing had dissolved: C_base + x (C_rock - C_base)
      loss       share of the feedstock's cation that has gone: (C_mix - C_post) / (x C_rock)

    Assumptions: the tracer does not move; the soil's own cations have not changed; the loss of rock
    mass on dissolving is small. A total analysis includes cations on exchange sites, so `loss` is
    what has left the sampled layer, which is the quantity Isometric counts at the Near-Field Zone.
    """
    tb, tp, tf = baseline[tracer], post[tracer], feedstock[tracer]
    if tf <= tb:
        raise ValueError(f"the feedstock must hold more {tracer} than the soil for the tracer to work")
    x = (tp - tb) / (tf - tb)
    out = dict(tracer=tracer, feedstock_fraction=x, loss={}, expected={}, eq_lost_per_kg_rock=0.0)
    if x <= 0:
        out["note"] = f"no feedstock detected: {tracer} in the post sample is not above the baseline"
        return out
    for c in cations:
        mix = baseline[c] + x * (feedstock[c] - baseline[c])
        out["expected"][c] = mix
        out["loss"][c] = (mix - post[c]) / (x * feedstock[c])
        out["eq_lost_per_kg_rock"] += out["loss"][c] * feedstock[c] / 1e6 / MOLAR[c] * 2 * 1e3     # mg/kg -> eq per kg rock
    if soil_mass_t_ha:
        rock_t_ha = x * soil_mass_t_ha / (1 - x)                               # rock found in the layer
        out["apparent_application_t_ha"] = rock_t_ha
        out["co2e_t_ha"] = out["eq_lost_per_kg_rock"] * rock_t_ha * 1e3 * 44.01 / 1e6
    return out


def bootstrap_mass_balance(baseline_samples: list[dict], post_samples: list[dict], feedstock: dict, tracer="Ti",
                           cations=("Ca", "Mg"), soil_mass_t_ha: float | None = None, n_boot: int = 2000,
                           seed: int = 1) -> dict:
    """The mass balance on the means of the composites, with a bootstrap over the composites for the
    spread. Returns the central estimate and P10, P30, P40, P50, P90 of the cation loss (and of the
    CO2e per hectare when the soil mass is given). P30 is the level Isometric issues at (P40 with a
    validation approach)."""
    keys = [tracer, *cations]
    mean = lambda rows: {k: float(np.mean([r[k] for r in rows])) for k in keys}
    central = tracer_mass_balance(mean(baseline_samples), mean(post_samples), feedstock, tracer, cations, soil_mass_t_ha)
    rng = np.random.default_rng(seed)
    vals = {"loss_" + c: [] for c in cations}
    vals["co2e_t_ha"] = []
    for _ in range(n_boot):
        b = [baseline_samples[i] for i in rng.integers(0, len(baseline_samples), len(baseline_samples))]
        p = [post_samples[i] for i in rng.integers(0, len(post_samples), len(post_samples))]
        try:
            r = tracer_mass_balance(mean(b), mean(p), feedstock, tracer, cations, soil_mass_t_ha)
        except ValueError:
            continue
        if not r["loss"]:
            continue
        for c in cations:
            vals["loss_" + c].append(r["loss"][c])
        if "co2e_t_ha" in r:
            vals["co2e_t_ha"].append(r["co2e_t_ha"])
    pct = {k: ({f"P{q}": float(np.percentile(v, q)) for q in (10, 30, 40, 50, 90)} if v else {})
           for k, v in vals.items()}
    return dict(central=central, percentiles=pct, n_boot=n_boot, n_valid=len(vals["loss_" + cations[0]]),
                n_baseline=len(baseline_samples), n_post=len(post_samples))
