"""
Trace metals added to the topsoil with the rock.

Inputs: the rock's metal contents in data/lab/06_rock_trace_metals.csv (mg/kg) and the soil limits and
baseline in data/lab/20_trace_metal_limits.csv. Both are blank until the lab results and the limits
arrive, and nothing is filled in: with blank inputs the functions return empty values and say why.

Two numbers per element:
  total added      everything in the rock, dissolved or not, mixed into the topsoil. This is what a
                   total (aqua regia) soil analysis would see, and what soil limits refer to.
  released         for Ni, Cr and Co only: the part set free as olivine and pyroxene dissolve, taken
                   as proportional to the share of those minerals dissolved in the model.

Assumptions, all on the cautious side for accumulation:
  - all of the rock's Ni, Cr and Co is taken to sit in olivine and pyroxene. Chromium held in chromite
    or spinel would not be released; the deportment test in table 06 will tell.
  - released metals stay in the topsoil (no leaching, no plant uptake).
"""
from __future__ import annotations

import pandas as pd

from . import data_io
from .data_io import num, parse_value, status_class

HOSTED_IN_OLIVINE_PYROXENE = ("Ni", "Cr", "Co")


def rock_metals(sample_id: str, lab_dir=data_io.LAB) -> dict[str, dict]:
    """{element: {mg_kg, below_detection, status}} for one rock sample. Blank rows are left out."""
    df = data_io.read_table("06", lab_dir)
    out = {}
    for _, r in df[df["sample_id"] == sample_id].iterrows():
        p = parse_value(r["mg_kg"])
        if p.is_number and len(r["element"]) <= 2:
            out[r["element"]] = dict(mg_kg=p.value, below_detection=p.below_detection,
                                     status=data_io.effective_status(sample_id, r["status"]))
    return out


def limits(lab_dir=data_io.LAB) -> dict[str, dict]:
    """{element: {limit, baseline, source}} from table 20. Values are None where the table is blank."""
    try:
        df = data_io.read_table("20", lab_dir)
    except KeyError:
        return {}
    return {r["element"]: dict(limit=num(r["soil_limit_mg_kg"]), baseline=num(r["baseline_soil_mg_kg"]),
                               source=r.get("source", ""), status=status_class(r["status"]))
            for _, r in df.iterrows()}


def topsoil_and_rock_t_ha(cfg) -> tuple[float, float]:
    """Dry mass of the topsoil the rock is mixed into, and of the rock, in t/ha."""
    soil = cfg["mix_depth_cm"] / 100.0 * cfg["bulk_density_kg_m3"] * 10.0
    return soil, float(cfg["application_t_ha"])


def accumulation(cfg, hist=None, lab_dir=data_io.LAB) -> tuple[pd.DataFrame, str]:
    """Table of metal additions to the topsoil for the rock of this run, and a note on what is missing.

    Columns: element, rock_mg_kg, added_total_mg_kg (rise in the topsoil's total content), released_mg_kg
    (Ni, Cr, Co at the end of the run; needs `hist` from a model run), baseline_mg_kg, after_mg_kg,
    limit_mg_kg, pct_of_limit, status.
    """
    rock = cfg.get("rock_key", "")
    metals = rock_metals(rock, lab_dir)
    lim = limits(lab_dir)
    if not metals:
        return (pd.DataFrame(), f"No trace-metal analysis for {rock} in data/lab/06_rock_trace_metals.csv yet, "
                                "so metal additions cannot be calculated. Nothing has been assumed.")
    soil_t, rock_t = topsoil_and_rock_t_ha(cfg)
    host = (hist or {}).get("host_frac_diss") or []
    rows = []
    for el, m in metals.items():
        base = lim.get(el, {}).get("baseline")
        limit = lim.get(el, {}).get("limit")
        mixed = (soil_t * (base or 0.0) + rock_t * m["mg_kg"]) / (soil_t + rock_t)      # content of the mixture
        added = mixed - (base or 0.0)
        released = None
        if el in HOSTED_IN_OLIVINE_PYROXENE and host:
            released = rock_t * m["mg_kg"] * host[-1] / (soil_t + rock_t)
        after = mixed if base is not None else None
        rows.append({"element": el, "rock_mg_kg": m["mg_kg"], "added_total_mg_kg": added,
                     "released_mg_kg": released, "baseline_mg_kg": base, "after_mg_kg": after, "limit_mg_kg": limit,
                     "pct_of_limit": (after / limit * 100) if (after is not None and limit) else None,
                     "status": ("upper bound: at the detection limit" if m["below_detection"] else m["status"])})
    missing = []
    if any(r["baseline_mg_kg"] is None for r in rows):
        missing.append("soil baseline")
    if any(r["limit_mg_kg"] is None for r in rows):
        missing.append("soil limits")
    note = ("Total added counts all the metal in the rock, dissolved or not. "
            + (f"Still blank in data/lab/20_trace_metal_limits.csv: {' and '.join(missing)}, so no comparison "
               "with a limit can be made for those elements." if missing else ""))
    return pd.DataFrame(rows), note.strip()


def release_series(cfg, hist, element: str, lab_dir=data_io.LAB) -> list[float]:
    """Released Ni, Cr or Co in the topsoil (mg/kg) at each year end."""
    metals = rock_metals(cfg.get("rock_key", ""), lab_dir)
    if element not in metals or element not in HOSTED_IN_OLIVINE_PYROXENE:
        return []
    soil_t, rock_t = topsoil_and_rock_t_ha(cfg)
    return [rock_t * metals[element]["mg_kg"] * f / (soil_t + rock_t) for f in hist["host_frac_diss"]]
