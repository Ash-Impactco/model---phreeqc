"""Soil (PHREEQC) -> aquifer cross-section (MODFLOW 6 + PHREEQC) for two scenarios.

    python scripts/run_groundwater.py

Writes outputs/groundwater.csv (one line per scenario) and outputs/groundwater_stream.csv (the
extra alkalinity arriving at the stream over time). Requires the MODFLOW 6 binaries in ./bin
(python scripts/get_modflow.py). All hydrogeology is assumed: data/lab/21_hydrogeology.csv.
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pandas as pd
from erw.model import default_config
from erw.groundwater import run_coupled

SCEN = {
    "S1 Base: soil A, 0/4 mm": ("A_Moreira", "as_is_0-4mm", {}),
    "S7 Optimistic: small pool, 100 um, 20% glass, ff 0.3": ("A_Moreira", "milled_p80_100um",
        {"acidity_mode": "effective", "glass_frac": 0.2, "field_factor": 0.3}),
}
pct = lambda x: "too small to resolve" if x is None else round(x * 100, 1)
rows, series = [], []
for name, (soil, psd, extra) in SCEN.items():
    cfg = default_config(soil, psd); cfg.update(extra)
    r = run_coupled(cfg, label=name.split(" ")[0])
    ok = r["resolvable"]
    b, s, d = r["balance_t_ha"], r["balance_share"], r["degassing"]
    rows.append({"scenario": name,
                 "extra_alk_into_aquifer_10y_t_ha": round(b["input"], 4),
                 "peak_extra_alk_in_drainage_ueq_L": round(max(r["d_alk_input_ueq_by_year"]), 1),
                 "extra_alk_at_stream_end_ueq_L": round(r["d_alk_ueq"][-1], 1),
                 "first_arrival_yr": round(r["first_arrival_yr"], 1) if r["first_arrival_yr"] else "not within the run",
                 "arrived_at_stream_pct": pct(s["stream"]), "still_in_aquifer_pct": pct(s["in_transit"]),
                 "lost_to_calcite_pct": pct(s["calcite"]), "balance_residual_pct": pct(s["residual"]),
                 "share_not_lost_on_the_way_pct": round(r["reach_fraction"] * 100, 1) if ok else "too small to resolve",
                 "carbon_kept_after_stream_degassing_pct": pct(d.get("carbon_kept_per_alk")),
                 "status": "model output; ASSUMED hydrogeology"})
    print(rows[-1])
    for t, v in zip(r["t_years"], r["d_alk_ueq"]):
        series.append({"scenario": name.split(" ")[0], "year": round(t, 4), "extra_alkalinity_at_stream_ueq_L": round(v, 3)})
out = ROOT / "outputs"; out.mkdir(exist_ok=True)
pd.DataFrame(rows).to_csv(out / "groundwater.csv", index=False)
pd.DataFrame(series).to_csv(out / "groundwater_stream.csv", index=False)
print("written:", out / "groundwater.csv", "and", out / "groundwater_stream.csv")
