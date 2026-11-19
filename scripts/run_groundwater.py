"""Soil (PHREEQC) -> aquifer (MODFLOW 6 + PHREEQC) for two scenarios; writes outputs/groundwater.csv.
Requires MODFLOW 6 binaries in ./bin (python scripts/get_modflow.py)."""
import sys, json
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
rows = []
for name, (soil, psd, extra) in SCEN.items():
    cfg = default_config(soil, psd); cfg.update(extra)
    r = run_coupled(cfg, {"years": cfg["years"]}, label=name.split(" ")[0])
    rows.append({"scenario": name, "soil_CO2_removed_10y_t_ha": round(r["soil_export_t_ha"], 4),
                 "extra_alk_in_drainage_ueq_L": round(r["d_alk_input_ueq"], 2),
                 "extra_alk_at_stream_end_ueq_L": round(r["d_alk_ueq"][-1], 2),
                 "share_reaching_stream_pct": round(r["reach_fraction"] * 100, 1) if r["d_alk_input_ueq"] >= 10 else "too small",
                 "status": "model output; ASSUMED hydrogeology"})
    print(rows[-1])
out = ROOT / "outputs"; out.mkdir(exist_ok=True)
pd.DataFrame(rows).to_csv(out / "groundwater.csv", index=False)
