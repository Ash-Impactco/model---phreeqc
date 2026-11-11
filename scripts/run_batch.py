"""Run every scenario in data/scenarios.json and write outputs/results.csv and outputs/scenarios.png.
Used locally (VS Code) and by the GitHub Actions workflow."""
import json, sys, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from erw.model import default_config, run, full_input
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

def one(item):
    name, spec = item
    cfg = default_config(spec.get("soil", "A_Moreira"), spec.get("psd", "as_is_0-4mm"))
    cfg.update(spec.get("set", {}))
    t = time.time(); res = run(cfg)
    return name, cfg, res, time.time() - t

if __name__ == "__main__":
    scen = json.loads((ROOT / "data" / "scenarios.json").read_text())
    out = ROOT / "outputs"; (out / "pqi").mkdir(parents=True, exist_ok=True)
    rows, fig, ax = [], *plt.subplots(figsize=(9, 5))
    with ProcessPoolExecutor() as ex:
        for name, cfg, res, dt in ex.map(one, scen.items()):
            h = res["hist"]
            print(f"{name:60s} released {h['released_t'][-1]:.3f}  removed {h['export_t'][-1]:.3f} tCO2/ha  ({dt:.0f}s)")
            rows.append({"scenario": name, "years": h["year"][-1], "rock_dissolved_pct": round(h["frac_diss"][-1]*100, 2),
                         "released_tCO2e_ha": round(h["released_t"][-1], 4), "co2_removed_t_ha": round(h["export_t"][-1], 4),
                         "on_exchange_tCO2e_ha": round(h["exchange_t"][-1], 4),
                         "topsoil_pH_change": round(h["ph_top_rock"][-1] - h["ph_top_ctrl"][-1], 3),
                         "iso_exported_from_NFZ_20cm_t_ha": round(h["iso"][-1]["exported_from_nfz"], 4),
                         "iso_net_sorption_t_ha": round(h["iso"][-1]["net_sorption"], 4),
                         "iso_stored_t_ha": round(h["iso"][-1]["stored"], 4),
                         "iso_emissions_t_ha": round(h["iso"][-1]["emissions"], 4),
                         "iso_net_removal_t_ha": round(h["iso"][-1]["net_removal"], 4)})
            l, = ax.plot(h["year"], h["released_t"], ":", lw=1.5)
            ax.plot(h["year"], h["export_t"], "-", lw=2.5, color=l.get_color(), label=name)
            safe = name.split(" ")[0]
            (out / "pqi" / f"{safe}_basalt.pqi").write_text(full_input(cfg, True))
    ax.set_xlabel("Years after spreading"); ax.set_ylabel("tCO2e/ha (solid: CO2 removed, dotted: released)")
    ax.legend(fontsize=7, frameon=False); fig.tight_layout(); fig.savefig(out / "scenarios.png", dpi=160)
    pd.DataFrame(rows).to_csv(out / "results.csv", index=False)
    print("written:", out / "results.csv")
