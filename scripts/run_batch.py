"""Run every scenario in data/scenarios.json and write outputs/results.csv and outputs/scenarios.png.
Used locally (VS Code) and by the GitHub Actions workflow.

    python scripts/run_batch.py                 # inputs from the lab files (data/lab/*.csv)
    python scripts/run_batch.py --source json   # the original inputs (data/inputs.json), Step 0 reference
    python scripts/run_batch.py --update-lab    # also rewrite data/lab/15_model_results_10yr.csv

The batch also writes the data validation report (outputs/validation.md) and the inputs of the base
case with their status (outputs/inputs_S1.csv).
"""
import argparse, json, sys, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from erw.model import default_config, run, full_input
from erw import data_io, validate
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SOURCE = "lab"


def one(item):
    name, spec, source = item
    cfg = default_config(spec.get("soil", "A_Moreira"), spec.get("psd", "as_is_0-4mm"), spec.get("rock"), source=source)
    cfg.update(spec.get("set", {}))
    t = time.time(); res = run(cfg)
    return name, cfg, res, time.time() - t


def data_status(cfg):
    """Words for the status column: what kind of data the forecast rests on."""
    prov = cfg.get("provenance", {})
    kinds = {prov.get(k, {}).get("status", "stand-in") for k in ("oxides", "soil")}
    data = "stand-in data" if kinds <= {"stand-in"} else "data: " + ", ".join(sorted(kinds))
    return f"model forecast ({data}, placeholder emission factors)"


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", choices=["lab", "json"], default="lab")
    ap.add_argument("--update-lab", action="store_true", help="rewrite data/lab/15_model_results_10yr.csv")
    a = ap.parse_args()
    scen = json.loads((ROOT / "data" / "scenarios.json").read_text())
    out = ROOT / "outputs"; (out / "pqi").mkdir(parents=True, exist_ok=True)
    rows, fig, ax = [], *plt.subplots(figsize=(9, 5))
    with ProcessPoolExecutor() as ex:
        for name, cfg, res, dt in ex.map(one, [(n, s, a.source) for n, s in scen.items()]):
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
                         "iso_net_removal_t_ha": round(h["iso"][-1]["net_removal"], 4),
                         "status": data_status(cfg)})
            l, = ax.plot(h["year"], h["released_t"], ":", lw=1.5)
            ax.plot(h["year"], h["export_t"], "-", lw=2.5, color=l.get_color(), label=name)
            safe = name.split(" ")[0]
            (out / "pqi" / f"{safe}_basalt.pqi").write_text(full_input(cfg, True))
            if safe == "S1" and cfg.get("provenance"):
                data_io.provenance_table(cfg).to_csv(out / "inputs_S1.csv", index=False)
    ax.set_xlabel("Years after spreading"); ax.set_ylabel("tCO2e/ha (solid: CO2 removed, dotted: released)")
    ax.legend(fontsize=7, frameon=False); fig.tight_layout(); fig.savefig(out / "scenarios.png", dpi=160)
    df = pd.DataFrame(rows)
    df.to_csv(out / "results.csv", index=False)
    print("written:", out / "results.csv")
    print("written:", validate.write_report(out / "validation.md"))
    if a.update_lab:
        if a.source != "lab":
            sys.exit("--update-lab only makes sense with --source lab")
        target = ROOT / "data" / "lab" / "15_model_results_10yr.csv"
        df.rename(columns={"rock_dissolved_pct": "rock_dissolved_10yr_pct"}).to_csv(target, index=False)
        print("updated:", target)
