"""Run the model once per soil composite, in parallel, plus a site-mean run.

    python scripts/run_samples.py                 # Lomba da Maia screening (data/lab/09), once filled
    python scripts/run_samples.py --table 08      # the Sanguinho stand-in soils, to see how it works
    python scripts/run_samples.py --years 3 --psd milled_p80_100um

Writes outputs/samples.csv. Rows that lack pH, CEC or exchangeable Ca, Mg, K are listed and skipped.
"""
import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pandas as pd

from erw.samples import run_many, sample_configs, summary_row

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--table", default="09", help="soil table: 09 (Lomba da Maia) or 08 (Sanguinho)")
    ap.add_argument("--psd", default="as_is_0-4mm")
    ap.add_argument("--rock", default="MB-BRITA2")
    ap.add_argument("--years", type=int, default=None)
    a = ap.parse_args()

    configs, skipped = sample_configs(a.table, a.psd, a.rock, {"years": a.years} if a.years else None)
    for sid, missing in skipped:
        print(f"skipped {sid}: no value yet for {', '.join(missing)}")
    if not configs:
        print(f"\nNo row in table {a.table} has pH, CEC and exchangeable Ca, Mg and K yet, so there is nothing "
              "to run. Fill the table and run this again.")
        sys.exit(0)
    t = time.time()
    rows = [summary_row(*r) for r in run_many(configs)]
    out = ROOT / "outputs"
    out.mkdir(exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "samples.csv", index=False)
    print(df.to_string(index=False))
    print(f"\n{len(rows)} runs in {time.time() - t:.0f} s. Written: {out / 'samples.csv'}")
