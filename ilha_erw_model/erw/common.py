"""Shared helpers: config access, input overrides, plotting style."""
import json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs")
INP = os.path.join(ROOT, "inputs")


class Config:
    """Wraps config.json. cfg.v('soil.pH_H2O') returns the value; cfg.status() lists every parameter."""

    def __init__(self, path=None):
        self.path = path or os.path.join(ROOT, "config.json")
        with open(self.path) as f:
            self.raw = json.load(f)
        self.overrides = {}

    def v(self, key):
        if key in self.overrides:
            return self.overrides[key]
        sec, name = key.split(".", 1)
        return self.raw[sec][name]["value"]

    def set(self, key, value, source="measured"):
        self.overrides[key] = value
        sec, name = key.split(".", 1)
        self.raw[sec][name]["status"] = source

    def status_table(self):
        rows = []
        for sec, d in self.raw.items():
            if sec.startswith("_"):
                continue
            for name, p in d.items():
                if name.startswith("_") and not isinstance(p, dict):
                    continue
                val = self.overrides.get(f"{sec}.{name}", p.get("value"))
                if isinstance(val, (dict, list)):
                    val = json.dumps(val)[:90] + ("..." if len(json.dumps(val)) > 90 else "")
                rows.append({"section": sec, "parameter": name, "value": val,
                             "status": p.get("status", ""), "note": p.get("note", "")})
        return rows


def mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlesize": 10, "legend.frameon": False,
                         "savefig.dpi": 170, "savefig.bbox": "tight"})
    return plt


GREYS = {0: "#ffffff", 20: "#cfcfcf", 40: "#8c8c8c", 80: "#3a3a3a", 160: "#000000"}
LINES = ["#000000", "#555555", "#8a8a8a", "#b5b5b5", "#d0d0d0"]
DASH = ["-", "--", "-.", ":", (0, (1, 1))]


def ensure_dirs():
    os.makedirs(OUT, exist_ok=True)


def save_csv(df, name):
    p = os.path.join(OUT, name)
    df.to_csv(p, index=False)
    return p
