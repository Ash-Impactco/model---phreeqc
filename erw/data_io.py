"""
Load the Ilha Carbon lab data pack (data/lab/*.csv) and turn it into model inputs.

Every value keeps its provenance. The `status` text of each row is kept as written, and is also
sorted into one of a few classes (see `status_class`):

    measured    a real lab or supplier value on the project's own rock or soil
    stand-in    a real value, but from another rock or another site (Marques Britas, Sanguinho)
    assumed     a modelling choice or literature value
    to measure  an empty row waiting for results
    calculated / plan / model output / unspecified   (bookkeeping classes)

Nothing in this module invents a number. Where the lab files have no value, the loader falls back
to data/inputs.json and says so in the provenance record.

The module does not import the model, so the model can import it.
"""
from __future__ import annotations

import copy
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
LAB = DATA / "lab"
INPUTS_JSON = DATA / "inputs.json"

MODEL_OXIDES = ["SiO2", "Al2O3", "Fe2O3", "CaO", "MgO", "TiO2", "Na2O", "K2O"]
MW_FEO, MW_FE2O3 = 71.844, 159.69

# Sample-id prefixes of the project's own materials. Anything else that is "measured" is a stand-in.
PROJECT_PREFIXES = ("TV-", "LDM-")

# Old names used in data/inputs.json and data/scenarios.json, kept so existing scenarios still run.
SOIL_ALIASES = {"A_Moreira": "SG-04", "B_TerraRiba": "SG-05"}
PSD_ALIASES = {"milled_p80_1mm": "MILL-P80-1MM", "milled_p80_100um": "MILL-P80-100UM"}
AS_IS_PSD = "as_is_0-4mm"
DEFAULT_ROCK = "MB-BRITA2"
DEFAULT_SOIL = "SG-04"


# ====================================================================== cell parsing
@dataclass(frozen=True)
class Parsed:
    """One table cell. `value` is the best single number, or None if the cell holds no number."""
    value: float | None = None
    low: float | None = None
    high: float | None = None
    qualifier: str = ""      # "", "<", ">", "~" or "range"
    text: str = ""           # the cell as written
    note: str = ""           # any words that came with the number

    @property
    def is_number(self) -> bool:
        return self.value is not None

    @property
    def below_detection(self) -> bool:
        return self.qualifier == "<"


_NUM = r"[-+]?\d+(?:[.,]\d+)?(?:[eE][-+]?\d+)?"


def _to_float(token: str) -> float:
    token = token.strip()
    if "," in token and "." not in token:
        head, tail = token.split(",", 1)
        # "4,000" is read as four thousand; "5,8" and "0,05" as decimal commas (Portuguese reports)
        token = head + tail if len(tail) == 3 and head.lstrip("+-") not in ("", "0") else head + "." + tail
    return float(token)


def parse_value(cell) -> Parsed:
    """Read a cell that may hold `12.3`, `<0.05`, `2850-3000`, `~11`, `0.09 (a note)` or nothing.

    - `<0.05`  -> value 0.05, low 0, high 0.05, qualifier "<" (a detection limit, not a measurement)
    - `2850-3000` or `-2.5 to -1.5` -> value = midpoint, low and high kept, qualifier "range"
    - `~11`    -> value 11, qualifier "~"
    - blank    -> value None
    - words    -> value None, note = the words
    """
    if cell is None or (isinstance(cell, float) and math.isnan(cell)):
        return Parsed()
    if isinstance(cell, (int, float)):
        v = float(cell)
        return Parsed(v, v, v, "", str(cell))
    text = str(cell).strip()
    if not text:
        return Parsed(text="")
    t = text.replace("−", "-").replace("–", "-").replace("≈", "~").replace("≤", "<").replace("≥", ">")

    m = re.fullmatch(rf"({_NUM})\s*(?:to|-)\s*({_NUM})", t)
    if m:
        lo, hi = _to_float(m.group(1)), _to_float(m.group(2))
        if lo <= hi:
            return Parsed((lo + hi) / 2, lo, hi, "range", text)
        return Parsed(text=text, note=text)          # e.g. a date such as 2026-11

    m = re.fullmatch(rf"([<>~])\s*=?\s*({_NUM})\s*(.*)", t)
    if not m:
        m2 = re.fullmatch(rf"(?:about|approx\.?|ca\.?|circa)\s*({_NUM})\s*(.*)", t, flags=re.I)
        if m2:
            v = _to_float(m2.group(1))
            return Parsed(v, v, v, "~", text, m2.group(2).strip())
    if m:
        q, v, note = m.group(1), _to_float(m.group(2)), m.group(3).strip()
        if q == "<":
            return Parsed(v, 0.0, v, "<", text, note)
        if q == ">":
            return Parsed(v, v, None, ">", text, note)
        return Parsed(v, v, v, "~", text, note)

    m = re.fullmatch(rf"({_NUM})\s*(\(.*\)|[A-Za-zµ%°].*)?", t)
    if m:
        v = _to_float(m.group(1))
        return Parsed(v, v, v, "", text, (m.group(2) or "").strip())
    return Parsed(text=text, note=text)


def num(cell, default=None):
    """The number in a cell, or `default` if there is none."""
    p = parse_value(cell)
    return p.value if p.is_number else default


def parse_latlon(cell):
    """`"37.7909, -25.5916"` -> (37.7909, -25.5916); anything else -> None."""
    m = re.fullmatch(rf"\s*({_NUM})\s*[,;]\s*({_NUM})\s*", str(cell or ""))
    if not m:
        return None
    return float(m.group(1)), float(m.group(2))


# ====================================================================== status handling
def status_class(status: str) -> str:
    """Sort a free-text status into one provenance class."""
    s = str(status or "").strip().lower()
    if not s:
        return "unspecified"
    if s.startswith("to ") or s in ("optional", "unknown", "planned") or "to measure" in s:
        return "to measure"
    if "forecast" in s or s.startswith("modelled") or "model output" in s:
        return "model output"
    if "stand-in" in s or "stand in" in s:
        return "stand-in"
    if "measured" in s and "to measure" not in s:
        return "measured"
    if any(w in s for w in ("assum", "placeholder", "estimate", "approximate", "model choice", "literature")):
        return "assumed"
    if "calculated" in s:
        return "calculated"
    if any(w in s for w in ("plan", "protocol", "isometric")):
        return "plan"
    return "unspecified"


def is_project_sample(sample_id: str) -> bool:
    return str(sample_id or "").upper().startswith(PROJECT_PREFIXES)


def effective_status(sample_id: str, status: str) -> str:
    """Provenance class as the model must report it: a measurement on a rock or soil that is not
    the project's own (not TECHNOVIA fines, not Lomba da Maia) is a stand-in, never 'measured'."""
    c = status_class(status)
    if c == "measured" and not is_project_sample(sample_id):
        return "stand-in"
    return c


# ====================================================================== tables
def table_files(lab_dir: Path | str = LAB) -> dict[str, Path]:
    """{'01_rock_major_oxides': Path, ...} for every CSV in the data pack, in file order."""
    return {p.stem: p for p in sorted(Path(lab_dir).glob("*.csv"))}


def read_table(key: str, lab_dir: Path | str = LAB) -> pd.DataFrame:
    """Read one table as text (so `<0.05` and blanks survive). `key` may be the number ("03"),
    the name ("rock_grain_size") or the full stem. Blank cells are empty strings."""
    files = table_files(lab_dir)
    hits = [s for s in files if s == key or s.startswith(f"{key}_") or s.split("_", 1)[-1] == key]
    if len(hits) != 1:
        raise KeyError(f"lab table {key!r}: {len(hits)} matches in {lab_dir}")
    path = files[hits[0]]
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = [r for r in csv.reader(fh) if any(c.strip() for c in r)]
    header = [c.strip() for c in rows[0]]
    body, ragged = [], []
    for i, r in enumerate(rows[1:], start=2):
        if len(r) != len(header):
            # Never let cells slide under the wrong heading without saying so (see erw.validate).
            ragged.append((i, len(r)))
            r = (r + [""] * len(header))[:len(header)]
        body.append([c.strip() for c in r])
    df = pd.DataFrame(body, columns=header, dtype=object)
    df.attrs["file"] = f"data/lab/{path.name}"
    df.attrs["key"] = hits[0]
    df.attrs["ragged_rows"] = ragged      # [(line number, number of cells)], empty when the file is sound
    return df


def load_all(lab_dir: Path | str = LAB) -> dict[str, pd.DataFrame]:
    """Every table in the data pack, keyed by file stem. A `status_class` column is added wherever
    there is a `status` column, and the original text is left untouched."""
    out = {}
    for stem in table_files(lab_dir):
        df = read_table(stem, lab_dir)
        if "status" in df.columns:
            ids = df["sample_id"] if "sample_id" in df.columns else [""] * len(df)
            df["status_class"] = [effective_status(i, s) for i, s in zip(ids, df["status"])]
        out[stem] = df
    return out


def _fallback_inputs() -> dict:
    return json.loads(INPUTS_JSON.read_text())


def _prov(value, status, source, unit=""):
    return {"value": value, "unit": unit, "status": status, "source": source}


# ====================================================================== rock
def rock_samples(lab_dir=LAB) -> dict[str, dict]:
    """Every row of 01 as a dict. `usable` is True when all eight model oxides are present."""
    df = read_table("01", lab_dir)
    out = {}
    for _, r in df.iterrows():
        sid = r["sample_id"]
        ox = {k: num(r.get(k)) for k in MODEL_OXIDES}
        feo = num(r.get("FeO"))
        usable = all(v is not None for v in ox.values())
        rec = dict(sample_id=sid, source=r.get("source", ""), material=r.get("material", ""),
                   status=r.get("status", ""), status_class=effective_status(sid, r.get("status", "")),
                   notes=r.get("notes", ""), usable=usable, feo=feo,
                   extras={k: parse_value(r.get(k)) for k in ("MnOx", "P2O5", "CO2", "LOI", "SO3")},
                   label=f"{sid}: {r.get('source', '')}, {r.get('material', '')}".strip(": ,"))
        if usable:
            rec["oxides"], rec["fe2_fraction"] = _iron(ox, feo)
        out[sid] = rec
    return out


def _iron(ox: dict, feo: float | None):
    """The model carries total iron as Fe2O3 plus a ferrous fraction. If FeO is reported, the
    Fe2O3 column is taken to be ferric iron only, and the ferrous fraction follows from the two.
    If FeO is blank, Fe2O3 is total iron and the fraction comes from the model settings (None here)."""
    ox = dict(ox)
    if feo is None:
        return ox, None
    fe2 = feo / MW_FEO
    fe3 = 2 * ox["Fe2O3"] / MW_FE2O3
    ox["Fe2O3"] = (fe2 + fe3) / 2 * MW_FE2O3
    return ox, (fe2 / (fe2 + fe3) if fe2 + fe3 > 0 else 0.0)


def rock_density(sample_id: str, lab_dir=LAB):
    """Particle density in kg/m3 from 02, with its provenance. Loose bulk density is not used.
    Falls back to data/inputs.json if the sample has no density row."""
    df = read_table("02", lab_dir)
    rows = df[(df["sample_id"] == sample_id) & df["property"].str.lower().str.contains("density")
              & ~df["property"].str.lower().str.contains("bulk|loose")]
    for _, r in rows.iterrows():
        p = parse_value(r["value"])
        if p.is_number:
            v = p.value * (1000.0 if r["unit"].strip().lower() in ("g/cm3", "t/m3", "g/cm³") else 1.0)
            detail = f" (midpoint of {p.text})" if p.qualifier == "range" else ""
            return v, _prov(v, effective_status(sample_id, r["status"]),
                            f"{df.attrs['file']}: {sample_id}, {r['property']}{detail}", "kg/m3")
    v = _fallback_inputs()["feedstock"]["density_g_m3"] / 1000.0
    return v, _prov(v, "assumed", f"fallback data/inputs.json (no density row for {sample_id} in 02)", "kg/m3")


# ---------------------------------------------------------------------- grain size
def psd_points(sample_id: str, lab_dir=LAB) -> list[dict]:
    """PSD points of one sample from 03, sorted by size: [{size_um, cum_pct, status_class}, ...]."""
    df = read_table("03", lab_dir)
    pts = []
    for _, r in df[df["sample_id"] == sample_id].iterrows():
        d, c = num(r["size_um"]), num(r["cumulative_pct_passing"])
        if d is not None and c is not None and d > 0:
            pts.append(dict(size_um=d, cum_pct=c, status_class=effective_status(sample_id, r["status"]),
                            status=r["status"]))
    return sorted(pts, key=lambda p: p["size_um"])


def psd_ids(lab_dir=LAB) -> list[str]:
    """Sample ids in 03 that have at least two usable points."""
    df = read_table("03", lab_dir)
    return [s for s in dict.fromkeys(df["sample_id"]) if len(psd_points(s, lab_dir)) >= 2]


def _surface_diameter(d1: float, d2: float) -> float:
    """Diameter with the same specific surface as mass spread evenly in log(size) between d1 and d2
    (a straight line on the usual semi-log grading curve): ln(d2/d1) / (1/d1 - 1/d2)."""
    return math.log(d2 / d1) / (1 / d1 - 1 / d2)


def build_bins(points: list[dict], n_min: int = 5, n_max: int = 8, fines_rep_um: float = 30.0,
               rep: str = "surface"):
    """Turn cumulative PSD points into 5-8 model bins.

    The grading curve is interpolated as a straight line in log(size) between points. Each bin is
    the mass between two sizes, with a representative diameter that keeps its surface area
    (`rep="surface"`) or the arithmetic midpoint (`rep="midpoint"`, the older hand-made bins).
    Mass finer than the smallest point has no lower size limit in a sieve analysis, so its
    diameter is an assumption: `fines_rep_um`, capped at half the smallest size.

    Returns (bins, detail): bins = [[mass_fraction, diameter_m], ...] finest first, for the model;
    detail = one dict per bin with its size range and whether it rests on assumed points.
    """
    pts = [p for p in sorted(points, key=lambda p: p["size_um"])]
    if len(pts) < 2:
        raise ValueError("a grain-size curve needs at least two points")
    if abs(pts[-1]["cum_pct"] - 100.0) > 0.5:
        raise ValueError("the coarsest PSD point must be 100 % passing")
    if any(b["cum_pct"] < a["cum_pct"] for a, b in zip(pts, pts[1:])):
        raise ValueError("cumulative % passing must not decrease with size")

    def assumed(a, b=None):
        return any(p is not None and p["status_class"] not in ("measured", "stand-in") for p in (a, b))

    edges = [(p["size_um"], p["cum_pct"], p) for p in pts]

    def interp(d):
        for (d1, c1, _), (d2, c2, _) in zip(edges, edges[1:]):
            if d1 <= d <= d2:
                return c1 + (c2 - c1) * math.log(d / d1) / math.log(d2 / d1)
        raise ValueError(d)

    have_fines = edges[0][1] > 0
    n_intervals_max = n_max - (1 if have_fines else 0)
    n_intervals_min = max(1, n_min - (1 if have_fines else 0))
    work = [(d, c, assumed(p)) for d, c, p in edges]
    if len(work) - 1 > n_intervals_max:                 # too many points (laser PSD): resample
        d0, dn = work[0][0], work[-1][0]
        any_assumed = any(w[2] for w in work)
        sizes = [d0 * (dn / d0) ** (i / n_intervals_max) for i in range(n_intervals_max + 1)]
        work = [(d, interp(d), any_assumed) for d in sizes]
    while len(work) - 1 < n_intervals_min:              # too few: split the widest interval in log size
        i = max(range(len(work) - 1), key=lambda j: math.log(work[j + 1][0] / work[j][0]))
        d = math.sqrt(work[i][0] * work[i + 1][0])
        work.insert(i + 1, (d, interp(d), work[i][2] or work[i + 1][2]))

    bins, detail = [], []
    if have_fines:
        d = min(fines_rep_um, 0.5 * work[0][0])
        bins.append([work[0][1] / 100.0, d * 1e-6])
        detail.append(dict(from_um=0.0, to_um=work[0][0], mass_pct=work[0][1], rep_um=d,
                           basis="assumed diameter for the fraction below the smallest size"
                                 + ("; amount from an assumed point" if work[0][2] else "")))
    for (d1, c1, a1), (d2, c2, a2) in zip(work, work[1:]):
        if c2 - c1 <= 0:
            continue
        d = _surface_diameter(d1, d2) if rep == "surface" else (d1 + d2) / 2
        bins.append([(c2 - c1) / 100.0, d * 1e-6])
        detail.append(dict(from_um=d1, to_um=d2, mass_pct=c2 - c1, rep_um=d,
                           basis="assumed split" if (a1 or a2) else "between measured points"))
    total = sum(b[0] for b in bins)
    bins = [[f / total, d] for f, d in bins]
    return bins, detail


# ---------------------------------------------------------------------- mineralogy
_PLAG_AN = {"albite": 0.05, "oligoclase": 0.2, "andesine": 0.4, "labradorite": 0.6, "bytownite": 0.8,
            "anorthite": 0.95}


def measured_mineralogy(sample_id: str, lab_dir=LAB) -> dict | None:
    """Measured (XRD) phases of a sample from 05, mapped onto the model's minerals. Returns None
    unless the sample has rows whose status is `measured` with numbers in them, so the normative
    estimate (status `assumed`) is never mistaken for a measurement.

    Returned dict: phases {Plag, Kfs, Cpx, Fo, Fa, Gls, Inert} in wt % (as reported, not
    renormalised), plag_an (anorthite fraction, or None if not stated), olivine_fo (or None),
    rows (how each reported phase was mapped), total_wt_pct, source.
    """
    df = read_table("05", lab_dir)
    rows = df[df["sample_id"] == sample_id]
    phases = dict(Plag=0.0, Kfs=0.0, Cpx=0.0, Fo=0.0, Fa=0.0, Gls=0.0, Inert=0.0)
    mapping, plag_an, ol_fo, ol_bulk, n = [], None, None, 0.0, 0
    for _, r in rows.iterrows():
        w = num(r["wt_pct"])
        if w is None or status_class(r["status"]) != "measured":
            continue
        n += 1
        name = r["phase"].lower()
        key = "Inert"
        if "glass" in name or "amorph" in name:
            key = "Gls"
        elif "plag" in name or any(k in name for k in _PLAG_AN):
            key = "Plag"
            m = re.search(r"an\s*(\d{1,3})", name)
            an = float(m.group(1)) / 100 if m else next((v for k, v in _PLAG_AN.items() if k in name), None)
            if an is not None:      # mass-weighted if several plagioclase rows
                plag_an = an if plag_an is None else (plag_an * phases["Plag"] + an * w) / (phases["Plag"] + w)
        elif any(k in name for k in ("k-feld", "kfeld", "k feld", "sanidine", "orthoclase", "anorthoclase",
                                     "alkali feld", "microcline")):
            key = "Kfs"
        elif any(k in name for k in ("pyroxene", "augite", "diopside", "pigeonite")):
            key = "Cpx"
        elif "forsterite" in name:
            key = "Fo"
        elif "fayalite" in name:
            key = "Fa"
        elif "olivine" in name:
            m = re.search(r"fo\s*(\d{1,3})", name)
            if m:
                x = float(m.group(1)) / 100
                wf = x * 140.69 / (x * 140.69 + (1 - x) * 203.77)      # mole fraction -> mass fraction
                phases["Fo"] += w * wf
                phases["Fa"] += w * (1 - wf)
                ol_fo = x
                mapping.append((r["phase"], w, f"olivine split Fo{x*100:.0f}"))
            else:
                ol_bulk += w
                mapping.append((r["phase"], w, "olivine, composition not stated"))
            continue
        phases[key] += w
        mapping.append((r["phase"], w, key))
    if n == 0:
        return None
    return dict(phases=phases, plag_an=plag_an, olivine_fo=ol_fo, olivine_unsplit_wt_pct=ol_bulk,
                rows=mapping, total_wt_pct=sum(phases.values()) + ol_bulk,
                source=f"{df.attrs['file']}: measured phases for {sample_id}")


# ====================================================================== soil
SOIL_REQUIRED = ("pH_H2O", "CEC_pH7_cmol_kg", "exch_Ca", "exch_Mg", "exch_K")


def soils(lab_dir=LAB, tables=("08", "09")) -> dict[str, dict]:
    """Every soil row in 08 and 09 that has pH, CEC and exchangeable Ca, Mg and K.

    Keys are sample ids. Units: cmol(+)/kg for exchange properties. Measured acidity is carried
    as `acidity_pH7` (exchange acidity at pH 7) and `acidity_kcl` (KCl-exchangeable acidity, from
    `exch_Al_KCl`, or ECEC minus the bases) when the lab reported them.
    """
    out = {}
    for key in tables:
        df = read_table(key, lab_dir)
        for _, r in df.iterrows():
            if any(num(r.get(c)) is None for c in SOIL_REQUIRED):
                continue
            sid = r["sample_id"]
            ca, mg, k, na = (num(r.get(c)) for c in ("exch_Ca", "exch_Mg", "exch_K", "exch_Na"))
            bases = ca + mg + k + (na or 0.0)
            s = dict(sample_id=sid, table=df.attrs["file"],
                     label=f"{sid}: {r.get('site', '')}, {r.get('plot', '')} ({r.get('land_use', '')})"
                           .replace(", ()", "").replace(",  (", " ("),
                     status=r.get("status", ""), status_class=effective_status(sid, r.get("status", "")),
                     site=r.get("site", ""), land_use=r.get("land_use", ""), notes=r.get("notes", ""),
                     pH=num(r["pH_H2O"]), cec=num(r["CEC_pH7_cmol_kg"]), ca=ca, mg=mg, k=k)
            if na is not None:
                s["na"] = na
            for src, dst in (("pH_KCl", "pH_kcl"), ("OM_pct", "om_pct"), ("base_saturation_pct", "bs_pct"),
                             ("exch_acidity_pH7", "acidity_pH7"), ("bulk_density_g_cm3", "bulk_density_g_cm3"),
                             ("latitude", "latitude"), ("longitude", "longitude"), ("ECEC", "ecec"),
                             ("exch_Al_KCl", "acidity_kcl")):
                v = num(r.get(src))
                if v is not None:
                    s[dst] = v
            if "acidity_kcl" not in s and "ecec" in s and s["ecec"] > bases:
                s["acidity_kcl"] = s["ecec"] - bases
            out[sid] = s
    return out


def resolve_soil(key: str) -> str:
    return SOIL_ALIASES.get(key, key)


def soil_gaps(lab_dir=LAB, table="09") -> list[tuple[str, list[str]]]:
    """Rows of a soil table that cannot be run yet, with the columns still missing."""
    df = read_table(table, lab_dir)
    return [(r["sample_id"], [c for c in SOIL_REQUIRED if num(r.get(c)) is None])
            for _, r in df.iterrows() if any(num(r.get(c)) is None for c in SOIL_REQUIRED)]


def in_site_mean(soil: dict) -> bool:
    """Field duplicates (QC) and the separate forest zone are run but left out of the site mean."""
    text = f"{soil.get('sample_id', '')} {soil.get('land_use', '')} {soil.get('notes', '')}".lower()
    return not any(w in text for w in ("dup", "forest"))


def site_mean_soil(soil_list: list[dict], label="Site mean") -> dict:
    """Arithmetic mean of the composites (pH averaged as pH, the usual soil-survey convention).
    Optional fields are averaged only if every composite has them."""
    if not soil_list:
        raise ValueError("no soils to average")
    keys = set(soil_list[0])
    for s in soil_list[1:]:
        keys &= set(s)
    out = dict(sample_id="SITE-MEAN", label=f"{label} of {len(soil_list)} composites",
               status="calculated", status_class=soil_list[0].get("status_class", "calculated"),
               members=[s["sample_id"] for s in soil_list])
    for k in sorted(keys):
        if all(isinstance(s[k], (int, float)) and not isinstance(s[k], bool) for s in soil_list):
            out[k] = sum(s[k] for s in soil_list) / len(soil_list)
    return out


# ====================================================================== settings, climate, site, rates
# (table, row label, cfg key, factor, unit). Row labels are matched without regard to case.
_SETTING_ROWS = [
    ("14", "setting", "Column depth", "column_depth_cm", 1.0, "cm"),
    ("14", "setting", "Cell size", "cell_cm", 1.0, "cm"),
    ("14", "setting", "Mixing depth", "mix_depth_cm", 1.0, "cm"),
    ("14", "setting", "Surface roughness", "roughness", 1.0, "x geometric"),
    ("14", "setting", "Field-to-lab rate factor", "field_factor", 1.0, ""),
    ("14", "setting", "Fe2+/Fe total", "fe2_fraction", 1.0, ""),
    ("14", "setting", "Soil-solution Ca used for calibration", "soil_solution_ca_mol", 1e-3, "mol/L"),
    ("14", "setting", "Glass fraction of the rock", "glass_frac", 1.0, ""),
    ("14", "setting", "Exchangeable (KCl) acidity when not measured", "effective_acidity_cmol_kg", 1.0, "cmol(+)/kg"),
    ("14", "setting", "Exchangeable Na share of CEC when not measured", "na_fraction_default", 1.0, ""),
    ("14", "setting", "Simulation length", "years", 1.0, "years"),
    ("14", "setting", "Diameter assumed for the fraction below the smallest sieve", "fines_rep_um", 1.0, "um"),
    ("14", "setting", "Half-point pH of the variable charge", "vc_half_ph", 1.0, "pH"),
    ("14", "setting", "Start month of the run", "start_month", 1.0, "month"),
    ("14", "setting", "Reactive Al hydroxide in the soil", "al_reservoir_mol_kgw", 1.0, "mol/kgw"),
    ("12", "parameter", "Mean soil temperature", "temp_c", 1.0, "°C"),
    ("12", "parameter", "Net drainage below root zone", "drainage_mm_yr", 1.0, "mm/yr"),
    ("12", "parameter", "Soil log pCO2", "log_pco2", 1.0, "log atm"),
    ("12", "parameter", "Soil water content", "water_content", 1.0, "m3/m3"),
    ("12", "parameter", "Soil bulk density", "bulk_density_kg_m3", 1000.0, "kg/m3"),
    ("11", "item", "Application rate", "application_t_ha", 1.0, "t/ha"),
]
_INT_KEYS = {"years", "start_month"}
MONTHS = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October",
          "November", "December"]

_ISO_ROWS = [
    ("14", "setting", "Near-Field Zone depth", "nfz_depth_cm", "cm"),
    ("14", "setting", "eta RiverOcean", "eta_river_ocean", ""),
    ("11", "item", "Untreated control", "control_area_pct", "% of area"),
]

_RATE_KEYS = [("plagioclase", "Plag"), ("k-feldspar", "Kfs"), ("diopside", "Cpx"), ("forsterite", "Fo"),
              ("fayalite", "Fa"), ("glass", "Gls")]

_EMISSION_ROWS = {
    "haul distance quarry to field": "haul_distance_km",
    "truck payload": "truck_payload_t",
    "diesel, loaded": "diesel_l_per_km_loaded",
    "diesel, empty return": "diesel_l_per_km_empty",
    "spreading diesel": "spreading_diesel_l_per_t",
    "diesel emission factor": "diesel_kgco2e_per_l",
    "grid emission factor, são miguel": "grid_kgco2e_per_kwh",
    "quarrying allocation to fines": "quarrying_kgco2e_per_t",
    "other project emissions": "other_project_t_co2e",
}
_MILLING_ROWS = {
    "milling energy, 0/4 mm as supplied": "as_is_0-4mm",
    "milling energy, p80 1 mm": "milled_p80_1mm",
    "milling energy, p80 100 um": "milled_p80_100um",
}


def _find_row(df: pd.DataFrame, col: str, label: str):
    hit = df[df[col].str.lower() == label.lower()]
    return hit.iloc[0] if len(hit) else None


def _value_col(df):
    return "value_used_in_model" if "value_used_in_model" in df.columns else "value"


def lab_defaults(lab_dir=LAB):
    """Model settings from 11, 12 and 14. Returns (defaults, isometric, provenance).
    Any setting missing from the lab files is taken from data/inputs.json and labelled as such."""
    fb = _fallback_inputs()
    tables = {k: read_table(k, lab_dir) for k in ("11", "12", "14")}
    d, prov = {}, {}
    for tkey, col, label, key, factor, unit in _SETTING_ROWS:
        df = tables[tkey]
        row = _find_row(df, col, label)
        v = num(row[_value_col(df)]) if row is not None else None
        if v is None:
            if key not in fb["defaults"]:
                continue
            d[key] = fb["defaults"][key]
            prov[key] = _prov(d[key], "assumed", f"fallback data/inputs.json (no value for '{label}' in {df.attrs['file']})", unit)
            continue
        v = v * factor
        d[key] = int(round(v)) if key in _INT_KEYS else v
        prov[key] = _prov(d[key], status_class(row["status"]) if row["status"] else "unspecified",
                          f"{df.attrs['file']}: {label}", unit)
    # the acidity pool is a text choice
    row = _find_row(tables["14"], "setting", "Soil acidity pool")
    txt = (row["value"] if row is not None else "").lower()
    d["acidity_mode"] = "andosol" if txt.startswith("andosol") else "effective" if txt.startswith("effective") else "pH7_CEC"
    prov["acidity_mode"] = _prov(d["acidity_mode"], status_class(row["status"]) if row is not None else "assumed",
                                 f"{tables['14'].attrs['file']}: Soil acidity pool" if row is not None
                                 else "fallback data/inputs.json", "")
    row = _find_row(tables["14"], "setting", "Seasonal drainage")
    d["seasonal_drainage"] = bool(row is not None and row["value"].strip().lower() in ("on", "yes", "true", "1"))
    prov["seasonal_drainage"] = _prov("on" if d["seasonal_drainage"] else "off",
                                      status_class(row["status"]) if row is not None else "assumed",
                                      f"{tables['14'].attrs['file']}: Seasonal drainage" if row is not None
                                      else "no row in table 14: steady drainage", "")
    shares, statuses = [], set()
    for m in MONTHS:
        r12 = _find_row(tables["12"], "parameter", f"Drainage share in {m}")
        v = num(r12[_value_col(tables["12"])]) if r12 is not None else None
        if v is None or v < 0:
            shares = []
            break
        shares.append(v)
        statuses.add(status_class(r12["status"]))
    if shares and sum(shares) > 0:
        d["drainage_monthly_share"] = [x / sum(shares) for x in shares]
        prov["drainage_monthly_share"] = _prov([round(x, 4) for x in d["drainage_monthly_share"]],
                                               "assumed" if "assumed" in statuses or len(statuses) != 1 else statuses.pop(),
                                               f"{tables['12'].attrs['file']}: Drainage share in January to December",
                                               "fraction of the year")
    row = _find_row(tables["14"], "setting", "Glass rate law")
    txt = (row["value"] if row is not None else "").lower()
    d["glass_rate_law"] = "gislason_oelkers_2003" if "gislason" in txt else "approximate"
    prov["glass_rate_law"] = _prov(d["glass_rate_law"], status_class(row["status"]) if row is not None else "assumed",
                                   f"{tables['14'].attrs['file']}: Glass rate law" if row is not None
                                   else "no row in table 14: older approximate rate", "")
    row = _find_row(tables["14"], "setting", "Aluminium phase controlling soil-water Al")
    txt = (row["value"] if row is not None else "").lower()
    d["al_phase"] = "gibbsite" if txt.startswith("gibbsite") else "amorphous"
    prov["al_phase"] = _prov(d["al_phase"], status_class(row["status"]) if row is not None else "assumed",
                             f"{tables['14'].attrs['file']}: Aluminium phase controlling soil-water Al"
                             if row is not None else "default in erw/data_io.py", "")
    for k, v in fb["defaults"].items():                      # anything still missing
        if k not in d:
            d[k] = v
            prov[k] = _prov(v, "assumed", "fallback data/inputs.json", "")
    iso = copy.deepcopy(fb.get("isometric", {}))
    for tkey, col, label, key, unit in _ISO_ROWS:
        df = tables[tkey]
        row = _find_row(df, col, label)
        v = num(row[_value_col(df)]) if row is not None else None
        if v is not None:
            iso[key] = v
            prov[f"isometric.{key}"] = _prov(v, status_class(row["status"]), f"{df.attrs['file']}: {label}", unit)
        elif key in iso:
            prov[f"isometric.{key}"] = _prov(iso[key], "assumed", "fallback data/inputs.json", unit)
    # Near-Field Zone: the deeper of the minimum depth and tillage depth + margin, on a cell boundary
    get14 = lambda label: (lambda r: num(r["value"]) if r is not None else None)(_find_row(tables["14"], "setting", label))
    till, margin = get14("Tillage depth"), get14("Near-Field Zone margin below tillage depth")
    nfz_min = float(iso.get("nfz_depth_cm", 20))
    if till is not None and margin is not None:
        nfz = max(nfz_min, till + margin)
        cell = float(d.get("cell_cm", 10))
        nfz_cells = math.ceil(nfz / cell - 1e-9)
        iso.update(nfz_depth_cm=nfz_cells * cell, tillage_depth_cm=till, nfz_margin_cm=margin, nfz_min_cm=nfz_min)
        prov["isometric.nfz_depth_cm"] = _prov(
            iso["nfz_depth_cm"], "plan",
            f"{tables['14'].attrs['file']}: deeper of {nfz_min:g} cm and tillage depth {till:g} cm + {margin:g} cm"
            + ("" if nfz_cells * cell == nfz else f", rounded up to the {cell:g} cm cells"), "cm")
    v = get14("Minimum control share of area")
    if v is not None:
        iso["control_min_pct"] = v
        prov["isometric.control_min_pct"] = _prov(v, "plan", f"{tables['14'].attrs['file']}: Minimum control share of area", "% of area")
    return d, iso, prov


def mc_parameters(lab_dir=LAB) -> list[dict]:
    """Inputs varied in the Monte Carlo and tornado runs, with their ranges, from tables 14 and 12.
    A row takes part when its `range` (table 14) or `range_tested` (table 12) cell holds a range."""
    spec = [("14", "setting", "Field-to-lab rate factor", "field_factor", "range", 1.0),
            ("14", "setting", "Surface roughness", "roughness", "range", 1.0),
            ("14", "setting", "Glass fraction of the rock", "glass_frac", "range", 1.0),
            ("12", "parameter", "Net drainage below root zone", "drainage_mm_yr", "range_tested", 1.0),
            ("12", "parameter", "Soil log pCO2", "log_pco2", "range_tested", 1.0),
            ("12", "parameter", "Mean soil temperature", "temp_c", "range_tested", 1.0)]
    out = []
    for tkey, col, label, key, rcol, factor in spec:
        df = read_table(tkey, lab_dir)
        row = _find_row(df, col, label)
        if row is None or rcol not in df.columns:
            continue
        p = parse_value(row[rcol])
        if p.qualifier != "range":
            continue
        dist = (row.get("distribution", "") or "uniform").strip().lower()
        out.append(dict(key=key, label=label, low=p.low * factor, high=p.high * factor, distribution=dist,
                        status=status_class(row["status"]), source=f"{df.attrs['file']}: {label}, range {row[rcol]}"))
    return out


def rate_parameters(lab_dir=LAB):
    """Kinetic constants from 13 as {Plag: [logk_acid, Ea_acid, n_acid, logk_neutral, Ea_neutral], ...}.
    A mineral with no acid term in the table (basaltic glass) gets a switched-off acid term."""
    df = read_table("13", lab_dir)
    fb = _fallback_inputs()["rates_palandri_kharaka_2004"]
    rates, prov, glass_go = {}, {}, None
    for _, r in df.iterrows():
        name = r["mineral"].lower()
        key = next((k for frag, k in _RATE_KEYS if frag in name), None)
        if key is None:
            continue
        vals = [num(r[c]) for c in ("logk_acid_25C", "Ea_acid_kJ_mol", "n_acid", "logk_neutral_25C",
                                    "Ea_neutral_kJ_mol")]
        if vals[3] is None or vals[4] is None:
            continue
        if vals[0] is None:
            vals[:3] = [-99.0, 0.0, 0.0]
        rates[key] = vals
        status = r["status"] if "status" in df.columns and r["status"] else "assumed"
        if key == "Gls":
            la, ea = num(r.get("logA_GO2003_mol_Si_m2_s")), num(r.get("Ea_GO2003_kJ_mol"))
            if la is not None and ea is not None:
                glass_go = dict(log_a=la, ea=ea)
                prov["glass_go"] = _prov(glass_go, status_class(status),
                                         f"{df.attrs['file']}: Basaltic glass, Gislason & Oelkers (2003) constants",
                                         "log mol Si/m2/s, kJ/mol")
        prov[f"rates.{key}"] = _prov(vals, status_class(status), f"{df.attrs['file']}: {r['mineral']} ({r['source']})",
                                     "log mol/m2/s, kJ/mol")
    for key, vals in fb.items():
        if not key.startswith("_") and key not in rates:
            rates[key] = vals
            prov[f"rates.{key}"] = _prov(vals, "assumed", "fallback data/inputs.json", "")
    if glass_go:
        rates["_glass_go"] = glass_go            # taken out again by build_config
    return rates, prov


def emission_factors(lab_dir=LAB):
    """Emission factors from 18, in the shape the model's emissions function expects."""
    df = read_table("18", lab_dir)
    e = copy.deepcopy(_fallback_inputs().get("emissions", {}))
    prov = {}
    for _, r in df.iterrows():
        label, v = r["item"].lower(), num(r["value"])
        if v is None:
            continue
        if label in _EMISSION_ROWS:
            e[_EMISSION_ROWS[label]] = v
            prov[f"emissions.{_EMISSION_ROWS[label]}"] = _prov(v, status_class(r["status"]),
                                                               f"{df.attrs['file']}: {r['item']}", r["unit"])
        elif label in _MILLING_ROWS:
            e.setdefault("milling_kwh_per_t", {})[_MILLING_ROWS[label]] = v
            prov[f"emissions.milling.{_MILLING_ROWS[label]}"] = _prov(v, status_class(r["status"]),
                                                                      f"{df.attrs['file']}: {r['item']}", r["unit"])
    return e, prov


def plant_uptake(lab_dir=LAB):
    """Net yearly removal of Ca, Mg and K from the field (kg/ha/yr) from table 19, by two routes:
    herbage carried off (yield x content x share leaving the field) and milk sold (litres x content).
    A blank share and a blank milk figure give zero: nothing is assumed. Returns (cfg entries, provenance)."""
    try:
        df = read_table("19", lab_dir)
    except KeyError:
        return {"uptake_kg_ha_yr": dict(Ca=0.0, Mg=0.0, K=0.0)}, {}
    get = lambda label: (lambda r: (num(r["value"]), status_class(r["status"])) if r is not None else (None, ""))(
        _find_row(df, "item", label))
    yield_t, y_status = get("Pasture dry-matter yield")
    share, s_status = get("Share of herbage cations that leave the field")
    milk, m_status = get("Milk sold per hectare")
    out, gross, statuses = {}, {}, set()
    for el in ("Ca", "Mg", "K"):
        herb, h_status = get(f"{el} in herbage")
        in_milk, k_status = get(f"{el} in milk")
        gross[el] = (yield_t or 0.0) * (herb or 0.0)                       # t DM/ha x g/kg = kg/ha
        v = 0.0
        if share is not None and yield_t is not None and herb is not None:
            v += share * gross[el]
            statuses |= {y_status, h_status, s_status}
        if milk is not None and in_milk is not None:
            v += milk * in_milk / 1000.0                                   # L/ha x g/L = g/ha -> kg/ha
            statuses |= {m_status, k_status}
        out[el] = v
    cfg = {"uptake_kg_ha_yr": out, "uptake_gross_kg_ha_yr": gross}
    depth, d_status = get("Rooting depth for uptake")
    if depth is not None:
        cfg["root_depth_cm"] = depth
    km, _ = get("Half-saturation of uptake as a share of the starting soil-water concentration")
    if km is not None:
        cfg["uptake_km_fraction"] = km
    active = any(v > 0 for v in out.values())
    status = ("assumed" if "assumed" in statuses or "to measure" in statuses else
              "stand-in" if "stand-in" in statuses else "measured") if active else "to measure"
    prov = {"uptake_kg_ha_yr": _prov({k: round(v, 2) for k, v in out.items()}, status,
                                     f"{df.attrs['file']}: " + ("yield x herbage content x share leaving, plus milk"
                                                                if active else "no removal: share leaving the field "
                                                                               "and milk sold are not filled in"),
                                     "kg/ha/yr")}
    return cfg, prov


# ====================================================================== the model configuration
def build_config(soil: str = DEFAULT_SOIL, psd: str = AS_IS_PSD, rock: str = DEFAULT_ROCK, lab_dir=LAB) -> dict:
    """A complete model configuration from the lab files, with a provenance record for every input.

    soil  a sample id from 08 or 09 (old names A_Moreira and B_TerraRiba still work)
    psd   "as_is_0-4mm" (the chosen rock's own grading curve), or a curve id from 03
          (old names milled_p80_1mm and milled_p80_100um still work)
    rock  a sample id from 01 with complete oxides
    """
    d, iso, prov = lab_defaults(lab_dir)
    rates, p2 = rate_parameters(lab_dir)
    emis, p3 = emission_factors(lab_dir)
    prov.update(p2)
    prov.update(p3)
    cfg = dict(d)
    cfg["isometric"] = iso
    cfg["emissions"] = emis
    up, p4 = plant_uptake(lab_dir)
    cfg.update(up)
    prov.update(p4)
    glass_go = rates.pop("_glass_go", None)
    if glass_go:
        cfg["glass_go"] = glass_go
    cfg["rates"] = rates
    cfg["data_source"] = "data/lab"

    # ---- rock
    rocks = rock_samples(lab_dir)
    if rock not in rocks or not rocks[rock]["usable"]:
        ok = [k for k, v in rocks.items() if v["usable"]]
        raise KeyError(f"rock sample {rock!r} has no complete oxide analysis in 01. Usable: {ok}")
    rk = rocks[rock]
    cfg["rock_key"] = rock
    cfg["rock_label"] = rk["label"]
    cfg["oxides"] = dict(rk["oxides"])
    prov["oxides"] = _prov(dict(rk["oxides"]), rk["status_class"], f"data/lab/01_rock_major_oxides.csv: {rock}", "wt %")
    if rk["fe2_fraction"] is not None:
        cfg["fe2_fraction"] = rk["fe2_fraction"]
        prov["fe2_fraction"] = _prov(rk["fe2_fraction"], rk["status_class"],
                                     f"data/lab/01_rock_major_oxides.csv: FeO reported for {rock}", "")
    dens, pd_ = rock_density(rock, lab_dir)
    cfg["density_g_m3"] = dens * 1000.0
    prov["density_g_m3"] = pd_

    # ---- grain size
    cfg["psd_key"] = psd
    psd_id = rock if psd == AS_IS_PSD else PSD_ALIASES.get(psd, psd)
    pts = psd_points(psd_id, lab_dir)
    borrowed = ""
    if len(pts) < 2 and psd == AS_IS_PSD and rock != DEFAULT_ROCK:
        pts = psd_points(DEFAULT_ROCK, lab_dir)
        borrowed = f" (no curve for {rock}; {DEFAULT_ROCK} used as a stand-in)"
        psd_id = DEFAULT_ROCK
    if len(pts) >= 2:
        cfg["psd"], cfg["psd_detail"] = build_bins(pts, fines_rep_um=cfg.get("fines_rep_um", 30.0))
        measured_only = all(p["status_class"] in ("measured", "stand-in") for p in pts)
        base = "stand-in" if not is_project_sample(psd_id) or borrowed else "measured"
        prov["psd"] = _prov([[round(f, 4), round(dm * 1e6, 1)] for f, dm in cfg["psd"]],
                            base if measured_only else "assumed",
                            f"data/lab/03_rock_grain_size.csv: {psd_id}, {len(pts)} points"
                            f"{'' if measured_only else ', some assumed'}{borrowed}", "fraction, um")
    else:
        fb = _fallback_inputs()["grain_size_bins"]
        if psd not in fb:
            raise KeyError(f"grain size {psd!r}: no curve in 03 and no entry in data/inputs.json")
        cfg["psd"] = copy.deepcopy(fb[psd]["bins"])
        prov["psd"] = _prov(cfg["psd"], "assumed", f"fallback data/inputs.json: {psd}", "fraction, m")

    # ---- mineralogy: measured XRD if 05 has it, otherwise the norm from the oxides
    xrd = measured_mineralogy(rock, lab_dir)
    if xrd:
        cfg["mineralogy"] = xrd
        prov["mineralogy"] = _prov("measured phases (XRD)", rk["status_class"], xrd["source"], "wt %")
    else:
        prov["mineralogy"] = _prov("normative estimate from the oxides", "assumed",
                                   "calculated in erw/model.py: no measured phases in data/lab/05_rock_mineralogy.csv "
                                   f"for {rock}", "wt %")

    # ---- soil
    all_soils = soils(lab_dir)
    sid = resolve_soil(soil)
    if sid not in all_soils:
        raise KeyError(f"soil {soil!r} is not a usable row in 08 or 09. Usable: {list(all_soils)}")
    s = copy.deepcopy(all_soils[sid])
    cfg["soil"] = s
    cfg["soil_key"] = sid
    prov["soil"] = _prov({k: s[k] for k in ("pH", "cec", "ca", "mg", "k")}, s["status_class"],
                         f"{s['table']}: {sid}", "pH; cmol(+)/kg")
    if "acidity_kcl" in s:
        cfg["effective_acidity_cmol_kg"] = s["acidity_kcl"]
        prov["effective_acidity_cmol_kg"] = _prov(s["acidity_kcl"], s["status_class"],
                                                  f"{s['table']}: {sid}, KCl-exchangeable acidity", "cmol(+)/kg")
    if "bulk_density_g_cm3" in s:
        cfg["bulk_density_kg_m3"] = s["bulk_density_g_cm3"] * 1000.0
        prov["bulk_density_kg_m3"] = _prov(cfg["bulk_density_kg_m3"], s["status_class"],
                                           f"{s['table']}: {sid}, bulk density", "kg/m3")
    cfg["provenance"] = prov
    return cfg


def provenance_table(cfg: dict) -> pd.DataFrame:
    """The inputs of a run as a table: parameter, value, unit, status, source."""
    rows = []
    for k, p in (cfg.get("provenance") or {}).items():
        v = p["value"]
        rows.append({"parameter": k, "value": json.dumps(v) if isinstance(v, (list, dict)) else str(v),
                     "unit": p.get("unit", ""), "status": p["status"], "source": p["source"]})
    return pd.DataFrame(rows)
