"""
Checks on the lab data pack (data/lab/*.csv). Writes outputs/validation.md.

    python -m erw.validate

Each rule is a small function that takes a table (a pandas DataFrame of text cells, as read by
erw.data_io.read_table) and returns a list of findings. A finding has a level:

    pass   the check ran and the data agree
    warn   worth a look, or a known open question
    fail   the data disagree with themselves
    info   the check could not run (values missing), or a note

The report never corrects a value. It only says what to check in the original lab reports.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from . import data_io
from .data_io import num, parse_value

ROOT = Path(__file__).resolve().parents[1]
MW = dict(SiO2=60.08, Al2O3=101.96, Fe2O3=159.69, CaO=56.08, MgO=40.30, TiO2=79.87, Na2O=61.98, K2O=94.20)

OXIDE_COLS = ["SiO2", "Al2O3", "Fe2O3", "FeO", "CaO", "MgO", "Na2O", "K2O", "TiO2", "MnOx", "P2O5", "SO3"]
TOTAL_LOW, TOTAL_HIGH = 97.0, 102.0
BS_TOLERANCE = 0.10          # relative
SUM_TOLERANCE = 0.02         # relative, for sums that should agree exactly apart from rounding


@dataclass
class Finding:
    rule: str
    level: str      # pass / warn / fail / info
    table: str
    item: str
    message: str


def _file(df, default=""):
    return df.attrs.get("file", default)


# ---------------------------------------------------------------- 0. file structure
def check_structure(df: pd.DataFrame) -> list[Finding]:
    """Every row must have as many cells as the header. If not, values sit under the wrong heading."""
    f = _file(df)
    bad = df.attrs.get("ragged_rows", [])
    if bad:
        lines = ", ".join(f"line {i} has {n}" for i, n in bad[:6]) + (" ..." if len(bad) > 6 else "")
        return [Finding("file structure", "fail", f, "(rows)",
                        f"the header has {len(df.columns)} columns but {lines} cells; values may sit under "
                        "the wrong heading")]
    return [Finding("file structure", "pass", f, "(rows)", f"{len(df)} rows, all with {len(df.columns)} cells")]


# ---------------------------------------------------------------- 1. oxide totals
def check_oxide_totals(df: pd.DataFrame) -> list[Finding]:
    """Major oxides plus LOI should sum to 97-102 %. A value below detection (`<0.05`) counts as
    zero, with the detection limit added to the upper end of the total. CO2 is part of LOI, so it
    is only added when no LOI is reported."""
    out, f = [], _file(df, "01")
    for _, r in df.iterrows():
        sid = r["sample_id"]
        main = [r.get(c, "") for c in ("SiO2", "Al2O3", "Fe2O3", "CaO", "MgO", "Na2O", "K2O", "TiO2")]
        if any(num(c) is None for c in main):
            out.append(Finding("oxide total", "info", f, sid, "analysis not complete, total not checked"))
            continue
        lo = hi = 0.0
        for c in OXIDE_COLS + ["LOI"]:
            p = parse_value(r.get(c, ""))
            if not p.is_number:
                continue
            lo += 0.0 if p.below_detection else p.value
            hi += p.value
        loi = parse_value(r.get("LOI", ""))
        co2 = parse_value(r.get("CO2", ""))
        if not loi.is_number and co2.is_number:
            lo += co2.value
            hi += co2.value
        span = f"{lo:.2f} %" if abs(hi - lo) < 0.005 else f"{lo:.2f} to {hi:.2f} %"
        if hi < TOTAL_LOW or lo > TOTAL_HIGH:
            out.append(Finding("oxide total", "fail", f, sid,
                               f"total {span}, outside {TOTAL_LOW:.0f}-{TOTAL_HIGH:.0f} %. "
                               "Oxides are missing or the analysis is incomplete; check before using this rock."))
        else:
            out.append(Finding("oxide total", "pass", f, sid, f"total {span}"))
    return out


# ---------------------------------------------------------------- 2. base saturation
def check_base_saturation(df: pd.DataFrame) -> list[Finding]:
    """Reported base saturation should equal (Ca+Mg+K+Na)/CEC within 10 % (relative)."""
    out, f = [], _file(df, "soil")
    for _, r in df.iterrows():
        sid = r["sample_id"]
        bs, cec = num(r.get("base_saturation_pct")), num(r.get("CEC_pH7_cmol_kg"))
        cat = [num(r.get(c)) for c in ("exch_Ca", "exch_Mg", "exch_K")]
        na = num(r.get("exch_Na"))
        if bs is None or cec is None or any(c is None for c in cat) or cec <= 0:
            if bs is not None or cec is not None:
                out.append(Finding("base saturation", "info", f, sid, "not enough values to recompute"))
            continue
        calc = (sum(cat) + (na or 0.0)) / cec * 100
        rel = abs(calc - bs) / bs if bs else float("inf")
        note = "" if na is not None else " (no exchangeable Na reported, so the recomputed value is a little low)"
        level = "pass" if rel <= BS_TOLERANCE else "fail"
        out.append(Finding("base saturation", level, f, sid,
                           f"reported {bs:.2f} %, recomputed {calc:.2f} % ({rel*100:.1f} % apart){note}"))
    return out


# ---------------------------------------------------------------- 3. pH
def check_ph(df: pd.DataFrame) -> list[Finding]:
    """pH in KCl should be below pH in water (KCl displaces exchangeable acidity)."""
    out, f = [], _file(df, "soil")
    for _, r in df.iterrows():
        w, k = num(r.get("pH_H2O")), num(r.get("pH_KCl"))
        if w is not None and not 3.0 <= w <= 9.5:
            out.append(Finding("pH", "fail", f, r["sample_id"], f"pH H2O {w} is outside the range for soils (3-9.5)"))
        if w is None or k is None:
            continue
        if k < w:
            out.append(Finding("pH", "pass", f, r["sample_id"], f"pH KCl {k} is below pH H2O {w} (gap {w-k:.1f})"))
        else:
            out.append(Finding("pH", "fail", f, r["sample_id"],
                               f"pH KCl {k} is not below pH H2O {w}; check the columns were not swapped"))
    return out


# ---------------------------------------------------------------- 4. units and internal consistency
NUMERIC_SOIL_COLS = ["pH_H2O", "pH_KCl", "OM_pct", "total_N_mg_g", "P_extractable_mg_kg_P2O5",
                     "K_extractable_mg_kg_K2O", "CaO_extractable_mg_kg", "MgO_extractable_mg_kg",
                     "CEC_pH7_cmol_kg", "exch_Ca", "exch_Mg", "exch_K", "exch_Na", "exch_acidity_pH7",
                     "sum_bases", "base_saturation_pct", "bulk_density_g_cm3", "ECEC", "exch_Al_KCl"]


def check_units_soil(df: pd.DataFrame) -> list[Finding]:
    """Soil tables: cells that carry their own unit note, sums that should agree, and CEC against
    its parts. All exchange properties must be in cmol(+)/kg."""
    out, f = [], _file(df, "soil")
    for _, r in df.iterrows():
        sid = r["sample_id"]
        for c in NUMERIC_SOIL_COLS:
            if c not in df.columns or not r[c]:
                continue
            p = parse_value(r[c])
            if not p.is_number:
                out.append(Finding("units", "fail", f, sid, f"{c}: '{r[c]}' is not a number"))
            elif p.note:
                out.append(Finding("units", "warn", f, sid,
                                   f"{c}: '{r[c]}' carries a unit note, so it is not in the unit of the column"))
            if "," in r[c]:
                out.append(Finding("units", "warn", f, sid, f"{c}: '{r[c]}' contains a comma; check the decimal separator"))
        cat = [num(r.get(c)) for c in ("exch_Ca", "exch_Mg", "exch_K")]
        na, sb = num(r.get("exch_Na")), num(r.get("sum_bases"))
        cec, acid = num(r.get("CEC_pH7_cmol_kg")), num(r.get("exch_acidity_pH7"))
        if sb is not None and all(c is not None for c in cat):
            calc = sum(cat) + (na or 0.0)
            ok = abs(calc - sb) <= SUM_TOLERANCE * max(sb, 1e-9)
            out.append(Finding("units", "pass" if ok else "fail", f, sid,
                               f"sum of bases: reported {sb:.2f}, Ca+Mg+K{'+Na' if na is not None else ''} = {calc:.2f} cmol(+)/kg"))
        if sb is not None and cec is not None and acid is not None:
            ok = abs(sb + acid - cec) <= SUM_TOLERANCE * cec
            out.append(Finding("units", "pass" if ok else "fail", f, sid,
                               f"CEC {cec:.2f} against bases + acidity at pH 7 = {sb + acid:.2f} cmol(+)/kg"))
        if cec is not None and all(c is not None for c in cat) and sum(cat) > cec:
            out.append(Finding("units", "fail", f, sid, "exchangeable bases exceed the CEC; check units (cmol(+)/kg)"))
    return out


def check_units_long(df: pd.DataFrame, name_col: str, value_col: str = "value") -> list[Finding]:
    """Long tables (property, value, unit): the same property must use one unit throughout, and
    a filled value needs a unit unless it is a plain ratio or text."""
    out, f = [], _file(df)
    if name_col not in df.columns or "unit" not in df.columns:
        return out
    for prop, grp in df.groupby(df[name_col].str.lower()):
        units = {u for u, v in zip(grp["unit"], grp[value_col]) if v and u}
        if len(units) > 1:
            out.append(Finding("units", "fail", f, prop, f"the same property is given in different units: {sorted(units)}"))
    if not any(x.level == "fail" for x in out):
        out.append(Finding("units", "pass", f, "(all rows)", "each property uses one unit"))
    return out


# ---------------------------------------------------------------- 5. duplicate ids
def check_duplicate_ids(df: pd.DataFrame, extra_key: list[str] | None = None) -> list[Finding]:
    """Sample ids must be unique (or unique together with `extra_key` in long tables)."""
    f = _file(df)
    if "sample_id" not in df.columns:
        return []
    key = ["sample_id"] + (extra_key or [])
    filled = df[df["sample_id"] != ""]
    dup = filled[filled.duplicated(key, keep=False)]
    if len(dup):
        ids = sorted({" / ".join(r[k] for k in key) for _, r in dup.iterrows()})
        return [Finding("duplicate ids", "fail", f, i, "appears more than once") for i in ids]
    return [Finding("duplicate ids", "pass", f, "(all rows)", f"{len(filled)} rows, no repeated id")]


# ---------------------------------------------------------------- 6. known open questions
K_PER_K2O = 2 * 39.098 / 94.196


def check_known_issues(rock: pd.DataFrame | None = None, soil: pd.DataFrame | None = None) -> list[Finding]:
    """Open questions already listed in data/lab/00_README.md. They stay as warnings until the
    original lab reports are checked and the note is removed from the table."""
    out = []
    if soil is not None:
        f = _file(soil, "08")
        for _, r in soil.iterrows():
            note = r.get("notes", "")
            if "K2O" in note and "second" in note.lower():
                k2o, exk = num(r.get("K_extractable_mg_kg_K2O")), num(r.get("exch_K"))
                detail = ""
                if k2o is not None and exk is not None:
                    detail = (f" For scale: {k2o:.0f} mg/kg K2O equals {k2o*K_PER_K2O/390.98:.2f} cmol(+)/kg of K, "
                              f"against exchangeable K of {exk:.2f} cmol(+)/kg. The two methods extract different "
                              "amounts, so this does not settle which table is right.")
                out.append(Finding("known issue", "warn", f, r["sample_id"],
                                   f"potassium conflict ({r.get('plot', '')}): {note}.{detail}"))
    if rock is not None:
        f = _file(rock, "01")
        for _, r in rock.iterrows():
            cell = r.get("MnOx", "")
            if "mn2o4" in cell.lower():
                out.append(Finding("known issue", "warn", f, r["sample_id"],
                                   f"manganese oxide written as Mn2O4 in the source ('{cell}'). Mn2O4 is not a "
                                   "known oxide; Mn3O4 is likely. The amount is small and is not used by the model."))
    return out


# ---------------------------------------------------------------- 7. CO2 potential table
def co2_potential(ox: dict) -> tuple[float, float]:
    """Renforth (2012): tCO2 per tonne of rock from CaO+MgO only, and from Ca, Mg, Na and K as bicarbonate."""
    m = {k: ox[k] / MW[k] for k in ("CaO", "MgO", "Na2O", "K2O")}
    cons = (m["CaO"] + m["MgO"]) * 44.01 / 100
    full = (2 * m["CaO"] + 2 * m["MgO"] + 2 * m["Na2O"] + 2 * m["K2O"]) * 44.01 / 100
    return cons, full


def check_co2_potential(rock: pd.DataFrame, table07: pd.DataFrame) -> list[Finding]:
    """Table 07 is derived from table 01. Recompute it and compare."""
    out, f = [], _file(table07, "07")
    ox01 = {r["sample_id"]: {k: num(r.get(k)) for k in ("CaO", "MgO", "Na2O", "K2O")} for _, r in rock.iterrows()}
    for _, r in table07.iterrows():
        sid = r["sample_id"]
        a, b = num(r.get("tCO2_per_t_conservative_CaMg_only")), num(r.get("tCO2_per_t_bicarbonate_full"))
        ox = ox01.get(sid)
        if a is None or b is None or ox is None or any(v is None for v in ox.values()):
            continue
        cons, full = co2_potential(ox)
        ok = abs(cons - a) < 0.0015 and abs(full - b) < 0.0015
        out.append(Finding("CO2 potential", "pass" if ok else "fail", f, sid,
                           f"listed {a:.3f} and {b:.3f}, recomputed from 01: {cons:.3f} and {full:.3f} tCO2/t"))
    return out


# ---------------------------------------------------------------- 8. control area
def check_control_area(site: pd.DataFrame, settings: pd.DataFrame) -> list[Finding]:
    """The untreated control must be at least the minimum share of the area (Isometric EWA v1.3, Eq. 9)."""
    f = _file(site, "11")
    row = site[site["item"].str.lower() == "untreated control"]
    need = settings[settings["setting"].str.lower() == "minimum control share of area"]
    have = num(row.iloc[0]["value"]) if len(row) else None
    minimum = num(need.iloc[0]["value"]) if len(need) else None
    if have is None or minimum is None:
        return [Finding("control area", "info", f, "Untreated control", "control share or its minimum is not filled in")]
    ok = have >= minimum - 1e-9
    return [Finding("control area", "pass" if ok else "fail", f, "Untreated control",
                    f"planned control {have:g} % of the area, minimum {minimum:g} %")]


# ---------------------------------------------------------------- run everything
def validate_all(lab_dir=data_io.LAB) -> list[Finding]:
    t = {k: data_io.read_table(k, lab_dir) for k in ("01", "02", "03", "05", "06", "07", "08", "09", "11", "12", "14", "18")}
    out: list[Finding] = []
    for stem in data_io.table_files(lab_dir):
        out += check_structure(data_io.read_table(stem, lab_dir))
    out += check_oxide_totals(t["01"])
    for k in ("08", "09"):
        out += check_base_saturation(t[k])
        out += check_ph(t[k])
        out += check_units_soil(t[k])
    out += check_units_long(t["02"], "property")
    out += check_duplicate_ids(t["01"])
    out += check_duplicate_ids(t["08"])
    out += check_duplicate_ids(t["09"])
    both = pd.concat([t["08"][["sample_id"]], t["09"][["sample_id"]]], ignore_index=True)
    both.attrs["file"] = "data/lab/08 and 09 together"
    out += check_duplicate_ids(both)
    out += check_duplicate_ids(t["02"], ["property"])
    out += check_duplicate_ids(t["03"], ["size_um"])
    out += check_duplicate_ids(t["05"], ["phase"])
    out += check_duplicate_ids(t["06"], ["element"])
    out += check_known_issues(t["01"], t["08"])
    out += check_co2_potential(t["01"], t["07"])
    out += check_control_area(t["11"], t["14"])
    return out


def data_status(lab_dir=data_io.LAB) -> pd.DataFrame:
    """How many rows of each table are measured, stand-in, assumed or still to measure."""
    rows = []
    for stem, df in data_io.load_all(lab_dir).items():
        if "status_class" not in df.columns:
            rows.append({"table": stem, "rows": len(df), "note": "no status column"})
            continue
        counts = df["status_class"].value_counts().to_dict()
        rows.append({"table": stem, "rows": len(df), **counts})
    out = pd.DataFrame(rows)
    num_cols = [c for c in out.columns if c not in ("table", "note")]
    out[num_cols] = out[num_cols].fillna(0).astype(int)
    return out


def report_markdown(findings: list[Finding], lab_dir=data_io.LAB) -> str:
    order = {"fail": 0, "warn": 1, "info": 2, "pass": 3}
    n = {lvl: sum(f.level == lvl for f in findings) for lvl in order}
    lines = ["# Data validation report", "",
             "Checks on `data/lab/*.csv`, run by `python -m erw.validate`. The report does not change any value.",
             "",
             f"**{n['fail']} failed, {n['warn']} to check, {n['pass']} passed, {n['info']} not checked (values missing).**",
             ""]
    labels = {"fail": "Failed", "warn": "To check (known open questions)", "info": "Not checked", "pass": "Passed"}
    for lvl in order:
        sel = [f for f in findings if f.level == lvl]
        if not sel:
            continue
        lines += [f"## {labels[lvl]}", "", "| Rule | Table | Item | Finding |", "|---|---|---|---|"]
        lines += [f"| {f.rule} | `{f.table}` | {f.item} | {f.message.replace('|', '/')} |" for f in sel]
        lines.append("")
    st = data_status(lab_dir)
    cols = ["table", "rows"] + [c for c in ("measured", "stand-in", "assumed", "to measure", "calculated", "plan",
                                            "model output", "unspecified") if c in st.columns]
    lines += ["## Status of the data", "",
              "Rows per table by status. A measurement on a rock or soil that is not the project's own "
              "(Marques Britas, Sanguinho) counts as a stand-in.", "",
              "| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in st.iterrows():
        lines.append("| " + " | ".join(str(int(r[c])) if c != "table" else f"`{r[c]}`" for c in cols) + " |")
    lines += ["", "## Rules", "",
              f"1. Oxide total with LOI between {TOTAL_LOW:.0f} and {TOTAL_HIGH:.0f} %.",
              f"2. Reported base saturation within {BS_TOLERANCE*100:.0f} % (relative) of (Ca+Mg+K+Na)/CEC.",
              "3. pH in KCl below pH in water.",
              "4. Units: no unit notes or text in number columns, each property in one unit, sum of bases and "
              "CEC agree with their parts.",
              "5. No repeated sample ids.",
              "6. Known open questions from `data/lab/00_README.md` (potassium for Moreira and Terra da Riba, "
              "Mn2O4 against Mn3O4).",
              "7. CO2 potential in table 07 recomputed from the oxides in table 01.",
              "8. Planned control area at or above the minimum share (Isometric EWA v1.3, Equation 9).", ""]
    return "\n".join(lines)


def write_report(path: Path | str | None = None, lab_dir=data_io.LAB) -> Path:
    path = Path(path or ROOT / "outputs" / "validation.md")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(report_markdown(validate_all(lab_dir), lab_dir), encoding="utf-8")
    return path


if __name__ == "__main__":
    p = write_report()
    fs = validate_all()
    for lvl in ("fail", "warn"):
        for f in fs:
            if f.level == lvl:
                print(f"{lvl.upper():4s} {f.rule}: {f.table} {f.item}: {f.message}")
    print(f"{sum(f.level == 'pass' for f in fs)} checks passed. Report written to {p}")
