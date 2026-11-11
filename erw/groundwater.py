"""
Groundwater step: MODFLOW 6 + PHREEQC (via mf6rtm / PhreeqcRM), database data/llnl.dat.

Takes the yearly drainage chemistry from the PHREEQC soil column (erw/model.py) and routes it through
a shallow basalt aquifer to a stream, as a vertical cross-section along the slope. It answers:
how long the extra alkalinity takes to reach the stream, how much of it arrives, how much is still
in transit, and how much is lost on the way (calcite in the aquifer, degassing at the stream).

What is modelled
  - Flow: MODFLOW 6, steady, a cross-section of `nlay` layers by `ncol` columns, 1 m wide. Recharge
    enters the top layer; the stream is a fixed head in the top layer at the downslope end.
  - Transport: MODFLOW 6 GWT, one model per chemical component (mf6rtm), with dispersion.
  - Chemistry: PHREEQC in every cell at every step. Calcite and other secondary phases can form.
    Optionally the aquifer basalt weathers slowly, so the untreated aquifer keeps its background
    alkalinity.
  - Recharge chemistry changes every year: each soil year is one stress period.
  - The unsaturated zone below the soil column is not simulated. Its travel time can be set as a
    delay in whole years (piston flow, no reaction).
  - It always runs as a pair: rock (field cells get the basalt-amended drainage) and control
    (every cell gets the control drainage). Every result is rock minus control.

All hydrogeology is ASSUMED until Azorean data are entered in data/lab/21_hydrogeology.csv.

Requires: flopy, mf6rtm, phreeqcrm, modflowapi, and the MODFLOW 6 binaries (mf6 + libmf6) in ./bin
(run scripts/get_modflow.py once).
"""
from __future__ import annotations

import copy
import csv
import json
import math
import platform
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from . import data_io
from .data_io import num, status_class

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
DB_PATH = ROOT / "data" / "llnl.dat"

RESOLVABLE_UEQ = 10.0       # below this, extra alkalinity cannot be told from numerical noise

# ---------------------------------------------------------------- assumed hydrogeology (fallback values)
AQUIFER_DEFAULTS = dict(
    length_m=1500.0,          # upslope edge of the field to the stream
    field_length_m=600.0,     # flow-path length under the treated field
    thickness_m=20.0,         # saturated thickness
    k_m_d=10.0,               # horizontal hydraulic conductivity
    anisotropy=10.0,          # horizontal / vertical conductivity
    porosity=0.15,
    dispersivity_m=10.0,      # longitudinal
    dispersivity_v_m=0.1,     # vertical transverse
    nlay=4, ncol=60, steps_per_year=24,
    top_upslope_m=300.0, top_stream_m=100.0,      # ground elevation, for the drawing only
    vadose_delay_yr=0,
    secondary_phases=["Calcite", "SiO2(am)"],
    aquifer_weathering=False,
    stream_log_pco2=-3.38,
    recharge_mm_yr=900.0, years=10, temp_c=17.0,
)

# Background groundwater (ASSUMED typical dilute basaltic groundwater, mol/kgw): replace with data
BACKGROUND_GW = {"pH": 7.0, "Ca": 3.0e-4, "Mg": 3.0e-4, "Na": 1.0e-3, "K": 1.0e-4,
                 "Cl": 8.0e-4, "Si": 5.0e-4, "Alkalinity": 1.2e-3}

_PARAM_ROWS = {
    "distance from the upslope edge of the field to the stream": ("length_m", float),
    "flow-path length under the treated field": ("field_length_m", float),
    "saturated thickness of the aquifer": ("thickness_m", float),
    "hydraulic conductivity": ("k_m_d", float),
    "ratio of horizontal to vertical conductivity": ("anisotropy", float),
    "effective porosity": ("porosity", float),
    "longitudinal dispersivity": ("dispersivity_m", float),
    "vertical transverse dispersivity": ("dispersivity_v_m", float),
    "number of model layers": ("nlay", int),
    "number of model columns": ("ncol", int),
    "transport steps per year": ("steps_per_year", int),
    "ground elevation at the upslope end": ("top_upslope_m", float),
    "ground elevation at the stream": ("top_stream_m", float),
    "travel time through the unsaturated zone below 50 cm": ("vadose_delay_yr", int),
    "stream co2 for degassing": ("stream_log_pco2", float),
}
_BG_ROWS = {"background groundwater ph": ("pH", 1.0), "background groundwater ca": ("Ca", 1e-3),
            "background groundwater mg": ("Mg", 1e-3), "background groundwater na": ("Na", 1e-3),
            "background groundwater k": ("K", 1e-3), "background groundwater cl": ("Cl", 1e-3),
            "background groundwater si": ("Si", 1e-3), "background groundwater alkalinity": ("Alkalinity", 1e-3)}

# Aquifer basalt for the optional slow weathering. The base cations and silica of the rock (the glass
# formula of the soil model without Al and Fe, which are taken to stay in secondary phases) enter the
# water as bicarbonate: weathering by CO2-bearing water (soil or volcanic CO2). Adding the cations
# without carbon would drive a closed aquifer to pH above 9, which the background water does not show.
AQUIFER_ROCK = ("Ca0.227Mg0.2175Na0.1378K0.0538SiC1.0806O4.7015", 1.0806)   # formula, eq of alkalinity per mole
AQUIFER_ROCK_DEFINITION = """
PHASES
AquiferBasalt
  Ca0.227Mg0.2175Na0.1378K0.0538SiC1.0806O4.7015 + 0.5403 H2O = 0.227 Ca+2 + 0.2175 Mg+2 + 0.1378 Na+ + 0.0538 K+ + SiO2 + 1.0806 HCO3-
  log_k 0
RATES
AquiferBasalt
  -start
  10 moles = parm(1) * time
  20 save moles
  -end
"""


def load_hydrogeology(lab_dir=data_io.LAB):
    """Aquifer parameters and background groundwater from data/lab/21_hydrogeology.csv.
    Returns (params, background, provenance). Missing rows fall back to the assumed defaults."""
    p, bg, prov = copy.deepcopy(AQUIFER_DEFAULTS), dict(BACKGROUND_GW), {}
    try:
        df = data_io.read_table("21", lab_dir)
    except KeyError:
        return p, bg, {k: dict(value=v, status="assumed", source="defaults in erw/groundwater.py", unit="")
                       for k, v in p.items()}
    src = df.attrs["file"]
    for _, r in df.iterrows():
        label, v = r["parameter"].lower(), num(r["value"])
        status = status_class(r["status"])
        if label in _PARAM_ROWS and v is not None:
            key, kind = _PARAM_ROWS[label]
            p[key] = kind(round(v)) if kind is int else float(v)
            prov[key] = dict(value=p[key], status=status, source=f"{src}: {r['parameter']}", unit=r["unit"])
        elif label in _BG_ROWS and v is not None:
            key, factor = _BG_ROWS[label]
            bg[key] = v * factor
            prov[f"background.{key}"] = dict(value=v, status=status, source=f"{src}: {r['parameter']}", unit=r["unit"])
        elif label.startswith("secondary phases") and r["value"]:
            p["secondary_phases"] = [x.strip() for x in r["value"].replace(",", ";").split(";") if x.strip()]
            prov["secondary_phases"] = dict(value="; ".join(p["secondary_phases"]), status=status,
                                            source=f"{src}: {r['parameter']}", unit="")
        elif label.startswith("weathering of the aquifer basalt"):
            p["aquifer_weathering"] = r["value"].strip().lower() in ("on", "yes", "true", "1")
            prov["aquifer_weathering"] = dict(value="on" if p["aquifer_weathering"] else "off", status=status,
                                              source=f"{src}: {r['parameter']}", unit="")
    return p, bg, prov


def provenance_table(lab_dir=data_io.LAB) -> pd.DataFrame:
    """The hydrogeology inputs as a table with their status, for the app and the reports."""
    _, _, prov = load_hydrogeology(lab_dir)
    return pd.DataFrame([{"parameter": k, "value": str(v["value"]), "unit": v.get("unit", ""),
                          "status": v["status"], "source": v["source"]} for k, v in prov.items()])


def load_transect(path) -> list[tuple[float, float]]:
    """Read a topographic transect from the field to the stream: [(distance_m, elevation_m), ...].

    Accepts a CSV with columns distance_m and elevation_m, or a GeoJSON LineString whose coordinates
    carry elevation as the third value ([lon, lat, z]); distance is then measured along the line.
    """
    path = Path(path)
    if path.suffix.lower() in (".geojson", ".json"):
        gj = json.loads(path.read_text())
        geom = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj.get("geometry", gj)
        coords = geom["coordinates"]
        if any(len(c) < 3 for c in coords):
            raise ValueError("the GeoJSON line needs elevation as a third coordinate")
        out, dist = [], 0.0
        for i, (lon, lat, z) in enumerate(c[:3] for c in coords):
            if i:
                lon0, lat0 = coords[i - 1][:2]
                a = (math.sin(math.radians(lat - lat0) / 2) ** 2 + math.cos(math.radians(lat0))
                     * math.cos(math.radians(lat)) * math.sin(math.radians(lon - lon0) / 2) ** 2)
                dist += 2 * 6371000.0 * math.asin(math.sqrt(a))
            out.append((dist, float(z)))
        return out
    with open(path, newline="") as fh:
        return sorted((float(r["distance_m"]), float(r["elevation_m"])) for r in csv.DictReader(fh))


def transect_params(points) -> dict:
    """Aquifer parameters taken from a transect: its length and the ground elevation at both ends."""
    return dict(length_m=points[-1][0] - points[0][0], top_upslope_m=points[0][1], top_stream_m=points[-1][1],
                transect=[(x - points[0][0], z) for x, z in points])


def _bin_files():
    sysname = platform.system().lower()
    lib = {"darwin": "libmf6.dylib", "linux": "libmf6.so", "windows": "libmf6.dll"}[sysname]
    exe = "mf6.exe" if sysname == "windows" else "mf6"
    return BIN / exe, BIN / lib


def check_binaries():
    exe, lib = _bin_files()
    return exe.exists() and lib.exists(), str(exe), str(lib)


# ---------------------------------------------------------------- recharge chemistry from the soil model
def _as_solution(c):
    return {"pH": c["pH"], "Ca": c["ca"], "Mg": c["mg"], "Na": c["na"], "K": c["k"], "Si": c["si"],
            "Cl": c.get("cl", 3.0e-4), "Alkalinity": max(c["alk"], 1e-7)}


def yearly_drainage(soil_result) -> dict[str, list[dict]]:
    """Rock and control drainage chemistry (mol/kgw) for every soil year, as PHREEQC solutions.
    Uses the yearly mean of the water leaving the column base. Older results without it fall back
    to the final-year chemistry repeated for every year."""
    drains = soil_result["hist"].get("drainage")
    if drains:
        return {kind: [_as_solution(d[kind]) for d in drains] for kind in ("rock", "ctrl")}
    ch = soil_result["effluent_chem"]
    n = len(soil_result["hist"]["year"])
    return {kind: [_as_solution(ch[kind])] * n for kind in ("rock", "ctrl")}


def soil_drainage_solutions(soil_result):
    """Final-year rock and control drainage chemistry (kept for older callers)."""
    d = yearly_drainage(soil_result)
    return {kind: d[kind][-1] for kind in ("rock", "ctrl")}


def _delayed(series, baseline, delay):
    """Shift a yearly series by `delay` whole years; before the front arrives the water is `baseline`."""
    delay = max(0, int(delay))
    return ([baseline] * delay + list(series))[:len(series)]


# ---------------------------------------------------------------- one aquifer run
def run_aquifer(recharge_field, recharge_outside, params=None, workdir=None, background=None, label="run"):
    """Build and run one MODFLOW 6 + PHREEQC cross-section.

    recharge_field, recharge_outside: one solution per year (lists of dicts), or a single dict to
    hold the chemistry constant. Returns the cell results over time (`sout`) and the grid."""
    import flopy
    from mf6rtm import mup3d

    p = copy.deepcopy(AQUIFER_DEFAULTS)
    p.update(params or {})
    bg = dict(background or BACKGROUND_GW)
    years = int(p["years"])
    if isinstance(recharge_field, dict):
        recharge_field = [recharge_field] * years
    if isinstance(recharge_outside, dict):
        recharge_outside = [recharge_outside] * years
    years = min(years, len(recharge_field), len(recharge_outside))
    ncol, nlay = int(p["ncol"]), int(p["nlay"])
    delr = p["length_m"] / ncol
    nfield = max(1, min(ncol - 1, int(round(p["field_length_m"] / delr))))
    dz = p["thickness_m"] / nlay
    top, botm = 0.0, [-dz * (k + 1) for k in range(nlay)]
    rch = p["recharge_mm_yr"] / 1000 / 365.25                 # m/d
    nstp = int(p["steps_per_year"])
    wd = Path(workdir or (ROOT / "outputs" / "gw" / label))
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True, exist_ok=True)

    sim = flopy.mf6.MFSimulation(sim_name=label, sim_ws=str(wd), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=years, perioddata=[(365.25, nstp, 1.0)] * years, time_units="days")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    ims = flopy.mf6.ModflowIms(sim, linear_acceleration="BICGSTAB" if nlay > 1 else "CG", filename="gwf.ims")
    sim.register_ims_package(ims, ["gwf"])
    flopy.mf6.ModflowGwfdis(gwf, length_units="meters", nlay=nlay, nrow=1, ncol=ncol, delr=delr, delc=1.0,
                            top=top, botm=botm, filename="gwf.dis")
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, icelltype=0, k=p["k_m_d"],
                            k33=p["k_m_d"] / p["anisotropy"], filename="gwf.npf")
    flopy.mf6.ModflowGwfic(gwf, strt=0.0, filename="gwf.ic")
    flopy.mf6.ModflowGwfsto(gwf, steady_state={0: True}, filename="gwf.sto")
    # the stream: fixed head in the top layer at the downslope end
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, ncol - 1), 0.0]], pname="chd", filename="gwf.chd")
    # recharge into the top layer (m3/d per 1 m wide strip), carrying a tracer aux that mf6rtm
    # replaces with the chemistry of the year
    rec = [[(0, 0, j), rch * delr * 1.0, 1.0] for j in range(ncol - 1)]
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={kper: rec for kper in range(years)}, auxiliary=["tracer"],
                            pname="rch", filename="gwf.rchwel")
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="gwf.hds", budget_filerecord="gwf.cbb",
                           saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")])

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt")
    imt = flopy.mf6.ModflowIms(sim, linear_acceleration="BICGSTAB", filename="gwt.ims")
    sim.register_ims_package(imt, ["gwt"])
    flopy.mf6.ModflowGwtdis(gwt, length_units="meters", nlay=nlay, nrow=1, ncol=ncol, delr=delr, delc=1.0,
                            top=top, botm=botm, filename="gwt.dis")
    flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename="gwt.ic")
    flopy.mf6.ModflowGwtssm(gwt, sources=[["rch", "aux", "tracer"]], filename="gwt.ssm")
    flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, alh=p["dispersivity_m"], ath1=0.0, atv=p["dispersivity_v_m"],
                            filename="gwt.dsp")
    flopy.mf6.ModflowGwtmst(gwt, porosity=p["porosity"], filename="gwt.mst")
    flopy.mf6.ModflowGwtoc(gwt, concentration_filerecord="gwt.ucn", saverecord=[("CONCENTRATION", "LAST")])
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt", filename="gwf-gwt.gwfgwt")

    # solutions: 1 background groundwater, then for each year: field drainage, drainage outside the field
    species = sorted(set(bg) | set(recharge_field[0]) | set(recharge_outside[0]))
    sols = [bg] + [s for y in range(years) for s in (recharge_field[y], recharge_outside[y])]
    solutions = mup3d.Solutions({sp: [s.get(sp, 0.0) for s in sols] for sp in species})
    solutions.set_ic(1)

    # database: llnl.dat itself, or a copy in the run folder with the aquifer basalt added (llnl.dat is not edited)
    database = DB_PATH
    weathering_rate = 0.0
    if p["aquifer_weathering"]:
        src = wd.parent / "_database"                      # mf6rtm copies the database into the run folder
        src.mkdir(parents=True, exist_ok=True)
        database = src / "llnl_with_aquifer_basalt.dat"
        text = DB_PATH.read_text().rstrip()
        if text.upper().endswith("END"):                   # PHREEQC stops reading a database at END
            text = text[:-3].rstrip()
        database.write_text(text + "\n" + AQUIFER_ROCK_DEFINITION + "END\n")
        # constant release that would raise the control recharge to the background alkalinity over one
        # mean residence time (pore water / recharge), so the untreated aquifer keeps its background
        residence_yr = p["thickness_m"] * p["porosity"] / (p["recharge_mm_yr"] / 1000.0)
        gap = max(bg["Alkalinity"] - float(np.mean([s["Alkalinity"] for s in recharge_outside])), 0.0)
        weathering_rate = gap / AQUIFER_ROCK[1] / (residence_yr * 365.25 * 86400)       # mol per kg water per s

    model = mup3d.Mup3d.from_mf6(sim, solutions, name=label, gwt_name="gwt")
    model.set_wd(str(wd))
    model.set_database(str(database))
    eq = mup3d.EquilibriumPhases({1: {ph: {"si": 0.0, "m0": 0.0} for ph in p["secondary_phases"]}})
    eq.set_ic(1)
    model.set_equilibrium_phases(eq)
    if weathering_rate > 0:
        kin = mup3d.KineticPhases({1: {"AquiferBasalt": {"m0": [1000.0], "parms": [weathering_rate],
                                                         "formula": f"{AQUIFER_ROCK[0]} 1"}}})
        kin.set_ic(1)
        model.set_phases(kin)
    model.set_initial_temp(p["temp_c"])
    model.initialize()
    rchchem = mup3d.ChemStress("rch", type="aux")
    rchchem.set_spd({kper: [2 + 2 * kper if j < nfield else 3 + 2 * kper for j in range(ncol - 1)]
                     for kper in range(years)})
    model.set_chem_stress(rchchem)
    model.write_simulation()
    exe, lib = _bin_files()
    for f in (exe, lib):
        shutil.copy2(f, wd / f.name)
    ok = model.run()
    sout = pd.read_csv(wd / "sout.csv")
    p["years"] = years
    return dict(ok=ok, sout=sout, params=p, nfield=nfield, delr=delr, dz=dz, wd=str(wd),
                weathering_rate=weathering_rate)


# ---------------------------------------------------------------- results
def _grid(run, column, t=None):
    """A (nlay, ncol) array of one output column at time t (default: the end)."""
    s = run["sout"]
    t = s["time_d"].max() if t is None else t
    sel = s[np.isclose(s["time_d"], t)]
    p = run["params"]
    out = np.zeros((int(p["nlay"]), int(p["ncol"])))
    out[sel["layer"].values - 1, sel["col"].values - 1] = sel[column].values if column in sel else 0.0
    return out


def stream_degassing(rock, ctrl, t=None):
    """What happens to the extra carbon when the groundwater reaches the open stream.

    The water of the stream cell (rock run and control run) is brought to equilibrium with the
    CO2 of the open air, and calcite may form if the water becomes supersaturated. Returns the
    extra dissolved inorganic carbon left after degassing, per unit of extra alkalinity arriving
    (rock minus control), and the calcite formed. The stream's own flow is ignored, so the water
    is undiluted groundwater: an upper bound for calcite.
    """
    from .model import _instance
    p = rock["params"]
    out = {}
    for kind, run in (("rock", rock), ("ctrl", ctrl)):
        s = run["sout"]
        tt = s["time_d"].max() if t is None else t
        r = s[np.isclose(s["time_d"], tt) & (s["layer"] == 1) & (s["col"] == s["col"].max())].iloc[0]
        tot = lambda el: float(r.get(f"solution_total_molality_{el}", 0.0))
        pp, ip = _instance()
        ip.run_string(f"""SOLUTION 1 stream cell
  temp {p['temp_c']}
  pH {r['solution_ph']}
  units mol/kgw
  Ca {tot('Ca'):.6e} ; Mg {tot('Mg'):.6e} ; Na {tot('Na'):.6e} ; K {tot('K'):.6e}
  Cl {tot('Cl'):.6e} ; Si {tot('Si'):.6e} ; C(4) {tot('C(4)'):.6e}
SELECTED_OUTPUT 1
  -reset false
USER_PUNCH 1
  -headings alk dic pH calcite
  10 PUNCH ALK, TOT("C"), -LA("H+"), EQUI("Calcite")
END
USE solution 1
EQUILIBRIUM_PHASES 1
  CO2(g) {p['stream_log_pco2']} 10
  Calcite 0 0
END""")
        arr = ip.get_selected_output_array()
        before, after = arr[1], arr[-1]
        out[kind] = dict(alk_in=float(before[0]), dic_in=float(before[1]), pH_in=float(before[2]),
                         alk_out=float(after[0]), dic_out=float(after[1]), pH_out=float(after[2]),
                         calcite=float(after[3]))
    d_alk = out["rock"]["alk_in"] - out["ctrl"]["alk_in"]
    res = dict(rock=out["rock"], ctrl=out["ctrl"], d_alk_in_ueq=d_alk * 1e6, resolvable=d_alk * 1e6 >= RESOLVABLE_UEQ)
    if res["resolvable"]:
        d_dic = out["rock"]["dic_out"] - out["ctrl"]["dic_out"]
        d_calc = out["rock"]["calcite"] - out["ctrl"]["calcite"]
        res.update(carbon_kept_per_alk=d_dic / d_alk, calcite_loss_per_alk=2 * d_calc / d_alk,
                   degassing_loss_per_alk=1 - d_dic / d_alk)
    return res


def summarise(rock, ctrl, soil, drainage):
    """Rock minus control: stream breakthrough, where the extra alkalinity is, and the balance."""
    p = rock["params"]
    ncol, nlay, years = int(p["ncol"]), int(p["nlay"]), int(p["years"])
    delr, dz = rock["delr"], rock["dz"]
    field_m2 = rock["nfield"] * delr * 1.0                               # field area per metre of width
    q_stream_m3_d = p["recharge_mm_yr"] / 1000 / 365.25 * delr * (ncol - 1)
    litres = p["porosity"] * delr * 1.0 * dz * 1000.0                    # pore water per cell
    to_t = lambda eq: eq / field_m2 * 44.01 / 100                        # eq per m width -> t CO2e per ha of field

    def stream_series(r):
        s = r["sout"]
        return s[(s["layer"] == 1) & (s["col"] == ncol)].sort_values("time_d")
    a, b = stream_series(rock), stream_series(ctrl)
    n = min(len(a), len(b))
    t_d = a["time_d"].values[:n]
    t = t_d / 365.25
    d_alk = a["solution_alkalinity"].values[:n] - b["solution_alkalinity"].values[:n]          # eq per kg water
    dt_d = np.diff(np.concatenate([[0.0], t_d]))
    cum_stream_eq = np.cumsum(d_alk * 1000.0 * q_stream_m3_d * dt_d)

    # what went in: extra alkalinity in the recharge under the field, year by year
    d_in = np.array([drainage["rock"][y]["Alkalinity"] - drainage["ctrl"][y]["Alkalinity"] for y in range(years)])
    d_in = np.array(_delayed(list(d_in), 0.0, p["vadose_delay_yr"]))
    in_eq_yr = d_in * 1000.0 * p["recharge_mm_yr"] / 1000.0 * field_m2
    cum_in_eq = np.cumsum(in_eq_yr)
    total_in = float(cum_in_eq[-1])

    alk_r, alk_c = _grid(rock, "solution_alkalinity"), _grid(ctrl, "solution_alkalinity")
    storage_eq = float(np.sum(alk_r - alk_c) * litres)
    phases = {}
    for ph in p["secondary_phases"]:
        col = f"equilibrium_phases_moles_{ph}"
        phases[ph] = float(np.sum(_grid(rock, col) - _grid(ctrl, col)) * litres)               # mol per m width
    calcite_eq = 2.0 * phases.get("Calcite", 0.0)
    stream_eq = float(cum_stream_eq[-1])
    residual_eq = total_in - stream_eq - storage_eq - calcite_eq

    resolvable = bool(d_in.max() * 1e6 >= RESOLVABLE_UEQ)
    above = np.where(d_alk * 1e6 >= RESOLVABLE_UEQ)[0]
    first_arrival = float(t[above[0]]) if len(above) else None
    half = np.where(d_alk >= 0.5 * d_alk[-1])[0] if d_alk[-1] > 0 else []
    share = lambda x: (x / total_in) if (resolvable and total_in > 0) else None
    degas = stream_degassing(rock, ctrl)

    frames = dict(t_years=list(range(1, years + 1)), alk_rock=[], alk_ctrl=[], ph_rock=[], ph_ctrl=[], calcite_rock=[])
    times = np.unique(rock["sout"]["time_d"].values)
    for y in range(1, years + 1):
        ty = times[np.argmin(np.abs(times - y * 365.25))]
        frames["alk_rock"].append(_grid(rock, "solution_alkalinity", ty) * 1e3)
        frames["alk_ctrl"].append(_grid(ctrl, "solution_alkalinity", ty) * 1e3)
        frames["ph_rock"].append(_grid(rock, "solution_ph", ty))
        frames["ph_ctrl"].append(_grid(ctrl, "solution_ph", ty))
        frames["calcite_rock"].append(_grid(rock, "equilibrium_phases_moles_Calcite", ty) * 1e3)
    x = (np.arange(ncol) + 0.5) * delr
    field_frac = rock["nfield"] / (ncol - 1)
    mean = lambda g: g.mean(axis=0)                                       # depth-averaged profile
    return dict(
        ok=bool(rock["ok"] and ctrl["ok"]), params=p, drainage=drainage, soil=soil,
        resolvable=resolvable, threshold_ueq=RESOLVABLE_UEQ,
        # time series at the stream
        t_years=t.tolist(), d_alk_ueq=(d_alk * 1e6).tolist(),
        alk_stream_rock_meq=(a["solution_alkalinity"].values[:n] * 1e3).tolist(),
        alk_stream_ctrl_meq=(b["solution_alkalinity"].values[:n] * 1e3).tolist(),
        cum_stream_t_ha=[to_t(v) for v in cum_stream_eq], cum_input_t_ha=[to_t(v) for v in cum_in_eq],
        d_alk_input_ueq_by_year=(d_in * 1e6).tolist(), d_alk_input_ueq=float(d_in[-1] * 1e6),
        first_arrival_yr=first_arrival, half_rise_yr=float(t[half[0]]) if len(half) else None,
        # where the extra alkalinity is at the end, t CO2e per ha of field and as shares of what went in
        balance_t_ha=dict(input=to_t(total_in), stream=to_t(stream_eq), in_transit=to_t(storage_eq),
                          calcite=to_t(calcite_eq), residual=to_t(residual_eq)),
        balance_share=dict(stream=share(stream_eq), in_transit=share(storage_eq), calcite=share(calcite_eq),
                           residual=share(residual_eq)),
        phases_mol_per_m=phases,
        # share that reaches the stream once it has passed through: everything not lost on the way
        reach_fraction=(1 - (calcite_eq + residual_eq) / total_in) if total_in > 0 else float("nan"),
        arrived_fraction=(stream_eq / total_in) if total_in > 0 else float("nan"),
        degassing=degas,
        soil_export_t_ha=soil["hist"]["export_t"][-1], field_end_m=field_m2, field_fraction=field_frac,
        frames=frames, x_m=x.tolist(), z_m=[-(k + 0.5) * dz for k in range(nlay)],
        profile=dict(x_m=x.tolist(), alk_rock=(mean(alk_r) * 1e3).tolist(), alk_ctrl=(mean(alk_c) * 1e3).tolist(),
                     ph_rock=mean(_grid(rock, "solution_ph")).tolist(), ph_ctrl=mean(_grid(ctrl, "solution_ph")).tolist(),
                     calcite_rock=mean(_grid(rock, "equilibrium_phases_moles_Calcite")).tolist()))


def run_coupled(soil_cfg, aquifer_params=None, soil_result=None, label="coupled", lab_dir=data_io.LAB):
    """Soil column (PHREEQC) -> aquifer cross-section (MODFLOW 6 + PHREEQC), rock against control.

    Each soil year's drainage chemistry is one stress period of recharge chemistry. Hydrogeology
    comes from data/lab/21_hydrogeology.csv; `aquifer_params` overrides single values."""
    from .model import run
    soil = soil_result or run(soil_cfg)
    drainage = yearly_drainage(soil)
    p, bg, _ = load_hydrogeology(lab_dir)
    p["years"] = len(drainage["rock"])
    p.update(aquifer_params or {})
    p["recharge_mm_yr"] = soil_cfg["drainage_mm_yr"]
    p["temp_c"] = soil_cfg["temp_c"]
    if "background" in p:
        bg = p.pop("background")
    first_ctrl = drainage["ctrl"][0]
    field_rock = _delayed(drainage["rock"], first_ctrl, p["vadose_delay_yr"])
    outside = _delayed(drainage["ctrl"], first_ctrl, p["vadose_delay_yr"])
    rock = run_aquifer(field_rock, outside, p, background=bg, label=f"{label}_rock")
    ctrl = run_aquifer(outside, outside, p, background=bg, label=f"{label}_ctrl")
    return summarise(rock, ctrl, soil, drainage)
