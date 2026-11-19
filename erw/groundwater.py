"""
Groundwater step: MODFLOW 6 + PHREEQC (via mf6rtm / PhreeqcRM).

Takes the soil drainage chemistry from the PHREEQC soil column (erw/model.py) and routes it
through a simple shallow basalt aquifer under the field to a stream, to see how much of the
extra alkalinity (the CO2 removed) arrives at the stream, how long it takes, and how much is
lost on the way (for example as calcite precipitation or by re-equilibration).

A deliberately simple, clearly labelled ASSUMED hydrogeology: one layer, a 1D flow line from the
upslope edge of the field to a stream. Real values (aquifer thickness, hydraulic conductivity,
porosity, distance to stream, groundwater chemistry) must come from Azorean hydrogeology data.

Requires: flopy, mf6rtm, phreeqcrm, modflowapi, and the MODFLOW 6 binaries (mf6 + libmf6) in ./bin
(run scripts/get_modflow.py once).
"""
from __future__ import annotations
import os, shutil, platform, copy
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
DB_PATH = ROOT / "data" / "llnl.dat"

# ---------------------------------------------------------------- assumed hydrogeology
AQUIFER_DEFAULTS = dict(
    length_m=1500.0,        # field upslope edge to stream (assumed)
    ncol=60,
    field_length_m=600.0,   # length of flow line under the treated field (assumed)
    thickness_m=20.0,       # saturated thickness (assumed)
    k_m_d=10.0,             # hydraulic conductivity, fractured basalt (assumed, range 1-100)
    porosity=0.15,          # effective porosity (assumed)
    dispersivity_m=10.0,    # longitudinal dispersivity (assumed)
    recharge_mm_yr=900.0,   # = drainage below the soil column
    years=10,
    steps_per_year=12,
    temp_c=17.0,
)

# Background groundwater (ASSUMED typical dilute basaltic groundwater, mol/kgw) -- replace with data
BACKGROUND_GW = {"pH": 7.0, "Ca": 3.0e-4, "Mg": 3.0e-4, "Na": 1.0e-3, "K": 1.0e-4,
                 "Cl": 8.0e-4, "Si": 5.0e-4, "Alkalinity": 1.2e-3}


def _bin_files():
    sysname = platform.system().lower()
    lib = {"darwin": "libmf6.dylib", "linux": "libmf6.so", "windows": "libmf6.dll"}[sysname]
    exe = "mf6.exe" if sysname == "windows" else "mf6"
    return BIN / exe, BIN / lib


def check_binaries():
    exe, lib = _bin_files()
    return exe.exists() and lib.exists(), str(exe), str(lib)


def soil_drainage_solutions(soil_result):
    """Extract rock and control drainage chemistry (mol/kgw) from the final year of a soil run.
    soil_result must come from erw.model.run_stream(..., keep_effluent_chem=True)."""
    ch = soil_result["effluent_chem"]
    out = {}
    for kind in ("rock", "ctrl"):
        c = ch[kind]
        out[kind] = {"pH": c["pH"], "Ca": c["ca"], "Mg": c["mg"], "Na": c["na"], "K": c["k"],
                     "Si": c["si"], "Cl": 3.0e-4, "Alkalinity": max(c["alk"], 1e-7)}
    return out


def run_aquifer(recharge_field, recharge_outside, params=None, workdir=None, background=None, label="run"):
    """Build and run one MODFLOW 6 + PHREEQC flow line. Returns a DataFrame of the stream-cell
    chemistry over time and the final profile."""
    import flopy
    from mf6rtm import utils, mup3d

    p = copy.deepcopy(AQUIFER_DEFAULTS); p.update(params or {})
    bg = background or BACKGROUND_GW
    ncol = int(p["ncol"]); delr = p["length_m"] / ncol
    nfield = max(1, int(round(p["field_length_m"] / delr)))
    top, bot = 0.0, -p["thickness_m"]
    rch = p["recharge_mm_yr"] / 1000 / 365.25          # m/d
    perlen = 365.25 * p["years"]; nstp = int(p["years"] * p["steps_per_year"])
    wd = Path(workdir or (ROOT / "outputs" / "gw" / label)); wd.mkdir(parents=True, exist_ok=True)

    sim = flopy.mf6.MFSimulation(sim_name=label, sim_ws=str(wd), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(perlen, nstp, 1.0)], time_units="days")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    ims = flopy.mf6.ModflowIms(sim, linear_acceleration="CG", filename="gwf.ims")
    sim.register_ims_package(ims, ["gwf"])
    flopy.mf6.ModflowGwfdis(gwf, length_units="meters", nlay=1, nrow=1, ncol=ncol, delr=delr, delc=1.0,
                            top=top, botm=bot, filename="gwf.dis")
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, icelltype=0, k=p["k_m_d"], filename="gwf.npf")
    flopy.mf6.ModflowGwfic(gwf, strt=0.0, filename="gwf.ic")
    flopy.mf6.ModflowGwfsto(gwf, steady_state={0: True}, filename="gwf.sto")
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, ncol - 1), 0.0]], pname="chd", filename="gwf.chd")
    # recharge per cell (m3/d per 1 m wide strip), carrying a tracer aux; field cells get solution 2, others 3
    rch_spd = [[(0, 0, j), rch * delr * 1.0, 1.0] for j in range(ncol - 1)]
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data=rch_spd, auxiliary=["tracer"], pname="rch",
                            filename="gwf.rchwel")
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="gwf.hds", budget_filerecord="gwf.cbb",
                           saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")])

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt")
    imt = flopy.mf6.ModflowIms(sim, linear_acceleration="BICGSTAB", filename="gwt.ims")
    sim.register_ims_package(imt, ["gwt"])
    flopy.mf6.ModflowGwtdis(gwt, length_units="meters", nlay=1, nrow=1, ncol=ncol, delr=delr, delc=1.0,
                            top=top, botm=bot, filename="gwt.dis")
    flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename="gwt.ic")
    flopy.mf6.ModflowGwtssm(gwt, sources=[["rch", "aux", "tracer"]], filename="gwt.ssm")
    flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, alh=p["dispersivity_m"], ath1=0.0, filename="gwt.dsp")
    flopy.mf6.ModflowGwtmst(gwt, porosity=p["porosity"], filename="gwt.mst")
    flopy.mf6.ModflowGwtoc(gwt, concentration_filerecord="gwt.ucn", saverecord=[("CONCENTRATION", "LAST")])
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt", filename="gwf-gwt.gwfgwt")

    # solutions: 1 background groundwater, 2 field drainage, 3 drainage outside the field
    species = sorted(set(bg) | set(recharge_field) | set(recharge_outside))
    sol = {s: [bg.get(s, 0.0), recharge_field.get(s, 0.0), recharge_outside.get(s, 0.0)] for s in species}
    solutions = mup3d.Solutions(sol); solutions.set_ic(1)
    eq = mup3d.EquilibriumPhases({1: {"Calcite": {"si": 0.0, "m0": 0.0}}}); eq.set_ic(1)

    model = mup3d.Mup3d.from_mf6(sim, solutions, name=label, gwt_name="gwt")
    model.set_wd(str(wd))
    model.set_database(str(DB_PATH))
    model.set_equilibrium_phases(eq)
    model.set_initial_temp(p["temp_c"])
    model.initialize()
    rchchem = mup3d.ChemStress("rch", type="aux")
    rchchem.set_spd([2 if j < nfield else 3 for j in range(ncol - 1)])
    model.set_chem_stress(rchchem)
    model.write_simulation()
    exe, lib = _bin_files()
    for f in (exe, lib):
        shutil.copy2(f, wd / f.name)
    ok = model.run()
    sout = pd.read_csv(wd / "sout.csv")
    return dict(ok=ok, sout=sout, params=p, nfield=nfield, delr=delr, wd=str(wd))


def stream_excess(rock_run, ctrl_run):
    """Extra alkalinity arriving at the stream (rock minus control), as tCO2/ha of treated field."""
    p = rock_run["params"]; ncol = int(p["ncol"])
    def at_outlet(r):
        s = r["sout"]; col = "col" if "col" in s.columns else None
        s = s[s[col] == s[col].max()] if col else s
        return s
    a, b = at_outlet(rock_run), at_outlet(ctrl_run)
    alk_col = [c for c in a.columns if c.lower().startswith("alk")]
    return a, b, alk_col


def run_coupled(soil_cfg, aquifer_params=None, soil_result=None, label="coupled"):
    """Soil column (PHREEQC) -> aquifer flow line (MODFLOW 6 + PHREEQC).
    The drainage chemistry from the final soil year is applied as constant recharge chemistry
    under the field (rock run) or everywhere (control run). This is a simplification: real
    drainage chemistry changes year by year."""
    from erw.model import run
    soil = soil_result or run(soil_cfg)
    sol = soil_drainage_solutions(soil)
    p = copy.deepcopy(AQUIFER_DEFAULTS); p.update(aquifer_params or {})
    p["recharge_mm_yr"] = soil_cfg["drainage_mm_yr"]; p["temp_c"] = soil_cfg["temp_c"]
    rock = run_aquifer(sol["rock"], sol["ctrl"], p, label=f"{label}_rock")
    ctrl = run_aquifer(sol["ctrl"], sol["ctrl"], p, label=f"{label}_ctrl")
    return summarise(rock, ctrl, soil, sol)


def summarise(rock, ctrl, soil, sol):
    p = rock["params"]; ncol = int(p["ncol"]); delr = rock["delr"]
    q_stream = p["recharge_mm_yr"] / 1000 * delr * (ncol - 1)          # m3/yr reaching stream per m width
    field_area_m2 = rock["nfield"] * delr * 1.0                          # per m width
    def series(r):
        s = r["sout"]
        out = s[s["col"] == s["col"].max() - 1].sort_values("time_d")  # last active cell before the stream
        return out
    a, b = series(rock), series(ctrl)
    n = min(len(a), len(b))
    t = a["time_d"].values[:n] / 365.25
    d_alk = a["solution_alkalinity"].values[:n] - b["solution_alkalinity"].values[:n]    # eq/kgw
    dt = np.diff(np.concatenate([[0.0], t]))
    flux_eq_m2 = d_alk * 1000 * q_stream * dt / field_area_m2                          # eq per m2 of field per step
    cum_t_ha = np.cumsum(flux_eq_m2) * 44.01 / 100
    soil_export = soil["hist"]["export_t"][-1]
    # profiles at end
    def prof(r):
        s = r["sout"]; last = s[s["time_d"] == s["time_d"].max()].sort_values("col")
        return last
    pa, pb = prof(rock), prof(ctrl)
    x = (pa["col"].values - 0.5) * delr
    calc_col = "equilibrium_phases_moles_Calcite"
    field_frac = rock["nfield"] / (ncol - 1)
    d_in = sol["rock"]["Alkalinity"] - sol["ctrl"]["Alkalinity"]
    reach_frac = float(d_alk[-1] / (d_in * field_frac)) if d_in > 0 else float("nan")
    return dict(
        reach_fraction=reach_frac, d_alk_input_ueq=d_in * 1e6, field_fraction=field_frac,
        t_years=t.tolist(), d_alk_ueq=(d_alk * 1e6).tolist(), cum_stream_t_ha=cum_t_ha.tolist(),
        soil_export_t_ha=soil_export, soil=soil, drainage=sol, params=p, field_end_m=field_area_m2,
        profile=dict(x_m=x.tolist(),
                     alk_rock=(pa["solution_alkalinity"].values * 1e3).tolist(),
                     alk_ctrl=(pb["solution_alkalinity"].values * 1e3).tolist(),
                     ph_rock=pa["solution_ph"].values.tolist(), ph_ctrl=pb["solution_ph"].values.tolist(),
                     calcite_rock=(pa[calc_col].values if calc_col in pa else np.zeros(len(pa))).tolist()),
        ok=bool(rock["ok"] and ctrl["ok"]))
