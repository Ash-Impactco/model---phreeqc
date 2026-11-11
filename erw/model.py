"""
Ilha Carbon ERW soil model (PHREEQC via phreeqpython).

A scenario model, not a validated prediction. 1D soil column with:
  - kinetic dissolution of a normative basalt mineral assemblage (Palandri & Kharaka 2004 rate law)
  - optional basaltic glass (approximate Gislason & Oelkers 2003 rate)
  - Gaines-Thomas cation exchange with a calibrated soil-acidity pool
  - fixed soil pCO2 and O2, secondary phases free to precipitate
  - downward advective transport (steady drainage)
CO2 removal = extra alkalinity (rock column minus identical no-rock control) leaving the column base.

run_stream(cfg) is a generator: it advances both columns one year at a time and yields results,
so a front end can update charts live while PHREEQC is running.
"""
from __future__ import annotations
import json, math, copy
from pathlib import Path
import numpy as np
import phreeqpython

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
INPUTS = json.loads((DATA / "inputs.json").read_text())
DB = (DATA / "llnl.dat").read_text()

MW = dict(SiO2=60.08, Al2O3=101.96, Fe2O3=159.69, CaO=56.08, MgO=40.30, TiO2=79.87, Na2O=61.98, K2O=94.20)
GLASS_FORMULA = ("Si1Al0.349Fe0.165Ca0.227Mg0.2175Na0.1378K0.0538O3.229", 118.05, 1.0806)


def default_config(soil="A_Moreira", psd="as_is_0-4mm"):
    d = copy.deepcopy(INPUTS["defaults"])
    d["isometric"] = copy.deepcopy(INPUTS.get("isometric", {}))
    d["emissions"] = copy.deepcopy(INPUTS.get("emissions", {}))
    d["soil"] = copy.deepcopy(INPUTS["soils"][soil]); d["soil_key"] = soil
    d["oxides"] = copy.deepcopy(INPUTS["feedstock"]["oxides_wt_pct"])
    d["fe2_fraction"] = INPUTS["feedstock"]["fe2_fraction_of_total_fe"]
    d["psd_key"] = psd
    d["psd"] = copy.deepcopy(INPUTS["grain_size_bins"][psd]["bins"])
    d["rates"] = {k: v for k, v in INPUTS["rates_palandri_kharaka_2004"].items() if not k.startswith("_")}
    return d


# ------------------------------------------------------------------ chemistry helpers
def co2_potential(ox):
    m = {k: ox[k] / MW[k] for k in MW}
    cons = (m["CaO"] + m["MgO"]) * 44.01 / 100
    full = (2 * m["CaO"] + 2 * m["MgO"] + 2 * m["Na2O"] + 2 * m["K2O"]) * 44.01 / 100
    return cons, full


def normative_minerals(ox, fe2_fraction=0.8):
    """Cation-conserving norm in g per 100 g rock. Ca, Mg, Na, K are conserved exactly."""
    m = {k: ox[k] / MW[k] for k in MW}
    K = 2 * m["K2O"]; Na = 2 * m["Na2O"]; Ca = m["CaO"]; Mg = m["MgO"]
    Al = 2 * m["Al2O3"]; Ti = m["TiO2"]; Fe = 2 * m["Fe2O3"]
    fe2 = Fe * fe2_fraction; fe3 = Fe - fe2
    kfs = K
    an = max((Al - kfs - Na) / 2, 0)
    plag = Na + an; xan = an / plag if plag else 0.0
    cpx = max(Ca - an, 0)
    fo = max(Mg - cpx, 0) / 2
    ilm = Ti; mag = fe3 / 2
    fa = max(fe2 - ilm - mag, 0) / 2
    mw_plag = 262.2 * (1 - xan) + 278.2 * xan
    g = dict(Plag=plag * mw_plag, Kfs=kfs * 278.3, Cpx=cpx * 216.55, Fo=fo * 140.69, Fa=fa * 203.77,
             Inert=ilm * 151.71 + mag * 231.53)
    s = 100 / sum(g.values())
    g = {k: v * s for k, v in g.items()}
    minerals = {
        "Plag": (f"Na{1-xan:.3f}Ca{xan:.3f}Al{1+xan:.3f}Si{3-xan:.3f}O8", mw_plag, (1 - xan) + 2 * xan),
        "Kfs": ("KAlSi3O8", 278.3, 1.0),
        "Cpx": ("CaMgSi2O6", 216.55, 4.0),
        "Fo": ("Mg2SiO4", 140.69, 4.0),
        "Fa": ("Fe2SiO4", 203.77, 0.0),
    }
    return g, minerals, xan


# ------------------------------------------------------------------ emissions (ASSUMED factors)
def project_emissions_t_ha(cfg):
    """One-off project emissions per hectare (t CO2e/ha) from haulage, milling, spreading and other.
    All factors are placeholders from inputs.json["emissions"] and must be replaced with measured data."""
    e = cfg.get("emissions") or INPUTS["emissions"]
    t_ha = cfg["application_t_ha"]
    trips_per_ha = t_ha / e["truck_payload_t"]
    diesel_haul = trips_per_ha * e["haul_distance_km"] * (e["diesel_l_per_km_loaded"] + e["diesel_l_per_km_empty"])
    diesel_spread = t_ha * e["spreading_diesel_l_per_t"]
    kwh = t_ha * e["milling_kwh_per_t"].get(cfg.get("psd_key", "as_is_0-4mm"), 0)
    out = dict(
        haulage=diesel_haul * e["diesel_kgco2e_per_l"] / 1000,
        spreading=diesel_spread * e["diesel_kgco2e_per_l"] / 1000,
        milling=kwh * e["grid_kgco2e_per_kwh"] / 1000,
        quarrying=t_ha * e["quarrying_kgco2e_per_t"] / 1000,
        other=e["other_project_t_co2e"] / 100.0,
    )
    out["total"] = sum(out.values())
    return out


# ------------------------------------------------------------------ PHREEQC input
def _geometry(cfg):
    ncells = int(round(cfg["column_depth_cm"] / cfg["cell_cm"]))
    rock_cells = max(1, int(round(cfg["mix_depth_cm"] / cfg["cell_cm"])))
    water = cfg["cell_cm"] / 100 * cfg["water_content"] * 1000           # kg water / m2 / cell
    soil_per_kgw = cfg["cell_cm"] / 100 * cfg["bulk_density_kg_m3"] / water
    rock_g_kgw = cfg["application_t_ha"] * 100 / rock_cells / water        # g rock / kg water (t/ha x 100 = g/m2)
    shift_s = water / cfg["drainage_mm_yr"] * 365.25 * 86400
    shifts_per_year = cfg["drainage_mm_yr"] / water
    return ncells, rock_cells, water, soil_per_kgw, rock_g_kgw, shift_s, shifts_per_year


def reactants(cfg):
    g100, minerals, xan = normative_minerals(cfg["oxides"], cfg["fe2_fraction"])
    _, _, _, _, rock_g_kgw, _, _ = _geometry(cfg)
    gf = cfg["glass_frac"]
    react_g100 = sum(v for k, v in g100.items() if k != "Inert")
    out = []
    for k, (formula, mw, chg) in minerals.items():
        out.append((k, formula, mw, chg, g100[k] / 100 * rock_g_kgw * (1 - gf), cfg["rates"][k]))
    if gf > 0:
        f, mw, chg = GLASS_FORMULA
        out.append(("Gls", f, mw, chg, react_g100 / 100 * rock_g_kgw * gf, cfg["rates"]["Gls"]))
    return out


def build_setup(cfg, rock=True):
    s = dict(cfg["soil"])
    ncells, rock_cells, water, soil_per_kgw, rock_g_kgw, shift_s, _ = _geometry(cfg)
    if cfg["acidity_mode"] == "effective":
        s["cec"] = s["ca"] + s["mg"] + s["k"] + cfg["effective_acidity_cmol_kg"]
    beta = dict(Ca=s["ca"] / s["cec"], Mg=s["mg"] / s["cec"], K=s["k"] / s["cec"], Na=0.01)
    beta["H"] = 1 - sum(beta.values())
    LOGK = dict(K=0.7, Na=0.0, Ca=0.8, Mg=0.6)
    aH = 10 ** (-s["pH"]); ca_m = cfg["soil_solution_ca_mol"]
    y = math.sqrt(beta["Ca"] / (10 ** LOGK["Ca"] * ca_m))
    mg_m = beta["Mg"] / (10 ** LOGK["Mg"] * y ** 2)
    k_m = beta["K"] / (10 ** LOGK["K"] * y)
    na_m = beta["Na"] / (10 ** LOGK["Na"] * y)
    logk_h = math.log10(beta["H"] / (aH * y))
    T = cfg["temp_c"]
    L = [f"""
EXCHANGE_SPECIES
  H+ + X- = HX
  log_k {logk_h:.4f}
  -gamma 9.0 0.0
SOLUTION 0 rain
  temp {T}
  pH 5.3
  units mmol/kgw
  Na 0.20 ; Cl 0.23 ; Mg 0.02 ; Ca 0.01 ; K 0.005
SOLUTION 1-{ncells} soil water
  temp {T}
  pH {s['pH']}
  units mmol/kgw
  Na {na_m*1e3:.6f} ; Cl 0.3 ; Ca {ca_m*1e3:.6f} ; Mg {mg_m*1e3:.6f} ; K {k_m*1e3:.6f}
EXCHANGE 1-{ncells}
  X {s['cec'] * soil_per_kgw / 100.0:.6f}
  -equilibrate with solution 1
EQUILIBRIUM_PHASES 1-{ncells}
  CO2(g) {cfg['log_pco2']} 10
  O2(g) -0.68 10
  Gibbsite 0 0
  Kaolinite 0 0
  SiO2(am) 0 0
  Calcite 0 0
  Fe(OH)3 0 0
"""]
    punch_bc = "0"
    bc0 = 0.0
    if rock:
        rates, kin, terms = [], [], []
        for name, formula, mw, chg, g_min, (lka, ea, na, lkn, en) in reactants(cfg):
            bc0 += g_min / mw * chg
            for b, (frac, d) in enumerate(cfg["psd"]):
                g = g_min * frac
                if g <= 0:
                    continue
                mol = g / mw
                area = g * 6.0 / (INPUTS["feedstock"]["density_g_m3"] * d) * cfg["roughness"]
                rn = f"{name}_{b}"
                rates.append(f"""{rn}
  -start
  10 if (m <= 0) then goto 200
  20 R = 8.314e-3
  30 ka = 10^({lka}) * exp(-{ea}/R*(1/TK - 1/298.15))
  40 kn = 10^({lkn}) * exp(-{en}/R*(1/TK - 1/298.15))
  50 rate = (ka*act("H+")^{na} + kn) * {area:.6e} * (m/m0)^0.67 * {cfg['field_factor']}
  60 moles = rate * time
  70 if (moles > m) then moles = m
  200 save moles
  -end""")
                kin.append(f"  {rn}\n    -formula {formula} 1\n    -m0 {mol:.6e}\n    -m {mol:.6e}")
                terms.append(f'KIN("{rn}")*{chg}')
        L.append("RATES\n" + "\n".join(rates))
        L.append(f"KINETICS 1-{rock_cells}\n" + "\n".join(kin) + "\n  -step_divide 10")
        punch_bc = " + ".join(terms)
    L.append(f"""
SELECTED_OUTPUT 1
  -reset false
USER_PUNCH 1
  -headings t_s cell pH alk ca mg na k si caX2 mgX2 kX naX hX bc_left calcite
  10 PUNCH TOTAL_TIME, CELL_NO, -LA("H+"), ALK, TOT("Ca"), TOT("Mg"), TOT("Na"), TOT("K"), TOT("Si")
  20 PUNCH MOL("CaX2"), MOL("MgX2"), MOL("KX"), MOL("NaX"), MOL("HX")
  30 PUNCH {punch_bc}
  40 PUNCH EQUI("Calcite")
END
""")
    return "\n".join(L), bc0


def transport_block(cfg, shifts, first):
    ncells, *_ , shift_s, _ = _geometry(cfg)
    if first:
        return f"""TRANSPORT
  -cells {ncells}
  -shifts {shifts}
  -time_step {shift_s:.6e}
  -flow_direction forward
  -boundary_conditions flux flux
  -lengths {cfg['cell_cm']/100}
  -dispersivities 0.01
  -punch_cells 1-{ncells}
  -punch_frequency 1
  -print_frequency 100000
END"""
    return f"""TRANSPORT
  -shifts {shifts}
END"""


def _instance():
    pp = phreeqpython.PhreeqPython()
    ip = pp.ip
    ip.load_database_string(DB)
    return pp, ip


def _rows(ip):
    arr = ip.get_selected_output_array()
    if not arr or len(arr) < 2:
        return np.zeros((0, 16))
    return np.array([[float(x) if isinstance(x, (int, float)) else np.nan for x in r] for r in arr[1:]])


COLS = ["t_s", "cell", "pH", "alk", "ca", "mg", "na", "k", "si", "caX2", "mgX2", "kX", "naX", "hX", "bc_left", "calcite"]
C = {c: i for i, c in enumerate(COLS)}


def run_stream(cfg):
    """Yield a results dict after each simulated year (rock and control advanced together)."""
    ncells, rock_cells, water, soil_per_kgw, rock_g_kgw, shift_s, spy = _geometry(cfg)
    cons, full = co2_potential(cfg["oxides"])
    sims = {}
    bc0 = 0.0
    for kind in ("rock", "ctrl"):
        pp, ip = _instance()
        setup, b = build_setup(cfg, rock=(kind == "rock"))
        ip.run_string(setup)
        if kind == "rock":
            bc0 = b * rock_cells
        sims[kind] = dict(pp=pp, ip=ip)
    years = int(cfg["years"])
    done_shifts = 0
    cum_export = 0.0
    hist = dict(year=[], export_t=[], released_t=[], exchange_t=[], frac_diss=[], ph_top_rock=[],
                ph_top_ctrl=[], profiles=[])
    effluent = dict(t=[], d_alk=[])
    iso_cfg = cfg.get("isometric") or INPUTS.get("isometric", {})
    nfz_cells = max(1, int(round(iso_cfg.get("nfz_depth_cm", 20) / cfg["cell_cm"])))
    eta = float(iso_cfg.get("eta_river_ocean", 1.0))
    emis = project_emissions_t_ha(cfg)
    cum_export_nfz = 0.0
    hist["iso"] = []
    for yr in range(1, years + 1):
        target = int(round(spy * yr))
        n = max(target - done_shifts, 1)
        last = {}
        for kind, s in sims.items():
            s["ip"].run_string(transport_block(cfg, n, first=(done_shifts == 0)))
            last[kind] = _rows(s["ip"])
        done_shifts += n
        r, c = last["rock"], last["ctrl"]
        # effluent: last cell, every shift
        er = r[r[:, C["cell"]] == ncells]; ec = c[c[:, C["cell"]] == ncells]
        k = min(len(er), len(ec))
        d_alk = er[:k, C["alk"]] - ec[:k, C["alk"]]
        cum_export += float(np.sum(d_alk) * water)                     # eq / m2
        effluent["t"] += list(er[:k, C["t_s"]] / (365.25 * 86400)); effluent["d_alk"] += list(d_alk)
        # export across the Near-Field Zone boundary (Isometric: removal counted when cations leave the NFZ)
        nr = r[r[:, C["cell"]] == nfz_cells]; nc = c[c[:, C["cell"]] == nfz_cells]
        kk = min(len(nr), len(nc))
        cum_export_nfz += float(np.sum(nr[:kk, C["alk"]] - nc[:kk, C["alk"]]) * water)
        # end-of-year state: last ncells rows
        R = r[-ncells:]; Cc = c[-ncells:]
        ex = lambda A: 2 * A[:, C["caX2"]] + 2 * A[:, C["mgX2"]] + A[:, C["kX"]] + A[:, C["naX"]]
        ex_gain = float((ex(R).sum() - ex(Cc).sum()) * water)
        bc_left = float(np.nansum(R[:rock_cells, C["bc_left"]]))
        released = (bc0 - bc_left) * water
        to_t = lambda eq: eq * 44.01 / 100
        hist["year"].append(yr)
        hist["export_t"].append(to_t(cum_export))
        hist["released_t"].append(to_t(released))
        hist["exchange_t"].append(to_t(ex_gain))
        hist["frac_diss"].append(1 - bc_left / bc0 if bc0 else 0.0)
        hist["ph_top_rock"].append(float(np.mean(R[:rock_cells, C["pH"]])))
        hist["ph_top_ctrl"].append(float(np.mean(Cc[:rock_cells, C["pH"]])))
        depth = (np.arange(ncells) + 0.5) * cfg["cell_cm"]
        cec_eq = lambda A: ex(A) + A[:, C["hX"]]
        hist["profiles"].append(dict(
            depth_cm=depth.tolist(),
            pH_rock=R[:, C["pH"]].tolist(), pH_ctrl=Cc[:, C["pH"]].tolist(),
            bs_rock=(ex(R) / cec_eq(R) * 100).tolist(), bs_ctrl=(ex(Cc) / cec_eq(Cc) * 100).tolist(),
            alk_rock=(R[:, C["alk"]] * 1e3).tolist(), alk_ctrl=(Cc[:, C["alk"]] * 1e3).tolist()))
        # ---- Isometric EWA v1.3 style accounting (Equation 4, NFZ), cumulative to this year, t CO2e/ha
        sorb_nfz = float((ex(R[:nfz_cells]).sum() - ex(Cc[:nfz_cells]).sum()) * water)          # eq/m2
        calc_nfz = float((np.nansum(R[:nfz_cells, C["calcite"]]) - np.nansum(Cc[:nfz_cells, C["calcite"]])) * water)  # mol/m2
        rel_nfz = released  # all rock sits inside the NFZ when mix depth <= NFZ depth
        noncarb = rel_nfz - cum_export_nfz - sorb_nfz - 2 * calc_nfz                             # residual, eq/m2
        bracket = to_t(rel_nfz) - 0.0 - to_t(calc_nfz) - to_t(sorb_nfz) - to_t(noncarb)
        stored = bracket * eta
        hist["iso"].append(dict(
            year=yr, weathering=to_t(rel_nfz), biomass_loss=0.0, net_new_carbonate=to_t(calc_nfz),
            net_new_silicate=0.0, net_sorption=to_t(sorb_nfz), non_carbonic_neut=to_t(noncarb),
            exported_from_nfz=to_t(cum_export_nfz), eta_river_ocean=eta, stored=stored,
            counterfactual="subtracted (rock minus no-rock control)", emissions=emis["total"],
            emissions_breakdown=emis, net_removal=stored - emis["total"]))
        chem = {}
        for kind, A in (("rock", er), ("ctrl", ec)):
            row = A[-1]
            chem[kind] = dict(pH=float(row[C["pH"]]), alk=float(row[C["alk"]]), ca=float(row[C["ca"]]),
                              mg=float(row[C["mg"]]), na=float(row[C["na"]]), k=float(row[C["k"]]),
                              si=float(row[C["si"]]))
        yield dict(year=yr, years=years, potential_full_t_ha=full * cfg["application_t_ha"],
                   effluent_chem=chem,
                   potential_cons_t_ha=cons * cfg["application_t_ha"],
                   bc0_t_ha=bc0 * water * 44.01 / 100, hist=hist, effluent=effluent)


def run(cfg):
    out = None
    for out in run_stream(cfg):
        pass
    return out


if __name__ == "__main__":
    import time, sys
    cfg = default_config()
    cfg["years"] = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    t = time.time()
    for res in run_stream(cfg):
        h = res["hist"]
        print(f"year {res['year']}: released {h['released_t'][-1]:.3f}  exported {h['export_t'][-1]:.4f}  "
              f"on exchange {h['exchange_t'][-1]:.3f} tCO2e/ha  dissolved {h['frac_diss'][-1]*100:.2f}%  "
              f"pH {h['ph_top_rock'][-1]:.2f}/{h['ph_top_ctrl'][-1]:.2f}  ({time.time()-t:.0f}s)")


def full_input(cfg, rock=True):
    """A standalone PHREEQC input file (.pqi) for the whole run, to open in the PHREEQC GUI or CLI
    with llnl.dat. Uses the same chemistry as the app; output goes to the selected-output file."""
    ncells, rock_cells, water, soil_per_kgw, rock_g_kgw, shift_s, spy = _geometry(cfg)
    setup, _ = build_setup(cfg, rock=rock)
    setup = setup.replace("SELECTED_OUTPUT 1\n  -reset false",
                          f"SELECTED_OUTPUT 1\n  -file erw_{'rock' if rock else 'control'}.tsv\n  -reset false")
    head = (f"# Ilha Carbon ERW scenario model ({'basalt' if rock else 'no-rock control'} column)\n"
            f"# Database: llnl.dat. Soil: {cfg['soil'].get('label', cfg.get('soil_key',''))}. "
            f"Grain size: {cfg['psd_key']}. {cfg['application_t_ha']} t/ha, {cfg['years']} years.\n")
    return head + setup + "\n" + transport_block(cfg, int(round(spy * cfg["years"])), first=True) + "\n"
