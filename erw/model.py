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

Inputs come from the lab data pack (data/lab/*.csv, read by erw/data_io.py), and every input carries
its status (measured, stand-in, assumed). data/inputs.json is kept only as a fallback and as the
frozen Step 0 reference: default_config(..., source="json") reproduces the original model exactly.
"""
from __future__ import annotations
import json, math, copy
from pathlib import Path
import numpy as np
import phreeqpython

try:
    from . import data_io
except ImportError:                      # run directly as a script: python erw/model.py
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from erw import data_io

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
INPUTS = json.loads((DATA / "inputs.json").read_text())
DB = (DATA / "llnl.dat").read_text()

MW = dict(SiO2=60.08, Al2O3=101.96, Fe2O3=159.69, CaO=56.08, MgO=40.30, TiO2=79.87, Na2O=61.98, K2O=94.20)
GLASS_FORMULA = ("Si1Al0.349Fe0.165Ca0.227Mg0.2175Na0.1378K0.0538O3.229", 118.05, 1.0806)
# Gislason & Oelkers (2003) Geochim. Cosmochim. Acta 67, 3817-3832: A = 10^-5.6 mol Si/cm2/s = 10^-1.6 mol Si/m2/s,
# Ea = 25.5 kJ/mol, per geometric surface area. Used only if data/lab/13 has no values.
GLASS_GO2003 = dict(log_a=-1.6, ea=25.5)


def default_config(soil="A_Moreira", psd="as_is_0-4mm", rock=None, source="lab"):
    """Model configuration for one scenario.

    source="lab"   (default) read everything from data/lab/*.csv, with provenance for each input.
                   `soil` is a sample id from tables 08/09, `rock` a sample id from table 01
                   (default MB-BRITA2), `psd` "as_is_0-4mm" or a grading-curve id from table 03.
                   The old names (A_Moreira, B_TerraRiba, milled_p80_1mm, milled_p80_100um) still work.
    source="json"  the original inputs from data/inputs.json (frozen reference).
    """
    if source == "json":
        return legacy_config(soil, psd)
    return data_io.build_config(soil, psd, rock or data_io.DEFAULT_ROCK)


def legacy_config(soil="A_Moreira", psd="as_is_0-4mm"):
    """The original configuration from data/inputs.json, kept as a fallback and for regression tests."""
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


def mineral_assemblage(cfg):
    """The minerals the model dissolves, in g per 100 g of rock.

    Returns (grams, minerals, xan, source). `grams` has the keys Plag, Kfs, Cpx, Fo, Fa, Gls, Inert.
    If cfg["mineralogy"] holds measured phases (XRD, from data/lab/05), those are used; otherwise
    the norm calculated from the oxides, with the assumed glass fraction taken out of the reactive
    minerals. `source` says which, so the front end can show it.
    """
    g_norm, minerals, xan = normative_minerals(cfg["oxides"], cfg["fe2_fraction"])
    xrd = cfg.get("mineralogy")
    if not xrd:
        gf = cfg["glass_frac"]
        react = sum(v for k, v in g_norm.items() if k != "Inert")
        g = {k: (v if k == "Inert" else v * (1 - gf)) for k, v in g_norm.items()}
        g["Gls"] = react * gf
        src = "normative estimate from the oxides (assumed; no measured XRD)"
        return g, minerals, xan, src + (f", with {gf:.0%} glass assumed" if gf else "")
    ph = dict(xrd["phases"])
    notes = []
    an = xrd.get("plag_an")
    if an is None:
        an = xan
        notes.append(f"plagioclase composition not reported, An{xan*100:.0f} taken from the norm")
    loose = xrd.get("olivine_unsplit_wt_pct", 0.0)
    if loose:
        ol = g_norm["Fo"] + g_norm["Fa"]
        share = g_norm["Fo"] / ol if ol else 1.0
        ph["Fo"] += loose * share
        ph["Fa"] += loose * (1 - share)
        notes.append("olivine composition not reported, split as in the norm")
    total = sum(ph.values())
    g = {k: v * 100.0 / total for k, v in ph.items()}
    mw_plag = 262.2 * (1 - an) + 278.2 * an
    minerals = dict(minerals)
    minerals["Plag"] = (f"Na{1-an:.3f}Ca{an:.3f}Al{1+an:.3f}Si{3-an:.3f}O8", mw_plag, (1 - an) + 2 * an)
    src = "measured phases (XRD), data/lab/05_rock_mineralogy.csv"
    return g, minerals, an, src + ("; " + "; ".join(notes) if notes else "")


def implied_oxides(cfg):
    """CaO, MgO, Na2O and K2O (wt %) carried by the dissolving minerals, to compare with the bulk
    analysis. With the norm they match by construction; with measured XRD the gap shows how well
    the idealised mineral formulas account for the rock chemistry."""
    g, minerals, _, _ = mineral_assemblage(cfg)
    f, mwg, _ = GLASS_FORMULA
    plag_f, mw_plag, _ = minerals["Plag"]
    an = float(plag_f.split("Ca")[1].split("Al")[0])
    mol = dict(Ca=g["Plag"] / mw_plag * an + g["Cpx"] / 216.55 + g["Gls"] / mwg * 0.227,
               Mg=g["Cpx"] / 216.55 + 2 * g["Fo"] / 140.69 + g["Gls"] / mwg * 0.2175,
               Na=g["Plag"] / mw_plag * (1 - an) + g["Gls"] / mwg * 0.1378,
               K=g["Kfs"] / 278.3 + g["Gls"] / mwg * 0.0538)
    return dict(CaO=mol["Ca"] * MW["CaO"], MgO=mol["Mg"] * MW["MgO"], Na2O=mol["Na"] / 2 * MW["Na2O"],
                K2O=mol["K"] / 2 * MW["K2O"])


def exchange_pool(cfg):
    """Size of the exchanger (cmol(+)/kg) and the equivalent fractions on it at the start.

    Two pool options:
      "pH7_CEC"    large pool: bases + exchange acidity at pH 7. Uses the measured acidity where the
                   lab reported it, otherwise CEC at pH 7 minus the bases.
      "effective"  small pool: bases + KCl-exchangeable acidity (measured if available, otherwise the
                   assumed value in the model settings).
    Exchangeable Na is used where measured; otherwise an assumed small share of the CEC.
    """
    s = cfg["soil"]
    na = s.get("na")
    bases = s["ca"] + s["mg"] + s["k"] + (na or 0.0)
    if cfg["acidity_mode"] == "effective":
        cec = bases + cfg["effective_acidity_cmol_kg"]
    elif "acidity_pH7" in s:
        cec = bases + s["acidity_pH7"]
    else:
        cec = s["cec"]
    beta = dict(Ca=s["ca"] / cec, Mg=s["mg"] / cec, K=s["k"] / cec,
                Na=(na / cec) if na is not None else cfg.get("na_fraction_default", 0.01))
    beta["H"] = 1 - sum(beta.values())
    if beta["H"] <= 0 or min(beta.values()) <= 0:
        raise ValueError(f"Soil {s.get('sample_id', '')}: the exchange pool cannot be set up "
                         f"(fractions {beta}). Check CEC, exchangeable cations and the acidity pool.")
    return cec, beta


AL_PHASES = {
    # name in the PHREEQC input, definition if it is not in llnl.dat, label
    "gibbsite": ("Gibbsite", "", "gibbsite (llnl.dat, log K 7.76 at 25 °C)"),
    # llnl.dat carries the Al13 polymer, Al13O4(OH)24+7. At amorphous Al(OH)3 saturation below pH 6 its
    # equilibrium concentration would be molar, which is not what soil solutions show (the polymer forms
    # only in partly neutralised concentrated Al solutions), and phreeqc.dat and wateq4f.dat leave it
    # out. It is switched off for this option only; nothing else in llnl.dat is changed.
    "amorphous": ("Al(OH)3(a)", """SOLUTION_SPECIES
28 H2O + 13 Al+3 = Al13O4(OH)24+7 + 32 H+
  log_k -300
  -llnl_gamma 6
PHASES
Al(OH)3(a)
  Al(OH)3 + 3 H+ = Al+3 + 3 H2O
  log_k 10.8
  delta_h -26.5 kcal
""", "amorphous Al(OH)3 (log K 10.8 at 25 °C, as in phreeqc.dat and wateq4f.dat; Al13 polymer of llnl.dat off)"),
}
LOGK_X = dict(K=0.7, Na=0.0, Ca=0.8, Mg=0.6)       # Gaines-Thomas constants of llnl.dat (and phreeqc.dat)
LOGK_ALX3_DB = 0.41                                 # Al+3 + 3X- = AlX3 in llnl.dat


def _free_al(cfg, s, ca_m, mg_m, k_m, na_m):
    """Free Al3+ (mol/kgw) in the starting soil water when it is in equilibrium with the chosen Al
    phase at the soil pH and temperature. A small PHREEQC speciation, so temperature and
    hydrolysis are handled by the database."""
    name, definition, _ = AL_PHASES[cfg.get("al_phase", "amorphous")]
    pp, ip = _instance()
    ip.run_string(f"""{definition}
SOLUTION 1
  temp {cfg['temp_c']}
  pH {s['pH']}
  units mmol/kgw
  Na {na_m*1e3:.6f} ; Cl 0.3 ; Ca {ca_m*1e3:.6f} ; Mg {mg_m*1e3:.6f} ; K {k_m*1e3:.6f}
  Al 0.001 {name} 0.0
SELECTED_OUTPUT 1
  -reset false
USER_PUNCH 1
  -headings al3 al_total
  10 PUNCH MOL("Al+3"), TOT("Al")
END""")
    row = ip.get_selected_output_array()[-1]
    return float(row[0]), float(row[1])


def andosol_pool(cfg):
    """The Andosol acidity model: exchangeable Al plus pH-dependent charge.

    Two exchangers:
      X  permanent charge, holding Ca, Mg, K, Na and exchangeable Al (AlX3). Soil-solution Al is
         held by an Al hydroxide phase (gibbsite or amorphous Al(OH)3), so when the pH rises the
         exchangeable Al comes off and precipitates, using up 3 OH- per Al.
      Xv variable charge (allophane, imogolite, organic matter): sites that are protonated (HXv,
         uncharged) at low pH and take up cations as the pH rises. One site type with a half-point pH.

    Calibration, all from the soil report:
      - charge at the soil pH  = effective CEC  = bases + KCl-exchangeable acidity
      - charge at pH 7         = CEC at pH 7    (bases + exchange acidity at pH 7 where reported)
      - exchangeable Al        = KCl-exchangeable acidity (measured, otherwise the assumed setting)
    The half-point pH of the variable charge is an assumption (default: midway between soil pH and 7).
    The Gaines-Thomas constant of AlX3 is fitted so that the exchangeable Al is in equilibrium with
    the Al phase at the soil pH; it is reported next to the database value (0.41) as a check.

    Returns a dict with the sizes (cmol(+)/kg), starting fractions, soil-water composition and
    fitted constants. Raises ValueError if the soil data cannot be fitted.
    """
    s = cfg["soil"]
    na = s.get("na")
    al = float(cfg["effective_acidity_cmol_kg"])
    cat = dict(Ca=s["ca"], Mg=s["mg"], K=s["k"])
    bases_no_na = sum(cat.values())
    ecec = bases_no_na + (na or 0.0) + al
    if na is None:                                    # assumed small Na share of the effective CEC
        na_amount = cfg.get("na_fraction_default", 0.01) * ecec
        al = al - na_amount                           # as in the older pools: taken out of the acidity
    else:
        na_amount = na
    cat["Na"] = na_amount
    if al <= 0:
        raise ValueError("Andosol acidity: the KCl-exchangeable acidity must be above zero")
    cec7 = (bases_no_na + (na or 0.0) + s["acidity_pH7"]) if "acidity_pH7" in s else s["cec"]
    vc = cec7 - ecec                                  # charge that appears between the soil pH and pH 7
    half = cfg.get("vc_half_ph") or (s["pH"] + 7.0) / 2
    ca_m = cfg["soil_solution_ca_mol"]
    cH = lambda ph: 10 ** (-ph)
    K = {k: 10 ** v for k, v in LOGK_X.items()}

    def site_y(k_h, conc, ph):
        """Solve the Gaines-Thomas sum of fractions = 1 on Y; return y and the protonated fraction."""
        a2 = K["Ca"] * conc["Ca"] + K["Mg"] * conc["Mg"]
        a1 = K["K"] * conc["K"] + K["Na"] * conc["Na"] + k_h * cH(ph)
        y = (-a1 + math.sqrt(a1 * a1 + 4 * a2)) / (2 * a2)
        return y, k_h * cH(ph) * y

    on_y = {k: 0.0 for k in cat}                      # cations held on Y at the soil pH (cmol/kg)
    s_y, k_hy, f0, f7 = 0.0, 0.0, 0.0, 0.0
    for _ in range(60):
        x_size = sum(cat[k] - on_y[k] for k in cat) + al
        if x_size <= al or any(cat[k] - on_y[k] <= 0 for k in cat):
            raise ValueError("Andosol acidity: the variable charge takes more cations than the soil holds. "
                             "Check ECEC, CEC at pH 7 and the half-point pH.")
        beta = {k: (cat[k] - on_y[k]) / x_size for k in cat}
        y = math.sqrt(beta["Ca"] / (K["Ca"] * ca_m))
        conc = dict(Ca=ca_m, Mg=beta["Mg"] / (K["Mg"] * y ** 2), K=beta["K"] / (K["K"] * y),
                    Na=beta["Na"] / (K["Na"] * y))
        if vc <= 0:
            break
        # half of Y protonated at the half-point pH, in this soil water
        a2 = K["Ca"] * conc["Ca"] + K["Mg"] * conc["Mg"]
        a1 = K["K"] * conc["K"] + K["Na"] * conc["Na"]
        y_half = (-a1 + math.sqrt(a1 * a1 + 2 * a2)) / (2 * a2)          # a2 y^2 + a1 y = 0.5
        k_hy = 0.5 / (cH(half) * y_half)
        y0, h0 = site_y(k_hy, conc, s["pH"])
        _, h7 = site_y(k_hy, conc, 7.0)
        f0, f7 = 1 - h0, 1 - h7
        s_y = vc / (f7 - f0)
        new = dict(Ca=s_y * K["Ca"] * conc["Ca"] * y0 ** 2, Mg=s_y * K["Mg"] * conc["Mg"] * y0 ** 2,
                   K=s_y * K["K"] * conc["K"] * y0, Na=s_y * K["Na"] * conc["Na"] * y0)
        change = max(abs(new[k] - on_y[k]) for k in cat)
        on_y = {k: 0.5 * (on_y[k] + new[k]) for k in cat}                # damped fixed point
        if change < 1e-10:
            break
    beta["Al"] = al / x_size
    al3, al_tot = _free_al(cfg, s, conc["Ca"], conc["Mg"], conc["K"], conc["Na"])
    logk_al = math.log10(beta["Al"] / (al3 * y ** 3))
    al_db = 10 ** LOGK_ALX3_DB * al3 * y ** 3 * x_size                  # what the database constant would hold
    return dict(x_cmol=x_size, y_cmol=s_y, beta_x=beta, conc=conc, logk_hy=math.log10(k_hy) if s_y else None,
                logk_alx3=logk_al, logk_alx3_database=LOGK_ALX3_DB, al_exch_cmol=al, al_exch_database_cmol=al_db,
                al3_mol=al3, al_total_mol=al_tot, ecec=ecec, cec7=cec7, variable_charge=vc, half_ph=half,
                frac_charged_at_soil_ph=f0, frac_charged_at_ph7=f7, bases_on_y=sum(on_y.values()),
                al_phase=AL_PHASES[cfg.get("al_phase", "amorphous")][2])


def reactants(cfg):
    _, _, _, _, rock_g_kgw, _, _ = _geometry(cfg)
    if cfg.get("mineralogy"):
        g, minerals, _, _ = mineral_assemblage(cfg)
        out = [(k, formula, mw, chg, g[k] / 100 * rock_g_kgw, cfg["rates"][k])
               for k, (formula, mw, chg) in minerals.items()]
        if g["Gls"] > 0:
            f, mw, chg = GLASS_FORMULA
            out.append(("Gls", f, mw, chg, g["Gls"] / 100 * rock_g_kgw, cfg["rates"]["Gls"]))
        return out
    g100, minerals, xan = normative_minerals(cfg["oxides"], cfg["fe2_fraction"])
    gf = cfg["glass_frac"]
    react_g100 = sum(v for k, v in g100.items() if k != "Inert")
    out = []
    for k, (formula, mw, chg) in minerals.items():
        out.append((k, formula, mw, chg, g100[k] / 100 * rock_g_kgw * (1 - gf), cfg["rates"][k]))
    if gf > 0:
        f, mw, chg = GLASS_FORMULA
        out.append(("Gls", f, mw, chg, react_g100 / 100 * rock_g_kgw * gf, cfg["rates"]["Gls"]))
    return out


def _classic_setup(cfg, s, ncells, soil_per_kgw, T):
    """The original acidity model: one exchanger X whose acidity is a proton pool (HX), with the
    constant fitted to the measured base saturation at the soil pH."""
    s["cec"], beta = exchange_pool(cfg)
    LOGK = dict(K=0.7, Na=0.0, Ca=0.8, Mg=0.6)
    aH = 10 ** (-s["pH"]); ca_m = cfg["soil_solution_ca_mol"]
    y = math.sqrt(beta["Ca"] / (10 ** LOGK["Ca"] * ca_m))
    mg_m = beta["Mg"] / (10 ** LOGK["Mg"] * y ** 2)
    k_m = beta["K"] / (10 ** LOGK["K"] * y)
    na_m = beta["Na"] / (10 ** LOGK["Na"] * y)
    logk_h = math.log10(beta["H"] / (aH * y))
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
    return L, 'MOL("CaX2"), MOL("MgX2"), MOL("KX"), MOL("NaX"), MOL("HX")', "0"


def build_setup(cfg, rock=True):
    s = dict(cfg["soil"])
    ncells, rock_cells, water, soil_per_kgw, rock_g_kgw, shift_s, _ = _geometry(cfg)
    T = cfg["temp_c"]
    if cfg["acidity_mode"] == "andosol":
        a = andosol_pool(cfg)
        c = a["conc"]
        phase, phase_def, _ = AL_PHASES[cfg.get("al_phase", "amorphous")]
        y_block = ""
        if a["y_cmol"] > 0:
            y_block = f"""EXCHANGE_MASTER_SPECIES
  Xv Xv-
EXCHANGE_SPECIES
  Xv- = Xv-
  log_k 0
  Na+ + Xv- = NaXv
  log_k {LOGK_X['Na']}
  -llnl_gamma 4
  K+ + Xv- = KXv
  log_k {LOGK_X['K']}
  -llnl_gamma 3
  Ca+2 + 2Xv- = CaXv2
  log_k {LOGK_X['Ca']}
  -llnl_gamma 6
  Mg+2 + 2Xv- = MgXv2
  log_k {LOGK_X['Mg']}
  -llnl_gamma 8
  H+ + Xv- = HXv
  log_k {a['logk_hy']:.4f}
  -gamma 9.0 0.0
"""
        L = [f"""
{phase_def}{y_block}EXCHANGE_SPECIES
  Al+3 + 3X- = AlX3
  log_k {a['logk_alx3']:.4f}
  -llnl_gamma 9
SOLUTION 0 rain
  temp {T}
  pH 5.3
  units mmol/kgw
  Na 0.20 ; Cl 0.23 ; Mg 0.02 ; Ca 0.01 ; K 0.005
SOLUTION 1-{ncells} soil water
  temp {T}
  pH {s['pH']}
  units mmol/kgw
  Na {c['Na']*1e3:.6f} ; Cl 0.3 ; Ca {c['Ca']*1e3:.6f} ; Mg {c['Mg']*1e3:.6f} ; K {c['K']*1e3:.6f}
  Al 0.001 {phase} 0.0
EXCHANGE 1-{ncells}
  X {a['x_cmol'] * soil_per_kgw / 100.0:.6f}
{f"  Xv {a['y_cmol'] * soil_per_kgw / 100.0:.6f}" if a['y_cmol'] > 0 else ""}
  -equilibrate with solution 1
EQUILIBRIUM_PHASES 1-{ncells}
  CO2(g) {cfg['log_pco2']} 10
  O2(g) -0.68 10
  {phase} 0 {cfg.get('al_reservoir_mol_kgw', 1.0)}
  Kaolinite 0 0
  SiO2(am) 0 0
  Calcite 0 0
  Fe(OH)3 0 0
"""]
        y_on = a["y_cmol"] > 0
        punch_ex = ('MOL("CaX2") + MOL("CaXv2"), MOL("MgX2") + MOL("MgXv2"), MOL("KX") + MOL("KXv"), '
                    'MOL("NaX") + MOL("NaXv"), MOL("HXv") + 3*MOL("AlX3") + 2*MOL("AlOHX2")') if y_on else \
                   'MOL("CaX2"), MOL("MgX2"), MOL("KX"), MOL("NaX"), 3*MOL("AlX3") + 2*MOL("AlOHX2")'
        punch_al = '3*MOL("AlX3") + 2*MOL("AlOHX2")'
    else:
        L, punch_ex, punch_al = _classic_setup(cfg, s, ncells, soil_per_kgw, T)
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
                geo_area = g * 6.0 / (cfg.get("density_g_m3", INPUTS["feedstock"]["density_g_m3"]) * d)
                area = geo_area * cfg["roughness"]
                rn = f"{name}_{b}"
                if name == "Gls" and cfg.get("glass_rate_law", "approximate") == "gislason_oelkers_2003":
                    go = cfg.get("glass_go") or GLASS_GO2003
                    al_phase = AL_PHASES[cfg.get("al_phase", "amorphous")][0] if cfg["acidity_mode"] == "andosol" else "Gibbsite"
                    rates.append(f"""{rn}
  -start
  # Gislason & Oelkers (2003): r = A exp(-Ea/RT) (aH^3/aAl)^(1/3), mol Si per m2 of geometric surface per s
  10 if (m <= 0) then goto 200
  20 R = 8.314e-3
  30 aal = act("Al+3")
  # no dissolved Al yet (first moment of the run): use the activity at saturation with the soil's Al phase
  40 if (aal < 1e-30) then aal = 10^(LK_PHASE("{al_phase}") + 3*LA("H+"))
  50 k = 10^({go['log_a']}) * exp(-{go['ea']}/(R*TK))
  60 rate = k * (act("H+")^3/aal)^(1/3) * {geo_area:.6e} * (m/m0)^0.67 * {cfg['field_factor']}
  70 moles = rate * time
  80 if (moles > m) then moles = m
  200 save moles
  -end""")
                    kin.append(f"  {rn}\n    -formula {formula} 1\n    -m0 {mol:.6e}\n    -m {mol:.6e}")
                    terms.append(f'KIN("{rn}")*{chg}')
                    continue
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
  -headings t_s cell pH alk ca mg na k si caX2 mgX2 kX naX hX bc_left calcite alX
  10 PUNCH TOTAL_TIME, CELL_NO, -LA("H+"), ALK, TOT("Ca"), TOT("Mg"), TOT("Na"), TOT("K"), TOT("Si")
  20 PUNCH {punch_ex}
  30 PUNCH {punch_bc}
  40 PUNCH EQUI("Calcite")
  50 PUNCH {punch_al}
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
        return np.zeros((0, len(COLS)))
    return np.array([[float(x) if isinstance(x, (int, float)) else np.nan for x in r] for r in arr[1:]])


COLS = ["t_s", "cell", "pH", "alk", "ca", "mg", "na", "k", "si", "caX2", "mgX2", "kX", "naX", "hX", "bc_left", "calcite", "alX"]
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
    pool = andosol_pool(cfg) if cfg["acidity_mode"] == "andosol" else None
    cec7_eq = pool["cec7"] * soil_per_kgw / 100.0 if pool else None      # eq per kg water, per cell
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
        # base saturation: against CEC at pH 7 in the Andosol mode, against the model exchanger otherwise
        cec_eq = (lambda A: np.full(len(A), cec7_eq)) if pool else (lambda A: ex(A) + A[:, C["hX"]])
        hist["profiles"].append(dict(
            depth_cm=depth.tolist(),
            pH_rock=R[:, C["pH"]].tolist(), pH_ctrl=Cc[:, C["pH"]].tolist(),
            bs_rock=(ex(R) / cec_eq(R) * 100).tolist(), bs_ctrl=(ex(Cc) / cec_eq(Cc) * 100).tolist(),
            alk_rock=(R[:, C["alk"]] * 1e3).tolist(), alk_ctrl=(Cc[:, C["alk"]] * 1e3).tolist(),
            exal_rock=(R[:, C["alX"]] / soil_per_kgw * 100).tolist(),          # exchangeable Al, cmol(+)/kg
            exal_ctrl=(Cc[:, C["alX"]] / soil_per_kgw * 100).tolist()))
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
                   effluent_chem=chem, acidity_pool=pool,
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
