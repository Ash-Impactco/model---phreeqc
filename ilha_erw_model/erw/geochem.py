"""Feedstock geochemistry, weathering kinetics and soil response.

1. Feedstock: cation charge from XRF oxides, CO2 removal potential (bicarbonate pathway,
   adjusted for inorganic C and S), weathering indices (CIA, MIA), mineral mode with
   glass composition solved by difference so the cations match the XRF.
2. Kinetics: shrinking-core dissolution of each mineral in each grain-size class,
   rate law k_acid*aH^n + k_neutral (Palandri and Kharaka 2004), Arrhenius to site
   temperature, reactive surface = geometric surface x roughness, x field factor.
3. Soil: monthly mass balance of released cations in the mixed layer. Cations leave as
   bicarbonate (carbon removal) in proportion to carbonate speciation at soil pH;
   the rest neutralise soil acidity (raising pH and base saturation) or leave with
   nitrate and other strong-acid anions (no removal). Background acid load from
   fertiliser and grazing acts on every plot, including controls.
4. Trace-metal mass balance against screening limits.
"""
import numpy as np
import pandas as pd
from .common import OUT, save_csv, mpl, LINES, DASH

R_GAS = 8.314
M_OX = {"CaO": 56.077, "MgO": 40.304, "Na2O": 61.979, "K2O": 94.196, "Al2O3": 101.961, "Fe2O3": 159.69,
        "SiO2": 60.084, "TiO2": 79.866}
CHG = {"CaO": 2, "MgO": 2, "Na2O": 2, "K2O": 2}       # equivalents of base cation per mol oxide
M_CO2, M_C, M_S = 44.01, 12.011, 32.06
REACTIVE = ["olivine", "augite", "labradorite", "glass"]


# ------------------------------------------------------------------ feedstock
def cation_eq_per_kg(ox):
    """Base-cation charge (eq per kg rock) from oxide wt%."""
    return {k: CHG[k] * ox.get(k, 0) * 10 / M_OX[k] for k in CHG}


def feedstock(cfg):
    ox = cfg.v("feedstock.oxides_pct")
    eq = cation_eq_per_kg(ox)
    tot = sum(eq.values())
    ic, s = cfg.v("feedstock.inorganic_C_pct"), cfg.v("feedstock.total_S_pct")
    ic_co2 = ic * 10 / M_C * M_CO2                           # kg CO2 per t already in the rock
    s_eq = 2 * s * 10 / M_S                                    # eq/kg acidity if sulfide oxidises
    pot_gross = tot * M_CO2                                     # kg CO2/t, 1 mol CO2 per eq
    pot_net = (tot - s_eq) * M_CO2 - ic_co2
    dr = cfg.v("kinetics.downstream_retention")
    # weathering indices (molar)
    mol = {k: ox.get(k, 0) / M_OX[k] for k in M_OX}
    cia = 100 * mol["Al2O3"] / (mol["Al2O3"] + mol["CaO"] + mol["Na2O"] + mol["K2O"])
    mia = 100 * (mol["Al2O3"] + mol["Fe2O3"]) / (mol["Al2O3"] + mol["Fe2O3"] + mol["MgO"] + mol["CaO"] + mol["Na2O"] + mol["K2O"])
    rows = [
        ("Ca charge", eq["CaO"], "eq/kg"), ("Mg charge", eq["MgO"], "eq/kg"), ("Na charge", eq["Na2O"], "eq/kg"),
        ("K charge", eq["K2O"], "eq/kg"), ("Total base-cation charge", tot, "eq/kg"),
        ("Ca + Mg charge only", eq["CaO"] + eq["MgO"], "eq/kg"),
        ("CO2 potential, gross (1 mol CO2 per eq)", pot_gross, "kg CO2/t"),
        ("Deduct: inorganic carbon in rock", ic_co2, "kg CO2/t"),
        ("Deduct: sulfide acidity (2 eq per mol S)", s_eq * M_CO2, "kg CO2/t"),
        ("CO2 potential, net of rock C and S", pot_net, "kg CO2/t"),
        (f"Net potential after downstream retention ({dr})", pot_net * dr, "kg CO2/t"),
        ("Ca + Mg only, net, after retention", ((eq["CaO"] + eq["MgO"] - s_eq) * M_CO2 - ic_co2) * dr, "kg CO2/t"),
        ("Chemical Index of Alteration (CIA)", cia, "fresh basalt ~35 to 45"),
        ("Mafic Index of Alteration (MIA)", mia, "fresh basalt ~35 to 45; higher = weathered"),
        ("Oxide total", sum(ox.values()), "wt%"),
    ]
    tab = pd.DataFrame(rows, columns=["quantity", "value", "unit"])
    tab["value"] = tab["value"].round(3)
    save_csv(tab, "feedstock_potential.csv")
    return tab, {"eq_kg": tot, "eq_kg_CaMg": eq["CaO"] + eq["MgO"], "pot_gross": pot_gross, "pot_net": pot_net,
                 "ic_co2": ic_co2, "s_eq": s_eq, "cia": cia, "mia": mia}


def mineral_table(cfg, fs):
    mode = cfg.v("feedstock.mineral_mode_wt_frac")
    mins = cfg.v("kinetics.minerals")
    cryst = 0.0
    rows = []
    for m in ["olivine", "augite", "labradorite"]:
        p = mins[m]
        eqkg = p["eq_per_mol"] / (p["M_g_mol"] / 1000)
        cryst += mode[m] * eqkg
        rows.append({"mineral": m, "formula": p["formula"], "mode_wt_frac": mode[m], "eq_per_kg_mineral": eqkg,
                     "eq_per_kg_rock": mode[m] * eqkg})
    rem = fs["eq_kg"] - cryst
    flag = ""
    if rem < 0:
        flag = "XRF cations are lower than the crystalline phases imply: check the mode"
        rem = 0.0
    g_eqkg = rem / mode["glass"] if mode["glass"] > 0 else 0.0
    rows.append({"mineral": "glass", "formula": mins["glass"]["formula"], "mode_wt_frac": mode["glass"],
                 "eq_per_kg_mineral": g_eqkg, "eq_per_kg_rock": rem})
    for m in ["oxides_inert", "other_inert"]:
        rows.append({"mineral": m, "formula": "Fe-Ti oxides, apatite etc. (not counted)", "mode_wt_frac": mode[m],
                     "eq_per_kg_mineral": 0.0, "eq_per_kg_rock": 0.0})
    tab = pd.DataFrame(rows)
    tab["share_of_cations"] = tab.eq_per_kg_rock / fs["eq_kg"]
    tab = tab.round(3)
    save_csv(tab, "feedstock_mineral_balance.csv")
    eqkg = {r.mineral: r.eq_per_kg_mineral for r in tab.itertuples()}
    return tab, eqkg, flag


# ------------------------------------------------------------------ kinetics
def rate_mol_m2_s(p, pH, T_C):
    """Far-from-equilibrium dissolution rate, mol per m2 per s."""
    T = T_C + 273.15
    arr = lambda Ea: np.exp(-Ea * 1000 / R_GAS * (1 / T - 1 / 298.15))
    aH = 10.0 ** (-np.asarray(pH, dtype=float))
    return 10 ** p["logk_acid"] * arr(p["Ea_acid"]) * aH ** p["n_acid"] + 10 ** p["logk_neut"] * arr(p["Ea_neut"])


def grain_classes(cfg, fines_scale=1.0, nsub=6):
    """Representative radii (m) and mass weights. fines_scale >1 shifts mass into the two finest bins."""
    bins = cfg.v("feedstock.psd_bins_mm")
    w = np.array(cfg.v("feedstock.psd_mass_frac"), dtype=float)
    w = w / w.sum()
    w2 = w.copy()
    w2[:2] *= fines_scale
    w2 = w2 / w2.sum()
    r, wt = [], []
    for (a, b), wi in zip(bins, w2):
        ds = np.geomspace(a, b, nsub + 1)
        mids = np.sqrt(ds[:-1] * ds[1:])
        r += list(mids / 2 / 1000)            # mm diameter -> m radius
        wt += [wi / nsub] * nsub
    return np.array(r), np.array(wt)


def simulate(cfg, rates, months, eqkg, S=1, P=None, T_C=None, acid_load=None, depth_cm=None, irrigation=None):
    """Monthly mass balance. P is a dict of per-sample parameter arrays (Monte Carlo) or None for central.
    Returns dict of arrays shaped (S, R, months+1)."""
    rates = np.asarray(rates, dtype=float)
    Rn = len(rates)
    mins = cfg.v("kinetics.minerals")
    mode = cfg.v("feedstock.mineral_mode_wt_frac")
    lam = cfg.v("feedstock.roughness")
    P = P or {}
    g = lambda k, default: np.broadcast_to(np.asarray(P.get(k, default), dtype=float), (S,)).copy()
    ff = g("field_factor", cfg.v("kinetics.field_factor")["central"])
    fines = g("fines_scale", 1.0)
    pH0 = g("pH0", cfg.v("soil.pH_H2O"))
    pHBC = g("pHBC", cfg.v("soil.pH_buffer_cmolc_kg_pH"))
    aload = g("acid_load", cfg.v("soil.acid_load_keq_ha_yr") if acid_load is None else acid_load) * 1000  # eq/ha/yr
    Tc = g("T_C", cfg.v("climate.mean_temp_C") if T_C is None else T_C)
    rex = g("exch_ret", cfg.v("soil.exchange_retention_frac"))
    sec = g("secondary", cfg.v("soil.secondary_loss_frac"))
    upt = cfg.v("soil.plant_uptake_frac")
    pKa = cfg.v("soil.pKa1_carbonic")
    bd = cfg.v("soil.bulk_density_g_cm3")
    dep = (cfg.v("soil.mixing_depth_cm") if depth_cm is None else depth_cm) / 100
    Msoil = dep * 1e4 * bd * 1000                                 # kg soil per ha in the mixed layer
    cec = cfg.v("soil.CEC_cmolc_kg")

    # grain classes per sample (fines scaling varies)
    ufines = np.unique(fines)
    gc = {f: grain_classes(cfg, f) for f in ufines}
    nG = len(next(iter(gc.values()))[0])
    r0 = np.stack([gc[f][0] for f in fines])                       # (S, G)
    w = np.stack([gc[f][1] for f in fines])
    dt = 365.25 * 86400 / 12
    retreat = {m: np.zeros((S, Rn, nG)) for m in REACTIVE}
    pH = np.repeat(pH0[:, None], Rn, axis=1)
    out = {k: np.zeros((S, Rn, months + 1)) for k in
           ["pH", "release", "bicarb", "retained", "leached_strong", "frac_dissolved", "base_sat"]}
    out["pH"][:, :, 0] = pH
    out["base_sat"][:, :, 0] = cfg.v("soil.base_saturation_frac")
    cum = {k: np.zeros((S, Rn)) for k in ["release", "bicarb", "retained", "leached_strong"]}
    tot_eq = sum(mode[m] * eqkg[m] for m in REACTIVE)
    prev_rel = np.zeros((S, Rn))
    for t in range(1, months + 1):
        rel_kg = np.zeros((S, Rn))                                  # eq per kg rock released so far
        for m in REACTIVE:
            p = mins[m]
            k = rate_mol_m2_s(p, pH, Tc[:, None])                   # (S, R)
            drdt = k * lam * p["Vm_cm3_mol"] * 1e-6 * ff[:, None]   # m/s
            retreat[m] += (drdt * dt)[:, :, None]
            rem = (np.clip(r0[:, None, :] - retreat[m], 0, None) / r0[:, None, :]) ** 3
            frac_left = (rem * w[:, None, :]).sum(axis=2)
            rel_kg += mode[m] * eqkg[m] * (1 - frac_left)
        d_rel = (rel_kg - prev_rel) * rates[None, :] * 1000         # eq/ha this month
        prev_rel = rel_kg
        avail = d_rel * (1 - sec[:, None] - upt)
        fb = 1 / (1 + 10 ** (pKa - pH))
        bic = avail * fb
        rest = avail - bic
        ret = rest * rex[:, None]
        cum["release"] += d_rel; cum["bicarb"] += bic; cum["retained"] += ret; cum["leached_strong"] += rest - ret
        dpH = (rest - aload[:, None] / 12) / (Msoil * pHBC[:, None] / 100)
        pH = np.clip(pH + dpH, 4.0, 7.8)
        out["pH"][:, :, t] = pH
        for kk in cum:
            out[kk][:, :, t] = cum[kk]
        out["frac_dissolved"][:, :, t] = rel_kg / tot_eq
        out["base_sat"][:, :, t] = np.clip(cfg.v("soil.base_saturation_frac") + cum["retained"] * 100 / Msoil / cec, 0, 1)
    out["Msoil_kg_ha"] = Msoil
    out["dCaMg_cmolc_kg"] = out["retained"] * 100 / Msoil
    return out


def accounting(cfg, rates, sim, fs, month):
    """Isometric-style carbon account per hectare at a given month (central or per sample)."""
    rates = np.asarray(rates, dtype=float)
    dr = cfg.v("kinetics.downstream_retention")
    e = cfg.raw["emissions"]
    emis_t = (e["haul_km"]["value"] * e["truck_kgCO2e_per_tkm"]["value"] + e["loading_kgCO2e_per_t"]["value"]
              + e["spreading_kgCO2e_per_t"]["value"] + e["quarry_allocation_kgCO2e_per_t"]["value"])
    gross = sim["bicarb"][..., month] * dr * M_CO2 / 1e6                     # t CO2 / ha
    ded_ic = rates * fs["ic_co2"] / 1000
    ded_s = rates * fs["s_eq"] * M_CO2 / 1000
    ded_em = rates * emis_t / 1000
    net = gross - ded_ic - ded_s - ded_em
    return {"gross": gross, "ic": ded_ic, "s": ded_s, "emissions": ded_em, "net": net, "emis_kg_per_t": emis_t}


# ------------------------------------------------------------------ Monte Carlo
MC_SPEC = [
    # name, distribution, parameters, meaning
    ("field_factor", "lognormal", (0.0, 0.8), "Lab-to-field rate scaling (median 1, 90% range ~0.27 to 3.7)"),
    ("fines_scale", "uniform", (0.6, 1.6), "Mass in the two finest PSD bins relative to the placeholder"),
    ("pH0", "uniform", (5.0, 5.9), "Starting soil pH"),
    ("pHBC", "uniform", (3.0, 8.0), "Soil pH buffer capacity, cmolc/kg/pH"),
    ("acid_load", "uniform", (0.5, 4.0), "Net background acid load, keq/ha/yr"),
    ("T_C", "uniform", (13.5, 16.5), "Mean temperature, C"),
    ("exch_ret", "uniform", (0.5, 0.9), "Share of non-bicarbonate cations held on exchange sites"),
    ("secondary", "uniform", (0.05, 0.25), "Cations lost to secondary minerals"),
    ("downstream", "uniform", (0.80, 0.92), "Downstream retention of exported alkalinity"),
]


def draw(S, seed=11):
    rng = np.random.default_rng(seed)
    P = {}
    for name, dist, par, _ in MC_SPEC:
        P[name] = rng.lognormal(*par, S) if dist == "lognormal" else rng.uniform(*par, S)
    return P


def monte_carlo(cfg, rates, eqkg, fs):
    S = cfg.v("uncertainty.n_monte_carlo")
    months = cfg.v("design.months")
    P = draw(S)
    sim = simulate(cfg, rates, months, eqkg, S=S, P=P)
    acc = accounting(cfg, rates, sim, fs, months)
    gross = sim["bicarb"][:, :, months] * P["downstream"][:, None] * M_CO2 / 1e6
    net = gross - acc["ic"] - acc["s"] - acc["emissions"]
    rows = []
    for j, r in enumerate(rates):
        if r == 0:
            continue
        for lab, arr in [("frac_of_cations_dissolved", sim["frac_dissolved"][:, j, months]),
                         ("d_pH_vs_control", sim["pH"][:, j, months] - sim["pH"][:, 0, months]),
                         ("d_exch_CaMg_cmolc_kg", sim["dCaMg_cmolc_kg"][:, j, months]),
                         ("gross_CDR_t_ha", gross[:, j]), ("net_CDR_t_ha", net[:, j]),
                         ("net_CDR_kg_per_t_rock", net[:, j] / r * 1000)]:
            q = np.percentile(arr, [10, 50, 90])
            rows.append({"rate_t_ha": r, "quantity": lab, "P10": q[0], "P50": q[1], "P90": q[2]})
    tab = pd.DataFrame(rows).round(4)
    save_csv(tab, "mc_uncertainty_30m.csv")
    # sensitivity: Spearman rank correlation with net CDR per tonne at 40 t/ha
    from scipy.stats import spearmanr
    j40 = list(rates).index(40) if 40 in rates else -1
    y = net[:, j40] / rates[j40] * 1000
    sens = pd.DataFrame([{"parameter": n, "meaning": m, "spearman_rho": spearmanr(P[n], y)[0]}
                         for n, _, _, m in MC_SPEC]).sort_values("spearman_rho", key=np.abs, ascending=False)
    sens["spearman_rho"] = sens.spearman_rho.round(3)
    save_csv(sens, "mc_sensitivity.csv")
    spec = pd.DataFrame([{"parameter": n, "distribution": d, "parameters": str(p), "meaning": m} for n, d, p, m in MC_SPEC])
    save_csv(spec, "mc_input_distributions.csv")
    return tab, sens, y


# ------------------------------------------------------------------ trace metals
def trace_metals(cfg, rates):
    soil = cfg.v("soil.elements")
    rock = cfg.v("feedstock.trace_ppm")
    pH = cfg.v("soil.pH_H2O")
    lim_all = cfg.v("limits.soil_limits_mg_kg")
    cls = "pH<=5.5" if pH <= 5.5 else ("5.5<pH<=7" if pH <= 7 else "pH>7")
    lim = lim_all[cls]
    bd, dep = cfg.v("soil.bulk_density_g_cm3"), cfg.v("soil.mixing_depth_cm") / 100
    Ms = dep * 1e4 * bd                                           # t soil per ha
    rows = []
    for el in ["Ni", "Cr", "Co", "Cu", "Zn", "Pb", "As", "Cd", "Hg", "V"]:
        cs, cr = soil.get(f"{el}_ppm"), rock.get(el)
        row = {"element": el, "soil_mg_kg": cs, "rock_mg_kg": cr, "limit_mg_kg": lim.get(el, np.nan), "limit_class": cls}
        for r in rates:
            if r == 0:
                continue
            row[f"after_{r}_t_ha"] = round((cs * Ms + cr * r) / (Ms + r), 2)
        row["increase_at_80_pct"] = round(100 * ((cs * Ms + cr * 80) / (Ms + 80) / cs - 1), 1)
        L = lim.get(el)
        if L is None:
            row["max_cumulative_t_ha"] = "no limit set"
        elif cs >= L:
            row["max_cumulative_t_ha"] = "background already above limit"
        elif cr <= L:
            row["max_cumulative_t_ha"] = "not limiting"
        else:
            row["max_cumulative_t_ha"] = round(Ms * (L - cs) / (cr - L), 0)
        rows.append(row)
    tab = pd.DataFrame(rows)
    save_csv(tab, "trace_metal_mass_balance.csv")
    return tab


# ------------------------------------------------------------------ figures
def fig_kinetics(cfg, rates, eqkg, fs, years):
    plt = mpl()
    months = years * 12
    ffs = cfg.v("kinetics.field_factor")
    t = np.arange(months + 1) / 12
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.5))
    res = {}
    for i, (lab, f) in enumerate(ffs.items()):
        sim = simulate(cfg, rates, months, eqkg, P={"field_factor": f})
        res[lab] = sim
        j = list(rates).index(40)
        axes[0].plot(t, 100 * sim["frac_dissolved"][0, j], color=LINES[i], ls=DASH[i], label=f"{lab} (field factor {f})")
    axes[0].set_title("Share of feedstock cations released (any rate)")
    axes[0].set_xlabel("Years after spreading"); axes[0].set_ylabel("%")
    axes[0].axvline(2.5, color="#999", lw=0.6); axes[0].legend(fontsize=7)
    sim = res["central"]
    acc_t = lambda m: accounting(cfg, rates, sim, fs, m)["net"][0]
    for j, r in enumerate(rates):
        if r == 0:
            continue
        net = np.array([acc_t(m)[j] for m in range(months + 1)])
        axes[1].plot(t, net, color=LINES[j - 1], ls=DASH[j - 1], label=f"{r} t/ha")
        axes[2].plot(t, sim["pH"][0, j], color=LINES[j - 1], ls=DASH[j - 1], label=f"{r} t/ha")
    axes[2].plot(t, sim["pH"][0, 0], color="#000", lw=0.6, ls=":", label="control")
    axes[1].axhline(0, color="#000", lw=0.5); axes[1].axvline(2.5, color="#999", lw=0.6)
    axes[1].set_title("Net CO2 removal, central case"); axes[1].set_ylabel("t CO2 per ha"); axes[1].set_xlabel("Years")
    axes[1].legend(fontsize=7)
    axes[2].set_title("Soil pH (0-10 cm), central case"); axes[2].set_xlabel("Years"); axes[2].legend(fontsize=7)
    axes[2].axvline(2.5, color="#999", lw=0.6)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig02_kinetics_cdr_pH.png")
    plt.close(fig)
    return res


def fig_minerals(cfg, eqkg):
    plt = mpl()
    mins = cfg.v("kinetics.minerals")
    pHs = np.linspace(4.5, 7.5, 61)
    T = cfg.v("climate.mean_temp_C")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    for i, m in enumerate(REACTIVE):
        ax[0].semilogy(pHs, rate_mol_m2_s(mins[m], pHs, T), color=LINES[i], ls=DASH[i], label=m)
    ax[0].set_xlabel("Soil pH"); ax[0].set_ylabel("mol per m2 per s"); ax[0].legend(fontsize=7)
    ax[0].set_title(f"Lab dissolution rate at {T} C")
    # time to fully dissolve a grain, central
    ff = cfg.v("kinetics.field_factor")["central"]
    lam = cfg.v("feedstock.roughness")
    d_um = np.geomspace(2, 4000, 80)
    for i, m in enumerate(REACTIVE):
        k = rate_mol_m2_s(mins[m], cfg.v("soil.pH_H2O"), T)
        v = k * lam * mins[m]["Vm_cm3_mol"] * 1e-6 * ff * 365.25 * 86400        # m/yr
        ax[1].loglog(d_um, (d_um / 2 * 1e-6) / v, color=LINES[i], ls=DASH[i], label=m)
    ax[1].axhline(2.5, color="#999", lw=0.6)
    ax[1].text(2.2, 2.9, "grant period", fontsize=7, color="#555")
    ax[1].set_xlabel("Grain diameter (micrometres)"); ax[1].set_ylabel("Years to dissolve completely")
    ax[1].set_title("Grain lifetime at starting soil pH, central field factor"); ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig03_mineral_rates.png")
    plt.close(fig)


def fig_mc(tab, sens, y):
    plt = mpl()
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    ax[0].hist(y, bins=40, color="#888", edgecolor="#000", lw=0.3)
    for q, ls in zip([10, 50, 90], [":", "-", ":"]):
        ax[0].axvline(np.percentile(y, q), color="#000", ls=ls, lw=0.8)
    ax[0].set_xlabel("Net CO2 removal at 30 months, kg CO2 per t rock (40 t/ha)")
    ax[0].set_ylabel("Monte Carlo runs"); ax[0].set_title("Uncertainty, 40 t/ha (lines: P10, P50, P90)")
    s = sens.iloc[::-1]
    ax[1].barh(s.parameter, s.spearman_rho, color=["#333" if v > 0 else "#bbb" for v in s.spearman_rho], edgecolor="#000", lw=0.3)
    ax[1].axvline(0, color="#000", lw=0.5)
    ax[1].set_xlabel("Rank correlation with net removal per tonne"); ax[1].set_title("What drives the uncertainty")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig04_uncertainty.png")
    plt.close(fig)
