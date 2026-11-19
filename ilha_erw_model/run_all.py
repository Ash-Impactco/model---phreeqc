"""Ilha Carbon, Lomba da Maia basalt rate trial: soil and geological model.

Run everything:            python3 run_all.py
Faster check run:          python3 run_all.py --quick
Analyse real plot results: python3 run_all.py analyse inputs/plot_results.csv

Inputs: config.json (every parameter, with status) and, when they exist, the CSVs in inputs/
(copy the matching file from templates/, fill it, drop it in inputs/). Soil and rock values start
empty in config.json; until every required file is filled, the run lays out the plots and then
stops with a list of what is missing (outputs/missing_inputs.csv). Filled values are marked measured.
Outputs: outputs/ (tables, figures, GeoJSON) and report/Ilha_ERW_Model_Report.pdf
"""
import os, sys, json, time
import numpy as np
import pandas as pd
import geopandas as gpd
from erw.common import Config, ensure_dirs, ROOT, OUT, INP, save_csv
from erw import site, geochem as G, mrv, stats as St


OXIDES = ["SiO2", "Al2O3", "Fe2O3", "CaO", "MgO", "TiO2", "Na2O", "K2O"]
TRACE = ["Zr", "Ni", "Cr", "Co", "Cu", "Zn", "Pb", "As", "Cd", "Hg", "V"]
SOIL_EL = ["CaO_pct", "MgO_pct", "Na2O_pct", "K2O_pct", "TiO2_pct", "Zr_ppm", "Ni_ppm", "Cr_ppm", "Co_ppm",
           "Cu_ppm", "Zn_ppm", "Pb_ppm", "As_ppm", "Cd_ppm", "Hg_ppm", "V_ppm"]
CORE_CV = {"CaO_pct": "CaO_pct", "MgO_pct": "MgO_pct", "TiO2_pct": "TiO2_pct", "Zr_ppm": "Zr_ppm", "pH": "pH_H2O",
           "exch_CaMg": "exch_CaMg_cmolc_kg"}
MINERALS = ["olivine", "augite", "labradorite", "glass", "oxides_inert", "other_inert"]
# site and lab data the model needs before it can run the chemistry: config key -> input file that fills it
REQUIRED = {
    "soil.pH_H2O": "soil_screening.csv", "soil.bulk_density_g_cm3": "soil_screening.csv",
    "soil.CEC_cmolc_kg": "soil_screening.csv", "soil.base_saturation_frac": "soil_screening.csv",
    "soil.pH_buffer_cmolc_kg_pH": "soil_screening.csv", "soil.elements": "soil_screening.csv",
    "soil.between_plot_sd": "soil_screening.csv (4 or more composites)",
    "soil.heterogeneity_cv": "soil_cores.csv", "soil.acid_load_keq_ha_yr": "site_other.csv",
    "soil.analytical_cv": "site_other.csv",
    "feedstock.oxides_pct": "feedstock_xrf.csv", "feedstock.trace_ppm": "feedstock_xrf.csv",
    "feedstock.inorganic_C_pct": "feedstock_other.csv", "feedstock.total_S_pct": "feedstock_other.csv",
    "feedstock.moisture_wet_frac": "feedstock_other.csv", "feedstock.mineral_mode_wt_frac": "feedstock_xrd.csv",
    "feedstock.psd_mass_frac": "feedstock_psd.csv",
}
EXPECTED_KEYS = {"soil.elements": SOIL_EL, "soil.heterogeneity_cv": list(CORE_CV), "soil.between_plot_sd": ["d_pH", "d_CaMg"],
                 "feedstock.oxides_pct": OXIDES, "feedstock.trace_ppm": TRACE, "feedstock.mineral_mode_wt_frac": MINERALS}


def _num(series):
    v = pd.to_numeric(series, errors="coerce")
    return v[v.notna()]


def load_inputs(cfg, log=print):
    """Read every filled CSV in inputs/ and put its values into cfg (marked measured)."""
    used = []
    p = lambda n: os.path.join(INP, n)
    if os.path.exists(p("feedstock_xrf.csv")):
        d = pd.read_csv(p("feedstock_xrf.csv"))
        m = d.assign(value=pd.to_numeric(d.value, errors="coerce")).groupby("analyte").value.mean()
        ox = {k: float(m[k]) for k in OXIDES if k in m and pd.notna(m[k])}
        tr = {k: float(m[k]) for k in TRACE if k in m and pd.notna(m[k])}
        if ox:
            cfg.set("feedstock.oxides_pct", ox)
        if tr:
            cfg.set("feedstock.trace_ppm", tr)
        used.append("feedstock_xrf.csv")
    if os.path.exists(p("feedstock_xrd.csv")):
        d = pd.read_csv(p("feedstock_xrd.csv")).set_index("mineral").wt_frac
        d = _num(d)
        if d.sum() > 0:
            cfg.set("feedstock.mineral_mode_wt_frac", {k: float(v) for k, v in (d / d.sum()).items()})
        used.append("feedstock_xrd.csv")
    if os.path.exists(p("feedstock_psd.csv")):
        d = pd.read_csv(p("feedstock_psd.csv"))
        d["mass_frac"] = pd.to_numeric(d.mass_frac, errors="coerce")
        d = d[d.mass_frac.notna()]
        if d.mass_frac.sum() > 0:
            cfg.set("feedstock.psd_bins_mm", d[["bin_min_mm", "bin_max_mm"]].values.tolist())
            cfg.set("feedstock.psd_mass_frac", (d.mass_frac / d.mass_frac.sum()).tolist())
        used.append("feedstock_psd.csv")
    if os.path.exists(p("feedstock_other.csv")):
        d = pd.read_csv(p("feedstock_other.csv")).set_index("parameter").value
        d = _num(d)
        for k in ["moisture_wet_frac", "inorganic_C_pct", "total_S_pct"]:
            if k in d:
                cfg.set(f"feedstock.{k}", float(d[k]))
        used.append("feedstock_other.csv")
    if os.path.exists(p("soil_screening.csv")):
        d = pd.read_csv(p("soil_screening.csv"))
        if "is_duplicate" in d:
            d = d[d.is_duplicate != 1]
        for col, key in [("pH_H2O", "soil.pH_H2O"), ("bulk_density_g_cm3", "soil.bulk_density_g_cm3"),
                         ("CEC_cmolc_kg", "soil.CEC_cmolc_kg"), ("pH_buffer_cmolc_kg_pH", "soil.pH_buffer_cmolc_kg_pH")]:
            if col in d and _num(d[col]).size:
                cfg.set(key, float(_num(d[col]).mean()))
        if "base_saturation_pct" in d and _num(d.base_saturation_pct).size:
            cfg.set("soil.base_saturation_frac", float(_num(d.base_saturation_pct).mean() / 100))
        el = {k: float(_num(d[k]).mean()) for k in SOIL_EL if k in d and _num(d[k]).size}
        if el:
            cfg.set("soil.elements", el)
        bp = {}
        if "pH_H2O" in d and _num(d.pH_H2O).size > 3:
            bp["d_pH"] = float(0.5 * _num(d.pH_H2O).std())
        if {"exch_Ca_cmolc_kg", "exch_Mg_cmolc_kg"} <= set(d.columns):
            camg = _num(d.exch_Ca_cmolc_kg + d.exch_Mg_cmolc_kg)
            if camg.size > 3:
                bp["d_CaMg"] = float(0.5 * camg.std())
        if bp:
            cfg.set("soil.between_plot_sd", bp, "derived (0.5 x between-composite SD)")
        used.append("soil_screening.csv")
    if os.path.exists(p("soil_cores.csv")):
        d = pd.read_csv(p("soil_cores.csv"))
        cv = {}
        for k, col in CORE_CV.items():
            if col in d and _num(d[col]).size:
                g = d.assign(**{col: pd.to_numeric(d[col], errors="coerce")}).groupby("group")[col]
                cv[k] = float((g.std() / g.mean()).mean())
        if cv:
            cfg.set("soil.heterogeneity_cv", cv)
        used.append("soil_cores.csv")
    if os.path.exists(p("site_other.csv")):
        d = _num(pd.read_csv(p("site_other.csv")).set_index("parameter").value)
        for k, key in [("acid_load_keq_ha_yr", "soil.acid_load_keq_ha_yr"), ("analytical_cv", "soil.analytical_cv"),
                       ("mean_temp_C", "climate.mean_temp_C"), ("rain_mm_yr", "climate.rain_mm_yr"),
                       ("aet_mm_yr", "climate.aet_mm_yr")]:
            if k in d:
                cfg.set(key, float(d[k]))
        used.append("site_other.csv")
    log("  inputs used: " + (", ".join(used) if used else "none yet (inputs/ is empty)"))
    return used


def missing_inputs(cfg):
    """Rows for every required value that is still empty or incomplete."""
    rows = []
    for key, src in REQUIRED.items():
        v = cfg.v(key)
        if v is None:
            rows.append({"parameter": key, "missing": "all", "fill_with": f"inputs/{src}"})
        elif key in EXPECTED_KEYS:
            gone = [k for k in EXPECTED_KEYS[key] if k not in v]
            if gone:
                rows.append({"parameter": key, "missing": ", ".join(gone), "fill_with": f"inputs/{src}"})
    return pd.DataFrame(rows, columns=["parameter", "missing", "fill_with"])


def main(quick=False):
    t0 = time.time()
    ensure_dirs()
    cfg = Config()
    if quick:
        cfg.raw["uncertainty"]["n_monte_carlo"]["value"] = 150
        cfg.raw["uncertainty"]["n_power_sims"]["value"] = 60
    print("Ilha Carbon ERW model")
    used = load_inputs(cfg)
    R = {"inputs_used": used}

    print("1. Site and trial layout")
    plots, R["layout"], _ = site.layout(cfg)
    rates = cfg.v("design.rates_t_ha")
    months = cfg.v("design.months")

    miss = missing_inputs(cfg)
    save_csv(miss, "missing_inputs.csv")
    save_csv(pd.DataFrame(cfg.status_table()), "parameter_status.csv")
    if len(miss):
        print(f"\nStopped after the layout: {len(miss)} soil and rock values are still missing.")
        print("Fill these templates (copy from templates/ to inputs/) and run again:")
        for src, g in miss.groupby("fill_with"):
            print(f"  {src}: {', '.join(g.parameter.str.split('.').str[1])}")
        print("Full list: outputs/missing_inputs.csv. Layout outputs are in outputs/.")
        return R

    print("2. Feedstock geochemistry")
    R["feed_tab"], fs = G.feedstock(cfg)
    R["fs"] = fs
    R["mineral_tab"], eqkg, R["mineral_flag"] = G.mineral_table(cfg, fs)

    print("3. Kinetics, soil response and carbon account")
    G.fig_minerals(cfg, eqkg)
    yrs = cfg.v("kinetics.projection_years")
    scen = G.fig_kinetics(cfg, rates, eqkg, fs, yrs)
    rows = []
    for lab, sim in scen.items():
        for m in sorted({months, 60, yrs * 12}):
            acc = G.accounting(cfg, rates, sim, fs, m)
            for j, r in enumerate(rates):
                rows.append({"scenario": lab, "month": m, "rate_t_ha": r,
                             "frac_cations_dissolved": sim["frac_dissolved"][0, j, m] if r > 0 else 0.0,
                             "pH": sim["pH"][0, j, m], "base_saturation": sim["base_sat"][0, j, m],
                             "d_exch_CaMg_cmolc_kg": sim["dCaMg_cmolc_kg"][0, j, m],
                             "released_keq_ha": sim["release"][0, j, m] / 1000,
                             "bicarbonate_keq_ha": sim["bicarb"][0, j, m] / 1000,
                             "retained_keq_ha": sim["retained"][0, j, m] / 1000,
                             "leached_with_strong_acid_keq_ha": sim["leached_strong"][0, j, m] / 1000,
                             "gross_CDR_t_ha": acc["gross"][0, j], "deduct_rock_C_t_ha": acc["ic"][j],
                             "deduct_S_t_ha": acc["s"][j], "deduct_emissions_t_ha": acc["emissions"][j],
                             "net_CDR_t_ha": acc["net"][0, j]})
    acct = pd.DataFrame(rows).round(4)
    save_csv(acct, "carbon_account_by_scenario.csv")
    R["acct"] = acct
    R["emis_kg_t"] = G.accounting(cfg, rates, scen["central"], fs, months)["emis_kg_per_t"]
    # whole trial totals
    trial = []
    for lab, sim in scen.items():
        acc = G.accounting(cfg, rates, sim, fs, months)
        tot = sum(plots[plots.rate_t_ha == r].gross_ha.sum() * acc["net"][0, j] for j, r in enumerate(rates))
        gro = sum(plots[plots.rate_t_ha == r].gross_ha.sum() * acc["gross"][0, j] for j, r in enumerate(rates))
        trial.append({"scenario": lab, "gross_t_CO2": round(gro, 2), "net_t_CO2": round(tot, 2)})
    R["trial_tot"] = pd.DataFrame(trial)
    save_csv(R["trial_tot"], "trial_total_CDR_30m.csv")

    print("4. Monte Carlo uncertainty")
    R["mc"], R["sens"], y = G.monte_carlo(cfg, rates, eqkg, fs)
    G.fig_mc(R["mc"], R["sens"], y)

    print("5. Trace metals")
    R["metals"] = G.trace_metals(cfg, rates)

    print("6. Measurement design")
    sim_c = scen["central"]
    tau40 = float(sim_c["frac_dissolved"][0, rates.index(40), months])
    R["contrast"] = mrv.tracer_contrast(cfg)
    R["tau"] = mrv.tau_detectability(cfg, tau40, S=300 if quick else 1500)
    mrv.fig_tau(R["tau"])
    water = gpd.read_file(os.path.join(ROOT, cfg.v("site.site_dir"), "water_points.geojson"))
    R["water"] = mrv.water_signal(cfg, sim_c, rates, plots, water)
    R["meso"], R["meso_meta"] = mrv.mesocosm(cfg, eqkg)
    R["sampling"] = mrv.sampling_plan(cfg, plots)

    print("7. Statistics")
    R["power"], R["effects"], R["noise"] = St.power(cfg, eqkg)
    if os.path.exists(os.path.join(INP, "plot_results.csv")):
        R["analysis"] = St.analyse(_prep_results(pd.read_csv(os.path.join(INP, "plot_results.csv"))), "analysis_field")
    ex = St.example_dataset(cfg, eqkg, plots)
    R["example_analysis"] = St.analyse(ex, "analysis_example_simulated")

    status = pd.DataFrame(cfg.status_table())
    save_csv(status, "parameter_status.csv")
    R["status"] = status
    R["plots"] = plots
    R["cfg"] = cfg
    print("8. Report")
    from erw import report
    pdf = report.build(R)
    print(f"Done in {time.time()-t0:.0f} s. Report: {pdf}")
    return R


def _prep_results(df):
    """Accept either changes (d_pH, d_CaMg, d_totCaMg) or raw values per round (pH_H2O, exch_CaMg_cmolc_kg, total_CaMg_eq_kg)."""
    raw = {"pH_H2O": "d_pH", "exch_CaMg_cmolc_kg": "d_CaMg", "total_CaMg_eq_kg": "d_totCaMg"}
    if any(c in df.columns for c in raw) and "month" in df.columns:
        base = df[df.month == df.month.min()].set_index("plot")
        for c, d in raw.items():
            if c in df.columns and d not in df.columns:
                df[d] = df[c] - df["plot"].map(base[c])
    return df


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "analyse":
        ensure_dirs()
        print(St.analyse(_prep_results(pd.read_csv(sys.argv[2])), "analysis_field"))
    else:
        main(quick="--quick" in sys.argv)
