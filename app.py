"""
Ilha Carbon: live PHREEQC enhanced rock weathering model.
Run:  streamlit run app.py
"""
import copy, io, json, time
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from erw.model import (default_config, run_stream, co2_potential, mineral_assemblage, full_input, andosol_pool)
from erw import data_io

st.set_page_config(page_title="Ilha Carbon ERW model", layout="wide")

TEAL, GREY, ORANGE, INK, RED = "#2F7D6B", "#9AA19D", "#C9733A", "#1E2326", "#B23A48"
PALETTE = ["#2F7D6B", "#C9733A", "#3B6EA8", "#B23A48", "#6B5CA5", "#8A918D"]

st.markdown("## Basalt weathering on Azorean pasture: live PHREEQC model")
st.caption("Scenario model for discussion, not a validated prediction. Rock and soil values are stand-ins until "
           "TECHNOVIA's rock analysis and the Lomba da Maia soil screening are available. Inputs are read from "
           "the lab files in data/lab, each with its status. Database: llnl.dat.")

if "runs" not in st.session_state:
    st.session_state.runs = {}

# ------------------------------------------------------------------ sidebar inputs
PSD_LABELS = {"as_is_0-4mm": "0/4 mm fines as supplied", "milled_p80_1mm": "Milled, P80 1 mm (assumed curve)",
              "milled_p80_100um": "Milled, P80 100 µm (assumed curve)"}


def with_value(options, value):
    """A slider's option list, with the lab value added if it is not one of the usual steps."""
    return sorted(set(options) | {value})


with st.sidebar:
    st.header("Scenario")
    rocks = data_io.rock_samples()
    usable_rocks = [k for k, v in rocks.items() if v["usable"]]
    rock_key = st.selectbox("Rock (table 01)", usable_rocks,
                            index=usable_rocks.index(data_io.DEFAULT_ROCK) if data_io.DEFAULT_ROCK in usable_rocks else 0,
                            format_func=lambda k: f"{rocks[k]['label']} [{rocks[k]['status_class']}]")
    waiting = [k for k, v in rocks.items() if not v["usable"]]
    if waiting:
        st.caption(f"Waiting for results: {', '.join(waiting)}. They can be chosen once the oxides are filled in.")
    all_soils = data_io.soils()
    soil_key = st.selectbox("Soil (tables 08 and 09)", list(all_soils),
                            format_func=lambda k: f"{all_soils[k]['label']} [{all_soils[k]['status_class']}]")
    extra_curves = [c for c in data_io.psd_ids() if c not in rocks and c not in data_io.PSD_ALIASES.values()]
    psd_key = st.selectbox("Grain size (table 03)", list(PSD_LABELS) + extra_curves,
                           format_func=lambda k: PSD_LABELS.get(k, k))
    cfg = default_config(soil_key, psd_key, rock_key)
    base = copy.deepcopy(cfg)

    with st.expander("Soil chemistry (0-20 cm)", expanded=False):
        s = cfg["soil"]
        s["pH"] = st.number_input("pH (H2O)", 3.5, 9.0, float(s["pH"]), 0.1)
        s["cec"] = st.number_input("CEC at pH 7, cmol(+)/kg", 1.0, 120.0, float(s["cec"]), 0.5)
        s["ca"] = st.number_input("Exchangeable Ca, cmol(+)/kg", 0.0, 60.0, float(s["ca"]), 0.1)
        s["mg"] = st.number_input("Exchangeable Mg, cmol(+)/kg", 0.0, 30.0, float(s["mg"]), 0.1)
        s["k"] = st.number_input("Exchangeable K, cmol(+)/kg", 0.0, 15.0, float(s["k"]), 0.1)
        modes = ["pH7_CEC", "effective", "andosol"]
        cfg["acidity_mode"] = st.radio(
            "Soil acidity pool", modes, index=modes.index(cfg["acidity_mode"]),
            format_func=lambda m: {"pH7_CEC": "Large: bases plus acidity at pH 7",
                                   "effective": "Small: bases plus KCl-exchangeable acidity",
                                   "andosol": "Andosol: exchangeable Al plus pH-dependent charge"}[m])
        if cfg["acidity_mode"] == "pH7_CEC":
            st.caption("Uses the measured exchange acidity at pH 7." if "acidity_pH7" in s
                       else "No measured acidity at pH 7 for this soil: CEC at pH 7 minus the bases is used.")
        else:
            cfg["effective_acidity_cmol_kg"] = st.number_input(
                "KCl-exchangeable acidity, cmol(+)/kg", 0.05, 20.0, float(cfg["effective_acidity_cmol_kg"]), 0.1)
            st.caption("Measured value for this soil." if "acidity_kcl" in s
                       else "Not measured for this soil: an assumed value. Ask the lab for 1 M KCl acidity.")
        if cfg["acidity_mode"] == "andosol":
            phases = ["amorphous", "gibbsite"]
            cfg["al_phase"] = st.selectbox("Al phase holding soil-water Al", phases, index=phases.index(cfg["al_phase"]),
                                           format_func=lambda k: {"amorphous": "Amorphous Al(OH)3",
                                                                  "gibbsite": "Gibbsite"}[k])
            half = st.number_input("Half-point pH of the variable charge", 4.5, 8.0,
                                   float(cfg.get("vc_half_ph") or (s["pH"] + 7) / 2), 0.1,
                                   help="Assumed. Default: midway between the soil pH and 7")
            cfg["vc_half_ph"] = half
        cfg["bulk_density_kg_m3"] = st.number_input("Bulk density, kg/m³", 300, 1600, int(cfg["bulk_density_kg_m3"]), 50)

    with st.expander("Rock chemistry (wt%)", expanded=False):
        for ox in cfg["oxides"]:
            cfg["oxides"][ox] = st.number_input(ox, 0.0, 80.0, float(cfg["oxides"][ox]), 0.01, format="%.2f")
        if cfg.get("mineralogy"):
            st.caption("Glass fraction comes from the measured mineralogy (table 05).")
        else:
            cfg["glass_frac"] = st.slider("Basaltic glass fraction", 0.0, 0.6, float(cfg["glass_frac"]), 0.05,
                                          help="Unknown until XRD with amorphous quantification")
        laws = ["gislason_oelkers_2003", "approximate"]
        cfg["glass_rate_law"] = st.selectbox(
            "Glass rate law", laws, index=laws.index(cfg["glass_rate_law"]),
            format_func=lambda k: {"gislason_oelkers_2003": "Gislason & Oelkers (2003), with Al3+ dependence",
                                   "approximate": "Older approximation (constant rate)"}[k],
            help="The full law slows down when dissolved Al is high, as it is where amorphous Al(OH)3 holds the Al")

    st.subheader("Field")
    cfg["application_t_ha"] = st.slider("Application rate, t/ha", 5, 100, int(cfg["application_t_ha"]), 5)
    cfg["drainage_mm_yr"] = st.slider("Drainage below roots, mm/yr", 300, 1800, int(cfg["drainage_mm_yr"]), 50)
    cfg["seasonal_drainage"] = st.checkbox(
        "Seasonal drainage (assumed monthly pattern)", value=bool(cfg["seasonal_drainage"]),
        help="Spreads the yearly drainage over the months using the shares in table 12 (75 % from October "
             "to March, an assumption until station data are entered). The run starts in the month set in table 14.")
    cfg["log_pco2"] = st.slider("Soil log pCO2", -3.0, -1.0, float(cfg["log_pco2"]), 0.1)
    cfg["temp_c"] = st.slider("Soil temperature, °C", 10.0, 22.0, float(cfg["temp_c"]), 0.5)

    st.subheader("Kinetics")
    cfg["field_factor"] = st.select_slider("Field-to-lab rate factor",
                                           with_value([0.01, 0.03, 0.1, 0.3, 1.0], cfg["field_factor"]),
                                           cfg["field_factor"],
                                           help="Field weathering is usually 10-1000x slower than lab rates")
    cfg["roughness"] = st.select_slider("Surface roughness (BET/geometric)",
                                        with_value([1.0, 2.0, 5.0, 10.0, 20.0], cfg["roughness"]), cfg["roughness"])
    cfg["years"] = st.slider("Years to simulate", 1, 20, int(cfg["years"]))
    name = st.text_input("Run name", f"{cfg['soil_key']} · {psd_key} · glass {cfg['glass_frac']:.0%} · "
                                     f"ff {cfg['field_factor']}")
    go_btn = st.button("Run PHREEQC", type="primary", use_container_width=True)

# Anything changed by hand in the sidebar is no longer the lab value: say so in the provenance.
EDITED = "assumed (edited in the app)"
for key in ("acidity_mode", "al_phase", "glass_rate_law", "seasonal_drainage", "effective_acidity_cmol_kg", "bulk_density_kg_m3", "glass_frac", "application_t_ha",
            "drainage_mm_yr", "log_pco2", "temp_c", "field_factor", "roughness", "years"):
    if cfg[key] != base[key] and not (isinstance(cfg[key], (int, float)) and abs(cfg[key] - base[key]) < 1e-9):
        cfg["provenance"][key] = dict(cfg["provenance"].get(key, {}), value=cfg[key], status=EDITED,
                                      source=f"sidebar (lab value {base[key]})")
if any(abs(cfg["soil"][k] - base["soil"][k]) > 1e-9 for k in ("pH", "cec", "ca", "mg", "k")):
    cfg["provenance"]["soil"] = dict(cfg["provenance"]["soil"], status=EDITED,
                                     value={k: cfg["soil"][k] for k in ("pH", "cec", "ca", "mg", "k")},
                                     source=f"sidebar, starting from {cfg['soil_key']}")
if any(abs(cfg["oxides"][k] - base["oxides"][k]) > 1e-9 for k in cfg["oxides"]):
    cfg["provenance"]["oxides"] = dict(cfg["provenance"]["oxides"], status=EDITED, value=dict(cfg["oxides"]),
                                       source=f"sidebar, starting from {cfg['rock_key']}")

# ------------------------------------------------------------------ static info
cons, full = co2_potential(cfg["oxides"])
g100, _, xan, mineral_source = mineral_assemblage(cfg)
c1, c2, c3 = st.columns([1, 1, 2])
c1.metric("CO2 potential (bicarbonate)", f"{full:.3f} t/t", f"{full*cfg['application_t_ha']:.1f} t/ha at "
          f"{cfg['application_t_ha']} t/ha", delta_color="off")
c2.metric("CO2 potential (divalent only)", f"{cons:.3f} t/t", f"{cons*cfg['application_t_ha']:.1f} t/ha",
          delta_color="off")
with c3:
    st.caption(f"**Mineralogy in use: {mineral_source}.** Minerals (g per 100 g): " + ", ".join(
        f"{ {'Plag': f'plagioclase An{xan*100:.0f}', 'Kfs':'K-feldspar','Cpx':'diopside','Fo':'forsterite','Fa':'fayalite','Inert':'Fe-Ti oxides and other (inert)','Gls':'basaltic glass'}[k]} {v:.1f}"
        for k, v in g100.items() if v > 0))

if cfg["acidity_mode"] == "andosol":
    try:
        ap = andosol_pool(cfg)
        st.caption(
            f"**Andosol acidity.** Permanent-charge sites {ap['x_cmol']:.1f} cmol(+)/kg, holding the bases and "
            f"{ap['al_exch_cmol']:.2f} cmol(+)/kg of exchangeable Al. Variable-charge sites {ap['y_cmol']:.1f} cmol(+)/kg, "
            f"{ap['frac_charged_at_soil_ph']:.0%} charged at the soil pH and {ap['frac_charged_at_ph7']:.0%} at pH 7 "
            f"(half-point pH {ap['half_ph']:.1f}, assumed). This reproduces an effective CEC of {ap['ecec']:.1f} at the "
            f"soil pH and a CEC of {ap['cec7']:.1f} at pH 7. Soil-water Al is held by {ap['al_phase']}. "
            f"Fitted exchange constant for Al: log K {ap['logk_alx3']:.2f}, against {ap['logk_alx3_database']:.2f} in "
            f"llnl.dat, which on its own would hold {ap['al_exch_database_cmol']:.2f} cmol(+)/kg of Al.")
    except ValueError as err:
        st.error(str(err))

with st.expander("Inputs of this scenario and their status (measured, stand-in, assumed)"):
    st.dataframe(data_io.provenance_table(cfg), use_container_width=True, hide_index=True)
    st.caption("Stand-in means a real measurement on another rock or at another site. Assumed means a modelling "
               "choice or a literature value. Nothing here is a measurement of TECHNOVIA's rock or of the "
               "Lomba da Maia soil until those results are entered in data/lab.")

# ------------------------------------------------------------------ live run
def fig_lines(h, title):
    f = go.Figure()
    f.add_trace(go.Scatter(x=h["year"], y=h["released_t"], name="Released from rock", line=dict(color=INK, width=3)))
    f.add_trace(go.Scatter(x=h["year"], y=h["exchange_t"], name="Held on soil exchange sites",
                           line=dict(color=GREY, width=3, dash="dot")))
    f.add_trace(go.Scatter(x=h["year"], y=h["export_t"], name="CO2 removed (leached as bicarbonate)",
                           line=dict(color=TEAL, width=4)))
    f.update_layout(title=title, xaxis_title="Years after spreading", yaxis_title="tCO2e per ha",
                    height=380, margin=dict(l=10, r=10, t=50, b=10), legend=dict(orientation="h", y=-0.25))
    return f


def fig_profile(p, title):
    f = go.Figure()
    f.add_trace(go.Scatter(x=p["pH_ctrl"], y=p["depth_cm"], name="pH, no rock", line=dict(color=GREY, dash="dot")))
    f.add_trace(go.Scatter(x=p["pH_rock"], y=p["depth_cm"], name="pH, with basalt", line=dict(color=TEAL, width=3)))
    f.add_trace(go.Scatter(x=p["bs_ctrl"], y=p["depth_cm"], name="Base saturation %, no rock", xaxis="x2",
                           line=dict(color="#D8B79A", dash="dot")))
    f.add_trace(go.Scatter(x=p["bs_rock"], y=p["depth_cm"], name="Base saturation %, with basalt", xaxis="x2",
                           line=dict(color=ORANGE, width=3)))
    f.update_layout(title=title, height=420, margin=dict(l=10, r=10, t=60, b=10),
                    yaxis=dict(title="Depth (cm)", autorange="reversed"),
                    xaxis=dict(title="pH", domain=[0, 1]),
                    xaxis2=dict(title="Base saturation (%)", overlaying="x", side="top"),
                    legend=dict(orientation="h", y=-0.2))
    return f


live = st.container()
if go_btn:
    with live:
        prog = st.progress(0.0, text="Starting PHREEQC…")
        k1, k2, k3, k4 = st.columns(4)
        m1, m2, m3, m4 = k1.empty(), k2.empty(), k3.empty(), k4.empty()
        cA, cB = st.columns(2)
        chartA, chartB = cA.empty(), cB.empty()
        t0 = time.time()
        res = None
        for res in run_stream(cfg):
            h = res["hist"]; y = res["year"]
            prog.progress(y / res["years"], text=f"Year {y} of {res['years']} · {time.time()-t0:.0f} s")
            m1.metric("Rock dissolved", f"{h['frac_diss'][-1]*100:.2f}%")
            m2.metric("Cations released (tCO2e/ha)", f"{h['released_t'][-1]:.3f}")
            m3.metric("CO2 removed (t/ha)", f"{max(h['export_t'][-1],0):.3f}")
            m4.metric("Topsoil pH change", f"{h['ph_top_rock'][-1]-h['ph_top_ctrl'][-1]:+.2f}")
            chartA.plotly_chart(fig_lines(h, "Where the weathered cations go"), use_container_width=True,
                                key=f"a{y}")
            chartB.plotly_chart(fig_profile(h["profiles"][-1], f"Soil profile after year {y}"),
                                use_container_width=True, key=f"b{y}")
        prog.progress(1.0, text=f"Done in {time.time()-t0:.0f} s")
        st.session_state.runs[name] = dict(cfg=copy.deepcopy(cfg), res=res)
        st.session_state.last = name

# ------------------------------------------------------------------ explore results
if st.session_state.runs:
    st.divider()
    last = st.session_state.get("last")
    tab1, tabI, tab2, tab3 = st.tabs(["Explore a run", "Isometric accounting", "Compare runs", "Download"])
    with tab1:
        pick = st.selectbox("Run", list(st.session_state.runs), index=list(st.session_state.runs).index(last))
        R = st.session_state.runs[pick]["res"]; h = R["hist"]
        yr = st.slider("Year", 1, R["years"], R["years"])
        a, b = st.columns(2)
        a.plotly_chart(fig_lines(h, "Where the weathered cations go"), use_container_width=True)
        b.plotly_chart(fig_profile(h["profiles"][yr - 1], f"Soil profile after year {yr}"), use_container_width=True)
        e = R["effluent"]
        f = go.Figure(go.Scatter(x=e["t"], y=np.array(e["d_alk"]) * 1e6, line=dict(color=TEAL)))
        f.update_layout(title="Extra alkalinity in drainage at the column base (µeq/L, vs control)",
                        xaxis_title="Years", height=300, margin=dict(l=10, r=10, t=50, b=10))
        st.plotly_chart(f, use_container_width=True)
    with tabI:
        st.caption("Laid out like Isometric Enhanced Weathering in Agriculture v1.3 (Equations 1, 2 and 4), "
                   "per hectare, cumulative to the chosen year. Removal is counted when base cations leave the "
                   "Near-Field Zone (0-20 cm). **This is a model forecast.** Under v1.3, credits must be quantified "
                   "from soil or porewater measurements (List 1); reactive transport models are not yet an accepted "
                   "quantification approach, but can support the sorption time-lag term and project design. "
                   "River/ocean factor and emission factors are placeholders.")
        Ri = st.session_state.runs[pick]["res"]["hist"].get("iso")
        if not Ri:
            st.info("Re-run this scenario to get the Isometric accounting.")
        else:
            yi = st.slider("Reporting period end (year)", 1, len(Ri), len(Ri), key="iso_year")
            I = Ri[yi - 1]
            gross = I["weathering"] - I["biomass_loss"] - I["net_new_carbonate"] - I["net_new_silicate"] - I["net_sorption"] - I["non_carbonic_neut"]
            eta_loss = gross * (1 - I["eta_river_ocean"])
            wf = go.Figure(go.Waterfall(
                orientation="v",
                measure=["absolute", "relative", "relative", "relative", "relative", "relative", "total", "relative", "relative", "total"],
                x=["Weathering (cations released)", "Biomass loss", "Net new carbonate", "Net new silicate",
                   "Net sorption (held on CEC)", "Non-carbonic neutralisation", "Exported from NFZ",
                   f"River/ocean (η = {I['eta_river_ocean']:.2f})", "Project emissions", "Net removal"],
                y=[I["weathering"], -I["biomass_loss"], -I["net_new_carbonate"], -I["net_new_silicate"],
                   -I["net_sorption"], -I["non_carbonic_neut"], 0, -eta_loss, -I["emissions"], 0],
                connector=dict(line=dict(color="#C9CFCB")),
                increasing=dict(marker=dict(color=TEAL)), decreasing=dict(marker=dict(color=ORANGE)),
                totals=dict(marker=dict(color=INK))))
            wf.update_layout(title=f"Net CO2e removal per hectare, years 1-{yi} (t CO2e/ha)", height=460,
                             margin=dict(l=10, r=10, t=50, b=10), yaxis_title="t CO2e/ha")
            st.plotly_chart(wf, use_container_width=True)
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("CO2e stored (t/ha)", f"{I['stored']:.3f}")
            m2.metric("Project emissions (t/ha)", f"{I['emissions']:.3f}")
            m3.metric("Net removal (t/ha)", f"{I['net_removal']:.3f}")
            m4.metric("Net removal, 100 ha (t)", f"{I['net_removal']*100:.1f}")
            st.dataframe(pd.DataFrame([{
                "Term (Isometric EWA v1.3)": k, "t CO2e/ha": f"{v:.4f}" if isinstance(v, float) else str(v),
                "How it is obtained here": h} for k, v, h in [
                ("CO2e Weathering", I["weathering"], "PHREEQC: base cations released from the rock"),
                ("CO2e BiomassLoss", I["biomass_loss"], "not modelled yet (0): measure Ca, Mg, K in harvested grass"),
                ("CO2e NetNewCarbonate", I["net_new_carbonate"], "PHREEQC: calcite formed in 0-20 cm, rock minus control"),
                ("CO2e NetNewSilicate", I["net_new_silicate"], "not resolved for base cations (kaolinite and gibbsite hold none)"),
                ("CO2e NetSorption", I["net_sorption"], "PHREEQC: extra exchangeable Ca, Mg, K, Na in 0-20 cm"),
                ("CO2e NonCarbonicNeut", I["non_carbonic_neut"], "residual of the charge balance"),
                ("Exported from NFZ", I["exported_from_nfz"], "extra alkalinity crossing 20 cm, rock minus control"),
                ("η RiverOcean", I["eta_river_ocean"], "PLACEHOLDER: set the value Isometric requires"),
                ("CO2e Stored", I["stored"], "Equation 4"),
                ("CO2e Counterfactual", "subtracted", "every term is rock minus a no-rock control column"),
                ("CO2e Emissions", I["emissions"], "haulage, spreading, milling, other: PLACEHOLDER factors"),
                ("CO2e Removal", I["net_removal"], "Equation 1: stored minus emissions")]]),
                use_container_width=True, hide_index=True)
            st.caption("Credits under Isometric are issued conservatively (around the 30th-40th percentile of the "
                       "measured CDR distribution), so a credited figure would sit below this central forecast.")
    with tab2:
        f = go.Figure()
        for i, (nm, run) in enumerate(st.session_state.runs.items()):
            h = run["res"]["hist"]; c = PALETTE[i % len(PALETTE)]
            f.add_trace(go.Scatter(x=h["year"], y=h["released_t"], name=f"{nm}: released", line=dict(color=c, dash="dot")))
            f.add_trace(go.Scatter(x=h["year"], y=h["export_t"], name=f"{nm}: CO2 removed", line=dict(color=c, width=3)))
        f.update_layout(height=450, xaxis_title="Years", yaxis_title="tCO2e per ha",
                        legend=dict(orientation="h", y=-0.25), margin=dict(l=10, r=10, t=30, b=10))
        st.plotly_chart(f, use_container_width=True)
        rows = []
        for nm, run in st.session_state.runs.items():
            h = run["res"]["hist"]; c = run["cfg"]
            rows.append({"run": nm, "soil": c["soil_key"], "grain size": c["psd_key"], "t/ha": c["application_t_ha"],
                         "glass": c["glass_frac"], "field factor": c["field_factor"], "acidity pool": c["acidity_mode"],
                         "years": len(h["year"]), "rock dissolved %": round(h["frac_diss"][-1] * 100, 2),
                         "released tCO2e/ha": round(h["released_t"][-1], 3),
                         "CO2 removed t/ha": round(h["export_t"][-1], 3),
                         "on exchange tCO2e/ha": round(h["exchange_t"][-1], 3)})
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        if st.button("Clear runs"):
            st.session_state.runs = {}; st.rerun()
    with tab3:
        run = st.session_state.runs[last]
        h = run["res"]["hist"]
        df = pd.DataFrame({k: h[k] for k in ["year", "frac_diss", "released_t", "exchange_t", "export_t",
                                             "ph_top_rock", "ph_top_ctrl"]})
        st.download_button("Yearly results (CSV)", df.to_csv(index=False), f"erw_results_{last}.csv")
        st.download_button("PHREEQC input, basalt column (.pqi)", full_input(run["cfg"], True), "erw_basalt.pqi")
        st.download_button("PHREEQC input, control column (.pqi)", full_input(run["cfg"], False), "erw_control.pqi")
        st.download_button("Scenario settings (JSON)", json.dumps(run["cfg"], indent=1), "scenario.json")
        st.caption("The .pqi files run in the PHREEQC GUI or command line with llnl.dat, and reproduce this run.")
else:
    st.info("Set a scenario in the sidebar and press **Run PHREEQC**. The charts update year by year while the model runs.")
# ------------------------------------------------------------------ groundwater (MODFLOW 6 + PHREEQC)
st.divider()
st.markdown("### Groundwater step: MODFLOW 6 + PHREEQC")
st.caption("Routes the soil drainage from the last run through a simple shallow basalt aquifer to a stream "
           "(MODFLOW 6 flow and transport, PHREEQC chemistry via mf6rtm). The hydrogeology below is ASSUMED "
           "until Azorean aquifer data is available. The final-year soil drainage is held constant as recharge.")
try:
    from erw.groundwater import AQUIFER_DEFAULTS, check_binaries, run_coupled
    have_bins, exe, lib = check_binaries()
except Exception as e:  # missing packages
    have_bins, exe, lib = False, "", ""
    st.warning(f"Groundwater packages not installed ({e}). Run: pip install -r requirements.txt")
if not have_bins:
    st.info("MODFLOW 6 binaries not found. Run once in the terminal: `python scripts/get_modflow.py`")
elif not st.session_state.get("runs"):
    st.info("Run the soil model first; the groundwater step uses its drainage.")
else:
    g1, g2, g3, g4 = st.columns(4)
    L = g1.number_input("Distance field edge to stream (m)", 300, 5000, int(AQUIFER_DEFAULTS["length_m"]), 100)
    F = g2.number_input("Flow-path length under the field (m)", 100, 3000, int(AQUIFER_DEFAULTS["field_length_m"]), 50)
    K = g3.select_slider("Hydraulic conductivity (m/day)", [1.0, 3.0, 10.0, 30.0, 100.0], AQUIFER_DEFAULTS["k_m_d"])
    n = g4.select_slider("Effective porosity", [0.05, 0.1, 0.15, 0.2, 0.3], AQUIFER_DEFAULTS["porosity"])
    run_name = st.session_state.get("last")
    if st.button(f"Run MODFLOW 6 + PHREEQC for: {run_name}", type="primary"):
        run = st.session_state.runs[run_name]
        with st.spinner("Running MODFLOW 6 with PHREEQC (rock and control aquifers)…"):
            res = run_coupled(run["cfg"], dict(length_m=float(L), field_length_m=float(min(F, L - 50)), k_m_d=K,
                                               porosity=n, years=run["cfg"]["years"]),
                              soil_result=run["res"], label="app")
        st.session_state.gw = res
    if "gw" in st.session_state:
        res = st.session_state.gw
        h1, h2, h3 = st.columns(3)
        h1.metric("Extra alkalinity in field drainage (µeq/L)", f"{res['d_alk_input_ueq']:.1f}")
        h2.metric("Extra alkalinity at stream, end (µeq/L)", f"{res['d_alk_ueq'][-1]:.1f}")
        resolvable = res["d_alk_input_ueq"] >= 10
        h3.metric("Share of field alkalinity reaching stream",
                  f"{res['reach_fraction']*100:.0f}%" if resolvable else "too small to resolve",
                  help="After dilution by recharge outside the field. Needs at least ~10 µeq/L of extra "
                       "alkalinity in the drainage; below that, numerical noise dominates.")
        a, b = st.columns(2)
        f = go.Figure(go.Scatter(x=res["t_years"], y=res["d_alk_ueq"], line=dict(color=TEAL, width=3)))
        f.update_layout(title="Extra alkalinity arriving at the stream (rock minus control)",
                        xaxis_title="Years", yaxis_title="µeq/L", height=360, margin=dict(l=10, r=10, t=50, b=10))
        a.plotly_chart(f, use_container_width=True)
        pr = res["profile"]
        f2 = go.Figure()
        f2.add_trace(go.Scatter(x=pr["x_m"], y=pr["alk_ctrl"], name="No rock", line=dict(color=GREY, dash="dot")))
        f2.add_trace(go.Scatter(x=pr["x_m"], y=pr["alk_rock"], name="With basalt", line=dict(color=TEAL, width=3)))
        f2.add_vrect(x0=0, x1=res["field_end_m"], fillcolor="#E8F2EE", opacity=0.5, line_width=0,
                     annotation_text="treated field", annotation_position="top left")
        f2.update_layout(title="Groundwater alkalinity along the flow path (end of run)", xaxis_title="Distance (m)",
                         yaxis_title="Alkalinity (meq/L)", height=360, margin=dict(l=10, r=10, t=50, b=10),
                         legend=dict(orientation="h", y=-0.25))
        b.plotly_chart(f2, use_container_width=True)

