"""Builds report/Ilha_ERW_Model_Report.html and .pdf from the run results."""
import os, subprocess, html
import numpy as np
import pandas as pd
from .common import ROOT, OUT

REP = os.path.join(ROOT, "report")


def f(x, d=2):
    if isinstance(x, str):
        return x
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    return f"{x:,.{d}f}"


def tbl(df, cols, heads=None, num=None, cls=""):
    heads = heads or cols
    num = num or {}
    h = "".join(f'<th class="{"n" if c in num else ""}">{html.escape(t)}</th>' for c, t in zip(cols, heads))
    rows = []
    for _, r in df.iterrows():
        tds = []
        for c in cols:
            v = r[c]
            if c in num and not isinstance(v, str):
                tds.append(f'<td class="n">{f(v, num[c])}</td>')
            else:
                tds.append(f"<td>{html.escape(str(v))}</td>")
        rows.append("<tr>" + "".join(tds) + "</tr>")
    return f'<table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{"".join(rows)}</tbody></table>'


def fig(name, cap, w="100%"):
    return (f'<figure><img src="../outputs/{name}" style="width:{w}"><figcaption>{cap}</figcaption></figure>')


def build(R):
    cfg = R["cfg"]
    L = R["layout"]
    fs = R["fs"]
    acct = R["acct"]
    c30 = acct[(acct.scenario == "central") & (acct.month == cfg.v("design.months"))].set_index("rate_t_ha")
    c120 = acct[(acct.scenario == "central") & (acct.month == acct.month.max())].set_index("rate_t_ha")
    tt = R["trial_tot"].set_index("scenario")
    mc = R["mc"]
    mc40 = mc[mc.rate_t_ha == 40].set_index("quantity")
    sens = R["sens"]
    tau = R["tau"]
    tauB80 = tau[(tau.estimator == "B_dose_trays") & (tau.rate_t_ha == 80)].set_index("cores_per_plot_per_round")
    tauA = tau[tau.estimator == "A_tracer_Ti"]
    water = R["water"]
    lys = water[water["where"] == "plot lysimeter"].set_index("rate_t_ha")
    wmax = water[water["where"] != "plot lysimeter"].alkalinity_increase_meq_L.max()
    meso = R["meso"].set_index("rate_t_ha")
    mm = R["meso_meta"]
    pw = R["power"]
    pcol = ["power_d_pH", "power_d_CaMg", "power_d_totCaMg", "power_MANOVA"]
    pD = pw[(pw.plots == 16) & (pw.noise == "medium")].iloc[0]
    pB = pw[(pw.plots == 12) & (pw.noise == "medium")].iloc[0]
    pDh = pw[(pw.plots == 16) & (pw.noise == "high")].iloc[0]
    metals = R["metals"]
    above = metals[metals.max_cumulative_t_ha == "background already above limit"].element.tolist()
    con = R["contrast"].set_index("element")
    eff = R["effects"].set_index("rate_t_ha")
    plots = R["plots"]
    nplace = (R["status"].status == "placeholder").sum()
    lo, hi = tt.loc["low", "net_t_CO2"], tt.loc["high", "net_t_CO2"]
    emis = R["emis_kg_t"]
    pot = fs["pot_net"] * cfg.v("kinetics.downstream_retention")
    loq = mm["loq_meq_L"]

    css = open(os.path.join(REP, "doc.css")).read()
    css += """
figure { margin: 6pt 0 10pt; page-break-inside: avoid; } figure img { display: block; }
figcaption { font-size: 8.6pt; color: #222; margin-top: 3pt; font-family: Inter, Arial, sans-serif; }
pre { font-size: 7.6pt; line-height: 1.3; white-space: pre-wrap; font-family: 'DejaVu Sans Mono', monospace; margin: 4pt 0 8pt; }
table.small td, table.small th { font-size: 8pt; } table.tiny td { font-size: 7.4pt; } table.tiny th { font-size: 6.6pt; }
.eq { font-family: 'DejaVu Serif', serif; font-size: 9.6pt; margin: 3pt 0 6pt 12pt; }
.pb { break-before: page; }
"""
    S = []
    A = S.append
    A(f"""<!doctype html><html><head><meta charset="utf-8"><title>Ilha ERW model</title><style>{css}</style></head><body>
<div class="titleblock"><div class="doclabel">Ilha Carbon · Technical note · Version 1</div>
<h1>Lomba da Maia basalt rate trial: soil and geological model</h1>
<div class="sub">Model description, results, and the data still needed</div>
<table class="meta"><tbody>
<tr><td>Date</td><td>{__import__("datetime").date.today().strftime("%d %B %Y")}</td></tr>
<tr><td>Prepared by</td><td>Ilha Carbon, for the science team (MRV lead and advisors) and the land partner</td></tr>
<tr><td>Status</td><td>Input files used: {', '.join(R['inputs_used']) or 'none'}. {nplace} parameters are still placeholders (Appendix B); results that depend on them show what the model does and where the design is sensitive, not what the site will deliver.</td></tr>
<tr><td>Run it</td><td><code>python3 run_all.py</code> (about 3 minutes). Drop filled templates into <code>inputs/</code> and run again.</td></tr>
</tbody></table></div>

<h2>1. Summary</h2>
<ol>
<li><b>Layout.</b> {L['plots']} plots of 0.64 ha in {L['blocks']} blocks fit in the open ground of units {', '.join(L['units_used'])}, avoiding trees plus 5 m, a {cfg.v('site.drainage_buffer_m')} m buffer on drainage lines and slopes over {cfg.v('site.max_plot_slope_pct')}% (plot means {L['plot_slope_mean_range_pct'][0]} to {L['plot_slope_mean_range_pct'][1]}%). Treated and control area {f(L['trial_gross_ha'])} ha, between {L['elev_range_m'][0]} and {L['elev_range_m'][1]} m.</li>
<li><b>Basalt.</b> {f(L['basalt_dry_t'],0)} t dry for the 16-plot design, {f(L['basalt_dry_t_with_contingency'],0)} t with 10% contingency, about {f(L['basalt_wet_t_with_contingency'],0)} t as delivered at {cfg.v('feedstock.moisture_wet_frac')*100:.0f}% moisture ({L['truckloads']} loads). If the screening allows 12 plots: {f(L['if_3_blocks_basalt_dry_t_with_contingency'],0)} t dry with contingency.</li>
<li><b>Carbon in the grant period is small and may be net negative.</b> The rock's long-run potential is about {f(pot,0)} kg CO<sub>2</sub> per t after downstream losses. In 30 months the central case dissolves about {c30.loc[40,'frac_cations_dissolved']*100:.0f}% of the cations, and in soil at pH {cfg.v('soil.pH_H2O')} most of them neutralise soil acidity rather than leave as bicarbonate. After haulage, spreading and rock carbon ({f(emis,1)} kg CO<sub>2</sub>e per t plus {f(fs['ic_co2'],1)} kg for rock carbon), the whole trial is {f(tt.loc['central','net_t_CO2'],1)} t CO<sub>2</sub> net at 30 months (range {f(lo,1)} to {f(hi,1)} t). The trial is a measurement and calibration trial; that is the right framing for an R&amp;D call, and no credits should be promised from it.</li>
<li><b>What drives the answer:</b> the lab-to-field rate factor, the starting soil pH and the share of fine grains (rank correlations {', '.join(f"{r.parameter} {r.spearman_rho:+.2f}" for r in sens.head(3).itertuples())}). All three are measurable this winter and next year.</li>
<li><b>What can be seen within 30 months:</b> soil pH, exchangeable Ca + Mg and total Ca + Mg respond clearly to rate. The dissolved fraction itself (soil mass balance) cannot be resolved in 30 months: even with 30 cores per plot at 80 t/ha the smallest detectable fraction is about {tauB80.loc[30,'min_detectable_tau_4_plots']:.2f} against about {tau.tau_true.iloc[0]:.2f} expected. Titanium gives too little contrast with these soils and zirconium is higher in the soil than in the rock. Lysimeters under the 40 and 80 t/ha plots (expected +{lys.loc[40,'alkalinity_increase_meq_L']:.2f} and +{lys.loc[80,'alkalinity_increase_meq_L']:.2f} meq/L alkalinity) and the mesocosms are where weathering will actually be measured. Stream points W1 to W6 will not see the trial (at most {wmax:.3f} meq/L).</li>
<li><b>Plots needed:</b> with medium variability, 16 plots detect the dose response with power {pD[pcol].min():.2f} to {pD[pcol].max():.2f}; 12 plots give {pB[pcol].min():.2f} to {pB[pcol].max():.2f}. Testing single rates against the control is weak in any design ({pD['power_pH_40_vs_control_Holm']:.2f} for 40 t/ha on pH with 16 plots), so the hypotheses in Part B should be written as dose responses.</li>
<li><b>Metals:</b> with placeholder soil values the background already exceeds the Portuguese sludge-application limits for {', '.join(above)} at pH 5.5 or below; basalt adds little on top (nickel +{metals.set_index('element').loc['Ni','increase_at_80_pct']:.0f}%, chromium +{metals.set_index('element').loc['Cr','increase_at_80_pct']:.0f}% at 80 t/ha). The real issue is which reference applies to volcanic soils; Juan and the regional authority need to settle that before the application, using the screening results.</li>
</ol>

<h2>2. How the model works</h2>
<p>The model is a set of Python modules run in sequence by <code>run_all.py</code>. Each stage passes its results to the next, so one change in the inputs (for example the real particle size distribution) flows through to the carbon account, the measurement design and the power analysis.</p>
<table><thead><tr><th>Stage</th><th>Module</th><th>What it does</th><th>Main outputs</th></tr></thead><tbody>
<tr><td>Site</td><td>site.py</td><td>Eligible ground from the QGIS layers; plots oriented to the slope; compact blocks along the elevation gradient; rates randomised within blocks (seed {cfg.v('site.seed')}); tonnage per plot</td><td>trial_plots.geojson / .gpkg / .csv, fig 1</td></tr>
<tr><td>Feedstock</td><td>geochem.py</td><td>Cation charge and CO<sub>2</sub> potential from XRF, deductions for rock carbon and sulfur, CIA and MIA, mineral cation budget with glass by difference</td><td>feedstock_potential.csv, feedstock_mineral_balance.csv</td></tr>
<tr><td>Kinetics</td><td>geochem.py</td><td>Shrinking-core dissolution of each mineral in each grain size, pH- and temperature-dependent rate laws</td><td>fig 2, fig 3</td></tr>
<tr><td>Soil</td><td>geochem.py</td><td>Monthly balance of released cations: bicarbonate export, exchange retention, leaching with strong acids, secondary minerals, plant uptake; pH and base saturation feedback; background acid load</td><td>carbon_account_by_scenario.csv</td></tr>
<tr><td>Uncertainty</td><td>geochem.py</td><td>Monte Carlo over nine uncertain inputs; rank sensitivity</td><td>mc_*.csv, fig 4</td></tr>
<tr><td>Metals</td><td>geochem.py</td><td>Soil mixing mass balance against screening limits</td><td>trace_metal_mass_balance.csv</td></tr>
<tr><td>Measurement</td><td>mrv.py</td><td>Tracer contrast, soil mass-balance detectability, lysimeter and stream signals, mesocosm design, sampling counts</td><td>tau_detectability.csv, water_signal.csv, mesocosm_design.csv, sampling_schedule_and_counts.csv, fig 5</td></tr>
<tr><td>Statistics</td><td>stats.py</td><td>Power by design using the model's effects and the screening noise; the full analysis for real data</td><td>power_by_design.csv, analysis_*.txt, figs 6 and 7</td></tr>
</tbody></table>
<p><b>Inputs.</b> <code>config.json</code> holds every parameter with its status (placeholder, literature, design, measured). The folder <code>templates/</code> has one CSV per data source. A filled CSV placed in <code>inputs/</code> overrides the matching placeholder and is marked measured in <code>outputs/parameter_status.csv</code>. Results of the trial itself go in <code>inputs/plot_results.csv</code>.</p>
""")

    # 3 layout
    pt = plots.copy()
    A(f"""<h2>3. Trial layout and basalt quantities</h2>
<p>Candidate ground is the open land of units {', '.join(cfg.v('site.candidate_units'))} (contiguous, rated easy or workable, 90th-percentile slope under 20%). From {f(L['trial_footprint_ha'],0)} ha of footprint the layout keeps the 16 most compact 80 × 80 m plots, separated by {cfg.v('site.alley_m')} m untreated alleys. Sampling is confined to the inner {f(plots.net_sampling_ha.iloc[0],2)} ha of each plot ({cfg.v('site.net_inset_m')} m inset) to stay clear of spreader edges. Blocks are groups of four neighbouring plots, numbered from the top of the slope down, so the block term in the analysis absorbs the elevation, aspect and drainage gradient. The 0 t/ha plot in each block is the control; the old 0.3 ha control squares from the 100 ha plan are no longer needed. Rates are randomised within blocks, redrawn if rate lines up with slope position or plot slope ({L['randomisation']}).</p>
{fig('fig01_trial_layout.png', 'Figure 1. Plot layout with rates (t/ha dry), blocks, no-spread areas and drainage lines. Files: outputs/trial_plots.geojson (WGS84) and trial_plots_EPSG5015.gpkg for QGIS.', '88%')}
{tbl(pt, ['plot','block','rate_t_ha','unit','gross_ha','basalt_dry_t','basalt_wet_t','slope_mean','elev_mean'],
     ['Plot','Block','Rate t/ha','Unit','Area ha','Dry t','Wet t','Slope %','Elev. m'],
     {'rate_t_ha':0,'gross_ha':2,'basalt_dry_t':1,'basalt_wet_t':1,'slope_mean':1,'elev_mean':0}, 'small')}
<p>Total {f(L['basalt_dry_t'],1)} t dry; with {cfg.v('design.contingency_frac')*100:.0f}% contingency {f(L['basalt_dry_t_with_contingency'],1)} t dry, or {f(L['basalt_wet_t_with_contingency'],1)} t delivered at the placeholder moisture. This is at the top of the 300 to 350 t used in the plan so far, before contingency, because 16 plots of 0.64 ha cover 10.24 ha at a mean rate of 35 t/ha. The 12-plot version (blocks B1 to B3) needs {f(L['if_3_blocks_basalt_dry_t_with_contingency'],0)} t on {f(L['if_3_blocks_gross_ha'],2)} ha. Application rate is verified with {cfg.v('design.catch_trays_per_plot')} catch trays per treated plot.</p>
""")

    # 4 feedstock
    ft = R["feed_tab"].copy()
    ft["value"] = ft["value"].map(lambda v: f(v, 2))
    mt = R["mineral_tab"].copy()
    A(f"""<h2>4. Feedstock geochemistry</h2>
<p>Values below come from the feedstock analyses in <code>inputs/</code> (fusion XRF, ICP-MS, XRD, PSD, moisture, carbon and sulfur).</p>
{tbl(ft, ['quantity','value','unit'], ['Quantity','Value','Unit'], cls='small')}
<p>The CO<sub>2</sub> potential counts one mole of CO<sub>2</sub> per equivalent of base cation (the bicarbonate pathway), subtracts carbon already in the rock and the acidity sulfide would generate, then applies downstream retention. Including sodium and potassium adds about 20%; the Ca + Mg line is the conservative figure. CIA and MIA near 35 indicate unweathered basalt, as expected for fresh quarry fines.</p>
<h3>Where the cations sit</h3>
{tbl(mt, ['mineral','formula','mode_wt_frac','eq_per_kg_mineral','share_of_cations'], ['Mineral','Assumed composition','Mass fraction','eq/kg mineral','Share of rock cations'], {'mode_wt_frac':2,'eq_per_kg_mineral':2,'share_of_cations':2}, 'small')}
<p>Crystalline phases are given ideal compositions; the glass is assigned whatever cations remain so that the total matches the XRF. {'<b>' + R['mineral_flag'] + '.</b>' if R['mineral_flag'] else 'The balance closes with a plausible glass composition.'} When the XRD arrives, check that the glass by difference is not negative and that its cation content is in a sensible range for basaltic glass (about 8 to 11 eq/kg); if not, the mode or the assumed compositions need revising with Juan.</p>
""")

    # 5 kinetics
    A(f"""<h2>5. Weathering kinetics</h2>
<p>Each reactive mineral dissolves at a rate per unit surface that depends on pH and temperature:</p>
<p class="eq">r = k<sub>acid</sub> · a<sub>H+</sub><sup>n</sup> · exp[−E<sub>a,acid</sub>/R (1/T − 1/298.15)] + k<sub>neutral</sub> · exp[−E<sub>a,neutral</sub>/R (1/T − 1/298.15)]</p>
<p>with constants from Palandri and Kharaka (2004) for olivine, augite and labradorite and an approximation for basaltic glass after Gislason and Oelkers (2003). Grains shrink from the surface inward: the radius falls by r · λ · V<sub>m</sub> · F per unit time, where λ is surface roughness ({cfg.v('feedstock.roughness')}), V<sub>m</sub> the molar volume and F the lab-to-field factor. Each sieve class is split into six sub-sizes, and every mineral is assumed present in every size class in proportion to the mode.</p>
{fig('fig03_mineral_rates.png', f"Figure 2. Left: rate per unit surface against soil pH at {cfg.v('climate.mean_temp_C')} C. Right: time for a grain to dissolve completely at the starting soil pH. Only grains below roughly 10 micrometres (olivine, glass) dissolve within the grant period; the 2 to 4 mm fraction is effectively inert on this time scale.")}
<p>The consequence for the feedstock: unmilled 0 to 4 mm fines react mostly through their finest tail. In the placeholder PSD {cfg.v('feedstock.psd_mass_frac')[0]*100:.0f}% of the mass is below 63 micrometres, and that fraction controls almost all weathering in the first years. The laser-diffraction curve of the fine tail is therefore one of the most valuable lab results, and the trial should report PSD and BET for the batch actually spread.</p>
""")

    # 6 soil and carbon
    ct = acct[(acct.month == cfg.v("design.months")) & (acct.rate_t_ha > 0)].copy()
    ct["dissolved_pct"] = ct.frac_cations_dissolved * 100
    A(f"""<h2 class="pb">6. Soil response and carbon account</h2>
<p>Released cations are tracked month by month in the top {cfg.v('soil.mixing_depth_cm')} cm. A share is lost to secondary minerals ({cfg.v('soil.secondary_loss_frac')}) and grass uptake ({cfg.v('soil.plant_uptake_frac')}). Of the rest, the share that leaves as bicarbonate follows carbonate speciation at the soil pH, 1/(1 + 10<sup>pKa1 − pH</sup>): about 10% at pH 5.4, 50% at 6.35, 80% at 7. The remainder neutralises acidity: {cfg.v('soil.exchange_retention_frac')*100:.0f}% is held on exchange sites (raising base saturation) and the rest leaves with nitrate, sulfate and chloride, which removes no carbon. Soil pH moves by the net acid neutralised divided by the soil's buffer capacity ({cfg.v('soil.pH_buffer_cmolc_kg_pH')} cmolc/kg per pH unit), and a background acid load of {cfg.v('soil.acid_load_keq_ha_yr')} keq/ha/yr from fertiliser and grazing acts on every plot. Rising pH speeds up bicarbonate export and slows acid-driven dissolution, so the response to rate is not linear.</p>
{fig('fig02_kinetics_cdr_pH.png', 'Figure 3. Left: share of rock cations released under the three field-factor scenarios. Centre: net CO2 removal per hectare after deductions (central). Right: soil pH, with the control drifting down under the background acid load. Vertical line: end of the 30-month grant period.')}
{tbl(ct, ['scenario','rate_t_ha','dissolved_pct','pH','d_exch_CaMg_cmolc_kg','bicarbonate_keq_ha','retained_keq_ha','gross_CDR_t_ha','deduct_emissions_t_ha','deduct_rock_C_t_ha','net_CDR_t_ha'],
     ['Scenario','Rate','Dissolved %','pH','Δ exch. Ca+Mg','HCO3 keq/ha','Retained keq/ha','Gross t/ha','Emissions t/ha','Rock C t/ha','Net t/ha'],
     {'rate_t_ha':0,'dissolved_pct':1,'pH':2,'d_exch_CaMg_cmolc_kg':2,'bicarbonate_keq_ha':2,'retained_keq_ha':2,'gross_CDR_t_ha':3,'deduct_emissions_t_ha':3,'deduct_rock_C_t_ha':3,'net_CDR_t_ha':3}, 'tiny')}
<p>Whole trial at 30 months: gross {f(tt.loc['central','gross_t_CO2'],1)} t CO<sub>2</sub> and net {f(tt.loc['central','net_t_CO2'],1)} t in the central case (low {f(tt.loc['low','net_t_CO2'],1)}, high {f(tt.loc['high','net_t_CO2'],1)}). Over ten years the central case turns positive at 40 and 80 t/ha (net {f(c120.loc[40,'net_CDR_t_ha'],2)} and {f(c120.loc[80,'net_CDR_t_ha'],2)} t/ha). Three points follow for the plan:</p>
<ul>
<li>Cations held on exchange sites are not lost. They are a store that can leave as bicarbonate later if pH rises, or be lost with nitrate if acidification continues. The model counts only what leaves as bicarbonate, which is the conservative reading; Isometric's treatment of non-carbonic-acid weathering should be checked against it.</li>
<li>Haulage and spreading ({f(emis,1)} kg CO<sub>2</sub>e/t) are small against the long-run potential but large against what dissolves in the first two years. Rock carbon matters as much: the placeholder 0.05% inorganic C costs {f(fs['ic_co2'],1)} kg CO<sub>2</sub>/t up front.</li>
<li>The farm's nitrogen and lime records set the background acid load and decide how much of the early weathering is spent neutralising fertiliser acidity. That data belongs in the land partner's work package.</li>
</ul>
""")

    # 7 uncertainty
    mt7 = mc[mc.rate_t_ha.isin([40, 80])].copy()
    A(f"""<h2>7. Uncertainty</h2>
<p>{cfg.v('uncertainty.n_monte_carlo')} runs draw nine uncertain inputs at once (distributions in outputs/mc_input_distributions.csv). The result is a range for each quantity at 30 months and a ranking of which inputs matter.</p>
{fig('fig04_uncertainty.png', 'Figure 4. Net CO2 removal per tonne of rock at 30 months for 40 t/ha, and rank correlation of each input with it.')}
{tbl(mt7, ['rate_t_ha','quantity','P10','P50','P90'], ['Rate','Quantity','P10','P50','P90'], {'rate_t_ha':0,'P10':3,'P50':3,'P90':3}, 'small')}
<p>The field factor dominates, followed by starting pH and the fines. The mesocosms and the first field rounds are designed to pin down the field factor; the screening gives the pH; the lab gives the fines. After those three, the range should narrow sharply.</p>
""")

    # 8 metals
    mtab = metals.copy()
    A(f"""<h2>8. Trace metals</h2>
<p>Concentrations after mixing each rate into the top {cfg.v('soil.mixing_depth_cm')} cm (no dissolution, no leaching: the worst case for accumulation), against the soil limits for sludge application in Decreto-Lei 276/2009 for the soil's pH class. These limits are a screening reference only and are not written for rock dust or volcanic soils.</p>
{tbl(mtab, ['element','soil_mg_kg','rock_mg_kg','limit_mg_kg','after_20_t_ha','after_40_t_ha','after_80_t_ha','increase_at_80_pct','max_cumulative_t_ha'],
     ['Element','Soil','Rock','Limit','20 t/ha','40 t/ha','80 t/ha','+% at 80','Max cumulative t/ha'],
     {'soil_mg_kg':2,'rock_mg_kg':2,'limit_mg_kg':1,'after_20_t_ha':2,'after_40_t_ha':2,'after_80_t_ha':2,'increase_at_80_pct':1}, 'tiny')}
<p>Two things to settle with the screening data: the real soil background (Azorean andosols are often naturally high in nickel and chromium), and the reference values that apply. The grant's DNSH answers should rest on the increment over background, the leaching tests on the fines (TCLP, SPLP, pH-dependent), and pasture tissue metals from the trial.</p>
""")

    # 9 measurement
    ctab = R["contrast"].reset_index(drop=True)
    tsub = tau[tau.cores_per_plot_per_round.isin([6, 30])].copy()
    wl = water.copy()
    wl["where"] = [f"Lysimeter, {r:.0f} t/ha" if w == "plot lysimeter" else f"{w.replace(' water sampling', '')} ({c:.0f} ha catchment)" for w, r, c in zip(wl["where"], wl.rate_t_ha, wl.get("catchment_ha", pd.Series([np.nan]*len(wl))))]
    A(f"""<h2 class="pb">9. What can be measured within 30 months</h2>
<h3>Tracer contrast</h3>
<p>Soil mass balance needs an element that is much richer in the rock than in the soil and does not move. The shift each rate causes, against the core-to-core noise with {cfg.v('design.individual_cores_per_plot')} cores before and after:</p>
{tbl(ctab, ['element','soil','rock','rock_to_soil','shift_at_40','signal_to_noise_40_n6','shift_at_80','signal_to_noise_80_n6'], ['Element','Soil','Rock','Rock/soil','Shift 40','S/N 40','Shift 80','S/N 80'], {'soil':2,'rock':2,'rock_to_soil':2,'shift_at_40':3,'signal_to_noise_40_n6':2,'shift_at_80':3,'signal_to_noise_80_n6':2}, 'small')}
<p>Calcium and magnesium move clearly, which confirms the rock is there. Titanium moves less than its own noise, and zirconium goes the wrong way because soils on trachytic pumice can carry more zirconium than basalt. Check both in the screening before relying on either.</p>
<h3>Dissolved fraction from the soil</h3>
{fig('fig05_tau_detectability.png', 'Figure 5. Smallest dissolved fraction that could be detected with four plots per rate, against the fraction expected after 30 months (dotted). Left: mixing estimated from titanium. Right: mixing estimated from the measured dose.')}
<p>Neither approach resolves 30-month dissolution: with the Ti tracer the estimate is far too noisy to use, and with the dose measured by catch trays the smallest detectable fraction at 80 t/ha is about {tauB80.loc[6,'min_detectable_tau_4_plots']:.2f} with 6 cores and {tauB80.loc[30,'min_detectable_tau_4_plots']:.2f} with 30, against {tau.tau_true.iloc[0]:.2f} expected. Soil mass balance stays in the protocol as a long-term measurement (it becomes useful as dissolution accumulates), but the grant-period evidence for weathering has to come from soil exchange chemistry, pore water and mesocosms.</p>
<h3>Water</h3>
{tbl(wl, ['where','bicarbonate_eq_ha_yr_in_year_3','alkalinity_increase_meq_L','base_cation_increase_meq_L'], ['Where','HCO3 export eq/yr','Alkalinity + meq/L','Base cations + meq/L'], {'bicarbonate_eq_ha_yr_in_year_3':0,'alkalinity_increase_meq_L':4,'base_cation_increase_meq_L':4}, 'small')}
<p>Bicarbonate export is per hectare for the lysimeter rows and the whole trial for the stream rows, both in the third year. Against a Gran-titration LOQ of about {loq} meq/L, suction-cup lysimeters at 30 cm should see the 40 and 80 t/ha plots, especially in base cations. The stream points drain 75 to 560 ha and dilute the trial below detection; keep them for context and for the DNSH record, not for carbon.</p>
""")

    # 10 mesocosm
    ms = R["meso"].copy()
    A(f"""<h2>10. Mesocosms</h2>
<p>{mm['columns']} columns of {cfg.v('mesocosm.column_diameter_cm')} cm diameter with {cfg.v('mesocosm.soil_height_cm')} cm of trial soil repacked by horizon: rates {', '.join(str(r) for r in cfg.v('mesocosm.rates_t_ha'))} t/ha with {cfg.v('mesocosm.replicates')} replicates, plus {cfg.v('mesocosm.sand_controls')} columns of basalt in acid-washed quartz sand to measure the rock's own dissolution rate. Artificial rain at the site drainage rate ({cfg.v('mesocosm.irrigation_mm_yr')} mm/yr) once a week for {cfg.v('mesocosm.months')} months, at about 20 C. They measure what the field cannot within the grant: the cation and alkalinity flux, the field factor, the split between bicarbonate and strong-acid export, and secondary mineral formation (destructive sampling at the end).</p>
{tbl(ms, ['rate_t_ha','rock_g_per_column','soil_kg_per_column','irrigation_L_per_week','alk_increase_meq_L_month3','alk_increase_meq_L_month12','base_cations_meq_L_month12','pH_0_10cm_at_12m'],
     ['Rate','Rock g','Soil kg','Water L/week','Alk. month 3','Alk. month 12','Cations month 12','pH at 12 m'],
     {'rate_t_ha':0,'rock_g_per_column':1,'soil_kg_per_column':2,'irrigation_L_per_week':3,'alk_increase_meq_L_month3':3,'alk_increase_meq_L_month12':3,'base_cations_meq_L_month12':3,'pH_0_10cm_at_12m':2}, 'small')}
<p>Every treated rate gives a leachate signal near or above the LOQ, and the 160 t/ha columns extend the curve beyond the field rates. About {mm['weekly_leachate_samples']} weekly leachate samples over the year; composite to fortnightly for the full cation panel if lab capacity is short.</p>
""")

    # 11 stats
    pwt = pw.copy()
    ex = R["example_analysis"]
    A(f"""<h2 class="pb">11. Statistics</h2>
<p>Power is simulated with the effects the model predicts at 30 months (table below) and plot-level noise built from the core-to-core variability, the number of cores per composite, lab precision and the between-plot variation that blocks do not remove. The three responses are change in pH, change in exchangeable Ca + Mg (composite of {cfg.v('design.cores_per_composite')} cores) and change in total Ca + Mg (mean of {cfg.v('design.individual_cores_per_plot')} individual cores).</p>
{tbl(R['effects'], ['rate_t_ha','d_pH','d_CaMg','d_totCaMg'], ['Rate t/ha','Δ pH','Δ exch. Ca+Mg cmolc/kg','Δ total Ca+Mg eq/kg'], {'rate_t_ha':0,'d_pH':3,'d_CaMg':3,'d_totCaMg':3}, 'small')}
{fig('fig06_power_by_design.png', 'Figure 6. Chance of detecting the dose response at 30 months by design and variability. Regression on rate with blocks for each response; MANOVA jointly; and 40 t/ha against control on pH with Holm correction.')}
{tbl(pwt, ['design','plots','noise','power_d_pH','power_d_CaMg','power_d_totCaMg','power_MANOVA','power_pH_40_vs_control_Holm'], ['Design','Plots','Noise','pH','Exch. Ca+Mg','Total Ca+Mg','MANOVA','40 vs 0 (pH)'], {'plots':0,'power_d_pH':2,'power_d_CaMg':2,'power_d_totCaMg':2,'power_MANOVA':2,'power_pH_40_vs_control_Holm':2}, 'tiny')}
<p><b>Recommendation.</b> Lay out 16 plots (4 rates × 4 blocks). Drop to 12 only if the screening shows low variability. Under high variability even 16 plots reach only {pDh[pcol].min():.2f} to {pDh[pcol].max():.2f}; in that case the trial should lean on the repeated rounds (mixed model) and the lysimeters.</p>
<h3>Statistical analysis plan (for Part B)</h3>
<ol>
<li>Primary: for each response, linear regression of the change from baseline on rate (continuous, t/ha) with block as a fixed effect, at 30 months. The slope and its 95% interval are the main results.</li>
<li>Shape: a quadratic term and a saturating curve a(1 − e<sup>−b·rate</sup>), compared by AIC, to find whether higher rates give diminishing returns (expected for pH).</li>
<li>Joint check: MANOVA (Pillai's trace) on the three responses with rate and block.</li>
<li>Over time: linear mixed model with rate × month, block fixed, plot as random intercept, using all six rounds.</li>
<li>Each rate against control: Holm-adjusted contrasts, reported as secondary because they are underpowered.</li>
<li>Covariates checked: baseline value of the response, mean plot slope and elevation. Missing plots are analysed as missing (no imputation), with a sensitivity run excluding any plot with a spreading failure (catch-tray rate off by more than 25%).</li>
<li>Success: a positive slope for pH and exchangeable Ca + Mg at p &lt; 0.05, alkalinity in the 80 t/ha lysimeters above control, and a mesocosm field factor estimated with its interval.</li>
</ol>
<h3>Example run of the analysis on simulated data</h3>
<p>The code that will analyse the real data was run on a simulated dataset for the actual plot list (outputs/example_simulated_field_results.csv, not data). Output:</p>
<pre>{html.escape(ex)}</pre>
{fig('fig07_analysis_example_simulated_dose_response.png', 'Figure 7. Example dose-response plot produced by the analysis (simulated data).', '92%')}
""")

    # 12 sampling
    smp = R["sampling"].copy()
    A(f"""<h2>12. Sampling schedule and analysis counts</h2>
<p>Counts only; prices come later from the lab quotes. These feed the budget workbook.</p>
{tbl(smp, ['month','item','samples'], ['When','What','Samples'], {'samples':0}, 'small')}
""")

    # 13 data needed
    st = R["status"]
    ph = st[st.status == "placeholder"]
    A(f"""<h2 class="pb">13. Data that replace the placeholders</h2>
<table class="small"><thead><tr><th>Priority</th><th>Data</th><th>Template</th><th>Source and timing</th><th>What it changes</th></tr></thead><tbody>
<tr><td>1</td><td>Soil pH, CEC, base saturation, bulk density, total elements, pH buffer capacity</td><td>soil_screening.csv</td><td>December screening (Juan's protocol); add a buffer titration on the composites</td><td>Carbon efficiency, pH response, metals, tracer choice, power</td></tr>
<tr><td>2</td><td>Core-to-core variability</td><td>soil_cores.csv</td><td>Analyse 5 to 6 subsamples of two composites individually</td><td>Number of plots and cores</td></tr>
<tr><td>3</td><td>Fines PSD including the fine tail; BET</td><td>feedstock_psd.csv</td><td>Tecnovia composites, January</td><td>Dissolution rate (second biggest driver)</td></tr>
<tr><td>4</td><td>XRF, ICP-MS, XRD, carbon, sulfur, moisture</td><td>feedstock_xrf.csv, _xrd.csv, _other.csv</td><td>Tecnovia composites, January</td><td>Potential, deductions, metals, tonnage wet/dry</td></tr>
<tr><td>5</td><td>Farm N and lime records, cut and grazing calendar</td><td>config.json (acid load)</td><td>Land partner</td><td>Background acid load, control drift</td></tr>
<tr><td>6</td><td>Temperature and rainfall at the site</td><td>config.json (climate)</td><td>IPMA or a farm logger</td><td>Rates (Arrhenius) and drainage</td></tr>
<tr><td>7</td><td>Field factor</td><td>config.json (kinetics)</td><td>Mesocosms and first field rounds</td><td>Main uncertainty</td></tr>
<tr><td>8</td><td>Haul distance, truck payload, fuel use; fate of fines</td><td>config.json (emissions)</td><td>Tecnovia and haulier</td><td>Emissions deduction and quarry allocation</td></tr>
</tbody></table>
<h2>14. Limits of the model</h2>
<ul>
<li>One mixed layer of {cfg.v('soil.mixing_depth_cm')} cm; no vertical transport of rock or solutes and no explicit pore-water chemistry. A reactive-transport model (for example PHREEQC or SCEPTER) is the next step once there are lysimeter data.</li>
<li>Dissolution rates are far from equilibrium; no saturation effects, no armouring, no secondary coatings other than the bulk secondary-loss share.</li>
<li>Minerals are assumed liberated in every size class. Real quarry fines carry composite grains; the BET and SEM work will show how far off that is.</li>
<li>The bicarbonate share follows carbonate speciation at bulk soil pH. Real soil CO<sub>2</sub>, organic acids and the Al chemistry of andosols will shift it; Juan should review this rule against the Isometric deduction for non-carbonic-acid weathering.</li>
<li>Exchange-retained cations are not re-released later. That understates long-run removal.</li>
<li>Metal limits are a screening reference; the applicable thresholds are to be confirmed.</li>
<li>Basaltic glass constants are approximate; glass carries a small share of cations in the placeholder mode but may carry more in the Tecnovia fines.</li>
</ul>

<h2>Appendix A. Main equations</h2>
<p class="eq">Cation charge (eq/kg) = Σ 2 · w<sub>oxide</sub> · 10 / M<sub>oxide</sub> for CaO, MgO, Na<sub>2</sub>O, K<sub>2</sub>O</p>
<p class="eq">CO<sub>2</sub> potential (kg/t) = [charge − 2 · S · 10/32.06] · 44.01 − C<sub>inorg</sub> · 10 · 44.01/12.011, then × downstream retention</p>
<p class="eq">Grain retreat: dr/dt = − r<sub>mineral</sub>(pH, T) · λ · V<sub>m</sub> · F; mass left in a class = (r/r<sub>0</sub>)<sup>3</sup></p>
<p class="eq">Bicarbonate share: f<sub>b</sub> = 1 / (1 + 10<sup>pKa1 − pH</sup>)</p>
<p class="eq">pH step: ΔpH = (non-bicarbonate cations − acid load) / (M<sub>soil</sub> · pHBC); M<sub>soil</sub> = depth · 10<sup>4</sup> m<sup>2</sup> · bulk density</p>
<p class="eq">Mixing: f = rate / (rate + M<sub>soil</sub>); C<sub>after</sub> = (C<sub>soil</sub> · M<sub>soil</sub> + C<sub>rock</sub> · rate) / (M<sub>soil</sub> + rate)</p>
<p class="eq">Dissolved fraction: τ = 1 − [C<sub>post</sub> − (1 − f) · C<sub>base</sub>] / (f · C<sub>rock</sub>), with f from Ti or from the measured dose</p>
<p class="eq">Net removal = exported bicarbonate · retention · 44.01 − rock carbon − sulfide acidity − (haul + loading + spreading + quarry) emissions</p>

<h2>Appendix B. Parameters still placeholders ({len(ph)})</h2>
{tbl(ph, ['section','parameter','value','note'], ['Section','Parameter','Value','Note'], cls='tiny')}
<p class="foot">Full parameter list with status: outputs/parameter_status.csv. Sources: Palandri, J.L. and Kharaka, Y.K. (2004) USGS Open File Report 2004-1068; Gislason, S.R. and Oelkers, E.H. (2003) Geochimica et Cosmochimica Acta 67; Isometric, Enhanced Weathering in Agriculture v1.3 and Rock and Mineral Feedstock Characterization v1.2; Decreto-Lei 276/2009, Annex I.</p>
</body></html>""")
    out_html = os.path.join(REP, "Ilha_ERW_Model_Report.html")
    with open(out_html, "w") as fh:
        fh.write("".join(S))
    out_pdf = os.path.join(REP, "Ilha_ERW_Model_Report.pdf")
    try:
        subprocess.run(["python3", os.path.join(REP, "topdf.py"), out_html, out_pdf,
                        "Ilha Carbon · Lomba da Maia ERW model · v1, placeholder run"], check=True)
    except Exception as e:
        print("  PDF step failed (HTML is still there):", e)
        return out_html
    return out_pdf
