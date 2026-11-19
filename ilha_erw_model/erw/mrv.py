"""Measurement design: can we see the weathering, and where.

1. Tracer contrast table (Ti, Zr against Ca, Mg).
2. Soil mass-balance detectability: Monte Carlo of individual cores before and after spreading,
   two estimators of the dissolved fraction tau:
     A  tracer-based mixing fraction (Ti), the classic immobile-element approach
     B  dose-based mixing fraction from the measured application rate (catch trays)
3. Water: alkalinity signal at the catchment sampling points W1/W2 and in plot-scale lysimeters.
4. Mesocosms: column masses, irrigation, expected leachate alkalinity against the lab LOQ.
5. Sampling schedule and analysis counts (no prices).
"""
import numpy as np
import pandas as pd
from .common import OUT, save_csv, mpl, LINES, DASH
from . import geochem as G


def _soil_eq(el):
    return 2 * el["CaO_pct"] * 10 / 56.077 + 2 * el["MgO_pct"] * 10 / 40.304


def tracer_contrast(cfg):
    soil, ox, tr = cfg.v("soil.elements"), cfg.v("feedstock.oxides_pct"), cfg.v("feedstock.trace_ppm")
    cv = cfg.v("soil.heterogeneity_cv")
    bd, dep = cfg.v("soil.bulk_density_g_cm3"), cfg.v("soil.mixing_depth_cm") / 100
    Ms = dep * 1e4 * bd
    n = cfg.v("design.individual_cores_per_plot")
    rows = []
    for name, s, r, c in [("TiO2 (%)", soil["TiO2_pct"], ox["TiO2"], cv["TiO2_pct"]),
                          ("Zr (ppm)", soil["Zr_ppm"], tr["Zr"], cv["Zr_ppm"]),
                          ("CaO (%)", soil["CaO_pct"], ox["CaO"], cv["CaO_pct"]),
                          ("MgO (%)", soil["MgO_pct"], ox["MgO"], cv["MgO_pct"])]:
        row = {"element": name, "soil": s, "rock": r, "rock_to_soil": round(r / s, 2), "soil_cv": c}
        for R in (20, 40, 80):
            f = R / (R + Ms)
            shift = f * (r - s)
            row[f"shift_at_{R}"] = round(shift, 4)
            row[f"signal_to_noise_{R}_n{n}"] = round(abs(shift) / (c * s / np.sqrt(n) * np.sqrt(2)), 2)
        rows.append(row)
    tab = pd.DataFrame(rows)
    save_csv(tab, "tracer_contrast.csv")
    return tab


def tau_detectability(cfg, tau_true, S=2000, seed=3):
    rng = np.random.default_rng(seed)
    soil, ox = cfg.v("soil.elements"), cfg.v("feedstock.oxides_pct")
    cv = cfg.v("soil.heterogeneity_cv")
    acv = cfg.v("soil.analytical_cv")
    bd, dep = cfg.v("soil.bulk_density_g_cm3"), cfg.v("soil.mixing_depth_cm") / 100
    Ms = dep * 1e4 * bd
    Ti_s, Ti_r = soil["TiO2_pct"], ox["TiO2"]
    C_s = _soil_eq(soil)
    C_r = 2 * ox["CaO"] * 10 / 56.077 + 2 * ox["MgO"] * 10 / 40.304
    cvC = np.sqrt((cv["CaO_pct"] ** 2 + cv["MgO_pct"] ** 2) / 2)
    spread_cv, tray_cv = 0.25, 0.08          # core-scale spreading unevenness; plot-mean dose error from trays
    rows = []
    for R in (20, 40, 80):
        f = R / (R + Ms)
        for n in (6, 15, 30):
            est = {"A_tracer_Ti": [], "B_dose_trays": []}
            for _ in range(S):
                def cores(nn):
                    z = rng.multivariate_normal([0, 0], [[1, 0.5], [0.5, 1]], nn)
                    return Ti_s * (1 + cv["TiO2_pct"] * z[:, 0]), C_s * (1 + cvC * z[:, 1])
                tb, cb = cores(n)
                tp, cp = cores(n)
                fi = np.clip(f * (1 + spread_cv * rng.standard_normal(n)), 0, None)
                Ti_post = fi * Ti_r + (1 - fi) * tp
                C_post = fi * C_r * (1 - tau_true) + (1 - fi) * cp
                noise = lambda x: x * (1 + acv * rng.standard_normal(x.shape))
                tb, cb, Ti_post, C_post = noise(tb), noise(cb), noise(Ti_post), noise(C_post)
                fA = (Ti_post.mean() - tb.mean()) / (Ti_r - tb.mean())
                fB = f * (1 + tray_cv * rng.standard_normal())
                for k, fh in (("A_tracer_Ti", fA), ("B_dose_trays", fB)):
                    est[k].append(1 - (C_post.mean() - (1 - fh) * cb.mean()) / (fh * C_r) if fh > 1e-4 else np.nan)
            for k, v in est.items():
                v = np.array(v)
                v = v[np.isfinite(v)]
                q1, q3 = np.percentile(v, [25, 75])
                sd = (q3 - q1) / 1.349          # robust spread: the tracer estimator has very heavy tails
                rows.append({"rate_t_ha": R, "cores_per_plot_per_round": n, "estimator": k, "tau_true": round(tau_true, 3),
                             "tau_median": round(float(np.median(v)), 3), "tau_robust_sd_one_plot": round(float(sd), 3),
                             "min_detectable_tau_one_plot": round(2.8 * sd, 3),
                             "min_detectable_tau_4_plots": round(2.8 * sd / 2, 3)})
    tab = pd.DataFrame(rows)
    save_csv(tab, "tau_detectability.csv")
    return tab


def water_signal(cfg, sim, rates, plots, water_pts):
    months = cfg.v("design.months")
    drain_m = (cfg.v("climate.rain_mm_yr") - cfg.v("climate.aet_mm_yr")) / 1000
    yrs = months / 12
    rows = []
    # plot-scale drainage below the mixed layer
    for j, r in enumerate(rates):
        bic_yr = (sim["bicarb"][0, j, months] - sim["bicarb"][0, j, months - 12]) if months >= 12 else sim["bicarb"][0, j, months] / yrs
        cat_yr = ((sim["bicarb"] + sim["leached_strong"])[0, j, months] - (sim["bicarb"] + sim["leached_strong"])[0, j, months - 12])
        L = 1e4 * drain_m * 1000
        rows.append({"where": "plot lysimeter", "rate_t_ha": r, "bicarbonate_eq_ha_yr_in_year_3": round(bic_yr),
                     "alkalinity_increase_meq_L": round(bic_yr / L * 1000, 4),
                     "base_cation_increase_meq_L": round(cat_yr / L * 1000, 4)})
    # whole trial at catchment outlets
    trial_bic = 0.0
    for j, r in enumerate(rates):
        ha = plots[plots.rate_t_ha == r].gross_ha.sum()
        trial_bic += ha * (sim["bicarb"][0, j, months] - sim["bicarb"][0, j, months - 12])
    for _, w in water_pts.iterrows():
        L = w.catchment_ha * 1e4 * drain_m * 1000
        rows.append({"where": w["name"], "catchment_ha": w.catchment_ha, "bicarbonate_eq_ha_yr_in_year_3": round(trial_bic),
                     "alkalinity_increase_meq_L": round(trial_bic / L * 1000, 5)})
    tab = pd.DataFrame(rows)
    save_csv(tab, "water_signal.csv")
    return tab


def mesocosm(cfg, eqkg):
    d = cfg.raw["mesocosm"]
    v = lambda k: d[k]["value"]
    A = np.pi * (v("column_diameter_cm") / 200) ** 2
    rates = v("rates_t_ha")
    months = v("months")
    sim = G.simulate(cfg, rates, months, eqkg, T_C=20.0, acid_load=0.0)
    irr_L_month = v("irrigation_mm_yr") / 1000 / 12 * A * 1000
    bd = cfg.v("soil.bulk_density_g_cm3")
    rows = []
    for j, r in enumerate(rates):
        bic = np.diff(sim["bicarb"][0, j]) / 1e4 * A          # eq per column per month
        cat = np.diff((sim["bicarb"] + sim["leached_strong"])[0, j]) / 1e4 * A
        rows.append({"rate_t_ha": r, "rock_g_per_column": round(r * 100 * A, 1),
                     "soil_kg_per_column": round(A * v("soil_height_cm") / 100 * bd * 1000, 2),
                     "irrigation_L_per_week": round(v("irrigation_mm_yr") / 1000 / 52 * A * 1000, 3),
                     "alk_increase_meq_L_month3": round(bic[2] / irr_L_month * 1000, 4),
                     "alk_increase_meq_L_month12": round(bic[-1] / irr_L_month * 1000, 4),
                     "base_cations_meq_L_month12": round(cat[-1] / irr_L_month * 1000, 4),
                     "frac_cations_dissolved_12m": round(sim["frac_dissolved"][0, j, -1], 4),
                     "pH_0_10cm_at_12m": round(sim["pH"][0, j, -1], 2)})
    tab = pd.DataFrame(rows)
    n_col = len(rates) * v("replicates") + v("sand_controls")
    meta = {"columns": n_col, "column_area_m2": round(A, 5), "loq_meq_L": v("alkalinity_loq_meq_L"),
            "weekly_leachate_samples": n_col * 4 * months // 1 // 1}
    meta["weekly_leachate_samples"] = n_col * 52 * months // 12
    save_csv(tab, "mesocosm_design.csv")
    return tab, meta


def sampling_plan(cfg, plots):
    rounds = cfg.v("design.sampling_months")
    nplots = len(plots)
    ind = cfg.v("design.individual_cores_per_plot")
    deep_rounds = [0, 12, rounds[-1]]
    rows = []
    for m in rounds:
        rows.append({"month": m, "item": "Composite 0-10 cm (15 cores): pH, exch. Ca Mg K Na Al, CEC, base sat., SOC, total N",
                     "samples": nplots})
        rows.append({"month": m, "item": f"Individual cores 0-10 cm ({ind}/plot): total elements by fusion XRF/ICP (Ca, Mg, Na, K, Ti, Zr, trace metals), inorganic C",
                     "samples": nplots * ind})
        if m in deep_rounds:
            rows.append({"month": m, "item": "Composites 10-20 and 20-30 cm: same as 0-10 composite plus total elements", "samples": nplots * 2})
        if m in (0, rounds[-1]):
            rows.append({"month": m, "item": "Bulk density rings 0-10 cm (2/plot)", "samples": nplots * 2})
    blocks = plots.block.nunique()
    rows.append({"month": "monthly", "item": f"Suction-cup lysimeters at 30 cm, 2 per plot in the 0 and 80 t/ha plots of each block ({blocks*2} plots): alkalinity, pH, DIC, Ca Mg Na K, Si, NO3, SO4, Cl, EC",
                 "samples": blocks * 2 * 2 * cfg.v("design.months")})
    rows.append({"month": "quarterly", "item": "Stream water W1 and W2 (context only, see water signal)", "samples": 2 * cfg.v("design.months") // 3})
    rows.append({"month": "each cut (about 6/yr)", "item": "Pasture yield (quadrats) and tissue Ca, Mg, K, P, Ni, Cr, Co (feed safety)",
                 "samples": nplots * 6 * cfg.v("design.months") // 12})
    rows.append({"month": "spreading", "item": "Catch trays for application rate", "samples": int(plots.catch_trays.sum())})
    rows.append({"month": "delivery", "item": "Feedstock moisture and XRF check per 50 t delivered", "samples": int(np.ceil(plots.basalt_wet_t.sum() * 1.1 / 50))})
    tab = pd.DataFrame(rows)
    save_csv(tab, "sampling_schedule_and_counts.csv")
    return tab


def fig_tau(tab):
    plt = mpl()
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.3), sharey=True)
    for a, est in zip(ax, ["A_tracer_Ti", "B_dose_trays"]):
        sub = tab[tab.estimator == est]
        for i, n in enumerate(sorted(sub.cores_per_plot_per_round.unique())):
            s = sub[sub.cores_per_plot_per_round == n]
            a.plot(s.rate_t_ha, s.min_detectable_tau_4_plots, marker="o", color=LINES[i], ls=DASH[i], label=f"{n} cores per plot")
        a.axhline(tab.tau_true.iloc[0], color="#000", lw=0.6, ls=":")
        a.text(21, tab.tau_true.iloc[0] * 1.08, "expected at 30 months (central)", fontsize=7)
        a.set_yscale("log"); a.set_xlabel("Application rate (t/ha)"); a.set_xticks([20, 40, 80])
        a.set_title({"A_tracer_Ti": "A. Mixing from Ti tracer", "B_dose_trays": "B. Mixing from measured dose"}[est])
        a.legend(fontsize=7)
    ax[0].set_ylabel("Smallest detectable dissolved fraction\n(4 plots per rate, 80% power)")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig05_tau_detectability.png")
    plt.close(fig)
