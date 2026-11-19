"""Statistics: power by design (effect sizes from the geochemical model, noise from screening)
and the analysis that will be run on the real plot data.

Responses (change from each plot's own baseline):
  d_pH        soil pH (H2O), 0-10 cm composite
  d_CaMg      exchangeable Ca + Mg, cmolc/kg, 0-10 cm composite
  d_totCaMg   total Ca + Mg (eq/kg) from fusion XRF, mean of individual cores: shows the rock is there
Main analysis: regression on rate as a continuous dose with blocks. Curvature tested with a
quadratic term and a saturating curve. MANOVA (Pillai) as a joint check. Each rate against the
control with Holm correction. Repeated rounds: linear mixed model with plot as random effect.
"""
import itertools
import warnings
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.multivariate.manova import MANOVA
from scipy.optimize import curve_fit
from .common import OUT, save_csv, mpl, LINES, DASH
from . import geochem as G

warnings.filterwarnings("ignore")
RESP = ["d_pH", "d_CaMg", "d_totCaMg"]
DESIGNS = {
    "A: 5 rates x 2 blocks": ([0, 10, 20, 40, 80], 2),
    "B: 4 rates x 3 blocks": ([0, 20, 40, 80], 3),
    "C: 5 rates x 3 blocks": ([0, 10, 20, 40, 80], 3),
    "D: 4 rates x 4 blocks": ([0, 20, 40, 80], 4),
}
NOISE_MULT = {"low": 0.6, "medium": 1.0, "high": 1.6}


def effects(cfg, eqkg, rates, months):
    """Expected change from baseline for each response at each rate, central parameters."""
    sim = G.simulate(cfg, rates, months, eqkg)
    soil, ox = cfg.v("soil.elements"), cfg.v("feedstock.oxides_pct")
    Cs = 2 * soil["CaO_pct"] * 10 / 56.077 + 2 * soil["MgO_pct"] * 10 / 40.304
    Cr = 2 * ox["CaO"] * 10 / 56.077 + 2 * ox["MgO"] * 10 / 40.304
    Ms = cfg.v("soil.mixing_depth_cm") / 100 * 1e4 * cfg.v("soil.bulk_density_g_cm3")
    out = {}
    for j, r in enumerate(rates):
        f = r / (r + Ms)
        tau = sim["frac_dissolved"][0, j, months]
        out[r] = {"d_pH": sim["pH"][0, j, months] - sim["pH"][0, j, 0],
                  "d_CaMg": sim["dCaMg_cmolc_kg"][0, j, months],
                  "d_totCaMg": f * (Cr * (1 - tau) - Cs)}
    return out, sim


def noise_sd(cfg):
    cv = cfg.v("soil.heterogeneity_cv")
    bp = cfg.v("soil.between_plot_sd")
    n_comp = cfg.v("design.cores_per_composite")
    n_ind = cfg.v("design.individual_cores_per_plot")
    pH = cfg.v("soil.pH_H2O")
    camg = cfg.v("soil.CEC_cmolc_kg") * cfg.v("soil.base_saturation_frac")
    soil = cfg.v("soil.elements")
    Cs = 2 * soil["CaO_pct"] * 10 / 56.077 + 2 * soil["MgO_pct"] * 10 / 40.304
    acv = cfg.v("soil.analytical_cv")
    meas = lambda sd_core, n, level: np.sqrt(2 * (sd_core ** 2 / n + (acv * level) ** 2))
    cvC = np.sqrt((cv["CaO_pct"] ** 2 + cv["MgO_pct"] ** 2) / 2)
    return {"d_pH": float(np.hypot(bp["d_pH"], meas(cv["pH"] * pH, n_comp, pH))),
            "d_CaMg": float(np.hypot(bp["d_CaMg"], meas(cv["exch_CaMg"] * camg, n_comp, camg))),
            "d_totCaMg": float(np.hypot(0.05, meas(cvC * Cs, n_ind, Cs)))}


def _sim_plots(rates, blocks, eff, sd, rng, corr=0.4):
    rows = [{"block": f"B{b+1}", "rate_t_ha": r} for b in range(blocks) for r in rates]
    df = pd.DataFrame(rows)
    C = np.full((3, 3), corr); np.fill_diagonal(C, 1.0)
    z = rng.multivariate_normal(np.zeros(3), C, size=len(df))
    bz = {b: rng.normal(0, 0.5, 3) for b in df.block.unique()}
    for i, k in enumerate(RESP):
        df[k] = [eff[r][k] for r in df.rate_t_ha] + (np.array([bz[b][i] for b in df.block]) + z[:, i]) * sd[k]
    return df


def _tests(df):
    p = {}
    for k in RESP:
        m = smf.ols(f"{k} ~ rate_t_ha + C(block)", data=df).fit()
        p[k] = m.pvalues["rate_t_ha"]
    mv = MANOVA.from_formula("d_pH + d_CaMg + d_totCaMg ~ rate_t_ha + C(block)", data=df).mv_test()
    p["MANOVA"] = mv.results["rate_t_ha"]["stat"].loc["Pillai's trace", "Pr > F"]
    # 40 t/ha vs control on pH, Holm across the non-zero rates
    m = smf.ols("d_pH ~ C(rate_t_ha) + C(block)", data=df).fit()
    names = [n for n in m.params.index if n.startswith("C(rate_t_ha)")]
    ps = sorted([(m.pvalues[n], n) for n in names])
    adj = {}
    for i, (pv, n) in enumerate(ps):
        adj[n] = min(1.0, max(pv * (len(ps) - i), max(adj.values()) if adj else 0))
    key = [n for n in names if "[T.40" in n]
    p["pH_40_vs_control_Holm"] = adj[key[0]] if key else np.nan
    return p


def power(cfg, eqkg):
    months = cfg.v("design.months")
    nsim = cfg.v("uncertainty.n_power_sims")
    eff, _ = effects(cfg, eqkg, [0, 10, 20, 40, 80], months)
    sd0 = noise_sd(cfg)
    rng = np.random.default_rng(7)
    rows = []
    for name, (rates, blocks) in DESIGNS.items():
        for nl, mult in NOISE_MULT.items():
            sd = {k: v * mult for k, v in sd0.items()}
            hits = {}
            for _ in range(nsim):
                p = _tests(_sim_plots(rates, blocks, eff, sd, rng))
                for k, v in p.items():
                    hits[k] = hits.get(k, 0) + (v < 0.05)
            rows.append({"design": name, "plots": len(rates) * blocks, "noise": nl,
                         **{f"power_{k}": hits[k] / nsim for k in hits}})
    tab = pd.DataFrame(rows).round(3)
    save_csv(tab, "power_by_design.csv")
    eff_tab = pd.DataFrame([{"rate_t_ha": r, **{k: round(v, 4) for k, v in e.items()}} for r, e in eff.items()])
    save_csv(eff_tab, "expected_effects_30m.csv")
    sd_tab = pd.DataFrame([{"response": k, "plot_level_sd_medium": round(v, 4)} for k, v in sd0.items()])
    save_csv(sd_tab, "plot_level_noise.csv")
    _fig_power(tab)
    return tab, eff_tab, sd_tab


def _fig_power(tab):
    plt = mpl()
    cols = [("power_d_pH", "pH"), ("power_d_CaMg", "Exch. Ca+Mg"), ("power_d_totCaMg", "Total Ca+Mg"),
            ("power_MANOVA", "MANOVA"), ("power_pH_40_vs_control_Holm", "pH: 40 t/ha vs control")]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)
    greys = ["#111", "#555", "#888", "#bbb", "#fff"]
    for ax, nl in zip(axes, NOISE_MULT):
        sub = tab[tab.noise == nl]
        x = np.arange(len(sub))
        w = 0.16
        for i, (c, lab) in enumerate(cols):
            ax.bar(x + (i - 2) * w, sub[c], w, label=lab, color=greys[i], edgecolor="#000", lw=0.4)
        ax.axhline(0.8, color="#000", lw=0.8, ls="--")
        ax.set_xticks(x); ax.set_xticklabels([f"{d.split(':')[0]}\n{p} plots" for d, p in zip(sub.design, sub.plots)])
        ax.set_title(f"{nl.capitalize()} variability", fontsize=9.5); ax.set_ylim(0, 1)
    axes[0].set_ylabel("Power at 30 months (dashed = 80%)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, fontsize=8)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(f"{OUT}/fig06_power_by_design.png")
    plt.close(fig)


# ------------------------------------------------------------------ example data and analysis
def example_dataset(cfg, eqkg, plots, seed=21):
    """Simulated field results for the actual plot list, every sampling round, so the analysis code
    can be run and checked before real data exist. NOT DATA."""
    rng = np.random.default_rng(seed)
    rounds = cfg.v("design.sampling_months")
    rates = sorted(plots.rate_t_ha.unique())
    sd = noise_sd(cfg)
    rows = []
    bz = {b: rng.normal(0, 0.5, 3) for b in plots.block.unique()}
    pz = {p: rng.normal(0, 1, 3) for p in plots["plot"]}
    for m in rounds:
        eff, _ = effects(cfg, eqkg, rates, max(m, 0)) if m > 0 else ({r: {k: 0.0 for k in RESP} for r in rates}, None)
        for _, p in plots.iterrows():
            row = {"plot": p["plot"], "block": p.block, "rate_t_ha": p.rate_t_ha, "month": m}
            for i, k in enumerate(RESP):
                if m == 0:
                    row[k] = 0.0
                else:
                    trend = (bz[p.block][i] + 0.6 * pz[p["plot"]][i]) * sd[k] * min(1, m / 30)
                    row[k] = eff[p.rate_t_ha][k] + trend + rng.normal(0, 0.6 * sd[k])
            rows.append(row)
    df = pd.DataFrame(rows).round(4)
    save_csv(df, "example_simulated_field_results.csv")
    return df


def _sat(x, a, b):
    return a * (1 - np.exp(-b * x))


def analyse(df, label="analysis", log=print):
    """Run the full analysis on a results table. Columns: plot, block, rate_t_ha, [month], and any of RESP."""
    lines = []
    P = lambda s="": lines.append(s)
    resp = [k for k in RESP if k in df.columns]
    has_time = "month" in df.columns and df.month.nunique() > 1
    final = df[df.month == df.month.max()] if has_time else df
    final = final.copy()
    P(f"Plots: {final['plot'].nunique()}, blocks: {final.block.nunique()}, rates: {[int(x) for x in sorted(final.rate_t_ha.unique())]}"
      + (f", rounds: {[int(x) for x in sorted(df.month.unique())]}" if has_time else ""))
    P("")
    P("1. Dose-response regression at the final round:  y = b0 + b1*rate + block")
    fits = {}
    for k in resp:
        m = smf.ols(f"{k} ~ rate_t_ha + C(block)", data=final).fit()
        q = smf.ols(f"{k} ~ rate_t_ha + I(rate_t_ha**2) + C(block)", data=final).fit()
        lo, hi = m.conf_int().loc["rate_t_ha"]
        fits[k] = m
        P(f"   {k:10s} slope per t/ha {m.params['rate_t_ha']:+.5f} (95% CI {lo:+.5f} to {hi:+.5f}), p = {m.pvalues['rate_t_ha']:.4f};"
          f" effect at 40 t/ha {40*m.params['rate_t_ha']:+.3f}; R2 {m.rsquared:.2f}")
        P(f"   {'':10s} curvature (quadratic term) p = {q.pvalues['I(rate_t_ha ** 2)']:.3f}; AIC linear {m.aic:.1f}, quadratic {q.aic:.1f}")
        # saturating fit on block-centred values
        y = final[k] - final.groupby("block")[k].transform("mean") + final[k].mean()
        try:
            popt, _ = curve_fit(_sat, final.rate_t_ha.values, y.values, p0=[y.max() * 1.5 + 1e-3, 0.02], maxfev=20000)
            res = y.values - _sat(final.rate_t_ha.values, *popt)
            n, kpar = len(y), 2
            aic = n * np.log(np.sum(res ** 2) / n) + 2 * kpar
            if popt[1] < 1e-3 or popt[1] > 0.2:
                P(f"   {'':10s} saturating fit: not identifiable from these data (no clear plateau within the rates tested)")
            else:
                P(f"   {'':10s} saturating fit a(1-exp(-b*rate)): a = {popt[0]:.3f}, b = {popt[1]:.4f} (half of the plateau at {np.log(2)/popt[1]:.0f} t/ha)")
        except Exception as e:
            P(f"   {'':10s} saturating fit did not converge ({e.__class__.__name__})")
    P("")
    if len(resp) >= 2:
        f = " + ".join(resp) + " ~ rate_t_ha + C(block)"
        mv = MANOVA.from_formula(f, data=final).mv_test()
        st = mv.results["rate_t_ha"]["stat"].loc["Pillai's trace"]
        P(f"2. MANOVA (joint check), Pillai's trace = {st['Value']:.3f}, F = {st['F Value']:.2f}, p = {st['Pr > F']:.4f}")
        P("")
    P("3. Each rate against the control (Holm-adjusted)")
    for k in resp:
        m = smf.ols(f"{k} ~ C(rate_t_ha) + C(block)", data=final).fit()
        names = [n for n in m.params.index if n.startswith("C(rate_t_ha)")]
        ps = sorted([(m.pvalues[n], n) for n in names])
        prev = 0
        for i, (pv, n) in enumerate(ps):
            adj = min(1.0, max(pv * (len(ps) - i), prev)); prev = adj
            P(f"   {k:10s} {n.replace('C(rate_t_ha)[T.', '').replace(']', '')} t/ha vs 0: diff {m.params[n]:+.3f}, Holm p = {adj:.4f}")
    P("")
    if has_time:
        P("4. Repeated rounds: linear mixed model  y ~ rate * month + block, random intercept per plot")
        d2 = df[df.month > 0]
        for k in resp:
            try:
                mm = smf.mixedlm(f"{k} ~ rate_t_ha * month + C(block)", d2, groups=d2["plot"]).fit(reml=True)
                b = mm.params.get("rate_t_ha:month", np.nan)
                pv = mm.pvalues.get("rate_t_ha:month", np.nan)
                P(f"   {k:10s} rate x month {b:+.6f} per t/ha per month, p = {pv:.4f} (does the rate effect grow over time)")
            except Exception as e:
                P(f"   {k:10s} mixed model failed: {e}")
    txt = "\n".join(lines)
    with open(f"{OUT}/{label}.txt", "w") as fh:
        fh.write(txt)
    _fig_dose(final, resp, fits, label)
    return txt


def _fig_dose(final, resp, fits, label):
    plt = mpl()
    titles = {"d_pH": "Change in soil pH", "d_CaMg": "Change in exch. Ca+Mg (cmolc/kg)", "d_totCaMg": "Change in total Ca+Mg (eq/kg)"}
    fig, axes = plt.subplots(1, len(resp), figsize=(3.7 * len(resp), 3.2))
    xs = np.linspace(0, final.rate_t_ha.max(), 50)
    for ax, k in zip(np.atleast_1d(axes), resp):
        for i, (b, g) in enumerate(final.groupby("block")):
            ax.scatter(g.rate_t_ha + (i - 1.5) * 1.2, g[k], s=16, facecolor=["#000", "#fff", "#777", "#bbb"][i % 4],
                       edgecolor="#000", lw=0.5, label=b)
        m = fits[k]
        ax.plot(xs, m.params["Intercept"] + m.params["rate_t_ha"] * xs + np.mean([m.params.get(n, 0) for n in m.params.index if "block" in n]) * 0,
                color="#000", lw=1)
        ax.set_title(titles.get(k, k)); ax.set_xlabel("Rate (t/ha)")
    np.atleast_1d(axes)[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig07_{label}_dose_response.png")
    plt.close(fig)
