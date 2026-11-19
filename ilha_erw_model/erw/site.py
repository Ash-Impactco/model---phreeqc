"""Site and trial layout from the QGIS layers.

Builds the eligible area (candidate units, minus trees + 5 m, minus drainage buffers,
minus steep ground), orients plots along the contour, lays out a randomised complete
block design with blocks along the elevation gradient, and computes tonnage per plot.
"""
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import geometry_mask
from shapely.geometry import Polygon
from shapely import affinity
from shapely.ops import unary_union
from .common import ROOT, OUT, save_csv, mpl, GREYS

CRS = 5015


def _load(cfg):
    d = os.path.join(ROOT, cfg.v("site.site_dir"))
    rd = lambda n: gpd.read_file(os.path.join(d, n)).to_crs(CRS)
    layers = {"units": rd("units_pilot_and_expansion.geojson"), "nospread": rd("no_spread_5m.geojson"),
              "drain": rd("drainage.geojson"), "water": rd("water_points.geojson"), "parcel": rd("parcel.geojson")}
    layers["dem"] = rasterio.open(os.path.join(d, "dem_farm.tif"))
    layers["slope"] = rasterio.open(os.path.join(d, "slope_percent.tif"))
    return layers


def _zonal(ras, geom):
    m = geometry_mask([geom], out_shape=ras.shape, transform=ras.transform, invert=True, all_touched=False)
    a = ras.read(1, masked=True)
    vals = a[m]
    vals = vals.compressed() if hasattr(vals, "compressed") else vals
    return vals[vals > -1000]


def _contour_angle(dem, poly):
    """Angle (deg) of the contour direction, from the mean DEM gradient over the area."""
    z = dem.read(1).astype(float)
    z[z < -1000] = np.nan
    gy, gx = np.gradient(z)
    m = geometry_mask([poly], out_shape=dem.shape, transform=dem.transform, invert=True)
    # raster rows run north to south, so flip gy to get a northward gradient
    gxm, gym = np.nanmean(gx[m]), -np.nanmean(gy[m])
    downslope = np.degrees(np.arctan2(-gym, -gxm))
    return downslope + 90.0, downslope


def _group(D, idx, k):
    """Greedy blocks of k nearest neighbours, then pairwise swaps until no gain. Returns groups and cost."""
    left, groups = list(idx), []
    while left:
        sub = D[np.ix_(left, left)]
        far = left[int(np.argmax(sub.sum(axis=1)))]
        d = D[far, left]
        g = [left[i] for i in np.argsort(d)[:k]]
        groups.append(g)
        left = [i for i in left if i not in g]
    within = lambda g: D[np.ix_(g, g)].sum()
    improved = True
    while improved:
        improved = False
        for a in range(len(groups)):
            for b in range(a + 1, len(groups)):
                for i in range(k):
                    for j in range(k):
                        ga, gb = groups[a][:], groups[b][:]
                        ga[i], gb[j] = gb[j], ga[i]
                        if within(ga) + within(gb) < within(groups[a]) + within(groups[b]) - 1e-6:
                            groups[a], groups[b] = ga, gb
                            improved = True
    return groups, sum(within(g) for g in groups) / 2


def layout(cfg, log=print):
    L = _load(cfg)
    units = L["units"]
    cand = units[units.unit.isin(cfg.v("site.candidate_units"))]
    area = unary_union(cand.geometry)
    elig = area.difference(unary_union(L["nospread"].geometry))
    elig = elig.difference(unary_union(L["drain"].geometry).buffer(cfg.v("site.drainage_buffer_m")))
    ang, downslope = _contour_angle(L["dem"], area)

    Lm, Wm, al = cfg.v("site.plot_length_m"), cfg.v("site.plot_width_m"), cfg.v("site.alley_m")
    maxs, maxc = cfg.v("site.max_plot_slope_pct"), cfg.v("site.max_cell_slope_pct")
    rates, nb = cfg.v("design.rates_t_ha"), cfg.v("design.blocks")
    n_need = len(rates) * nb

    # grid in a frame rotated so x runs along the contour
    c = area.centroid
    best, best_ang = None, ang
    for da in (-20, -10, 0, 10, 20):
      rot = affinity.rotate(elig, -(ang + da), origin=c)
      minx, miny, maxx, maxy = rot.bounds
      for ox in np.linspace(0, Lm + al, 5, endpoint=False):
        for oy in np.linspace(0, Wm + al, 5, endpoint=False):
            cells = []
            x = minx + ox
            while x + Lm < maxx:
                y = miny + oy
                while y + Wm < maxy:
                    r = Polygon([(x, y), (x + Lm, y), (x + Lm, y + Wm), (x, y + Wm)])
                    if r.intersection(rot).area / r.area > 0.98:
                        cells.append(affinity.rotate(r, ang + da, origin=c))
                    y += Wm + al
                x += Lm + al
            if best is None or len(cells) > len(best):
                best, best_ang = cells, ang + da
    rows = []
    for g in best:
        s = _zonal(L["slope"], g)
        e = _zonal(L["dem"], g)
        if len(s) == 0:
            continue
        rows.append({"geometry": g, "slope_mean": float(np.mean(s)), "slope_p90": float(np.percentile(s, 90)),
                     "steep_share": float(np.mean(s > maxc)), "elev_mean": float(np.mean(e)),
                     "elev_min": float(np.min(e)), "elev_max": float(np.max(e))})
    cells = gpd.GeoDataFrame(rows, crs=CRS)
    ok = cells[(cells.slope_mean <= maxs) & (cells.steep_share <= 0.10)].reset_index(drop=True)
    log(f"  eligible area {elig.area/1e4:.1f} ha in {', '.join(cand.unit)}; {len(cells)} candidate cells, {len(ok)} pass slope tests")
    if len(ok) < n_need:
        raise SystemExit(f"Only {len(ok)} valid plots for {n_need} needed: add units or shrink plots in config.json")

    # choose n_need of the valid cells and group them into blocks of k neighbours so that
    # blocks are as compact as possible (sum of within-block distances), with a mild
    # preference for gentle slopes
    k = len(rates)
    xy_all = np.array([[g.centroid.x, g.centroid.y] for g in ok.geometry])
    D_all = np.hypot(*(xy_all[:, None, :] - xy_all[None, :, :]).transpose(2, 0, 1))
    slope_all = ok.slope_mean.values
    import itertools
    n_drop = len(ok) - n_need
    combos = list(itertools.combinations(range(len(ok)), n_drop)) if n_drop <= 4 else None
    rng0 = np.random.default_rng(1)
    if combos is None:
        combos = [tuple(rng0.choice(len(ok), n_drop, replace=False)) for _ in range(3000)]
    best = (None, None, None)
    for drop in combos:
        idx = [i for i in range(len(ok)) if i not in drop]
        groups, c = _group(D_all, idx, k)
        score = c + 10 * slope_all[idx].sum()
        if best[0] is None or score < best[0]:
            best = (score, idx, groups)
    _, sel, groups_abs = best
    tp = ok.iloc[sel].copy().reset_index(drop=True)
    pos = {a: i for i, a in enumerate(sel)}
    groups = [[pos[a] for a in g] for g in groups_abs]
    order = sorted(groups, key=lambda g: -tp.elev_mean.values[g].mean())
    tp["block"] = ""
    for bi, g in enumerate(order):
        tp.loc[g, "block"] = f"B{bi+1}"
    # constrained randomisation: redraw while rate is correlated with position on the slope
    # within blocks (|r| >= 0.3) or with plot slope, so slope position cannot mimic a rate effect
    rng = np.random.default_rng(cfg.v("site.seed"))
    tp = tp.sort_values(["block", "elev_mean"], ascending=[True, False]).reset_index(drop=True)
    erank = tp.groupby("block").elev_mean.rank().values
    for draw in range(1, 1001):
        rr = np.concatenate([rng.permutation(rates) for _ in range(nb)])
        c1 = abs(np.corrcoef(rr, erank)[0, 1])
        c2 = abs(np.corrcoef(rr, tp.slope_mean.values)[0, 1])
        if c1 < 0.3 and c2 < 0.3:
            break
    tp["rate_t_ha"] = rr
    rand_note = f"draw {draw} from seed {cfg.v('site.seed')}; r(rate, elevation rank in block) = {np.corrcoef(rr, erank)[0,1]:+.2f}, r(rate, slope) = {np.corrcoef(rr, tp.slope_mean.values)[0,1]:+.2f}"
    tp["plot"] = [f"P{i+1:02d}" for i in range(len(tp))]
    tp["unit"] = [", ".join(cand[cand.intersects(g)].unit) for g in tp.geometry]
    tp["gross_ha"] = tp.geometry.area / 1e4
    inset = cfg.v("site.net_inset_m")
    tp["net_sampling_ha"] = tp.geometry.buffer(-inset, join_style=2).area / 1e4

    mo = cfg.v("feedstock.moisture_wet_frac")
    con = cfg.v("design.contingency_frac")
    tp["basalt_dry_t"] = tp.rate_t_ha * tp.gross_ha
    tp["basalt_wet_t"] = tp.basalt_dry_t / (1 - mo) if mo is not None else np.nan   # needs feedstock moisture
    tp["catch_trays"] = np.where(tp.rate_t_ha > 0, cfg.v("design.catch_trays_per_plot"), 0)
    cent = tp.geometry.centroid.to_crs(4326)
    tp["lon"], tp["lat"] = cent.x.round(6), cent.y.round(6)

    summary = {
        "plots": len(tp), "blocks": nb, "rates": rates, "plot_grid_rotation_deg": round(best_ang % 180, 1),
        "trial_gross_ha": round(tp.gross_ha.sum(), 2),
        "trial_footprint_ha": round(unary_union(tp.geometry).convex_hull.area / 1e4, 1),
        "basalt_dry_t": round(tp.basalt_dry_t.sum(), 1),
        "basalt_dry_t_with_contingency": round(tp.basalt_dry_t.sum() * (1 + con), 1),
        "basalt_wet_t_with_contingency": round(tp.basalt_wet_t.sum() * (1 + con), 1) if mo is not None else None,
        "truckloads": int(np.ceil(tp.basalt_wet_t.sum() * (1 + con) / cfg.v("design.truck_payload_t_wet"))) if mo is not None else None,
        "elev_range_m": [round(tp.elev_min.min()), round(tp.elev_max.max())],
        "plot_slope_mean_range_pct": [round(tp.slope_mean.min(), 1), round(tp.slope_mean.max(), 1)],
        "units_used": sorted(set(", ".join(tp.unit).split(", "))),
        "randomisation": rand_note,
    }
    # 12-plot fallback: drop the last block
    t12 = tp[tp.block != f"B{nb}"]
    summary["if_3_blocks_basalt_dry_t_with_contingency"] = round(t12.basalt_dry_t.sum() * (1 + con), 1)
    summary["if_3_blocks_gross_ha"] = round(t12.gross_ha.sum(), 2)

    cols = ["plot", "block", "rate_t_ha", "unit", "gross_ha", "net_sampling_ha", "basalt_dry_t", "basalt_wet_t",
            "catch_trays", "slope_mean", "slope_p90", "elev_mean", "elev_min", "elev_max", "lon", "lat"]
    tab = tp[cols].round({"gross_ha": 3, "net_sampling_ha": 3, "basalt_dry_t": 1, "basalt_wet_t": 1,
                          "slope_mean": 1, "slope_p90": 1, "elev_mean": 0, "elev_min": 0, "elev_max": 0})
    save_csv(tab, "trial_plots.csv")
    tp[cols + ["geometry"]].to_crs(4326).to_file(os.path.join(OUT, "trial_plots.geojson"), driver="GeoJSON")
    tp[cols + ["geometry"]].to_file(os.path.join(OUT, "trial_plots_EPSG5015.gpkg"), driver="GPKG")
    _map(L, cand, elig, tp, summary)
    return tab, summary, tp


def _map(L, cand, elig, tp, s):
    plt = mpl()
    from matplotlib.patches import Patch
    dem = L["dem"]
    z = dem.read(1).astype(float)
    z[z < -1000] = np.nan
    gy, gx = np.gradient(z, 15)
    hs = np.cos(np.arctan(np.hypot(gx, gy))) * 0.8 + 0.2 * (-gx - gy) / (np.hypot(gx, gy) + 1e-6) * np.sin(np.arctan(np.hypot(gx, gy)))
    b = dem.bounds
    xmin, ymin, xmax, ymax = cand.total_bounds
    pad = 120
    fig, ax = plt.subplots(figsize=(7.2, 8.6))
    ax.imshow(hs, extent=[b.left, b.right, b.bottom, b.top], cmap="gray", vmin=0.2, vmax=1.15, alpha=0.55)
    ax.contour(np.linspace(b.left, b.right, z.shape[1]), np.linspace(b.top, b.bottom, z.shape[0]), z,
               levels=np.arange(300, 700, 10), colors="#777", linewidths=0.35)
    L["units"].boundary.plot(ax=ax, color="#999", lw=0.6)
    cand.boundary.plot(ax=ax, color="#000", lw=1.0)
    L["nospread"].plot(ax=ax, facecolor="#e6e6e6", edgecolor="#aaa", hatch="....", lw=0.3)
    L["drain"].plot(ax=ax, color="#000", lw=0.8, ls="--")
    for _, r in cand.iterrows():
        p = r.geometry.representative_point()
        ax.text(p.x, p.y, r.unit, fontsize=8, color="#444", style="italic", zorder=1)
    for r, col in GREYS.items():
        sub = tp[tp.rate_t_ha == r]
        if len(sub):
            sub.plot(ax=ax, facecolor=col, edgecolor="#000", lw=0.9, zorder=3)
    for _, r in tp.iterrows():
        p = r.geometry.centroid
        ax.text(p.x, p.y, zorder=4, s=f"{r['plot']}\n{r['rate_t_ha']:g}", ha="center", va="center", fontsize=6.3,
                color="#fff" if r["rate_t_ha"] >= 40 else "#000")
    for blk, g in tp.groupby("block"):
        hull = unary_union(g.geometry).convex_hull.buffer(4)
        ax.plot(*hull.exterior.xy, color="#000", lw=0.6, ls=":")
        tc = unary_union(tp.geometry).centroid
        hc = hull.centroid
        v = np.array([hc.x - tc.x, hc.y - tc.y]); v = v / (np.hypot(*v) + 1e-9)
        ext = np.array(hull.exterior.coords)
        far = ext[np.argmax(ext @ v)]
        ax.text(far[0] + 22 * v[0], far[1] + 22 * v[1], blk, ha="center", va="center", fontsize=9, weight="bold",
                bbox=dict(facecolor="#fff", edgecolor="#000", lw=0.5, pad=1.5))
    ax.set_xlim(xmin - pad, xmax + pad); ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color("#000"); sp.set_linewidth(0.6)
    # scale bar and north arrow
    sx, sy = xmin - pad + 30, ymin - pad + 30
    ax.plot([sx, sx + 200], [sy, sy], color="#000", lw=2); ax.text(sx + 100, sy + 12, "200 m", ha="center", fontsize=7)
    ax.annotate("N", xy=(xmax + pad - 40, ymax + pad - 30), xytext=(xmax + pad - 40, ymax + pad - 110),
                arrowprops=dict(arrowstyle="-|>", color="#000"), ha="center", fontsize=9)
    leg = [Patch(facecolor=GREYS[r], edgecolor="#000", label=f"{r} t/ha") for r in sorted(tp.rate_t_ha.unique())]
    leg += [Patch(facecolor="#e6e6e6", edgecolor="#aaa", hatch="....", label="No-spread: trees + 5 m")]
    ax.legend(handles=leg, loc="lower right", fontsize=7, frameon=True, facecolor="#fff", edgecolor="#000")
    ax.set_title(f"Basalt rate trial layout: {s['plots']} plots, {s['blocks']} blocks, {s['trial_gross_ha']} ha "
                 f"(EPSG:5015, 10 m contours, drainage dashed)", fontsize=8.5)
    fig.savefig(os.path.join(OUT, "fig01_trial_layout.png"))
    plt.close(fig)
