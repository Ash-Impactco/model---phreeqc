"""
Uncertainty of the forecast: Monte Carlo over the model inputs, and a one-at-a-time tornado.

Inputs and ranges come from data/lab (erw.data_io.mc_parameters): the field-to-lab rate factor,
surface roughness, glass fraction, drainage, soil pCO2 and soil temperature. All the ranges are
assumptions, so the spread shows how much the forecast depends on what has not been measured. It is
not a statistical confidence interval.

    summary = monte_carlo(cfg, n=50)            # runs in parallel, cached on disk by input hash
    bars = tornado(cfg)                         # low and high value of each input, one at a time

Percentiles reported: P10, P30, P40, P50, P90. Under Isometric EWA v1.3 credits are issued at a
conservative percentile of the MEASURED removal (P30; P40 with a validation approach). The same
percentiles are shown here for the forecast as an indication only. They are not creditable tonnes.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import multiprocessing
import os
import pickle
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from . import data_io

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "outputs" / "cache"
PERCENTILES = (10, 30, 40, 50, 90)
METRICS = ("released", "exported_from_nfz", "net_sorption", "stored", "net_removal")
METRIC_LABELS = {"released": "Cations released", "exported_from_nfz": "Leaving the Near-Field Zone",
                 "net_sorption": "Held on exchange sites", "stored": "CO2e stored (Eq. 4)",
                 "net_removal": "Net removal (Eq. 1)"}


def _clean(cfg) -> dict:
    """The part of a configuration that decides the result (no provenance text)."""
    return {k: v for k, v in cfg.items() if k not in ("provenance", "psd_detail")}


def input_hash(cfg, extra=None) -> str:
    text = json.dumps([_clean(cfg), extra], sort_keys=True, default=str)
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def sample_inputs(params: list[dict], n: int, seed: int = 1) -> list[dict]:
    """Latin hypercube sample: n sets of input values, each range cut into n equal slices."""
    rng = np.random.default_rng(seed)
    out = [dict() for _ in range(n)]
    for p in params:
        u = (rng.permutation(n) + rng.random(n)) / n
        lo, hi = p["low"], p["high"]
        if p["distribution"].startswith("log") and lo > 0:
            vals = np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo)))
        else:
            vals = lo + u * (hi - lo)
        for i, v in enumerate(vals):
            out[i][p["key"]] = float(v)
    return out


def _series(hist) -> dict:
    iso = hist["iso"]
    return dict(released=[i["weathering"] for i in iso], exported_from_nfz=[i["exported_from_nfz"] for i in iso],
                net_sorption=[i["net_sorption"] for i in iso], stored=[i["stored"] for i in iso],
                net_removal=[i["net_removal"] for i in iso])


def _run_one(item):
    from .model import run
    idx, cfg = item
    try:
        return idx, _series(run(cfg)["hist"]), None
    except Exception as err:                          # a sample PHREEQC cannot solve: reported, not hidden
        return idx, None, str(err)[:200]


def _pool(items, workers, progress):
    done = 0
    results = {}
    workers = workers or min(len(items), os.cpu_count() or 1)
    if workers <= 1:
        for it in items:
            idx, res, err = _run_one(it)
            results[idx] = (res, err)
            done += 1
            if progress:
                progress(done, len(items))
        return results
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as ex:
        for fut in as_completed([ex.submit(_run_one, it) for it in items]):
            idx, res, err = fut.result()
            results[idx] = (res, err)
            done += 1
            if progress:
                progress(done, len(items))
    return results


def monte_carlo(cfg, n: int = 50, seed: int = 1, params: list[dict] | None = None, workers: int | None = None,
                progress=None, use_cache: bool = True) -> dict:
    """Run the model n times with sampled inputs. Returns percentiles per year for each metric.

    Result keys: years, params, samples (the sampled inputs), failed (count), percentiles
    {metric: {"P10": [...], ...}}, final {metric: [value per run]}, base {metric: [...]} is not
    included (run the base case separately). Cached in outputs/cache by a hash of the inputs.
    """
    params = params if params is not None else data_io.mc_parameters()
    key = input_hash(cfg, dict(n=n, seed=seed, params=[(p["key"], p["low"], p["high"], p["distribution"]) for p in params]))
    path = CACHE / f"mc_{key}.pkl"
    if use_cache and path.exists():
        out = pickle.loads(path.read_bytes())
        out["from_cache"] = True
        if progress:
            progress(n, n)
        return out
    samples = sample_inputs(params, n, seed)
    items = []
    for i, s in enumerate(samples):
        c = copy.deepcopy(_clean(cfg))
        c.update(s)
        items.append((i, c))
    results = _pool(items, workers, progress)
    good = [i for i in range(n) if results[i][0] is not None]
    years = list(range(1, int(cfg["years"]) + 1))
    pct, final = {}, {}
    for m in METRICS:
        arr = np.array([results[i][0][m] for i in good]) if good else np.zeros((0, len(years)))
        pct[m] = {f"P{q}": (np.percentile(arr, q, axis=0).tolist() if len(arr) else []) for q in PERCENTILES}
        final[m] = arr[:, -1].tolist() if len(arr) else []
    out = dict(years=years, n=n, seed=seed, params=params, samples=[samples[i] for i in good], failed=n - len(good),
               errors=[results[i][1] for i in range(n) if results[i][1]][:5], percentiles=pct, final=final,
               hash=key, from_cache=False)
    if use_cache:
        CACHE.mkdir(parents=True, exist_ok=True)
        path.write_bytes(pickle.dumps(out))
    return out


def tornado(cfg, metric: str = "exported_from_nfz", params: list[dict] | None = None, workers: int | None = None,
            progress=None, use_cache: bool = True) -> dict:
    """Each input set to the low and the high end of its range in turn, everything else at the base
    value. Returns the base value and one bar per input, widest first."""
    params = params if params is not None else data_io.mc_parameters()
    key = input_hash(cfg, dict(tornado=[(p["key"], p["low"], p["high"]) for p in params]))
    path = CACHE / f"tornado_{key}.pkl"
    if use_cache and path.exists():
        res = pickle.loads(path.read_bytes())
        if progress:
            progress(1, 1)
    else:
        items = [(0, copy.deepcopy(_clean(cfg)))]
        for j, p in enumerate(params):
            for k, end in enumerate(("low", "high")):
                c = copy.deepcopy(_clean(cfg))
                c[p["key"]] = p[end]
                items.append((1 + 2 * j + k, c))
        res = _pool(items, workers, progress)
        if use_cache:
            CACHE.mkdir(parents=True, exist_ok=True)
            path.write_bytes(pickle.dumps(res))
    val = lambda i: res[i][0][metric][-1] if res[i][0] is not None else float("nan")
    bars = []
    for j, p in enumerate(params):
        lo, hi = val(1 + 2 * j), val(2 + 2 * j)
        bars.append(dict(key=p["key"], label=p["label"], low_input=p["low"], high_input=p["high"],
                         at_low=lo, at_high=hi, base_input=cfg.get(p["key"]), spread=abs(hi - lo)))
    bars.sort(key=lambda b: -(b["spread"] if not math.isnan(b["spread"]) else -1))
    return dict(metric=metric, label=METRIC_LABELS[metric], base=val(0), bars=bars)


def percentile_table(mc: dict, year: int | None = None) -> list[dict]:
    """P10 to P90 of every metric at one year (default: the last), for a table."""
    yi = (year or mc["years"][-1]) - 1
    rows = []
    for m in METRICS:
        row = {"Term": METRIC_LABELS[m]}
        for q in PERCENTILES:
            v = mc["percentiles"][m][f"P{q}"]
            row[f"P{q}"] = round(v[yi], 4) if v else None
        rows.append(row)
    return rows


ISSUANCE_NOTE = ("P30 is the level at which Isometric issues credits from MEASURED removal (P40 with a validation "
                 "approach). It is shown here for the forecast as an indication only. These are not creditable tonnes.")
