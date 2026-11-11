"""
Per-sample runs: one model run for each soil composite in a table (the Lomba da Maia screening,
data/lab/09, once it is filled), in parallel, plus one run on the site mean.

    from erw.samples import sample_configs, run_many
    cfgs, skipped = sample_configs("09")
    results = run_many(cfgs)

Rows without pH, CEC and exchangeable Ca, Mg and K are skipped and listed; nothing is filled in.
"""
from __future__ import annotations

import copy
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor

from . import data_io


def sample_configs(table: str = "09", psd: str = data_io.AS_IS_PSD, rock: str = data_io.DEFAULT_ROCK,
                   overrides: dict | None = None, site_mean: bool = True, lab_dir=data_io.LAB):
    """Model configurations for every usable soil row of `table`, plus the site mean.

    Returns (configs, skipped): configs = [(name, cfg), ...]; skipped = [(sample_id, [missing columns]), ...].
    Field duplicates and the forest zone get their own run but are left out of the mean.
    """
    usable = data_io.soils(lab_dir, tables=(table,))
    skipped = data_io.soil_gaps(lab_dir, table)
    configs = []
    for sid in usable:
        cfg = data_io.build_config(sid, psd, rock, lab_dir)
        cfg.update(overrides or {})
        configs.append((sid, cfg))
    members = [s for s in usable.values() if data_io.in_site_mean(s)]
    if site_mean and len(members) >= 2:
        mean = data_io.site_mean_soil(members)
        mean["table"] = members[0]["table"]
        cfg = data_io.build_config(members[0]["sample_id"], psd, rock, lab_dir)
        cfg.update(overrides or {})
        cfg["soil"], cfg["soil_key"] = mean, "SITE-MEAN"
        cfg["provenance"]["soil"] = data_io._prov(
            {k: round(mean[k], 4) for k in ("pH", "cec", "ca", "mg", "k")}, mean["status_class"],
            f"{mean['table']}: arithmetic mean of {', '.join(mean['members'])}", "pH; cmol(+)/kg")
        for key in ("effective_acidity_cmol_kg", "bulk_density_kg_m3"):      # keep these consistent with the mean
            if key == "effective_acidity_cmol_kg" and "acidity_kcl" in mean:
                cfg[key] = mean["acidity_kcl"]
            if key == "bulk_density_kg_m3" and "bulk_density_g_cm3" in mean:
                cfg[key] = mean["bulk_density_g_cm3"] * 1000.0
        configs.append(("SITE-MEAN", cfg))
    return configs, skipped


def _run_one(item):
    from .model import run            # imported in the worker process
    name, cfg = item
    return name, cfg, run(cfg)


def run_many(configs, max_workers: int | None = None):
    """Run [(name, cfg), ...] in a process pool. Returns [(name, cfg, result), ...] in the same order."""
    if not configs:
        return []
    workers = max_workers or min(len(configs), os.cpu_count() or 1)
    if workers <= 1:
        return [_run_one(c) for c in configs]
    # "spawn" starts clean worker processes on every platform (the default on macOS and Windows)
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as ex:
        return list(ex.map(_run_one, [(n, copy.deepcopy(c)) for n, c in configs]))


def summary_row(name, cfg, res) -> dict:
    """One line of results for a run, with the status of the soil and rock data behind it."""
    h = res["hist"]
    i = h["iso"][-1]
    s = cfg["soil"]
    prov = cfg.get("provenance", {})
    return {"sample": name, "soil_status": s.get("status_class", ""),
            "rock": cfg.get("rock_key", ""), "rock_status": prov.get("oxides", {}).get("status", ""),
            "pH": round(s["pH"], 2), "CEC_cmol_kg": round(s["cec"], 2),
            "base_saturation_pct": round((s["ca"] + s["mg"] + s["k"] + s.get("na", 0.0)) / s["cec"] * 100, 1),
            "years": h["year"][-1], "rock_dissolved_pct": round(h["frac_diss"][-1] * 100, 2),
            "released_tCO2e_ha": round(h["released_t"][-1], 4),
            "co2_removed_50cm_t_ha": round(h["export_t"][-1], 4),
            "iso_exported_from_NFZ_t_ha": round(i["exported_from_nfz"], 4),
            "iso_net_sorption_t_ha": round(i["net_sorption"], 4),
            "iso_net_removal_t_ha": round(i["net_removal"], 4),
            "topsoil_pH_change": round(h["ph_top_rock"][-1] - h["ph_top_ctrl"][-1], 3),
            "status": "model forecast, not creditable"}
