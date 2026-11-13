"""
Project emissions (life cycle) per hectare, for the Emissions term of Isometric EWA v1.3, Equation 1.

Five items, each from a factor in data/lab/18_emission_factors.csv with its status and source:
quarrying (the share allocated to the fines), milling, haulage, spreading, and sampling and lab
logistics. Every factor is a placeholder until fuel, electricity and weighbridge records exist.

    table, total = life_cycle(cfg)        # one row per item, t CO2e per hectare
    year = payback_year(hist)             # first year in which the stored CO2 exceeds the emissions

The arithmetic is the same as erw.model.project_emissions_t_ha, which the model uses for the net
removal. This module adds the activity data, the factor, its status and its source for each line.
"""
from __future__ import annotations

import pandas as pd


def _p(cfg, key):
    """(status, source) of one emission factor from the provenance of the run."""
    p = (cfg.get("provenance") or {}).get(f"emissions.{key}", {})
    return p.get("status", "assumed"), p.get("source", "data/inputs.json (fallback)")


def life_cycle(cfg) -> tuple[pd.DataFrame, float]:
    """Emissions per hectare by item, with the activity, factor, status and source of each."""
    from .model import project_emissions_t_ha
    e = cfg["emissions"]
    t_ha = float(cfg["application_t_ha"])
    psd = cfg.get("psd_key", "as_is_0-4mm")
    kwh_t = e["milling_kwh_per_t"].get(psd, 0)
    trips = t_ha / e["truck_payload_t"]
    litres_haul = trips * e["haul_distance_km"] * (e["diesel_l_per_km_loaded"] + e["diesel_l_per_km_empty"])
    out = project_emissions_t_ha(cfg)
    worst = lambda *keys: ("assumed" if any(_p(cfg, k)[0] == "assumed" for k in keys) else _p(cfg, keys[0])[0])
    rows = [
        dict(item="Quarrying (share allocated to the fines)",
             activity=f"{t_ha:g} t rock/ha", factor=f"{e['quarrying_kgco2e_per_t']:g} kg CO2e per t",
             t_co2e_ha=out["quarrying"], status=_p(cfg, "quarrying_kgco2e_per_t")[0],
             source=_p(cfg, "quarrying_kgco2e_per_t")[1],
             note="0 treats the fines as a by-product with no share of the quarry's emissions. This must be justified."),
        dict(item="Milling", activity=f"{t_ha:g} t/ha x {kwh_t:g} kWh/t = {t_ha * kwh_t:g} kWh/ha",
             factor=f"{e['grid_kgco2e_per_kwh']:g} kg CO2e per kWh",
             t_co2e_ha=out["milling"], status=worst("grid_kgco2e_per_kwh", f"milling.{psd}"),
             source=_p(cfg, "grid_kgco2e_per_kwh")[1],
             note="No milling for the fines as supplied." if kwh_t == 0 else "Milling energy is an assumed figure, not a milling trial."),
        dict(item="Haulage, quarry to field and back",
             activity=f"{trips:g} loads/ha x {e['haul_distance_km']:g} km x "
                      f"({e['diesel_l_per_km_loaded']:g} + {e['diesel_l_per_km_empty']:g}) L/km = {litres_haul:.1f} L diesel/ha",
             factor=f"{e['diesel_kgco2e_per_l']:g} kg CO2e per L",
             t_co2e_ha=out["haulage"],
             status=worst("haul_distance_km", "truck_payload_t", "diesel_l_per_km_loaded", "diesel_kgco2e_per_l"),
             source=_p(cfg, "haul_distance_km")[1],
             note="Diesel factor is combustion only. Isometric may require a well-to-wheel factor."),
        dict(item="Spreading", activity=f"{t_ha:g} t/ha x {e['spreading_diesel_l_per_t']:g} L/t = "
                                        f"{t_ha * e['spreading_diesel_l_per_t']:g} L diesel/ha",
             factor=f"{e['diesel_kgco2e_per_l']:g} kg CO2e per L",
             t_co2e_ha=out["spreading"], status=worst("spreading_diesel_l_per_t", "diesel_kgco2e_per_l"),
             source=_p(cfg, "spreading_diesel_l_per_t")[1], note=""),
        dict(item="Sampling and lab logistics", activity=f"{e['other_project_t_co2e']:g} t CO2e for the whole project",
             factor="spread over 100 ha",
             t_co2e_ha=out["other"], status=_p(cfg, "other_project_t_co2e")[0],
             source=_p(cfg, "other_project_t_co2e")[1],
             note="Sampling trips, lab shipping and staff travel: a rough allowance."),
    ]
    df = pd.DataFrame(rows)
    return df, float(out["total"])


def payback_year(hist):
    """First year in which the CO2e stored (Isometric Eq. 4) reaches the project emissions.
    Returns (year or None, emissions, stored at the end). None means it is not reached within the run."""
    iso = hist["iso"]
    for i in iso:
        if i["stored"] >= i["emissions"]:
            return i["year"], i["emissions"], iso[-1]["stored"]
    return None, iso[-1]["emissions"], iso[-1]["stored"]


def payback_text(hist) -> str:
    year, emis, stored = payback_year(hist)
    n = hist["iso"][-1]["year"]
    if year is not None:
        return f"In this forecast the stored CO2e passes the project emissions ({emis:.3f} t CO2e/ha) in year {year}."
    share = stored / emis * 100 if emis else 0.0
    return (f"In this forecast the project emissions ({emis:.3f} t CO2e/ha) are not paid back within {n} years: "
            f"the stored CO2e reaches {stored:.3f} t/ha, {share:.0f} % of the emissions.")
