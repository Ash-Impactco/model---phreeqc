"""Emissions life cycle, Monte Carlo and tornado, detectability, tracer mass balance, database comparison,
and the Isometric settings (Near-Field Zone depth, control area)."""
import numpy as np
import pandas as pd
import pytest

from erw import data_io, dbcompare, detect, emissions, uncertainty, validate
from erw.model import default_config, project_emissions_t_ha, run


@pytest.fixture(scope="module")
def short_run():
    cfg = default_config("SG-04", "milled_p80_100um")
    cfg.update(acidity_mode="effective", glass_frac=0.2, field_factor=0.3, years=3)
    return cfg, run(cfg)


# ---------------------------------------------------------------- emissions
def test_life_cycle_lines_add_up_and_are_all_placeholders():
    cfg = default_config()
    table, total = emissions.life_cycle(cfg)
    assert list(table["item"].str.split(" ").str[0]) == ["Quarrying", "Milling", "Haulage,", "Spreading", "Sampling"]
    assert table["t_co2e_ha"].sum() == pytest.approx(total) == pytest.approx(project_emissions_t_ha(cfg)["total"])
    assert total == pytest.approx(0.1704, abs=1e-4)                       # 0.097 haulage + 0.054 spreading + 0.02 other
    assert set(table["status"]) == {"assumed"} and all(table["source"].str.contains("18_emission_factors"))
    milled, total_milled = emissions.life_cycle(default_config(psd="milled_p80_100um"))
    assert milled.set_index("item").loc["Milling", "t_co2e_ha"] == pytest.approx(40 * 30 * 0.45 / 1000)
    assert total_milled == pytest.approx(0.7104, abs=1e-4)


def test_payback_year(short_run):
    cfg, res = short_run
    year, emis, stored = emissions.payback_year(res["hist"])
    assert year is None and stored < emis                                 # not paid back in 3 years
    assert "not paid back within 3 years" in emissions.payback_text(res["hist"])
    fake = dict(iso=[dict(year=1, stored=0.1, emissions=0.5), dict(year=2, stored=0.6, emissions=0.5)])
    assert emissions.payback_year(fake)[0] == 2 and "in year 2" in emissions.payback_text(fake)


# ---------------------------------------------------------------- Isometric settings
def test_nfz_is_the_deeper_of_20_cm_and_tillage_plus_margin(lab_copy):
    assert data_io.build_config()["isometric"]["nfz_depth_cm"] == 20          # no tillage on pasture
    p = lab_copy / "14_model_settings.csv"
    p.write_text(p.read_text().replace("Tillage depth,0,cm", "Tillage depth,20,cm"))
    cfg = data_io.build_config(lab_dir=lab_copy)
    assert cfg["isometric"]["nfz_depth_cm"] == 30                             # 20 cm tillage + 10 cm margin
    assert "tillage depth 20 cm + 10 cm" in cfg["provenance"]["isometric.nfz_depth_cm"]["source"]
    p.write_text(p.read_text().replace("Near-Field Zone margin below tillage depth,10,cm",
                                       "Near-Field Zone margin below tillage depth,5,cm"))
    cfg = data_io.build_config(lab_dir=lab_copy)
    assert cfg["isometric"]["nfz_depth_cm"] == 30                             # 25 cm rounded up to the 10 cm cells
    assert "rounded up" in cfg["provenance"]["isometric.nfz_depth_cm"]["source"]


def test_control_area_rule():
    site = pd.DataFrame([dict(item="Untreated control", value="2.5")])
    need = pd.DataFrame([dict(setting="Minimum control share of area", value="2.5")])
    assert validate.check_control_area(site, need)[0].level == "pass"
    small = pd.DataFrame([dict(item="Untreated control", value="1.0")])
    assert validate.check_control_area(small, need)[0].level == "fail"
    assert validate.check_control_area(pd.DataFrame([dict(item="x", value="1")]), need)[0].level == "info"
    assert default_config()["isometric"]["control_min_pct"] == 2.5


# ---------------------------------------------------------------- Monte Carlo and tornado
def test_ranges_come_from_the_lab_tables_and_are_assumed():
    params = {p["key"]: p for p in data_io.mc_parameters()}
    assert {"field_factor", "roughness", "glass_frac", "drainage_mm_yr", "log_pco2"} <= set(params)
    assert (params["field_factor"]["low"], params["field_factor"]["high"]) == (0.01, 1.0)
    assert params["field_factor"]["distribution"] == "log-uniform" and params["drainage_mm_yr"]["low"] == 500
    assert all(p["status"] == "assumed" for p in params.values())


def test_latin_hypercube_covers_each_range():
    params = data_io.mc_parameters()
    s = uncertainty.sample_inputs(params, 40, seed=3)
    ff = np.array([x["field_factor"] for x in s])
    assert ff.min() >= 0.01 and ff.max() <= 1.0
    assert np.sum(ff < 0.1) == 20                                         # log-uniform: half below the log midpoint
    dr = np.sort([x["drainage_mm_yr"] for x in s])
    assert all(500 + 22.5 * i <= d <= 500 + 22.5 * (i + 1) for i, d in enumerate(dr))   # one sample per slice
    assert s == uncertainty.sample_inputs(params, 40, seed=3)             # reproducible


def test_monte_carlo_percentiles_and_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(uncertainty, "CACHE", tmp_path)
    cfg = default_config()
    cfg["years"] = 1
    params = [p for p in data_io.mc_parameters() if p["key"] in ("field_factor", "roughness")]
    seen = []
    mc = uncertainty.monte_carlo(cfg, n=6, params=params, workers=3, progress=lambda a, b: seen.append((a, b)))
    assert mc["failed"] == 0 and not mc["from_cache"] and seen[-1] == (6, 6)
    p = mc["percentiles"]["released"]
    assert p["P10"][0] <= p["P30"][0] <= p["P40"][0] <= p["P50"][0] <= p["P90"][0]
    assert p["P90"][0] > 5 * p["P10"][0]                                  # the rate factor spans two orders of magnitude
    assert len(mc["final"]["net_removal"]) == 6
    again = uncertainty.monte_carlo(cfg, n=6, params=params, workers=3)
    assert again["from_cache"] and again["percentiles"] == mc["percentiles"]
    rows = uncertainty.percentile_table(mc)
    assert rows[0]["Term"] == "Cations released" and set(rows[0]) == {"Term", "P10", "P30", "P40", "P50", "P90"}
    assert "not creditable" in uncertainty.ISSUANCE_NOTE


def test_tornado_orders_inputs_by_their_effect(tmp_path, monkeypatch):
    monkeypatch.setattr(uncertainty, "CACHE", tmp_path)
    cfg = default_config()
    cfg["years"] = 1
    params = [p for p in data_io.mc_parameters() if p["key"] in ("field_factor", "log_pco2")]
    t = uncertainty.tornado(cfg, metric="released", params=params, workers=3)
    assert [b["key"] for b in t["bars"]][0] == "field_factor"            # the largest lever
    ff = t["bars"][0]
    assert ff["at_low"] < t["base"] < ff["at_high"] and ff["at_high"] / ff["at_low"] > 20


# ---------------------------------------------------------------- minimum detectable change
def test_minimum_detectable_change():
    mdc = detect.minimum_detectable_change
    assert mdc(30, 10) == pytest.approx((1.734 + 0.862) * 30 * (2 / 10) ** 0.5, rel=0.01)    # t(0.95, 18), t(0.80, 18)
    assert mdc(30, 40) < mdc(30, 10) < mdc(60, 10)
    assert mdc(30, 10, paired=True, correlation=0.8) < mdc(30, 10)
    assert mdc(30, 10, one_sided=False) > mdc(30, 10)
    assert mdc(30, 1) == float("inf")
    n = detect.composites_needed(30, 5)
    assert mdc(30, n) <= 5 < mdc(30, n - 1)
    assert detect.composites_needed(30, 0.001, n_max=200) is None


def test_modelled_signals(short_run):
    cfg, res = short_run
    ex = detect.exchangeable_signal(res["hist"])
    assert len(ex["relative_pct"]) == 3 and 0 < ex["relative_pct"][0] < ex["relative_pct"][-1]
    rock_only = detect.ratio_signal(res["hist"], cfg)
    assert "rock alone" in rock_only["basis"] and rock_only["relative_pct"][-1] < 0
    mixed = detect.ratio_signal(res["hist"], cfg, soil_camg_mg_kg=20000, soil_ti_mg_kg=8000)     # test fixture values
    assert abs(mixed["relative_pct"][-1]) < abs(rock_only["relative_pct"][-1])               # the soil dilutes the signal


# ---------------------------------------------------------------- tracer mass balance
ROCK = dict(Ca=71800.0, Mg=41700.0, Ti=23200.0)          # MB-BRITA2 from its oxides, mg/kg
SOIL = dict(Ca=9000.0, Mg=6000.0, Ti=7000.0)             # made-up baseline soil (test fixture)


def mixture(x, loss):
    """A post-application sample: fraction x of rock in the soil, a share `loss` of its Ca and Mg gone."""
    return dict(Ti=SOIL["Ti"] + x * (ROCK["Ti"] - SOIL["Ti"]),
                Ca=SOIL["Ca"] + x * (ROCK["Ca"] - SOIL["Ca"]) - loss * x * ROCK["Ca"],
                Mg=SOIL["Mg"] + x * (ROCK["Mg"] - SOIL["Mg"]) - loss * x * ROCK["Mg"])


def test_feedstock_elements_from_oxides():
    f = detect.feedstock_camg_ti(default_config()["oxides"])
    assert f["Ca"] == pytest.approx(ROCK["Ca"], rel=0.005) and f["Mg"] == pytest.approx(ROCK["Mg"], rel=0.005)
    assert f["Ti"] == pytest.approx(ROCK["Ti"], rel=0.005)


def test_mass_balance_recovers_a_known_loss():
    x = 40 / 1640                                         # 40 t/ha of rock in 1600 t/ha of soil
    r = detect.tracer_mass_balance(SOIL, mixture(x, 0.15), ROCK, soil_mass_t_ha=1600)
    assert r["feedstock_fraction"] == pytest.approx(x)
    assert r["loss"]["Ca"] == pytest.approx(0.15) and r["loss"]["Mg"] == pytest.approx(0.15)
    assert r["apparent_application_t_ha"] == pytest.approx(40.0)
    eq_per_kg = 0.15 * (2 * 71.8 / 40.078 + 2 * 41.7 / 24.305)
    assert r["co2e_t_ha"] == pytest.approx(eq_per_kg * 40e3 * 44.01 / 1e6, rel=1e-6)
    none = detect.tracer_mass_balance(SOIL, SOIL, ROCK)
    assert none["feedstock_fraction"] == 0 and "no feedstock detected" in none["note"]
    with pytest.raises(ValueError):
        detect.tracer_mass_balance(dict(SOIL, Ti=30000.0), SOIL, ROCK)


def test_bootstrap_gives_the_issuance_percentiles():
    rng = np.random.default_rng(0)
    x = 40 / 1640
    noisy = lambda base: {k: v * (1 + rng.normal(0, 0.02)) for k, v in base.items()}        # 2 % scatter between composites
    baseline = [noisy(SOIL) for _ in range(20)]
    post = [noisy(mixture(x, 0.3)) for _ in range(20)]
    b = detect.bootstrap_mass_balance(baseline, post, ROCK, soil_mass_t_ha=1600, n_boot=500)
    p = b["percentiles"]["co2e_t_ha"]
    assert p["P10"] <= p["P30"] <= p["P40"] <= p["P50"] <= p["P90"]
    assert b["n_baseline"] == 20 and b["n_valid"] > 400
    # the rock is only 2.4 % of the soil, so 2 % scatter between composites makes the estimate noisy:
    # the central value lands within about 0.2 of the true loss (0.3) and the P10-P90 band is wide
    lo, hi = b["percentiles"]["loss_Ca"]["P10"], b["percentiles"]["loss_Ca"]["P90"]
    assert abs(b["central"]["loss"]["Ca"] - 0.3) < 0.2 and hi - lo > 0.1
    assert lo < b["central"]["loss"]["Ca"] < hi


# ---------------------------------------------------------------- database comparison
def test_databases_agree_closely_on_the_carbonate_system(short_run):
    cfg, res = short_run
    out = dbcompare.compare_run(res, cfg)
    t = out["rock"]
    assert list(t.index[:2]) == ["llnl.dat", "phreeqc.dat"]
    diff = t.loc["difference (phreeqc.dat minus llnl.dat)"]
    assert abs(diff["pH at soil CO2"]) < 0.05 and abs(diff["DIC after degassing (mmol/kg)"]) < 0.01
    assert t.loc["llnl.dat", "pH after degassing to open air"] > t.loc["llnl.dat", "pH at soil CO2"]
    kept = out["carbon_kept_per_alkalinity"]
    assert 0.9 < kept["llnl.dat"] <= 1.0 and abs(kept["llnl.dat"] - kept["phreeqc.dat"]) < 0.02
