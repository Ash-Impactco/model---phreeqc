"""The PHREEQC model: CO2 potential, the exchange pool, mineralogy sources, and regressions on S1.

Two regressions:
  - the original inputs (data/inputs.json) must still give the Step 0 numbers exactly
    (tests/reference/step0_results.csv), so the chemistry core is unchanged;
  - the lab-fed base case must match data/lab/15_model_results_10yr.csv, the table the batch
    rewrites (`python scripts/run_batch.py --update-lab`). If this fails after a deliberate model
    change, regenerate the table and record the change in docs/changelog.md.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from conftest import append_rows
from erw import data_io
from erw.model import (andosol_pool, build_setup, co2_potential, default_config, exchange_pool, implied_oxides,
                       mineral_assemblage, normative_minerals, run, full_input, shift_times, _geometry)

ROOT = Path(__file__).resolve().parents[1]
COLUMNS = {"rock_dissolved_10yr_pct": lambda h: h["frac_diss"][-1] * 100,
           "released_tCO2e_ha": lambda h: h["released_t"][-1],
           "co2_removed_t_ha": lambda h: h["export_t"][-1],
           "on_exchange_tCO2e_ha": lambda h: h["exchange_t"][-1],
           "iso_exported_from_NFZ_20cm_t_ha": lambda h: h["iso"][-1]["exported_from_nfz"],
           "iso_net_sorption_t_ha": lambda h: h["iso"][-1]["net_sorption"],
           "iso_stored_t_ha": lambda h: h["iso"][-1]["stored"],
           "iso_emissions_t_ha": lambda h: h["iso"][-1]["emissions"],
           "iso_net_removal_t_ha": lambda h: h["iso"][-1]["net_removal"]}


def s1_row(path):
    df = pd.read_csv(path)
    return df[df["scenario"].str.startswith("S1")].iloc[0]


def assert_matches(hist, ref, places):
    for col, get in COLUMNS.items():
        tol = 0.5 * 10 ** -(2 if col == "rock_dissolved_10yr_pct" else places) + 1e-9
        assert abs(get(hist) - ref[col]) <= tol, f"{col}: model {get(hist):.5f}, reference {ref[col]}"


# ---------------------------------------------------------------- quick checks (no PHREEQC run)
def test_co2_potential_mb_brita2():
    cons, full = co2_potential(default_config()["oxides"])
    assert round(cons, 3) == 0.154
    assert round(full, 3) == 0.375


def test_norm_conserves_base_cations():
    """The norm puts all Ca, Mg, Na and K into dissolving minerals: the minerals carry the bulk
    oxides in the same proportions (one common factor, from scaling the assemblage to 100 g)."""
    cfg = default_config()
    implied = implied_oxides(cfg)
    ratios = [implied[ox] / cfg["oxides"][ox] for ox in ("CaO", "MgO", "Na2O", "K2O")]
    assert max(ratios) - min(ratios) < 2e-3
    assert 0.95 < ratios[0] < 1.05


def test_mineralogy_source_is_reported():
    cfg = default_config()
    g, _, _, source = mineral_assemblage(cfg)
    assert "normative" in source and "assumed" in source and g["Gls"] == 0
    cfg["glass_frac"] = 0.2
    g, _, _, source = mineral_assemblage(cfg)
    assert "20% glass" in source and g["Gls"] > 0
    assert sum(g.values()) == pytest.approx(100.0)


def test_exchange_pool_options():
    cfg = default_config("SG-04")
    cec, beta = exchange_pool(cfg)
    assert cec == 23.61 and beta["Na"] == 0.01 and sum(beta.values()) == pytest.approx(1.0)
    cfg["acidity_mode"] = "effective"
    cec, beta = exchange_pool(cfg)
    assert cec == pytest.approx(3.54 + 1.74 + 1.2 + 1.0)                  # assumed 1.0 cmol/kg KCl acidity
    full = default_config("SG-01")                                         # report with Na and acidity at pH 7
    cec, beta = exchange_pool(full)
    assert cec == pytest.approx(5.76 + 1.59 + 2.25 + 0.17 + 8.23) and beta["Na"] == pytest.approx(0.17 / cec)
    assert beta["H"] == pytest.approx(8.23 / cec)


def test_exchange_pool_refuses_impossible_soils():
    cfg = default_config()
    cfg["soil"]["ca"] = 40.0                                               # more bases than CEC
    with pytest.raises(ValueError):
        exchange_pool(cfg)


def test_pqi_export_names_the_database_and_soil():
    text = full_input(default_config(), rock=True)
    assert "llnl.dat" in text and "SG-04" in text and "KINETICS" in text and "TRANSPORT" in text
    assert "KINETICS" not in full_input(default_config(), rock=False)


# ---------------------------------------------------------------- regressions on S1 (two 10-year runs)
@pytest.fixture(scope="module")
def s1_original():
    return run(default_config(source="json"))["hist"]


@pytest.fixture(scope="module")
def s1_lab():
    return run(default_config())["hist"]


def test_original_inputs_still_give_the_step0_numbers(s1_original):
    assert_matches(s1_original, s1_row(ROOT / "tests" / "reference" / "step0_results.csv"), places=4)
    assert round(s1_original["released_t"][-1], 3) == 0.096
    assert round(s1_original["export_t"][-1], 4) == 0.0027
    assert round(s1_original["iso"][-1]["exported_from_nfz"], 4) == 0.0081
    assert round(s1_original["iso"][-1]["net_removal"], 2) == -0.16


def test_lab_fed_base_case_matches_the_results_table(s1_lab):
    assert_matches(s1_lab, s1_row(ROOT / "data" / "lab" / "15_model_results_10yr.csv"), places=4)


def test_cation_budget_closes(s1_lab):
    """Released cations = held on exchange + leaving the Near-Field Zone + calcite, within a few per cent.
    (The flux is summed over the water that leaves at each shift; the year-end state is not counted twice.)"""
    i = s1_lab["iso"][-1]
    assert abs(i["non_carbonic_neut"]) < 0.03 * i["weathering"]
    assert i["weathering"] == pytest.approx(i["exported_from_nfz"] + i["net_sorption"] + i["net_new_carbonate"]
                                            + i["non_carbonic_neut"], rel=1e-9)


def test_lab_fed_base_case_is_physically_sensible(s1_lab, s1_original):
    h = s1_lab
    assert 0 < h["export_t"][-1] < h["iso"][-1]["exported_from_nfz"] < h["released_t"][-1]
    assert all(b >= a for a, b in zip(h["released_t"], h["released_t"][1:]))            # release only grows
    # finer representative diameters from the interpolated curve: a little more surface, same order of magnitude
    assert 1.0 < h["released_t"][-1] / s1_original["released_t"][-1] < 1.2
    assert h["iso"][-1]["net_removal"] < 0                                              # still net negative


# ---------------------------------------------------------------- measured mineralogy drives the model
def test_measured_xrd_replaces_the_norm_and_runs(lab_copy):
    # test fixture: oxides and XRD for a made-up TECHNOVIA sample (copied from MB-BRITA2 chemistry)
    append_rows(lab_copy / "01_rock_major_oxides.csv", [
        "TV-TEST,test,fixture,measured,test,47.44,14.05,13.02,,10.05,6.92,3.37,2.0,3.87,,,,,,,,"])
    append_rows(lab_copy / "05_rock_mineralogy.csv", [
        "TV-TEST,Plagioclase (An55),40,measured,XRD Rietveld,test fixture",
        "TV-TEST,Clinopyroxene,22,measured,XRD Rietveld,test fixture",
        "TV-TEST,Olivine,8,measured,XRD Rietveld,test fixture",
        "TV-TEST,K-feldspar,5,measured,XRD Rietveld,test fixture",
        "TV-TEST,Ilmenite,7,measured,XRD Rietveld,test fixture",
        "TV-TEST,Basaltic glass (amorphous),18,measured,XRD Rietveld,test fixture",
    ])
    cfg = data_io.build_config(rock="TV-TEST", lab_dir=lab_copy)
    assert cfg["provenance"]["mineralogy"]["status"] == "measured"
    assert cfg["provenance"]["oxides"]["status"] == "measured"
    assert "stand-in" in cfg["provenance"]["psd"]["source"]               # no grading curve for it: borrowed, and said so
    g, minerals, an, source = mineral_assemblage(cfg)
    assert "measured phases (XRD)" in source and "olivine composition not reported" in source
    assert g["Gls"] == pytest.approx(18.0) and g["Plag"] == pytest.approx(40.0) and an == pytest.approx(0.55)
    assert g["Fo"] + g["Fa"] == pytest.approx(8.0)
    assert minerals["Plag"][0].startswith("Na0.450Ca0.550")
    cfg["years"] = 1
    h = run(cfg)["hist"]
    norm = data_io.build_config(lab_dir=lab_copy)
    norm["years"] = 1
    assert h["released_t"][-1] > run(norm)["hist"]["released_t"][-1] > 0  # 18 % glass dissolves faster than none


# ---------------------------------------------------------------- Andosol acidity: exchangeable Al + variable charge
@pytest.mark.parametrize("soil", ["SG-04", "SG-05", "SG-01", "SG-02"])
@pytest.mark.parametrize("phase", ["amorphous", "gibbsite"])
def test_andosol_pool_reproduces_the_two_measured_charges(soil, phase):
    cfg = default_config(soil)
    cfg.update(acidity_mode="andosol", al_phase=phase)
    a = andosol_pool(cfg)
    s = cfg["soil"]
    bases = s["ca"] + s["mg"] + s["k"] + s.get("na", 0.0)
    charge_soil_ph = a["x_cmol"] + a["frac_charged_at_soil_ph"] * a["y_cmol"]
    charge_ph7 = a["x_cmol"] + a["frac_charged_at_ph7"] * a["y_cmol"]
    assert charge_soil_ph == pytest.approx(a["ecec"], rel=1e-6)            # effective CEC at the soil pH
    assert charge_ph7 == pytest.approx(a["cec7"], rel=1e-6)                # CEC at pH 7
    assert a["ecec"] == pytest.approx(bases + cfg["effective_acidity_cmol_kg"])
    assert a["x_cmol"] - a["al_exch_cmol"] + a["bases_on_y"] == pytest.approx(
        bases if "na" in s else bases + 0.01 * a["ecec"], rel=1e-6)        # every measured cation is placed once
    assert s["pH"] < a["half_ph"] < 7.0
    assert sum(a["beta_x"].values()) == pytest.approx(1.0)


def test_amorphous_control_is_closer_to_the_database_al_constant_than_gibbsite():
    for soil in ("SG-04", "SG-05", "SG-01", "SG-02"):
        fit = {}
        for phase in ("amorphous", "gibbsite"):
            cfg = default_config(soil)
            cfg.update(acidity_mode="andosol", al_phase=phase)
            fit[phase] = andosol_pool(cfg)["logk_alx3"]
        assert abs(fit["amorphous"] - 0.41) < abs(fit["gibbsite"] - 0.41), soil
        assert fit["gibbsite"] - fit["amorphous"] == pytest.approx(3.1, abs=0.3)   # the gap between the two phases


def test_andosol_pool_refuses_zero_acidity():
    cfg = default_config()
    cfg.update(acidity_mode="andosol", effective_acidity_cmol_kg=0.0)
    with pytest.raises(ValueError):
        andosol_pool(cfg)


def test_andosol_run_starts_at_the_measured_exchange_state_and_buffers():
    cfg = default_config()
    cfg.update(acidity_mode="andosol", years=1)
    res = run(cfg)
    a, p = res["acidity_pool"], res["hist"]["profiles"][0]
    # deepest cell is still close to the starting state: exchangeable Al near the KCl acidity, base
    # saturation near bases / CEC at pH 7
    assert p["exal_ctrl"][-1] == pytest.approx(a["al_exch_cmol"], rel=0.35)
    assert p["bs_ctrl"][-1] == pytest.approx((3.54 + 1.74 + 1.2) / 23.61 * 100, rel=0.08)
    assert res["hist"]["released_t"][-1] > 0
    text = full_input(cfg, rock=True)
    assert "Al(OH)3(a)" in text and "HXv" in text and "AlX3" in text and "Al13O4(OH)24+7" in text
    cfg["al_phase"] = "gibbsite"
    assert "Al13O4(OH)24+7" not in full_input(cfg, rock=True)             # llnl.dat left untouched for gibbsite


# ---------------------------------------------------------------- glass: Gislason & Oelkers (2003)
def glass_cfg(**kw):
    cfg = default_config()
    cfg.update(glass_frac=0.5, years=1)
    cfg.update(kw)
    return cfg


def test_glass_rate_law_constants_and_area():
    cfg = glass_cfg()
    assert cfg["glass_rate_law"] == "gislason_oelkers_2003"
    assert cfg["glass_go"] == dict(log_a=-1.6, ea=25.5)                   # 10^-5.6 mol Si/cm2/s = 10^-1.6 /m2/s
    assert cfg["provenance"]["glass_go"]["status"] == "assumed"
    text, _ = build_setup(cfg, rock=True)
    block = text[text.index("Gls_0"):text.index("-end", text.index("Gls_0"))]
    assert "10^(-1.6) * exp(-25.5/(R*TK))" in block
    assert '(act("H+")^3/aal)^(1/3)' in block and 'LK_PHASE("Gibbsite")' in block
    rougher = glass_cfg(roughness=20.0)
    text2, _ = build_setup(rougher, rock=True)
    block2 = text2[text2.index("Gls_0"):text2.index("-end", text2.index("Gls_0"))]
    assert block == block2                                                # per geometric area: no roughness factor
    old, _ = build_setup(glass_cfg(glass_rate_law="approximate"), rock=True)
    assert "aal" not in old


def test_glass_law_matches_the_old_rate_with_gibbsite_and_slows_with_amorphous_al():
    released = lambda cfg: run(cfg)["hist"]["released_t"][-1]
    full = released(glass_cfg())                                          # Al held by gibbsite (large pool)
    approx = released(glass_cfg(glass_rate_law="approximate"))
    assert full == pytest.approx(approx, rel=0.15)
    # amorphous Al(OH)3 holds about 1000 times more Al3+ in solution: the law gives a rate about 10 times lower
    slow = released(glass_cfg(acidity_mode="andosol", al_phase="amorphous"))
    fast = released(glass_cfg(acidity_mode="andosol", al_phase="gibbsite"))
    assert slow < 0.5 * fast


# ---------------------------------------------------------------- seasonal drainage
def test_shift_times_follow_the_monthly_shares():
    cfg = default_config()
    per_year = int(round(_geometry(cfg)[-1]))
    t = np.array([0.0] + shift_times(cfg, 3 * per_year)) / 86400.0
    assert t[per_year] == pytest.approx(365.25) and t[3 * per_year] == pytest.approx(3 * 365.25)
    steps = np.diff(t)[:per_year]
    assert steps.max() / steps.min() == pytest.approx(3.0, rel=0.15)      # dry months drain a third as fast
    # the run starts in October: the first six months take 75 % of the year's shifts
    first_half = np.searchsorted(t[1:per_year + 1], 31 + 30 + 31 + 31 + 28.25 + 31 + 1.0)      # within a day of 1 April
    assert first_half == round(0.75 * per_year)
    # equal shares give steady drainage
    cfg["drainage_monthly_share"] = [1 / 12] * 12
    even = np.diff([0.0] + shift_times(cfg, per_year))
    assert even.max() / even.min() < 1.15                                 # only the unequal month lengths remain


def test_seasonal_run_keeps_the_yearly_totals():
    steady = default_config(); steady["years"] = 2
    seasonal = default_config(); seasonal.update(years=2, seasonal_drainage=True)
    a, b = run(steady), run(seasonal)
    assert b["effluent"]["t"][-1] < 2.0 and len(a["effluent"]["t"]) == len(b["effluent"]["t"])   # same number of shifts
    assert b["hist"]["released_t"][-1] == pytest.approx(a["hist"]["released_t"][-1], rel=0.05)
    assert full_input(seasonal).count("TRANSPORT") > 4 and full_input(steady).count("TRANSPORT") == 1


# ---------------------------------------------------------------- plant uptake
def test_uptake_is_off_by_default_and_leaves_the_input_unchanged():
    text, _ = build_setup(default_config(), rock=True)
    assert "Uptake_" not in text
    ctrl, _ = build_setup(default_config(), rock=False)
    assert "KINETICS" not in ctrl


def test_uptake_removes_close_to_the_stated_amount_in_both_columns():
    cfg = default_config()
    cfg.update(years=2, uptake_share_removed=0.1)                 # 2.4 Ca, 1.2 Mg, 18.3 K kg/ha/yr
    text, _ = build_setup(cfg, rock=False)
    assert "Uptake_K" in text and "-formula K 1 H -1" in text and "-formula Ca 1 H -2" in text
    i = run(cfg)["hist"]["iso"][-1]
    stated = (2 * 2.4 / 40.078 + 2 * 1.2 / 24.305 + 18.3 / 39.098) * 0.1 * 44.01 / 100 * 2      # t CO2e/ha in 2 years
    assert i["uptake_ctrl"] == pytest.approx(stated, rel=0.15)
    assert i["uptake_rock"] == pytest.approx(i["uptake_ctrl"], rel=0.01)   # demand-driven: nearly the same with rock
    assert i["biomass_loss"] == pytest.approx(i["uptake_rock"] - i["uptake_ctrl"], abs=1e-9)
    assert abs(i["non_carbonic_neut"]) < 0.05 * i["weathering"]   # the budget still closes with the sink in it
