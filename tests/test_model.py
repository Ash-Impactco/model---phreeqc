"""The PHREEQC model: CO2 potential, the exchange pool, mineralogy sources, and regressions on S1.

Two regressions:
  - the original inputs (data/inputs.json) must still give the Step 0 numbers exactly
    (tests/reference/step0_results.csv), so the chemistry core is unchanged;
  - the lab-fed base case must match data/lab/15_model_results_10yr.csv, the table the batch
    rewrites (`python scripts/run_batch.py --update-lab`). If this fails after a deliberate model
    change, regenerate the table and record the change in docs/changelog.md.
"""
from pathlib import Path

import pandas as pd
import pytest

from conftest import append_rows
from erw import data_io
from erw.model import (co2_potential, default_config, exchange_pool, implied_oxides, mineral_assemblage,
                       normative_minerals, run, full_input)

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
