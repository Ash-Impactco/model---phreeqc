"""Trace metals: blank until table 06 is filled; made-up contents (test fixtures) show the arithmetic."""
import pytest

from conftest import append_rows
from erw import data_io, metals
from erw.model import run


def test_no_numbers_without_an_analysis():
    cfg = data_io.build_config()
    table, note = metals.accumulation(cfg)
    assert table.empty and "No trace-metal analysis for MB-BRITA2" in note and "Nothing has been assumed" in note
    assert metals.rock_metals("TV-FINES-01") == {}                      # rows exist but are blank
    assert all(v["limit"] is None and v["baseline"] is None for v in metals.limits().values())
    assert set(metals.HOSTED_IN_OLIVINE_PYROXENE) <= set(metals.limits())


def test_accumulation_with_fixture_values(lab_copy):
    append_rows(lab_copy / "06_rock_trace_metals.csv", [
        "MB-BRITA2,Ni,100,1,measured,test fixture,",
        "MB-BRITA2,Cr,250,1,measured,test fixture,",
        "MB-BRITA2,Cd,<0.2,0.2,measured,test fixture,"])
    p = lab_copy / "20_trace_metal_limits.csv"
    p.write_text(p.read_text().replace("Ni,,,TO PROVIDE", "Ni,50,20,TO PROVIDE"))      # limit 50, baseline 20: fixture
    cfg = data_io.build_config(lab_dir=lab_copy)
    cfg["years"] = 2
    hist = run(cfg)["hist"]
    table, note = metals.accumulation(cfg, hist, lab_copy)
    t = table.set_index("element")
    soil_t, rock_t = metals.topsoil_and_rock_t_ha(cfg)
    assert (soil_t, rock_t) == (1600.0, 40.0)                           # 20 cm at 800 kg/m3; 40 t/ha
    assert t.loc["Cr", "added_total_mg_kg"] == pytest.approx(40 * 250 / 1640)           # no baseline: plain dilution
    assert t.loc["Ni", "after_mg_kg"] == pytest.approx((1600 * 20 + 40 * 100) / 1640)   # mixing of soil and rock
    assert t.loc["Ni", "added_total_mg_kg"] == pytest.approx(t.loc["Ni", "after_mg_kg"] - 20)
    assert t.loc["Ni", "pct_of_limit"] == pytest.approx(t.loc["Ni", "after_mg_kg"] / 50 * 100)
    assert t.loc["Cr", "pct_of_limit"] is None or t.loc["Cr", "pct_of_limit"] != t.loc["Cr", "pct_of_limit"]
    assert t.loc["Ni", "status"] == "stand-in" and "detection limit" in t.loc["Cd", "status"]
    # released Ni follows the share of olivine + pyroxene dissolved
    host = hist["host_frac_diss"]
    assert 0 < host[0] < host[1] < 0.01
    assert t.loc["Ni", "released_mg_kg"] == pytest.approx(40 * 100 * host[-1] / 1640)
    assert metals.release_series(cfg, hist, "Ni", lab_copy)[-1] == pytest.approx(t.loc["Ni", "released_mg_kg"])
    assert metals.release_series(cfg, hist, "Cd", lab_copy) == []       # not tied to olivine and pyroxene
    assert "soil baseline" in note and "soil limits" in note            # still blank for Cr and Cd


def test_mass_by_mineral_and_bin_adds_up():
    cfg = data_io.build_config()
    cfg.update(years=1, glass_frac=0.2)
    h = run(cfg)["hist"]
    total0 = sum(h["mineral_init_g_m2"].values())
    assert total0 == pytest.approx(sum(h["bin_init_g_m2"]))
    assert total0 == pytest.approx(4000 * (1 - 0.106), rel=0.01)        # 40 t/ha less the inert oxides
    left_m = sum(v[-1] for v in h["mineral_left_g_m2"].values())
    left_b = sum(v[-1] for v in h["bin_left_g_m2"])
    assert left_m == pytest.approx(left_b, rel=1e-6) and left_m < total0
    finest, coarsest = h["bin_left_g_m2"][0][-1] / h["bin_init_g_m2"][0], h["bin_left_g_m2"][-1][-1] / h["bin_init_g_m2"][-1]
    assert finest < coarsest                                            # fine grains dissolve faster
