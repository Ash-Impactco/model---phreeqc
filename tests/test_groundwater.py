"""Groundwater step: inputs, transects, and the coupled soil -> aquifer -> stream run.

The coupled tests need the MODFLOW 6 binaries (python scripts/get_modflow.py) and are skipped without them.
"""
import json

import numpy as np
import pytest

from erw import groundwater as gw
from erw.model import default_config, run

needs_mf6 = pytest.mark.skipif(not gw.check_binaries()[0], reason="MODFLOW 6 binaries not in ./bin")


# ---------------------------------------------------------------- inputs
def test_hydrogeology_is_read_from_table_21_and_all_assumed():
    p, bg, prov = gw.load_hydrogeology()
    assert (p["length_m"], p["field_length_m"], p["thickness_m"], p["k_m_d"], p["porosity"]) == (1500, 600, 20, 10, 0.15)
    assert p["nlay"] == 4 and p["ncol"] == 60 and p["secondary_phases"] == ["Calcite", "SiO2(am)"]
    assert p["aquifer_weathering"] is False and p["vadose_delay_yr"] == 0
    assert bg["Alkalinity"] == pytest.approx(1.2e-3) and bg["Ca"] == pytest.approx(3e-4) and bg["pH"] == 7.0
    assert {v["status"] for v in prov.values()} <= {"assumed", "plan"}              # nothing measured yet
    assert not any(v["status"] == "measured" for v in prov.values())
    assert len(gw.provenance_table()) == len(prov)


def test_transect_from_csv_and_geojson(tmp_path):
    csv_file = tmp_path / "t.csv"
    csv_file.write_text("distance_m,elevation_m\n0,310\n700,220\n1400,95\n")       # test fixture
    pts = gw.load_transect(csv_file)
    assert pts == [(0.0, 310.0), (700.0, 220.0), (1400.0, 95.0)]
    assert gw.transect_params(pts) == dict(length_m=1400.0, top_upslope_m=310.0, top_stream_m=95.0, transect=pts)
    gj = tmp_path / "t.geojson"
    gj.write_text(json.dumps({"type": "Feature", "geometry": {"type": "LineString", "coordinates": [
        [-25.3534, 37.8296, 300.0], [-25.3534, 37.8386, 120.0]]}}))                 # 0.009 degrees of latitude north
    pts = gw.load_transect(gj)
    assert pts[0] == (0.0, 300.0) and pts[1][1] == 120.0
    assert pts[1][0] == pytest.approx(1000.0, rel=0.01)                             # about 1 km
    flat = tmp_path / "flat.geojson"
    flat.write_text(json.dumps({"type": "LineString", "coordinates": [[-25.35, 37.83], [-25.35, 37.84]]}))
    with pytest.raises(ValueError):
        gw.load_transect(flat)


def test_delay_shifts_the_yearly_series():
    assert gw._delayed([1, 2, 3, 4], 0, 0) == [1, 2, 3, 4]
    assert gw._delayed([1, 2, 3, 4], 0, 2) == [0, 0, 1, 2]


# ---------------------------------------------------------------- coupled run: S7, regression and conservation
@pytest.fixture(scope="module")
def s7():
    cfg = default_config("SG-04", "milled_p80_100um")
    cfg.update(acidity_mode="effective", glass_frac=0.2, field_factor=0.3)
    return cfg, run(cfg)


@pytest.fixture(scope="module")
def s7_coupled(s7):
    cfg, soil = s7
    return gw.run_coupled(cfg, soil_result=soil, label="test_s7")


@needs_mf6
def test_s7_share_reaching_the_stream(s7_coupled):
    r = s7_coupled
    assert r["ok"] and r["resolvable"]
    assert 0.90 <= r["reach_fraction"] <= 1.0001               # nothing is lost on the way in this scenario
    s = r["balance_share"]
    assert s["calcite"] == pytest.approx(0.0, abs=1e-6)
    assert 0.1 < s["stream"] < 0.5 and 0.5 < s["in_transit"] < 0.9       # most is still in the aquifer after 10 years
    assert 4.0 < r["first_arrival_yr"] < 9.0
    assert r["d_alk_ueq"][-1] > gw.RESOLVABLE_UEQ and all(np.diff(r["d_alk_ueq"][-24:]) > 0)   # still rising at the end


@needs_mf6
def test_alkalinity_is_conserved(s7_coupled, s7):
    r = s7_coupled
    b = r["balance_t_ha"]
    assert abs(b["residual"]) < 0.01 * b["input"]                                   # in = stream + in transit + calcite
    assert b["input"] == pytest.approx(b["stream"] + b["in_transit"] + b["calcite"] + b["residual"], rel=1e-9)
    # what enters the aquifer is what the soil model sent out at 50 cm
    assert b["input"] == pytest.approx(s7[1]["hist"]["export_t"][-1], rel=0.01)
    assert r["cum_stream_t_ha"][-1] == pytest.approx(b["stream"], rel=1e-9)


@needs_mf6
def test_stream_degassing_keeps_most_of_the_carbon(s7_coupled):
    d = s7_coupled["degassing"]
    assert d["resolvable"] and 0.9 < d["carbon_kept_per_alk"] <= 1.0
    assert d["calcite_loss_per_alk"] == pytest.approx(0.0, abs=1e-6)
    assert d["rock"]["pH_out"] > d["rock"]["pH_in"]                                 # CO2 leaves, pH rises
    assert d["rock"]["alk_out"] == pytest.approx(d["rock"]["alk_in"], rel=1e-6)     # degassing does not change alkalinity


@needs_mf6
def test_frames_for_the_cross_section(s7_coupled):
    f, p = s7_coupled["frames"], s7_coupled["params"]
    assert len(f["alk_rock"]) == p["years"] and f["alk_rock"][0].shape == (p["nlay"], p["ncol"])
    extra = f["alk_rock"][-1] - f["alk_ctrl"][-1]
    assert extra.max() > 0.05 and extra[:, :5].max() == pytest.approx(extra.max(), rel=0.5)   # plume starts under the field
    assert len(s7_coupled["x_m"]) == p["ncol"] and len(s7_coupled["z_m"]) == p["nlay"]


@needs_mf6
def test_small_signals_are_reported_as_not_resolvable():
    cfg = default_config()                                     # S1: about 1 ueq/L of extra alkalinity in the drainage
    cfg["years"] = 2
    r = gw.run_coupled(cfg, dict(ncol=20, nlay=1, steps_per_year=12), label="test_s1")
    assert not r["resolvable"] and r["balance_share"]["stream"] is None
    assert r["first_arrival_yr"] is None and not r["degassing"]["resolvable"]


@needs_mf6
def test_aquifer_weathering_keeps_the_background_alkalinity(s7):
    cfg, soil = s7
    small = dict(ncol=20, nlay=2, steps_per_year=12, years=5)
    off = gw.run_coupled(cfg, dict(small, aquifer_weathering=False), soil_result=soil, label="test_w_off")
    on = gw.run_coupled(cfg, dict(small, aquifer_weathering=True), soil_result=soil, label="test_w_on")
    assert off["alk_stream_ctrl_meq"][-1] < 0.8                # flushed towards the dilute soil drainage
    assert on["alk_stream_ctrl_meq"][-1] == pytest.approx(1.2, rel=0.2)             # stays near the background
    assert on["degassing"]["rock"]["pH_in"] < 8.0              # and at a pH a basalt aquifer would show
    assert abs(on["balance_t_ha"]["residual"]) < 0.02 * on["balance_t_ha"]["input"]
