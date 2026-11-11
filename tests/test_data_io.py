"""The lab-data loader: cell parsing, status classes, tables, and the model configuration built from them.

Numbers written into a temporary copy of the data pack in these tests are made-up test fixtures.
They never touch data/lab.
"""
import math

import pytest

from conftest import append_rows
from erw import data_io
from erw.data_io import parse_value, status_class, effective_status


# ---------------------------------------------------------------- cell parsing
@pytest.mark.parametrize("cell, value, low, high, qualifier", [
    ("47.44", 47.44, 47.44, 47.44, ""),
    ("<0.05", 0.05, 0.0, 0.05, "<"),
    ("< 0.01", 0.01, 0.0, 0.01, "<"),
    ("<0.00", 0.0, 0.0, 0.0, "<"),
    ("2850-3000", 2925.0, 2850.0, 3000.0, "range"),
    ("-2.5 to -1.5", -2.0, -2.5, -1.5, "range"),
    ("500–1400", 950.0, 500.0, 1400.0, "range"),
    ("~11", 11.0, 11.0, 11.0, "~"),
    ("about 300", 300.0, 300.0, 300.0, "~"),
    ("-2.0", -2.0, -2.0, -2.0, ""),
    ("5,8", 5.8, 5.8, 5.8, ""),               # decimal comma, as in Portuguese lab reports
    ("1e-4", 1e-4, 1e-4, 1e-4, ""),
])
def test_parse_numbers(cell, value, low, high, qualifier):
    p = parse_value(cell)
    assert p.is_number and p.qualifier == qualifier
    assert p.value == pytest.approx(value) and p.low == pytest.approx(low) and p.high == pytest.approx(high)
    assert p.text == cell


def test_parse_number_with_note():
    p = parse_value("0.09 (reported as Mn2O4, likely Mn3O4)")
    assert p.value == 0.09 and "Mn2O4" in p.note
    assert parse_value("53 (as P)").note == "(as P)"


@pytest.mark.parametrize("cell", ["", "   ", None, float("nan")])
def test_parse_blank(cell):
    p = parse_value(cell)
    assert not p.is_number and p.value is None


@pytest.mark.parametrize("cell", ["non-reactive", "0/4", "2026-11", "required", "Sept to early Oct 2027"])
def test_parse_text_is_not_a_number(cell):
    assert not parse_value(cell).is_number


def test_below_detection_is_flagged_not_treated_as_measured():
    p = parse_value("<0.05")
    assert p.below_detection and p.low == 0.0


def test_latlon():
    assert data_io.parse_latlon("37.8296, -25.3534") == (37.8296, -25.3534)
    assert data_io.parse_latlon("approximate") is None


# ---------------------------------------------------------------- status
@pytest.mark.parametrize("status, cls", [
    ("measured", "measured"), ("measured, stand-in", "stand-in"), ("measured, reference", "measured"),
    ("measured, stand-in for Lomba da Maia", "stand-in"), ("assumed", "assumed"), ("TO MEASURE", "to measure"),
    ("TO FIND", "to measure"), ("OPTIONAL", "to measure"), ("unknown", "to measure"),
    ("assumed (combustion only)", "assumed"), ("model choice", "assumed"), ("calculated (OSRM)", "calculated"),
    ("plan", "plan"), ("model forecast (stand-in data, placeholder emission factors)", "model output"),
    ("", "unspecified"),
])
def test_status_class(status, cls):
    assert status_class(status) == cls


def test_measurement_on_another_rock_or_site_is_a_stand_in():
    assert effective_status("MB-BLOCK-OL", "measured") == "stand-in"
    assert effective_status("SG-01", "measured") == "stand-in"
    assert effective_status("TV-FINES-01", "measured") == "measured"
    assert effective_status("LDM-P01", "measured") == "measured"
    assert effective_status("LDM-P01", "TO MEASURE") == "to measure"
    assert effective_status("MB-BRITA2", "assumed") == "assumed"


# ---------------------------------------------------------------- tables
def test_every_table_loads_with_status():
    tables = data_io.load_all()
    assert len(tables) == 20
    for stem, df in tables.items():
        assert len(df) > 0, stem
        assert df.attrs["ragged_rows"] == [], f"{stem}: rows do not match the header"
        if stem.startswith("16"):
            continue                                     # references: no status
        assert "status" in df.columns and "status_class" in df.columns, stem


def test_read_table_by_number_name_or_stem():
    a = data_io.read_table("03")
    b = data_io.read_table("rock_grain_size")
    c = data_io.read_table("03_rock_grain_size")
    assert list(a.columns) == list(b.columns) == list(c.columns) and len(a) == len(b) == len(c)
    with pytest.raises(KeyError):
        data_io.read_table("99")


def test_blank_cells_stay_blank():
    df = data_io.read_table("09")
    assert (df["pH_H2O"] == "").all()
    assert df.loc[df["sample_id"] == "LDM-P01", "n_subsamples"].item() == "15"
    assert df.loc[df["sample_id"] == "LDM-P05-DUP", "notes"].item() == "QC duplicate"


def test_misaligned_rows_are_reported(lab_copy):
    append_rows(lab_copy / "03_rock_grain_size.csv", ["TEST-X,10,5,assumed,test fixture,note,one cell too many"])
    df = data_io.read_table("03", lab_copy)
    assert len(df.attrs["ragged_rows"]) == 1 and df.attrs["ragged_rows"][0][1] == 7


# ---------------------------------------------------------------- rock
def test_rock_samples():
    rocks = data_io.rock_samples()
    assert rocks["MB-BRITA2"]["usable"] and rocks["MB-BRITA2"]["status_class"] == "stand-in"
    assert rocks["MB-BRITA2"]["oxides"]["CaO"] == 10.05 and rocks["MB-BRITA2"]["fe2_fraction"] is None
    assert not rocks["TV-FINES-01"]["usable"] and rocks["TV-FINES-01"]["status_class"] == "to measure"
    assert rocks["MB-BRITA2"]["extras"]["LOI"].below_detection


def test_feo_sets_the_ferrous_fraction(lab_copy):
    # test fixture: 9.00 % FeO and 3.00 % Fe2O3
    append_rows(lab_copy / "01_rock_major_oxides.csv",
                ["TEST-FEO,test,fixture,measured,test,47,14,3.00,9.00,10,7,3,2,3,,,,,,,,"])
    r = data_io.rock_samples(lab_copy)["TEST-FEO"]
    fe2, fe3 = 9.00 / 71.844, 2 * 3.00 / 159.69
    assert r["fe2_fraction"] == pytest.approx(fe2 / (fe2 + fe3))
    assert r["oxides"]["Fe2O3"] == pytest.approx((fe2 + fe3) / 2 * 159.69)      # total iron as Fe2O3


def test_density_comes_from_table_02():
    v, prov = data_io.rock_density("MB-BRITA2")
    assert v == 2900.0 and prov["status"] == "assumed" and "02_rock_physical_properties" in prov["source"]
    v, prov = data_io.rock_density("MB-BLOCK-OL")
    assert v == 2925.0 and prov["status"] == "stand-in" and "midpoint" in prov["source"]
    v, prov = data_io.rock_density("TV-FINES-01")                # no particle density yet
    assert "fallback" in prov["source"] and prov["status"] == "assumed"


# ---------------------------------------------------------------- grain size
def specific_surface(bins):
    return sum(f / d for f, d in bins)


def test_bins_from_the_feedstock_curve():
    pts = data_io.psd_points("MB-BRITA2")
    assert [p["size_um"] for p in pts] == [63, 250, 1000, 2000, 4000]
    assert [p["status_class"] for p in pts] == ["stand-in", "assumed", "assumed", "assumed", "stand-in"]
    bins, detail = data_io.build_bins(pts)
    assert 5 <= len(bins) <= 8
    assert sum(f for f, _ in bins) == pytest.approx(1.0)
    assert bins[0] == [pytest.approx(0.11), pytest.approx(30e-6)]            # fines: measured amount, assumed size
    assert "assumed diameter" in detail[0]["basis"]
    assert [d["basis"] for d in detail[1:]] == ["assumed split"] * 4          # assumed points stay flagged
    for (f, d), det in zip(bins[1:], detail[1:]):                             # diameter lies inside its bin
        assert det["from_um"] < d * 1e6 < det["to_um"]


def test_surface_diameter_formula():
    # mass spread evenly in log(size) between d1 and d2: <1/d> = (1/d1 - 1/d2) / ln(d2/d1)
    d = data_io._surface_diameter(63.0, 250.0)
    assert 1 / d == pytest.approx((1 / 63 - 1 / 250) / math.log(250 / 63))


def test_few_points_are_split_and_many_are_merged_without_changing_surface_area():
    two = [dict(size_um=63, cum_pct=50, status_class="assumed"), dict(size_um=136, cum_pct=100, status_class="assumed")]
    coarse, _ = data_io.build_bins(two, n_min=2)
    fine, _ = data_io.build_bins(two, n_min=5)
    assert len(coarse) == 2 and len(fine) == 5
    assert specific_surface(fine) == pytest.approx(specific_surface(coarse), rel=1e-9)
    many = [dict(size_um=2.0 * 1.5 ** i, cum_pct=min(100.0, 5.0 * (i + 1)), status_class="measured") for i in range(20)]
    bins, _ = data_io.build_bins(many)
    assert len(bins) == 8 and sum(f for f, _ in bins) == pytest.approx(1.0)


def test_bad_curves_are_refused():
    with pytest.raises(ValueError):
        data_io.build_bins([dict(size_um=63, cum_pct=11, status_class="measured")])
    with pytest.raises(ValueError):                              # does not reach 100 %
        data_io.build_bins([dict(size_um=63, cum_pct=11, status_class="measured"),
                            dict(size_um=250, cum_pct=60, status_class="measured")])
    with pytest.raises(ValueError):                              # goes down
        data_io.build_bins([dict(size_um=63, cum_pct=30, status_class="measured"),
                            dict(size_um=250, cum_pct=20, status_class="measured"),
                            dict(size_um=500, cum_pct=100, status_class="measured")])


def test_milled_curves_have_the_p80_they_claim():
    for psd_id, p80 in (("MILL-P80-1MM", 1000.0), ("MILL-P80-100UM", 100.0)):
        pts = data_io.psd_points(psd_id)
        assert all(p["status_class"] == "assumed" for p in pts)
        (d1, c1), (d2, c2) = [(p["size_um"], p["cum_pct"]) for p in pts[-2:]]
        d80 = d1 * (d2 / d1) ** ((80 - c1) / (c2 - c1))
        assert d80 == pytest.approx(p80, rel=0.01)


# ---------------------------------------------------------------- mineralogy
def test_normative_rows_are_not_taken_for_measurements():
    assert data_io.measured_mineralogy("MB-BRITA2") is None
    assert data_io.build_config()["provenance"]["mineralogy"]["status"] == "assumed"


def test_measured_xrd_is_used_when_present(lab_copy):
    # test fixture: an XRD result for a made-up TECHNOVIA sample
    append_rows(lab_copy / "05_rock_mineralogy.csv", [
        "TV-FINES-01,Plagioclase (An55),40,measured,XRD Rietveld,test fixture",
        "TV-FINES-01,Augite,22,measured,XRD Rietveld,test fixture",
        "TV-FINES-01,Olivine (Fo80),8,measured,XRD Rietveld,test fixture",
        "TV-FINES-01,Sanidine,5,measured,XRD Rietveld,test fixture",
        "TV-FINES-01,Titanomagnetite,7,measured,XRD Rietveld,test fixture",
        "TV-FINES-01,Amorphous (glass),18,measured,XRD Rietveld,test fixture",
    ])
    m = data_io.measured_mineralogy("TV-FINES-01", lab_copy)
    ph = m["phases"]
    assert ph["Plag"] == 40 and ph["Cpx"] == 22 and ph["Kfs"] == 5 and ph["Inert"] == 7 and ph["Gls"] == 18
    assert m["plag_an"] == pytest.approx(0.55) and m["olivine_fo"] == pytest.approx(0.80)
    assert ph["Fo"] + ph["Fa"] == pytest.approx(8.0) and ph["Fo"] > ph["Fa"]
    assert m["total_wt_pct"] == pytest.approx(100.0)


# ---------------------------------------------------------------- soil
def test_soils_with_enough_data_are_selectable():
    soils = data_io.soils()
    assert list(soils) == ["SG-04", "SG-05", "SG-01", "SG-02"]           # nothing from Lomba da Maia yet
    assert all(s["status_class"] == "stand-in" for s in soils.values())
    a = soils["SG-04"]
    assert (a["pH"], a["cec"], a["ca"], a["mg"], a["k"]) == (5.8, 23.61, 3.54, 1.74, 1.2)
    assert "na" not in a and "acidity_pH7" not in a
    full = soils["SG-01"]
    assert full["na"] == 0.17 and full["acidity_pH7"] == 8.23 and full["pH_kcl"] == 4.8


def test_filled_screening_rows_become_selectable(lab_copy):
    p = lab_copy / "09_soil_lomba_da_maia_screening.csv"
    head = p.read_text().splitlines()[0].split(",")
    row = {c: "" for c in head}
    row.update(sample_id="LDM-TEST", site="test fixture", land_use="Pasture", pH_H2O="5.6", CEC_pH7_cmol_kg="20",
               exch_Ca="3", exch_Mg="1.5", exch_K="0.8", exch_Na="0.2", exch_Al_KCl="1.4", bulk_density_g_cm3="0.7",
               status="measured")
    append_rows(p, [",".join(row[c] for c in head)])
    s = data_io.soils(lab_copy)["LDM-TEST"]
    assert s["status_class"] == "measured" and s["acidity_kcl"] == 1.4
    cfg = data_io.build_config("LDM-TEST", lab_dir=lab_copy)
    assert cfg["effective_acidity_cmol_kg"] == 1.4                      # measured KCl acidity replaces the assumed 1.0
    assert cfg["provenance"]["effective_acidity_cmol_kg"]["status"] == "measured"
    assert cfg["bulk_density_kg_m3"] == pytest.approx(700.0)
    assert ("LDM-TEST", []) not in data_io.soil_gaps(lab_copy)


def test_site_mean():
    soils = data_io.soils()
    mean = data_io.site_mean_soil([soils["SG-04"], soils["SG-05"]])
    assert mean["pH"] == pytest.approx(5.5) and mean["cec"] == pytest.approx((23.61 + 9.16) / 2)
    assert mean["members"] == ["SG-04", "SG-05"]
    assert not data_io.in_site_mean(dict(sample_id="LDM-P05-DUP", land_use="Pasture (field duplicate of P05)"))
    assert not data_io.in_site_mean(dict(sample_id="LDM-F01", land_use="Forest (separate zone)"))
    assert data_io.in_site_mean(dict(sample_id="LDM-P01", land_use="Pasture"))


# ---------------------------------------------------------------- configuration
def test_lab_config_matches_the_original_inputs_except_grain_size():
    """The lab files carry the same values as data/inputs.json. Only the grain-size bins differ,
    because they are now interpolated from table 03 instead of typed in by hand."""
    from erw.model import default_config
    lab, old = default_config(), default_config(source="json")
    for k, v in old.items():
        if k in ("psd", "soil", "soil_key", "isometric", "emissions", "export_accounting"):
            continue
        assert lab[k] == pytest.approx(v) if isinstance(v, (int, float)) else lab[k] == v, k
    for k in ("pH", "cec", "ca", "mg", "k"):
        assert lab["soil"][k] == old["soil"][k]
    assert lab["density_g_m3"] == 2.9e6
    assert lab["isometric"]["nfz_depth_cm"] == 20 and lab["isometric"]["eta_river_ocean"] == 0.85
    for k, v in old["emissions"].items():
        if not k.startswith("_") and not k.endswith("status"):
            assert lab["emissions"][k] == v, k
    assert [f for f, _ in lab["psd"]] == pytest.approx([f for f, _ in old["psd"]])     # same mass split


def test_every_input_has_a_status_and_a_source():
    cfg = data_io.build_config()
    prov = cfg["provenance"]
    for key in ("oxides", "soil", "psd", "mineralogy", "density_g_m3", "field_factor", "roughness", "drainage_mm_yr",
                "log_pco2", "temp_c", "application_t_ha", "rates.Plag", "rates.Gls", "emissions.haul_distance_km"):
        assert prov[key]["status"] in ("measured", "stand-in", "assumed", "calculated", "plan"), key
        assert prov[key]["source"], key
    assert prov["oxides"]["status"] == "stand-in" and prov["soil"]["status"] == "stand-in"
    assert not any(p["status"] == "measured" for p in prov.values())      # nothing is the project's own yet
    table = data_io.provenance_table(cfg)
    assert {"parameter", "value", "unit", "status", "source"} <= set(table.columns) and len(table) == len(prov)


def test_old_scenario_names_still_work():
    assert data_io.build_config("B_TerraRiba", "milled_p80_100um")["soil_key"] == "SG-05"
    with pytest.raises(KeyError):
        data_io.build_config("no-such-soil")
    with pytest.raises(KeyError):
        data_io.build_config(rock="TV-FINES-01")                         # oxides not filled in yet


# ---------------------------------------------------------------- plant uptake (table 19)
def test_plant_uptake_is_zero_until_the_share_leaving_the_field_is_known():
    cfg = data_io.build_config()
    assert cfg["uptake_kg_ha_yr"] == dict(Ca=0.0, Mg=0.0, K=0.0)
    assert cfg["uptake_gross_kg_ha_yr"] == dict(Ca=24.0, Mg=12.0, K=183.0)      # 6 t DM/ha x 4, 2 and 30.5 g/kg
    assert cfg["provenance"]["uptake_kg_ha_yr"]["status"] == "to measure"


def test_plant_uptake_from_share_and_milk(lab_copy):
    p = lab_copy / "19_plant_uptake_and_grazing.csv"
    text = p.read_text()
    text = text.replace("Share of herbage cations that leave the field,,", "Share of herbage cations that leave the field,0.5,")
    text = text.replace("Milk sold per hectare,,", "Milk sold per hectare,8000,")      # test fixture
    p.write_text(text)
    cfg = data_io.build_config(lab_dir=lab_copy)
    assert cfg["uptake_kg_ha_yr"]["Ca"] == pytest.approx(0.5 * 24 + 8000 * 1.2 / 1000)
    assert cfg["uptake_kg_ha_yr"]["K"] == pytest.approx(0.5 * 183 + 8000 * 1.3 / 1000)
    assert cfg["provenance"]["uptake_kg_ha_yr"]["status"] == "assumed"
