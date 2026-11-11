"""Each validation rule, on small made-up tables (test fixtures), then the real data pack."""
import pandas as pd

from conftest import append_rows
from erw import data_io, validate


def levels(findings, rule=None):
    return [f.level for f in findings if rule is None or f.rule == rule]


def rock_row(**kw):
    base = dict(sample_id="R1", SiO2="47.4", Al2O3="14.0", Fe2O3="13.0", FeO="", CaO="10.0", MgO="6.9", Na2O="3.4",
                K2O="2.0", TiO2="3.9", MnOx="0.1", P2O5="", CO2="", LOI="<0.05", SO3="")
    base.update(kw)
    return base


def soil_row(**kw):
    base = dict(sample_id="S1", plot="", notes="", pH_H2O="6.0", pH_KCl="4.8", CEC_pH7_cmol_kg="18.0", exch_Ca="5.76",
                exch_Mg="1.59", exch_K="2.25", exch_Na="0.17", exch_acidity_pH7="8.23", sum_bases="9.77",
                base_saturation_pct="54.28", K_extractable_mg_kg_K2O="902", P_extractable_mg_kg_P2O5="332")
    base.update(kw)
    return base


# ---------------------------------------------------------------- rule 0: structure
def test_structure_rule(lab_copy):
    assert levels(validate.check_structure(data_io.read_table("03", lab_copy))) == ["pass"]
    append_rows(lab_copy / "03_rock_grain_size.csv", ["TEST,10,5,assumed,fixture,note,extra cell"])
    assert levels(validate.check_structure(data_io.read_table("03", lab_copy))) == ["fail"]


# ---------------------------------------------------------------- rule 1: oxide totals
def test_oxide_total_passes_between_97_and_102():
    f = validate.check_oxide_totals(pd.DataFrame([rock_row()]))
    assert levels(f) == ["pass"] and "100.70 to 100.75" in f[0].message      # LOI below detection widens the total


def test_oxide_total_fails_when_low_or_high():
    low = validate.check_oxide_totals(pd.DataFrame([rock_row(SiO2="40.0")]))
    high = validate.check_oxide_totals(pd.DataFrame([rock_row(SiO2="52.0")]))
    assert levels(low) == ["fail"] and levels(high) == ["fail"]


def test_oxide_total_uses_feo_loi_and_co2():
    with_feo = validate.check_oxide_totals(pd.DataFrame([rock_row(Fe2O3="3.0", FeO="9.0", LOI="0.5")]))
    assert levels(with_feo) == ["pass"] and "100.20" in with_feo[0].message
    # CO2 is part of LOI: counted only when LOI is missing
    a = validate.check_oxide_totals(pd.DataFrame([rock_row(LOI="1.0", CO2="0.8")]))[0].message
    b = validate.check_oxide_totals(pd.DataFrame([rock_row(LOI="", CO2="0.8")]))[0].message
    assert "101.70" in a and "101.50" in b


def test_oxide_total_not_checked_when_incomplete():
    assert levels(validate.check_oxide_totals(pd.DataFrame([rock_row(CaO="")]))) == ["info"]


# ---------------------------------------------------------------- rule 2: base saturation
def test_base_saturation_within_ten_percent():
    assert levels(validate.check_base_saturation(pd.DataFrame([soil_row()]))) == ["pass"]
    assert levels(validate.check_base_saturation(pd.DataFrame([soil_row(base_saturation_pct="57.0")]))) == ["pass"]


def test_base_saturation_fails_beyond_ten_percent():
    assert levels(validate.check_base_saturation(pd.DataFrame([soil_row(base_saturation_pct="70")]))) == ["fail"]
    assert levels(validate.check_base_saturation(pd.DataFrame([soil_row(base_saturation_pct="40")]))) == ["fail"]


def test_base_saturation_not_checked_without_cations():
    assert levels(validate.check_base_saturation(pd.DataFrame([soil_row(exch_Ca="")]))) == ["info"]


# ---------------------------------------------------------------- rule 3: pH
def test_ph_kcl_must_be_below_ph_water():
    assert levels(validate.check_ph(pd.DataFrame([soil_row()]))) == ["pass"]
    assert levels(validate.check_ph(pd.DataFrame([soil_row(pH_KCl="6.3")]))) == ["fail"]
    assert levels(validate.check_ph(pd.DataFrame([soil_row(pH_KCl="")]))) == []
    assert levels(validate.check_ph(pd.DataFrame([soil_row(pH_H2O="58", pH_KCl="")]))) == ["fail"]     # typing slip


# ---------------------------------------------------------------- rule 4: units
def test_units_soil_consistent():
    assert set(levels(validate.check_units_soil(pd.DataFrame([soil_row()])))) == {"pass"}


def test_units_soil_flags_unit_notes_text_and_bad_sums():
    note = validate.check_units_soil(pd.DataFrame([soil_row(K_extractable_mg_kg_K2O="281 (as K)")]))
    assert "warn" in levels(note)
    text = validate.check_units_soil(pd.DataFrame([soil_row(exch_Ca="trace")]))
    assert "fail" in levels(text)
    wrong_sum = validate.check_units_soil(pd.DataFrame([soil_row(sum_bases="12.0")]))
    assert "fail" in levels(wrong_sum)
    # cations typed in mg/kg instead of cmol(+)/kg would exceed the CEC
    wrong_unit = validate.check_units_soil(pd.DataFrame([soil_row(exch_Ca="1150", sum_bases="")]))
    assert any("exceed the CEC" in f.message for f in wrong_unit)


def test_units_long_table_one_unit_per_property():
    ok = pd.DataFrame([dict(sample_id="A", property="Density", value="2900", unit="kg/m3"),
                       dict(sample_id="B", property="Density", value="2400", unit="kg/m3")])
    bad = pd.DataFrame([dict(sample_id="A", property="Density", value="2900", unit="kg/m3"),
                        dict(sample_id="B", property="Density", value="2.4", unit="g/cm3")])
    assert levels(validate.check_units_long(ok, "property")) == ["pass"]
    assert levels(validate.check_units_long(bad, "property")) == ["fail"]


# ---------------------------------------------------------------- rule 5: duplicate ids
def test_duplicate_ids():
    ok = pd.DataFrame([dict(sample_id="A"), dict(sample_id="B")])
    bad = pd.DataFrame([dict(sample_id="A"), dict(sample_id="B"), dict(sample_id="A")])
    assert levels(validate.check_duplicate_ids(ok)) == ["pass"]
    f = validate.check_duplicate_ids(bad)
    assert levels(f) == ["fail"] and f[0].item == "A"
    long = pd.DataFrame([dict(sample_id="A", property="x"), dict(sample_id="A", property="y")])
    assert levels(validate.check_duplicate_ids(long, ["property"])) == ["pass"]


# ---------------------------------------------------------------- rule 6: known issues
def test_known_issues_are_flagged_from_the_notes():
    soil = pd.DataFrame([soil_row(sample_id="SG-04", plot="Moreira", K_extractable_mg_kg_K2O="1019", exch_K="1.2",
                                  notes="A second Notion table lists K2O 286: check original report"),
                         soil_row(sample_id="SG-01")])
    rock = pd.DataFrame([rock_row(MnOx="0.09 (reported as Mn2O4, likely Mn3O4)"), rock_row(sample_id="R2")])
    f = validate.check_known_issues(rock, soil)
    assert levels(f) == ["warn", "warn"]
    assert f[0].item == "SG-04" and "2.16 cmol" in f[0].message and f[1].item == "R1"


# ---------------------------------------------------------------- rule 7: CO2 potential
def test_co2_potential_of_the_feedstock():
    ox = data_io.rock_samples()["MB-BRITA2"]["oxides"]
    cons, full = validate.co2_potential(ox)
    assert round(cons, 3) == 0.154 and round(full, 3) == 0.375


def test_co2_table_is_recomputed():
    rock = pd.DataFrame([rock_row(sample_id="MB-BRITA2", CaO="10.05", MgO="6.92", Na2O="3.37", K2O="2.0")])
    good = pd.DataFrame([dict(sample_id="MB-BRITA2", tCO2_per_t_conservative_CaMg_only="0.154",
                              tCO2_per_t_bicarbonate_full="0.375")])
    bad = pd.DataFrame([dict(sample_id="MB-BRITA2", tCO2_per_t_conservative_CaMg_only="0.154",
                             tCO2_per_t_bicarbonate_full="0.400")])
    assert levels(validate.check_co2_potential(rock, good)) == ["pass"]
    assert levels(validate.check_co2_potential(rock, bad)) == ["fail"]


# ---------------------------------------------------------------- the real data pack
def test_the_data_pack_as_it_stands(tmp_path):
    f = validate.validate_all()
    fails = [(x.rule, x.item) for x in f if x.level == "fail"]
    assert fails == [("oxide total", "MB-POLISHED")]                      # 93.2 %: an incomplete analysis
    warns = [(x.rule, x.item) for x in f if x.level == "warn"]
    assert ("known issue", "SG-04") in warns and ("known issue", "SG-05") in warns     # potassium conflict
    assert ("known issue", "MB-BRITA2") in warns                                         # Mn2O4 or Mn3O4
    assert ("units", "SG-TEA") in warns                                                  # P and K not as oxides
    report = validate.write_report(tmp_path / "validation.md").read_text()
    assert "# Data validation report" in report and "MB-POLISHED" in report and "## Status of the data" in report
