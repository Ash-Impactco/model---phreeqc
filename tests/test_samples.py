"""Per-sample runs in a process pool, with a site mean."""
from erw.samples import run_many, sample_configs, summary_row


def test_screening_table_is_empty_so_nothing_runs():
    configs, skipped = sample_configs("09")
    assert configs == []
    assert len(skipped) == 12 and all("pH_H2O" in missing for _, missing in skipped)


def test_one_run_per_composite_plus_site_mean_in_parallel():
    # the Sanguinho stand-ins (table 08) show the mechanism; one simulated year keeps the test short
    configs, skipped = sample_configs("08", overrides={"years": 1})
    names = [n for n, _ in configs]
    assert names == ["SG-04", "SG-05", "SG-01", "SG-02", "SITE-MEAN"]
    assert {s for s, _ in skipped} == {"SG-H1", "SG-H2", "SG-H3", "SG-HEMP3", "SG-TEA"}
    mean_cfg = dict(configs)["SITE-MEAN"]
    assert mean_cfg["soil"]["members"] == ["SG-04", "SG-05", "SG-01", "SG-02"]
    assert "arithmetic mean" in mean_cfg["provenance"]["soil"]["source"]

    results = run_many(configs, max_workers=4)
    assert [n for n, _, _ in results] == names
    rows = [summary_row(*r) for r in results]
    assert all(r["released_tCO2e_ha"] > 0 and r["years"] == 1 for r in rows)
    assert all(r["soil_status"] == "stand-in" and r["status"] == "model forecast, not creditable" for r in rows)
    released = {r["sample"]: r["released_tCO2e_ha"] for r in rows}
    assert released["SG-05"] > released["SG-01"]                # the more acid soil (pH 5.2) dissolves more than pH 6.0
