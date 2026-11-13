"""
Does the choice of thermodynamic database change the carbonate-system results?

The silicate kinetics need llnl.dat (phreeqc.dat lacks the basalt silicates). The carbonate system
(pH, dissolved inorganic carbon, calcite saturation, what is kept after degassing) can be computed
with any database, and Isometric notes that the choice affects such results. This module takes the
drainage water of a run and repeats the carbonate calculation with llnl.dat and with phreeqc.dat
(shipped with phreeqpython), and reports the difference.

    table = compare(water, log_pco2_soil=-2.0, temp_c=17.0)

`water` is a dict of total concentrations in mol/kgw with the keys alk, ca, mg, na, k, cl, si
(the yearly drainage chemistry of the soil model).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import phreeqpython

ROOT = Path(__file__).resolve().parents[1]
DATABASES = {
    "llnl.dat": ROOT / "data" / "llnl.dat",
    "phreeqc.dat": Path(phreeqpython.__file__).parent / "database" / "phreeqc.dat",
}
ATMOSPHERE_LOG_PCO2 = -3.38        # about 420 ppm


def _carbonate(db_path: Path, water: dict, log_pco2_soil: float, temp_c: float, log_pco2_air: float) -> dict:
    pp = phreeqpython.PhreeqPython()
    pp.ip.load_database_string(db_path.read_text(encoding="latin-1"))
    w = {k: max(float(water.get(k, 0.0)), 0.0) for k in ("alk", "ca", "mg", "na", "k", "cl", "si")}
    pp.ip.run_string(f"""SOLUTION 1 drainage water at the soil CO2
  temp {temp_c}
  units mol/kgw
  pH 7 CO2(g) {log_pco2_soil}
  Alkalinity {w['alk']:.6e}
  Ca {w['ca']:.6e} ; Mg {w['mg']:.6e} ; Na {w['na']:.6e} ; K {w['k']:.6e} ; Cl {w['cl']:.6e} ; Si {w['si']:.6e}
SELECTED_OUTPUT 1
  -reset false
USER_PUNCH 1
  -headings pH dic alk si_calcite calcite
  10 PUNCH -LA("H+"), TOT("C"), ALK, SI("Calcite"), EQUI("Calcite")
END
USE solution 1
EQUILIBRIUM_PHASES 1
  CO2(g) {log_pco2_air} 10
  Calcite 0 0
END""")
    arr = pp.ip.get_selected_output_array()
    soil, air = arr[1], arr[-1]
    return dict(pH_soil=float(soil[0]), dic_soil_mmol=float(soil[1]) * 1e3, si_calcite_soil=float(soil[3]),
                pH_air=float(air[0]), dic_air_mmol=float(air[1]) * 1e3, alk_air_meq=float(air[2]) * 1e3,
                calcite_formed_mmol=float(air[4]) * 1e3)


def compare(water: dict, log_pco2_soil: float = -2.0, temp_c: float = 17.0,
            log_pco2_air: float = ATMOSPHERE_LOG_PCO2) -> pd.DataFrame:
    """One row per database, plus a row with the difference (phreeqc.dat minus llnl.dat)."""
    rows = {name: _carbonate(path, water, log_pco2_soil, temp_c, log_pco2_air) for name, path in DATABASES.items()}
    df = pd.DataFrame(rows).T
    df.loc["difference (phreeqc.dat minus llnl.dat)"] = df.loc["phreeqc.dat"] - df.loc["llnl.dat"]
    return df.rename(columns={
        "pH_soil": "pH at soil CO2", "dic_soil_mmol": "DIC at soil CO2 (mmol/kg)",
        "si_calcite_soil": "Calcite saturation index at soil CO2", "pH_air": "pH after degassing to open air",
        "dic_air_mmol": "DIC after degassing (mmol/kg)", "alk_air_meq": "Alkalinity after degassing (meq/kg)",
        "calcite_formed_mmol": "Calcite formed on degassing (mmol/kg)"})


def compare_run(result, cfg, year: int | None = None) -> dict:
    """The comparison for the rock and the control drainage of a run, and for the extra carbon
    (rock minus control) kept after degassing, which is the number that matters for removal."""
    d = result["hist"]["drainage"][(year or len(result["hist"]["drainage"])) - 1]
    out = {kind: compare(d[kind], cfg["log_pco2"], cfg["temp_c"]) for kind in ("rock", "ctrl")}
    col = "DIC after degassing (mmol/kg)"
    extra_alk = (d["rock"]["alk"] - d["ctrl"]["alk"]) * 1e3
    kept = {}
    for db in DATABASES:
        extra_dic = out["rock"].loc[db, col] - out["ctrl"].loc[db, col]
        kept[db] = extra_dic / extra_alk if extra_alk > 0 else float("nan")
    out["extra_alkalinity_meq"] = extra_alk
    out["carbon_kept_per_alkalinity"] = kept
    return out
