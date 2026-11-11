"""Download the MODFLOW 6 executable and shared library (mf6 + libmf6) into ./bin for this computer.
Run once:  python scripts/get_modflow.py"""
from pathlib import Path
from flopy.utils import get_modflow
bindir = Path(__file__).resolve().parents[1] / "bin"
bindir.mkdir(exist_ok=True)
get_modflow(str(bindir), subset="mf6,libmf6")
print("MODFLOW 6 installed in", bindir)
