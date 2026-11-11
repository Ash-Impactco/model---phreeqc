"""Shared test set-up: make the project importable and give tests a private copy of the data pack."""
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


@pytest.fixture
def lab_copy(tmp_path):
    """A throw-away copy of data/lab that a test may edit (to stand in for results not yet received)."""
    dst = tmp_path / "lab"
    shutil.copytree(ROOT / "data" / "lab", dst)
    return dst


def append_rows(path: Path, lines: list[str]):
    raw = path.read_bytes()
    nl = b"\r\n" if b"\r\n" in raw else b"\n"
    if not raw.endswith(nl):
        raw += nl
    path.write_bytes(raw + nl.join(line.encode() for line in lines) + nl)
