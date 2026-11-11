"""Download the MODFLOW 6 executable and shared library (mf6 + libmf6) into ./bin for this computer.
Run once:  python scripts/get_modflow.py

First choice is FloPy's own downloader. It asks the GitHub API which release is the latest, and that
call can be refused (rate limit without a GITHUB_TOKEN, or a restricted network). In that case the
script falls back to the direct download link of the same official release
(github.com/MODFLOW-ORG/executables), which needs no API call.
"""
import io
import stat
import sys
import urllib.request
import zipfile
from pathlib import Path

from flopy.utils import get_modflow
from flopy.utils.get_modflow import get_ostag

bindir = Path(__file__).resolve().parents[1] / "bin"
bindir.mkdir(exist_ok=True)

WANTED = {"mf6", "mf6.exe", "libmf6.so", "libmf6.dylib", "libmf6.dll"}


def direct_download():
    url = f"https://github.com/MODFLOW-ORG/executables/releases/latest/download/{get_ostag()}.zip"
    print("Falling back to direct download:", url)
    with urllib.request.urlopen(url, timeout=300) as resp:
        archive = zipfile.ZipFile(io.BytesIO(resp.read()))
    found = []
    for info in archive.infolist():
        name = Path(info.filename).name
        if name in WANTED:
            target = bindir / name
            target.write_bytes(archive.read(info))
            target.chmod(target.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            found.append(name)
    if len(found) < 2:
        sys.exit(f"Release archive did not contain mf6 and libmf6 (found: {found})")


try:
    get_modflow(str(bindir), subset="mf6,libmf6")
except Exception as err:  # GitHub API refused or unreachable
    print(f"FloPy downloader failed ({err}).")
    direct_download()
print("MODFLOW 6 installed in", bindir, "->", sorted(p.name for p in bindir.iterdir()))
