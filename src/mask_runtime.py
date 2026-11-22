"""Ask cpptraj which atoms each AMBER mask selects. Check those atom numbers
against the QM and movable MM selections before starting the calculation.
"""

from pathlib import Path
import re
import shutil
import subprocess
from common import text


def check_native_mask(prmtop, mask, expected, log=None):
    executable = shutil.which("cpptraj")
    if not executable:
        raise ValueError("Installed cpptraj is required for the native moving-mask check")
    argv = [executable, "-p", str(Path(prmtop).resolve()), "-ms", mask]
    process = subprocess.run(argv, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if log is not None:
        text(log, process.stdout)
    found = re.findall(r"^Selected=([^\n]*)$", process.stdout, re.M)
    if process.returncode or len(found) != 1:
        raise ValueError(f"CPPTRAJ could not verify mask {mask}: " + process.stdout[-4000:])
    try:
        selected = [int(v) - 1 for v in found[0].split()]
    except ValueError as exc:
        raise ValueError("Unexpected CPPTRAJ selection output") from exc
    if selected != sorted(set(expected)):
        raise ValueError(f"CPPTRAJ mask {mask} selects different atoms from the input")
    return dict(
        passed=True,
        executable=executable,
        argv=argv,
        selected_indices_0based=selected,
        selected_count=len(selected),
        scope="Native static mask evaluation only; no SCF, forces, minimization or MD",
    )
