"""Give sander a short path to the ORCA launcher so it can start ORCA even when
the run directory has a long name.
"""

from contextlib import contextmanager
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from common import sha


@contextmanager
def short_orca_launcher(run, base_dir="/tmp", path_limit=64):
    """Keep one copied launcher alive for the complete serial Sander call.

    Production calls use /tmp and a full filename limit of 64 characters.
    No quantum-chemistry inputs/checkpoints/outputs are moved to local scratch.
    """
    run = Path(run).resolve()
    source = run / "bin/orca"
    if not source.is_file():
        raise ValueError("Frozen ORCA launcher is missing")
    with tempfile.TemporaryDirectory(prefix="aoc-", dir=base_dir) as temporary:
        folder = Path(temporary).resolve()
        launcher = folder / "orca"
        if len(str(launcher)) > path_limit or " " in str(launcher):
            raise ValueError("Local ORCA launcher path is not short enough: " + str(launcher))
        folder.chmod(0o700)
        shutil.copyfile(source, launcher)
        launcher.chmod(0o700)
        if sha(launcher) != sha(source):
            raise ValueError("Local ORCA launcher copy differs from frozen source")
        env = os.environ.copy()
        env["AMBER_ORCA_SHIM_DIR"] = str(folder)
        env["PATH"] = str(folder) + os.pathsep + env.get("PATH", "")
        discovered = shutil.which("orca", path=env["PATH"])
        if discovered is None or Path(discovered).resolve() != launcher:
            raise ValueError("ORCA discovery does not select the private short launcher")
        # Execute directly so a noexec filesystem is detected, rather than
        # bypassing the kernel's execute check with "bash launcher".
        check = subprocess.run(
            [str(launcher), "--check-launcher"],
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if check.returncode or check.stdout.strip() != "QMMM_LAUNCHER_OK":
            raise ValueError("Short ORCA launcher execution check failed: " + check.stdout[-2000:])
        yield env, dict(
            launcher=str(launcher),
            directory=str(folder),
            path_length=len(str(launcher)),
            path_limit=path_limit,
            source=str(source),
            source_sha256=sha(source),
            discovery_path=discovered,
            direct_execution_check_passed=True,
            scope="Only the launcher is node-local; all calculation data stay in the run directory",
        )
