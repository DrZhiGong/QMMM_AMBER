"""Use qmmm to check your input files, submit a job, and read its results."""

import argparse
from datetime import datetime
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

if sys.version_info < (3, 11):
    raise SystemExit("Run this program with Python 3.11 or newer.")

ROOT = Path(__file__).resolve().parent
DEFAULT_SYSTEM = ROOT / "examples/water_pair/system.json"
# An editable install keeps runs in the checkout. A wheel uses the caller's working directory.
RUNS = (ROOT if (ROOT / "requirements.txt").is_file() else Path.cwd()) / "runs"
sys.path.insert(0, str(ROOT / "src"))
from common import read, seal, write
from custom_model import describe_system, load_system, read_system_spec
from endpoint import classify_mm_gradient


def load_settings():
    try:
        settings = read(ROOT / "settings.json")
    except json.JSONDecodeError as exc:
        raise ValueError(f"settings.json, line {exc.lineno}, column {exc.colno}: {exc.msg}") from exc
    if not isinstance(settings, dict):
        raise ValueError("settings.json must contain a JSON object with the calculation settings.")
    target = settings.get("drms_kcal_mol_angstrom")
    if type(target) not in (int, float) or not 0 < target <= sys.float_info.max:
        raise ValueError("settings.json requires a positive, finite 'drms_kcal_mol_angstrom'.")
    for key in ("orca_ranks", "maxcore_mb_per_rank", "max_cycles"):
        value = settings.get(key)
        if type(value) is not int or value < 1:
            raise ValueError(f"settings.json requires a positive integer for '{key}'.")
    return settings


def copy_inputs(system, inputs):
    spec = read_system_spec(system)
    inputs.mkdir()
    shutil.copyfile(spec["topology"], inputs / "system.prmtop")
    shutil.copyfile(spec["coordinates"], inputs / "initial.rst7")
    write(inputs / "system.json", {**spec, "topology": "system.prmtop", "coordinates": "initial.rst7"})
    describe_system(load_system(inputs / "system.json"), stream=sys.stderr)


def snapshot(settings, system):
    parent = RUNS
    parent.mkdir(exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    run = Path(tempfile.mkdtemp(dir=parent, prefix=f"qm-mm-{stamp}-"))
    try:
        for name in ("tools", "logs", "bin"):
            (run / name).mkdir()
        copy_inputs(system, run / "inputs")
        shutil.copytree(ROOT / "dat", run / "dat")
        for src in sorted((ROOT / "src").glob("*.py")):
            shutil.copyfile(src, run / "tools" / src.name)
        for name in ("worker.sh", "runtime.sh"):
            shutil.copyfile(ROOT / "src/runtime" / name, run / "tools" / name)
        launcher = run / "bin/orca"
        shutil.copyfile(ROOT / "src/runtime/orca_wrapper.sh", launcher)
        launcher.chmod(0o700)
        write(run / "config.json", settings)
        seal(run, [p.relative_to(run).as_posix() for p in run.rglob("*") if p.is_file()])
    except Exception:
        shutil.rmtree(run)
        raise
    return run


def show_result(run):
    path = run / "result.json"
    if not path.is_file():
        return False
    status_path = run / "setup_status.json"
    worker = read(status_path) if status_path.is_file() else {}
    if worker.get("stage") != "finished":
        print(f"No finished result (worker stage: {worker.get('stage', 'unknown')}).")
        return False
    try:
        result = read(path)
        energy = result["full_qmmm_potential_energy_kcal_mol"]
        if type(energy) not in (int, float) or not math.isfinite(energy):
            raise ValueError("Invalid full QM/MM potential")
        target = read(run / "config.json")["drms_kcal_mol_angstrom"]
        rms = result["relaxation_rms_gradient"]
        converged, label = classify_mm_gradient(rms, target)
        if (result["drms_kcal_mol_angstrom"] != target
                or result["mm_geometry_converged"] is not converged or result["status"] != label
                or worker["mm_geometry_converged"] is not converged or worker["status"] != label):
            raise ValueError("MM convergence disagrees with the gradient, run settings or worker status")
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid result in {path}: {exc}") from exc
    print(f"Full QM/MM potential: {energy:.6f} kcal/mol")
    if converged is None:
        print("MM convergence unverified: no verified relaxation gradient for the final single point.")
    else:
        print(f"MM converged: {converged} (RMS gradient {rms:.6g}, target {target:.6g} kcal/mol/angstrom)")
    print(f"Results: {path}")
    return converged is True


def status(run):
    failure_log = run / "logs/submission_error.txt"
    if failure_log.is_file():
        print("Job startup issue: " + failure_log.read_text(encoding="utf-8", errors="replace").strip())
    if (run / "jobs.tsv").is_file():
        lines = (run / "jobs.tsv").read_text().splitlines()
        jobs = dict(line.split("\t", 1) for line in lines if line)
    elif (run / "jobs.json").is_file():  # Runs submitted with the previous driver.
        jobs = read(run / "jobs.json")
    else:
        jobs = {}
    if jobs:
        print("Slurm jobs: " + ", ".join(f"{phase} {job}" for phase, job in jobs.items()), flush=True)
        if shutil.which("sacct"):
            subprocess.run(["sacct", "-X", "-j", ",".join(map(str, jobs.values())),
                            "--format=JobID%22,State%22,ExitCode,Elapsed"])
    worker_status = run / "setup_status.json"
    if worker_status.is_file():
        worker = read(worker_status)
        print("Stage:", worker.get("stage", "unknown"))
        if worker.get("error"):
            print(worker["error"])
        # Older multi-state runs keep separate status files under each state.
        for path in sorted(run.glob("states/*/status.json")):
            item = read(path)
            print(f"{path.parent.name}: {item.get('error') or item.get('stage', 'unknown')}")
    elif not failure_log.is_file():
        print("Inputs prepared. Waiting for the worker." if jobs else "Inputs prepared.")
    if (run / "logs/modules.log").is_file():
        print(f"Software check: {run / 'logs/modules.log'}")
    show_result(run)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("check", help="verify the inputs locally")
    submit = sub.add_parser("submit", help="submit one calculation to Slurm")
    submit.add_argument("--account", help="Slurm allocation")
    submit.add_argument("--partition", help="Slurm partition")
    prepare = sub.add_parser("prepare", help="check the inputs and print the new run directory")
    for command in (check, prepare, submit):
        command.add_argument("system", nargs="?", type=Path, default=DEFAULT_SYSTEM, metavar="system.json",
                             help="AMBER system description (default: the water pair example)")
    for command in ("status", "collect"):
        sub.add_parser(command).add_argument(
            "run", nargs="?", type=Path, help="default: last prepared run"
        )
    args = parser.parse_args()
    if args.command == "check":
        load_settings()
        describe_system(load_system(args.system))
        return
    if args.command == "prepare":
        run = snapshot(load_settings(), args.system)
        write(RUNS / "latest.json", {"run": str(run)})
        print(run)
        return
    if args.command == "submit":
        settings = load_settings()
        system = args.system.resolve(strict=True)
        describe_system(load_system(system))
        cmd = ["sbatch", f"--ntasks={settings['orca_ranks']}"]
        if args.account:
            cmd.append(f"--account={args.account}")
        if args.partition:
            cmd.append(f"--partition={args.partition}")
        cmd.extend([str(ROOT / "run.sh"), str(system)])
        env = {**os.environ, "QMMM_PROJECT_DIR": str(ROOT), "QMMM_PYTHON": sys.executable}
        subprocess.run(cmd, check=True, env=env)
        return
    run = args.run
    if run is None:
        latest = RUNS / "latest.json"
        if not latest.is_file():
            parser.error("No prepared run yet. Check squeue while your Slurm job is pending. "
                         "Once it starts, slurm-<job ID>.out gives the run directory.")
        run = Path(read(latest)["run"])
    run = run.resolve()
    if not (run / "config.json").is_file():
        parser.error(f"Could not find config.json in {run}. Pass the run directory printed "
                     "during preparation or submission.")
    print(f"Run: {run}", flush=True)
    if args.command == "status":
        status(run)
    elif not (run / "inputs/system.json").is_file() and (run / "tools/report.py").is_file():
        # An archived run can still collect its results using its own frozen worker.
        subprocess.run([sys.executable, "-B", str(run / "tools/report.py"), str(run)], check=True)
    elif (run / "result.json").is_file():
        raise SystemExit(0 if show_result(run) else 2)
    else:
        print("No result yet. Use status to follow this run.")
        raise SystemExit(2)


def cli():
    try:
        main()
    except subprocess.CalledProcessError as exc:
        raise SystemExit(exc.returncode)
    except (OSError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc))


if __name__ == "__main__":
    cli()
