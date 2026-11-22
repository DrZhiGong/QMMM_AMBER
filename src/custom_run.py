"""Run one calculation on the compute node. Relax the selected MM atoms around the
fixed QM region, run a final single point and save the results in result.json.
"""

import argparse
import os
from pathlib import Path
import traceback

from common import read, read_rst, sha, verify, write
from custom_model import load_system
from endpoint import relaxation_status
from mask_runtime import check_native_mask
from sander import frozen_check, model_context, run_sander, stage


def calculate(run):
    stage(run, "checking_inputs")
    verify(run)
    inputs = run / "inputs"
    model = load_system(inputs / "system.json")
    config = read(run / "config.json")
    topology = inputs / "system.prmtop"
    initial = inputs / "initial.rst7"
    coords = read_rst(initial, model["natom"])
    qm, moving = model["qm_indices_0based"], model["movable_mm_indices_0based"]
    charge, mult = model["qm_charge"], model["multiplicity"]
    context = model_context(config, model, coords)
    check_native_mask(topology, context["qm_mask"], qm, run / "logs/qm_mask.cpptraj.log")
    check_native_mask(topology, context["movable_mm_mask"], moving, run / "logs/moving_mask.cpptraj.log")

    stage(run, "coupled_mm_relaxation")
    target = config["drms_kcal_mol_angstrom"]
    relaxation = run_sander(
        run, run / "relax", topology, initial, context,
        cycles=config["max_cycles"], moving=True, drms=target,
    )
    restart = Path(relaxation["restart_file"])
    frozen = frozen_check(coords, read_rst(restart, model["natom"]), moving, qm)

    stage(run, "full_qmmm_final_single_point")
    final = run_sander(run, run / "final", topology, restart, context, cycles=0, seed_result=relaxation["bridge"])
    bridge = final["bridge"]
    expected_spin = (mult, model["nalpha"], model["nbeta"])
    endpoint = relaxation_status(run / "relax", relaxation, final, target, expected_spin=expected_spin)
    gbw = Path(bridge["gbw_file"])
    if gbw.resolve() != run / "final/orc_job.gbw" or sha(gbw) != bridge["gbw_sha256"]:
        raise ValueError("The final native GBW differs from the checked single point")
    write(run / "result.json", dict(
        **endpoint, qm_charge=charge, multiplicity=mult,
        nalpha=bridge["nalpha"], nbeta=bridge["nbeta"], s_squared=bridge["s_squared"],
        atom_mulliken_populations=bridge["atom_mulliken_populations"],
        embedded_qm_energy_hartree=bridge["energy_hartree"],
        full_qmmm_potential_energy_kcal_mol=final["sander_energy_kcal_mol"],
        full_qmmm_energy_audit=final["energy_audit"], geometry_audit=frozen, drms_kcal_mol_angstrom=target,
        relaxed_structure=str(restart), restart_sha256=sha(restart),
        gbw_file=str(gbw), gbw_sha256=bridge["gbw_sha256"],
        pointcharges_file=bridge["pointcharges_file"], orbital_stability_established=False,
    ))
    converged = endpoint["mm_geometry_converged"]
    stage(run, "finished", status=endpoint["status"], mm_geometry_converged=converged)
    print(f"Results in {run / 'result.json'}", flush=True)
    return converged is True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    job = os.environ.get("SLURM_JOB_ID")
    if not job:
        parser.error("Run AMBER/ORCA with the submitted compute-node worker")
    run = args.run.resolve(strict=True)
    with (run / "setup.claim").open("x", encoding="utf-8") as claim:
        claim.write(job + "\n")
    try:
        return 0 if calculate(run) else 2
    except Exception as exc:
        stage(run, "FAILED", error=str(exc), traceback=traceback.format_exc())
        raise


if __name__ == "__main__":
    raise SystemExit(main())
