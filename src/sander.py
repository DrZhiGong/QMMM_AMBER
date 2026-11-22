"""Write the inputs for one sander job, run it with ORCA and read back the energy and gradient.
The same steps handle both the MM relaxation and the final single point.
"""

from datetime import datetime, timezone
import math
import os
from pathlib import Path
import re
import shutil
from common import ELEMENTS, read, run_command, sha, text, write
from inputs_control import sander_input, orca_template
from short_launcher import short_orca_launcher
from energy_precision import energy_matches


def stage(run, stage_name, **extra):
    write(
        run / "setup_status.json",
        dict(
            stage=stage_name,
            updated_at=datetime.now(timezone.utc).isoformat(),
            job_id=os.environ.get("SLURM_JOB_ID"),
            **extra,
        ),
    )
    message = {
        "checking_inputs": "Checking the topology and atom selections",
        "coupled_mm_relaxation": "Relaxing the movable MM atoms",
        "full_qmmm_final_single_point": "Calculating the final QM/MM single point",
        "RELAXED_UNREVIEWED": "The final MM gradient met the target",
        "MM_NOT_CONVERGED_UNREVIEWED": "The final MM gradient is above the target",
        "FINAL_SINGLE_POINT_MM_UNVERIFIED": "The final geometry could not be matched to an MM gradient",
    }.get(extra.get("status", stage_name), stage_name.capitalize())
    if extra.get("error"):
        message += f": {extra['error']}"
    print(message, flush=True)


def model_context(config, model, coords):
    """Describe the fixed QM atoms and the QM parent of each boundary link."""
    qm = model["qm_indices_0based"]
    qm_position = {atom: i for i, atom in enumerate(qm)}
    return dict(
        charge=model["qm_charge"],
        input_multiplicity=model["multiplicity"],
        nprocs=config["orca_ranks"],
        maxcore_mb=config["maxcore_mb_per_rank"],
        qm_mask=model["qm_mask"],
        movable_mm_mask=model["movable_mm_mask"],
        real_qm_atoms=[
            {"element": ELEMENTS[model["atomic_numbers"][i] - 1], "coord": coords[i]} for i in qm
        ],
        cap_parent_real_indices=[
            qm_position[a] if a in qm_position else qm_position[b] for a, b in model["link_bonds_0based"]
        ],
        expected_qm_atoms_with_links=model["expected_qm_atoms_with_links"],
        expected_mm_charge_sum_e=model["expected_mm_charge_sum_e"],
    )


def _write_coupled_inputs(folder, context, cycles, moving, drms, seed_result):
    """Set the QM and movable-MM masks, and stage the orbitals for the first ORCA call."""
    if seed_result:
        source = Path(seed_result["gbw_file"])
        if sha(source) != seed_result["gbw_sha256"]:
            raise ValueError("Seed native GBW changed")
        shutil.copyfile(source, folder / "orbital_seed.gbw")
        context["native_real_qm_indices"] = seed_result["native_real_qm_indices"]
        context["native_link_parents"] = seed_result["native_link_parents"]
    write(folder / "bridge_context.json", context)
    text(
        folder / "sander.in",
        sander_input(
            qm_mask=context["qm_mask"],
            multiplicity=context["input_multiplicity"],
            cycles=cycles,
            drms=drms,
            moving_mask=context["movable_mm_mask"] if moving else None,
            charge=context["charge"],
        ),
    )
    text(
        folder / "orc_job.tpl",
        orca_template(context, has_seed=bool(seed_result)),
    )


def _read_coupled_result(folder, moving):
    out = (folder / "sander.out").read_text(errors="replace")
    if re.search(r"\b(?:SANDER BOMB|QM calculation failed|Error in QM calculation)\b", out, re.I):
        raise RuntimeError("Sander reported an external-QM failure")
    if re.search(r"\b(?:NaN|Infinity)\b", out, re.I):
        raise ValueError("Sander reported a nonfinite energy/gradient")
    if not (folder / "bridge_last.json").is_file():
        raise ValueError("No checked ORCA evaluation returned through Sander")
    bridge = read(folder / "bridge_last.json")
    # A minimization can return to an earlier step after a line search.
    for _, parsed in energy_matches(out, bridge, all_events=moving):
        return bridge, parsed
    raise ValueError("No archived ORCA gradient reproduces the final Sander EXTERNESCF")


def run_sander(
    run, folder, prmtop, coordinates, context, cycles=0, moving=False, drms=0.01, seed_result=None
):
    """Run one sander job and return its final energy, gradient norm and ORCA orbitals."""
    run = Path(run).resolve()
    folder = Path(folder).resolve()
    folder.mkdir(parents=True)
    context = context.copy()
    _write_coupled_inputs(folder, context, cycles, moving, drms, seed_result)
    command = [
        os.environ["AMBER_SANDER"], "-O",
        "-i", "sander.in",
        "-o", "sander.out",
        "-p", str(Path(prmtop).resolve()),
        "-c", str(Path(coordinates).resolve()),
        "-r", "relaxed.rst7",
        "-inf", "mdinfo",
    ]
    write(folder / "sander_command.json", dict(argv=command, serial_sander=True, context=context))
    # AMBER truncates long executable paths when launching ORCA.
    with short_orca_launcher(run, base_dir="/tmp", path_limit=64) as (env, launcher):
        env.update(AMBER_ORCA_RUN=str(run), AMBER_ORCA_CONTEXT=str(folder / "bridge_context.json"))
        write(folder / "launcher.json", launcher)
        run_command(command, folder, folder / "sander.launch.log", env)

    bridge, parsed = _read_coupled_result(folder, moving)
    rms = parsed["sander_rms_gradient"]
    result = dict(
        bridge=bridge,
        **parsed,
        requested_cycles=cycles,
        geometry_converged=(cycles > 0 and rms <= drms),
        folder=str(folder),
        restart_file=str(folder / "relaxed.rst7"),
    )
    write(folder / "evaluation_summary.json", result)
    return result


def frozen_check(before, after, moving, qm):
    """Measure displacements in angstroms against the starting coordinates."""
    moving_indices = set(moving)
    displacements = [math.dist(start, end) for start, end in zip(before, after)]
    max_frozen = max(d for i, d in enumerate(displacements) if i not in moving_indices)
    max_qm = max(displacements[i] for i in qm)
    if max_frozen > 2e-6 or max_qm > 2e-6:
        raise ValueError("Frozen or real-QM coordinates moved")
    return dict(
        max_frozen_displacement_angstrom=max_frozen,
        max_real_qm_displacement_angstrom=max_qm,
        max_movable_mm_displacement_angstrom=max(displacements[i] for i in moving),
    )
