"""When sander calls ORCA, pass along the QM geometry and surrounding MM point charges.
Check the energy and gradients and keep the orbitals for the next call.
"""

import collections
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from common import Z, read, write, text, sha, xyz, write_xyz
from inputs_control import orca_template

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
# ORCA 6.1 manual, §10.3. The .engrad coordinates are printed in bohr.
ORCA_ANGSTROM_PER_BOHR = 0.5291772083


def final_number(s, pattern, label):
    found = list(re.finditer(pattern, s, re.M | re.I))
    if not found:
        raise ValueError("ORCA missing " + label)
    value = float(found[-1].group(1).replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError("Nonfinite " + label)
    return value


def _spin_electron_counts(final_scf):
    """Sum the occupations in the last pair of spin-up and spin-down tables."""
    heading = r"^\s*SPIN\s+(UP|DOWN)\s+ORBITALS\s*$"
    orbital_row = re.compile(rf"^\s*(\d+)\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})(?:\s+\S+)?\s*$")
    tables = list(re.finditer(heading, final_scf, re.M | re.I))[-2:]
    if [table[1].upper() for table in tables] != ["UP", "DOWN"]:
        raise ValueError("Expected final spin-up and spin-down orbital tables, in that order")
    counts = []
    for table, end in zip(tables, (tables[1].start(), len(final_scf))):
        occupations = []
        orbital_indices = []
        for line in final_scf[table.end():end].splitlines():
            match = orbital_row.match(line)
            if match:
                orbital_indices.append(int(match[1]))
                occupations.append(float(match[2]))
            elif occupations:
                break
        if not occupations or occupations[-1] != 0 or orbital_indices != list(range(len(occupations))):
            raise ValueError("Missing/incomplete final orbital occupations")
        if any(abs(occ - round(occ)) > 1e-5 or round(occ) not in (0, 1) for occ in occupations):
            raise ValueError("Unexpected fractional/restricted occupations")
        counts.append(round(sum(occupations)))
    return counts


def _mulliken_populations(final_scf, atoms, charge, spin):
    heading = r"^\s*MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS\s*$"
    tables = list(re.finditer(heading, final_scf, re.M))
    if not tables:
        raise ValueError("Final Mulliken population table missing")
    atom_row = re.compile(rf"^\s*(\d+)\s+([A-Za-z]+)\s*:\s*({NUMBER})\s+({NUMBER})\s*$")
    populations = []
    for line in final_scf[tables[-1].end():].splitlines():
        match = atom_row.match(line)
        if match:
            populations.append({
                "index": int(match[1]),
                "element": match[2],
                "charge": float(match[3]),
                "spin": float(match[4]),
            })
        elif populations:
            break
    if len(populations) != len(atoms) or any(
        row["index"] != i or row["element"].lower() != atoms[i][0].lower()
        for i, row in enumerate(populations)
    ):
        raise ValueError("Population atom count/element order differs")
    total_charge = sum(row["charge"] for row in populations)
    total_spin = sum(row["spin"] for row in populations)
    if abs(total_charge - charge) > 0.003 or abs(total_spin - spin) > 0.003:
        raise ValueError("Final Mulliken charge/spin sum differs")
    return populations


def parse_output(content, atoms, charge, mult):
    # ORCA echoes the input near the top of its output.
    output = "\n".join(line for line in content.splitlines() if not re.match(r"^\s*\|\s*\d+>", line))
    failure = (
        r"error termination|WAVEFUNCTION IS NOT YET CONVERGED|"
        r"SCF NOT CONVERGED|SCF DID NOT CONVERGE"
    )
    if output.count("ORCA TERMINATED NORMALLY") != 1 or re.search(failure, output, re.I):
        raise ValueError("ORCA did not terminate with a converged wavefunction")
    converged = list(re.finditer(r"^.*\bSCF\s+CONVERGED\b.*$", output, re.M))
    if not converged:
        raise ValueError("Missing SCF convergence marker")

    final_scf = output[converged[-1].end():]
    energy = final_number(
        final_scf, rf"^\s*FINAL SINGLE POINT ENERGY\s+({NUMBER})\s*$", "total energy"
    )
    dispersion = final_number(
        final_scf, rf"^\s*Dispersion correction\s*[:=]?\s*({NUMBER})\s*$", "dispersion"
    )
    s2 = final_number(
        final_scf, rf"^\s*Expectation value of\s*<S\*\*2>\s*[:=]?\s*({NUMBER})\s*$", "S squared"
    )
    if s2 < -1e-6:
        raise ValueError("Unphysical negative S squared")

    electrons = sum(Z[element] for element, _ in atoms) - charge
    if (electrons + mult - 1) % 2:
        raise ValueError("Invalid electron/spin parity")
    alpha = (electrons + mult - 1) // 2
    beta = electrons - alpha
    counted = _spin_electron_counts(final_scf)
    if counted != [alpha, beta]:
        raise ValueError(f"Final occupations differ: {counted} expected {[alpha, beta]}")
    populations = _mulliken_populations(final_scf, atoms, charge, alpha - beta)

    removed = re.findall(r"^\s*Number of eigenvalues below threshold\s+\.{2,}\s+(\d+)\s*$", output, re.M)
    if not removed or any(int(count) != 0 for count in removed):
        raise ValueError("Overlap directions removed or overlap evidence missing")
    return dict(
        converged=True,
        energy_hartree=energy,
        e_dispersion_hartree=dispersion,
        e_scf_without_dispersion_hartree=energy - dispersion,
        s_squared=s2,
        nalpha=alpha,
        nbeta=beta,
        atom_mulliken_populations=populations,
        state_selection="UNREVIEWED",
        stability_test_performed=False,
    )


def generated_input(path):
    """Read the QM coordinates and the point-charge file named by Sander."""
    path = Path(path).resolve()
    contents = path.read_text()
    point_charge_files = re.findall(r"^\s*%pointcharges\s+(.+?)\s*$", contents, re.M | re.I)
    if len(point_charge_files) != 1:
        raise ValueError("Expected one AMBER-generated pointcharge file")
    filenames = shlex.split(point_charge_files[0])
    if len(filenames) != 1:
        raise ValueError("Ambiguous pointcharge path")
    point_charge_file = (path.parent / filenames[0]).resolve(strict=True)

    external_xyz = re.search(
        r"^\s*\*\s*xyzfile\s+(-?\d+)\s+(\d+)\s+(.+?)\s*$", contents, re.M | re.I
    )
    if external_xyz:
        filenames = shlex.split(external_xyz[3])
        if len(filenames) != 1:
            raise ValueError("Ambiguous AMBER XYZ path")
        atoms = xyz(path.parent / filenames[0])
        charge, multiplicity = int(external_xyz[1]), int(external_xyz[2])
    else:
        inline_xyz = re.search(
            r"^\s*\*\s*xyz\s+(-?\d+)\s+(\d+)\s*\n(.*?)^\s*\*\s*$", contents, re.M | re.S | re.I
        )
        if not inline_xyz:
            raise ValueError("Unsupported AMBER ORCA coordinate block")
        charge, multiplicity = int(inline_xyz[1]), int(inline_xyz[2])
        atoms = []
        for line in inline_xyz[3].splitlines():
            fields = line.split()
            if len(fields) != 4 or fields[0] not in Z:
                raise ValueError("Unexpected generated atom")
            atoms.append((fields[0], list(map(float, fields[1:]))))

    lines = point_charge_file.read_text().splitlines()
    count = int(lines[0])
    rows = [list(map(float, line.split())) for line in lines[1:] if line.strip()]
    if len(rows) != count or any(len(row) != 4 or not all(map(math.isfinite, row)) for row in rows):
        raise ValueError("Malformed point-charge file")
    return atoms, charge, multiplicity, point_charge_file, count, math.fsum(row[0] for row in rows)


def native_atom_order(atoms, context):
    """Map the fixed QM atoms to Sander's order, then identify the parent of each link H."""
    fixed_atoms = context["real_qm_atoms"]
    expected = context["expected_qm_atoms_with_links"]
    if len(atoms) != expected:
        raise ValueError("Unexpected native QM/link count: " + str(len(atoms)))
    real_qm_indices = []
    matched_atoms = set()
    for source in fixed_atoms:
        matches = [
            i
            for i, (element, coord) in enumerate(atoms)
            if element == source["element"] and math.dist(coord, source["coord"]) < 2e-6
        ]
        if len(matches) != 1 or matches[0] in matched_atoms:
            raise ValueError("Fixed QM geometry/order could not be verified")
        real_qm_indices.append(matches[0])
        matched_atoms.add(matches[0])

    parents = []
    for i, (element, coord) in enumerate(atoms):
        if i in matched_atoms:
            continue
        if element != "H":
            raise ValueError("Unexpected non-H link atom")
        nearby_carbons = [
            k
            for k, atom in enumerate(fixed_atoms)
            if atom["element"] == "C" and abs(math.dist(coord, atom["coord"]) - 1.09) < 2e-5
        ]
        if len(nearby_carbons) != 1:
            raise ValueError("Cannot identify native link-H parent")
        parents.append(nearby_carbons[0])
    if collections.Counter(parents) != collections.Counter(context["cap_parent_real_indices"]):
        raise ValueError("AMBER native link assignments differ from the audited cuts")
    return real_qm_indices, parents


def check_gradient(path, atoms, energy, npc):
    stem = Path(path).with_suffix("")
    values = []
    for line in stem.with_suffix(".engrad").read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            values.extend(map(float, line.split()))
    n = len(atoms)
    if (
        len(values) != 2 + 7 * n
        or values[0] != n
        or abs(values[1] - energy) > 1e-7
        or not all(map(math.isfinite, values))
    ):
        raise ValueError("ORCA energy/gradient file inconsistent with converged result")
    tail = values[2 + 3 * n :]
    for i, (element, position) in enumerate(atoms):
        number, *coords_bohr = tail[4 * i : 4 * i + 4]
        if number != Z[element]:
            raise ValueError("ORCA gradient atom order differs")
        # Allow rounding of printed coordinates and small differences in Bohr conversion constants.
        if not all(math.isclose(x * ORCA_ANGSTROM_PER_BOHR, ref, abs_tol=5e-7, rel_tol=1e-8)
                   for x, ref in zip(coords_bohr, position)):
            raise ValueError(f"ORCA gradient geometry differs from the input at atom {i + 1}")
    pcgrad = stem.with_suffix(".pcgrad").read_text().split()
    if (
        not pcgrad
        or int(pcgrad[0]) != npc
        or len(pcgrad) != 1 + 3 * npc
        or not all(map(math.isfinite, map(float, pcgrad[1:])))
    ):
        raise ValueError("ORCA external-charge gradients missing/incomplete")


def callback_index(history):
    """Number only complete, contiguous callback directories; never hide a gap."""
    entries = [p for p in history.iterdir() if p.is_dir()]
    names = sorted(p.name for p in entries)
    expected = [f"{i:05d}" for i in range(1, len(names) + 1)]
    if names != expected:
        raise ValueError("ORCA callback history is not contiguous")
    return len(names) + 1


def _write_orca_input(path, atoms, charge, multiplicity, point_charges, context, has_seed):
    template = orca_template(context, has_seed=has_seed)
    geometry = [f"* xyz {charge} {multiplicity}"]
    for element, coords in atoms:
        geometry.append(element + " " + " ".join(f"{v:.12f}" for v in coords))
    geometry.append("*")
    charge_directive = f'\n%pointcharges "{point_charges}"\n! EnGrad Angs NoUseSym\n'
    text(path, template + charge_directive + "\n".join(geometry) + "\n")


def _run_orca(command, output_file):
    # Sander reads ORCA's stdout. Keep the same bytes in the per-step log.
    with output_file.open("xb") as log:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while chunk := process.stdout.read(65536):
            log.write(chunk)
            log.flush()
            sys.stdout.buffer.write(chunk)
            sys.stdout.buffer.flush()
        return_code = process.wait()
    if return_code:
        raise RuntimeError(f"ORCA subprocess failed: {return_code}")


def run(input_file):
    """Evaluate one Sander step and save the orbitals for its next call."""
    context = read(os.environ["AMBER_ORCA_CONTEXT"])
    work = Path.cwd().resolve()
    path = Path(input_file).resolve()
    if path.parent != work or path.suffix.lower() != ".inp":
        raise ValueError("Unexpected ORCA driver input")
    atoms, charge, sander_multiplicity, point_charge_file, charge_count, charge_sum = generated_input(path)
    if charge != context["charge"] or sander_multiplicity != context["input_multiplicity"]:
        raise ValueError("AMBER-generated QM charge/multiplicity differs from configured target")
    real_qm_indices, link_parents = native_atom_order(atoms, context)
    expected_charge_sum = context.get("expected_mm_charge_sum_e")
    if expected_charge_sum is not None and abs(charge_sum - expected_charge_sum) > 0.02:
        raise ValueError(
            "Native embedding charge sum violates configured full-system charge conservation"
        )
    seed_atom_order = context.get("native_real_qm_indices")
    if seed_atom_order is not None and real_qm_indices != seed_atom_order:
        raise ValueError("Native AO atom order differs from the seeded calculation")
    if "native_link_parents" in context and link_parents != context["native_link_parents"]:
        raise ValueError("Native link-H order differs from the seeded calculation")

    history = work / "evaluations"
    history.mkdir(exist_ok=True)
    index = callback_index(history)
    seed = work / "orbital_seed.gbw"
    gbw = path.with_suffix(".gbw")
    if gbw == seed:
        raise ValueError("The ORCA output basename would overwrite the orbital seed")
    if index > 1:
        previous = read(history / f"{index - 1:05d}" / "result.json")
        if (previous != read(work / "bridge_last.json")
                or previous["native_real_qm_indices"] != real_qm_indices
                or previous["native_link_parents"] != link_parents
                or Path(previous["gbw_file"]) != gbw
                or sha(gbw) != previous["gbw_sha256"]):
            raise ValueError("Native GBW or QM/link order changed since the previous callback")
    elif gbw.exists():
        raise ValueError("Unexpected ORCA orbitals before the first callback")
    event = history / f"{index:05d}"
    event.mkdir()
    shutil.copyfile(path, event / "sander_generated.inp")
    embedding = event / "embedding.pc"
    shutil.copyfile(point_charge_file, embedding)
    embedding_sha = sha(embedding)
    write_xyz(event / "native_qm.xyz", atoms, "AMBER native ordered QM atoms and link H")
    if gbw.is_file():
        gbw.replace(seed)
    _write_orca_input(
        path, atoms, charge, sander_multiplicity, embedding.relative_to(work), context, seed.is_file()
    )
    shutil.copyfile(path, event / "executed.inp")
    command = [os.environ["ORCA_BIN"], str(path)]
    if context["nprocs"] > 1 and os.environ.get("ORCA_MPI_FLAGS"):
        command.append(os.environ["ORCA_MPI_FLAGS"])
    write(
        event / "command.json",
        dict(
            argv=command,
            callback_index=index,
            coordinate_multiplicity=sander_multiplicity,
            pointcharge_count=charge_count,
            pointcharge_sha256=embedding_sha,
            native_real_qm_indices=real_qm_indices,
            native_link_parents=link_parents,
        ),
    )
    # Sander reuses this basename. Require fresh output even at a repeated geometry.
    path.with_suffix(".engrad").unlink(missing_ok=True)
    path.with_suffix(".pcgrad").unlink(missing_ok=True)
    _run_orca(command, event / "orca.out")
    result = parse_output(
        (event / "orca.out").read_text(errors="replace"), atoms, charge, sander_multiplicity
    )
    check_gradient(path, atoms, result["energy_hartree"], charge_count)
    if not gbw.is_file() or gbw.stat().st_size == 0:
        raise ValueError("Missing native GBW")
    shutil.copyfile(path.with_suffix(".engrad"), event / "orca.engrad")
    shutil.copyfile(path.with_suffix(".pcgrad"), event / "orca.pcgrad")
    result.update(
        event=str(event),
        gbw_file=str(gbw),
        gbw_sha256=sha(gbw),
        executed_input_sha256=sha(path),
        pointcharges_file=str(embedding),
        pointcharge_sha256=embedding_sha,
        native_qm_atoms=atoms,
        native_real_qm_indices=real_qm_indices,
        native_link_parents=link_parents,
        callback_index=index,
        coordinate_multiplicity=sander_multiplicity,
        output_sha256=sha(event / "orca.out"),
        native_qm_sha256=sha(event / "native_qm.xyz"),
        gradient_file=str(event / "orca.engrad"),
        gradient_sha256=sha(event / "orca.engrad"),
        pointcharge_gradient_file=str(event / "orca.pcgrad"),
        pointcharge_gradient_sha256=sha(event / "orca.pcgrad"),
        pointcharge_count=charge_count,
        pointcharge_sum_e=charge_sum,
        expected_mm_charge_sum_e=expected_charge_sum,
    )
    write(event / "result.json", result)
    write(work / "bridge_last.json", result)


if __name__ == "__main__":
    try:
        if len(sys.argv) != 2:
            raise ValueError("Expected one generated ORCA input filename")
        run(sys.argv[1])
    except BaseException as exc:
        print("AMBER_ORCA_BRIDGE_FAILED: " + str(exc), file=sys.stderr, flush=True)
        raise
