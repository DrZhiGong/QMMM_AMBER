"""Read your system.json and check the topology, coordinates and atom selections.
Use the printed summary to see which atoms are in the QM region and which MM atoms can move.
"""

import json
import math
from pathlib import Path
import re

from common import CHARGE_SCALE, Prmtop, mask, read
from moving_selection import LIMIT, build_moving_selection


def _selection(spec, field, natom):
    value = spec[field]
    selected = set()
    for part in value.split(","):
        match = re.fullmatch(r"\s*([0-9]+)\s*(?:-\s*([0-9]+)\s*)?", part)
        if not match:
            raise ValueError(f"{field}: invalid range {part!r}. Example: '1-3,5'.")
        first, last = int(match[1]), int(match[2] or match[1])
        if not 1 <= first <= last <= natom:
            raise ValueError(f"{field}: range {part!r} must lie between 1 and {natom}.")
        selected.update(range(first - 1, last))
    return sorted(selected)


def _validate_coordinates(path, natom):
    try:
        lines = path.read_text(encoding="ascii").splitlines()
    except UnicodeDecodeError as exc:
        raise ValueError("Export coordinates as an ASCII AMBER restart (inpcrd/rst7).") from exc
    while lines and not lines[-1].strip():
        lines.pop()
    header = lines[1].split() if len(lines) > 1 else []
    if len(header) not in (1, 2) or not header[0].isdigit() or int(header[0]) != natom:
        raise ValueError(f"Expected an ASCII AMBER restart with {natom} atoms on line 2.")
    if len(header) == 2 and not math.isfinite(float(header[1])):
        raise ValueError("Invalid restart time on line 2.")
    if len(lines) != 2 + (natom + 1) // 2:
        raise ValueError("Supply one coordinate-only ASCII AMBER restart without box or velocity records.")
    for i, line in enumerate(lines[2:]):
        fields = min(6, 3 * natom - 6 * i)
        if len(line) < fields * 12 or line[fields * 12:].strip():
            raise ValueError(f"Line {i + 3} must contain {fields} AMBER F12 coordinate fields.")
        for j in range(0, fields * 12, 12):
            try:
                value = float(line[j:j + 12])
            except ValueError as exc:
                raise ValueError(f"Invalid fixed-width coordinate on line {i + 3}.") from exc
            if not math.isfinite(value):
                raise ValueError(f"Nonfinite coordinate on line {i + 3}.")


def _validate_topology(parm):
    pointers = parm.get("POINTERS")
    if len(pointers) < 31 or any(type(p) is not int or p < 0 for p in pointers):
        raise ValueError("Invalid AMBER POINTERS section.")
    natom = pointers[0]
    if natom < 2:
        raise ValueError("The topology needs at least one QM and one MM atom.")
    if pointers[27] or "BOX_DIMENSIONS" in parm.sections:
        raise ValueError("Periodic topology detected. Supply a nonperiodic model.")
    if (pointers[20] or pointers[29] or pointers[30] or "CTITLE" in parm.sections
            or any(k.startswith(("CHARMM_", "AMOEBA_")) for k in parm.sections)
            or ("IPOL" in parm.sections and parm.get("IPOL") != [0])
            or ("POLARIZABILITY" in parm.sections and "IPOL" not in parm.sections)):
        raise ValueError("Supply a standard fixed-charge AMBER topology without CHAMBER, "
                         "polarization, extra points, solvent caps or perturbations.")
    numbers = parm.get("ATOMIC_NUMBER")
    names = parm.get("ATOM_NAME")
    charges = parm.get("CHARGE")
    if (len(names) != natom or not all(names) or len(numbers) != natom
            or any(type(z) is not int or not 1 <= z <= 118 for z in numbers)):
        raise ValueError("ATOM_NAME and ATOMIC_NUMBER must identify every atom.")
    if len(charges) != natom or not all(math.isfinite(q) for q in charges):
        raise ValueError("CHARGE must contain one finite charge per atom.")
    masses = parm.get("MASS")
    if len(masses) != natom or any(not math.isfinite(m) or m <= 0 for m in masses):
        raise ValueError("MASS must contain one finite positive mass per atom.")
    _validate_nonbonded(parm, pointers)
    starts = parm.get("RESIDUE_POINTER")
    if (not starts or len(starts) != pointers[11] or any(type(i) is not int for i in starts)
            or starts[0] != 1 or starts[-1] > natom or starts != sorted(set(starts))):
        raise ValueError("Invalid AMBER residue pointers or residue count.")
    ntypes = pointers[15]
    force = parm.get("BOND_FORCE_CONSTANT")
    length = parm.get("BOND_EQUIL_VALUE")
    if (len(force) != ntypes or len(length) != ntypes
            or any(not math.isfinite(k) or k < 0 for k in force)
            or any(not math.isfinite(r) or r <= 0 for r in length)):
        raise ValueError("Invalid bond parameter arrays in topology.")
    for name, count in (("BONDS_INC_HYDROGEN", pointers[2]), ("BONDS_WITHOUT_HYDROGEN", pointers[12])):
        records = parm.get(name)
        if len(records) != count * 3 or any(type(x) is not int for x in records):
            raise ValueError(f"Invalid bond records in {name}.")
        for a, b, kind in zip(records[::3], records[1::3], records[2::3]):
            if (a == b or min(a, b) < 0 or max(a, b) >= 3 * natom or a % 3 or b % 3
                    or not 1 <= kind <= ntypes):
                raise ValueError(f"Invalid atom offset or parameter index in {name}.")
    if len(parm.bonds()) != pointers[2] + pointers[12]:
        raise ValueError("Duplicate bonds in topology. Check the source topology.")


def _validate_nonbonded(parm, pointers):
    """Check the dimensions and references used for MM nonbonded interactions."""
    natom, ntypes, nexcluded, nhbond = pointers[0], pointers[1], pointers[10], pointers[19]
    nlj = ntypes * (ntypes + 1) // 2
    types = parm.get("ATOM_TYPE_INDEX")
    if len(types) != natom or any(type(i) is not int or not 1 <= i <= ntypes for i in types):
        raise ValueError("Invalid ATOM_TYPE_INDEX or nonbonded atom type count.")
    pairs = parm.get("NONBONDED_PARM_INDEX")
    if (len(pairs) != ntypes * ntypes
            or any(type(i) is not int or i < -nhbond or i > nlj for i in pairs)):
        raise ValueError("Invalid NONBONDED_PARM_INDEX table or coefficient index.")
    for name, size in (("LENNARD_JONES_ACOEF", nlj), ("LENNARD_JONES_BCOEF", nlj),
                       ("HBOND_ACOEF", nhbond), ("HBOND_BCOEF", nhbond)):
        if size or name in parm.sections:
            coeffs = parm.get(name)
            if len(coeffs) != size or not all(math.isfinite(c) for c in coeffs):
                raise ValueError(f"{name} must contain {size} finite coefficients.")
    counts = parm.get("NUMBER_EXCLUDED_ATOMS")
    if (len(counts) != natom or any(type(i) is not int or i < 0 for i in counts)
            or sum(counts) != nexcluded):
        raise ValueError("NUMBER_EXCLUDED_ATOMS does not match the topology exclusion count.")
    excluded = parm.get("EXCLUDED_ATOMS_LIST")
    if len(excluded) != nexcluded or any(type(i) is not int or not 0 <= i <= natom for i in excluded):
        raise ValueError("Invalid atom index or count in EXCLUDED_ATOMS_LIST.")


def read_system_spec(path):
    """Read the six-field JSON and resolve both file names against its directory."""
    path = Path(path).resolve()
    try:
        spec = read(path)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{path}, line {exc.lineno}, column {exc.colno}: {exc.msg}") from exc
    required = {"topology", "coordinates", "qm_atoms", "movable_mm_atoms", "qm_charge", "multiplicity"}
    if not isinstance(spec, dict):
        raise ValueError(f"{path}: expected a JSON object.")
    if set(spec) != required:
        missing = ", ".join(sorted(required - spec.keys())) or "none"
        unknown = ", ".join(sorted(spec.keys() - required)) or "none"
        raise ValueError(f"{path}: missing fields: {missing}. Unknown fields: {unknown}.")
    for field in ("qm_atoms", "movable_mm_atoms"):
        if not isinstance(spec[field], str) or not spec[field].strip():
            raise ValueError(f"{field}: give 1-based atom numbers or ranges, for example '1-3,5'.")
    if type(spec["qm_charge"]) is not int or type(spec["multiplicity"]) is not int or spec["multiplicity"] < 1:
        raise ValueError("qm_charge must be an integer and multiplicity must be a positive integer.")
    for key in ("topology", "coordinates"):
        if not isinstance(spec[key], str) or not spec[key].strip():
            raise ValueError(f"{key}: supply an AMBER file path.")
        file = (path.parent / spec[key]).resolve()
        if not file.is_file():
            raise ValueError(f"{key}: file does not exist: {file}")
        spec[key] = str(file)
    return spec


def _short_mask(parm, indices):
    numeric = mask(indices)
    if len(numeric) <= LIMIT:
        return numeric
    names = parm.get("ATOM_NAME")
    chosen = {names[i] for i in indices}
    named_mask = "@" + ",".join(sorted(chosen))
    if (len(named_mask) <= LIMIT and all(re.fullmatch(r"[A-Za-z0-9_+]+", n) and not n.isdigit() for n in chosen)
            and sum(name in chosen for name in names) == len(indices)):
        return named_mask
    return build_moving_selection(parm, indices)["mask"]


def _amber_mask(parm, indices, field):
    # Try both equivalent descriptions. All atoms except a short exclusion list
    # can be much more compact than enumerating the selected atoms.
    try:
        return _short_mask(parm, indices)
    except ValueError as selection_error:
        excluded = set(range(parm.get("POINTERS")[0])).difference(indices)
        try:
            excluded_mask = _short_mask(parm, excluded)
            if re.fullmatch(r"[@:][A-Za-z0-9_+,\-]+", excluded_mask):
                inverted = "!" + excluded_mask
            else:
                inverted = f"!({excluded_mask})"
            if len(inverted) <= LIMIT:
                return inverted
        except ValueError:
            pass
        raise ValueError(f"{field}: cannot compact the selection into an exact AMBER mask "
                         f"of at most {LIMIT} characters. {selection_error}") from selection_error


def load_system(path):
    """Validate the topology and coordinates. All returned atom indices are zero-based."""
    spec = read_system_spec(path)
    try:
        parm = Prmtop(spec["topology"])
        _validate_topology(parm)
    except KeyError as exc:
        raise ValueError(f"{spec['topology']}: missing %FLAG {exc.args[0]}. Supply a modern AMBER prmtop.") from exc
    except (IndexError, TypeError, ValueError) as exc:
        raise ValueError(f"{spec['topology']}: invalid or unsupported AMBER topology: {exc}") from exc
    natom = parm.get("POINTERS")[0]
    try:
        _validate_coordinates(Path(spec["coordinates"]), natom)
    except ValueError as exc:
        raise ValueError(f"{spec['coordinates']}: {exc}") from exc
    qm = _selection(spec, "qm_atoms", natom)
    moving = _selection(spec, "movable_mm_atoms", natom)
    qm_set = set(qm)
    if qm_set.intersection(moving):
        raise ValueError("qm_atoms and movable_mm_atoms overlap. Only MM atoms can move.")
    numbers = parm.get("ATOMIC_NUMBER")
    if any(numbers[i] > 36 for i in qm):
        raise ValueError("QM atoms must be H through Kr for the all-electron def2-SVP calculation.")
    cuts = sorted((a, b) for a, b in parm.bonds() if (a in qm_set) != (b in qm_set))
    mm_anchors = set()
    for a, b in cuts:
        if numbers[a] != 6 or numbers[b] != 6:
            raise ValueError(f"QM/MM boundary bond {a + 1}-{b + 1}: only C-C cuts with link H are supported.")
        anchor = b if a in qm_set else a
        if anchor in mm_anchors:
            raise ValueError(f"MM atom {anchor + 1} connects to multiple QM atoms. "
                             "Choose a boundary with one link per MM anchor.")
        mm_anchors.add(anchor)
    # AMBER 22 Reference Manual §§10.1.5, 10.1.7: one link per MM anchor;
    # adjust_q=2 requires NATOM - NQUANT - NLINK >= NLINK.
    remaining_mm = natom - len(qm) - len(cuts)
    if remaining_mm < len(cuts):
        raise ValueError(f"adjust_q=2 needs at least {len(cuts)} MM atoms outside the link pairs; "
                         f"the selection leaves {remaining_mm}.")
    charge, mult = spec["qm_charge"], spec["multiplicity"]
    electrons = sum(numbers[i] for i in qm) + len(cuts) - charge
    if electrons <= 0 or electrons < mult - 1 or (electrons + mult - 1) % 2:
        raise ValueError(f"QM charge {charge} and multiplicity {mult} are inconsistent with "
                         f"{electrons} electrons including {len(cuts)} link H.")
    total_charge = math.fsum(q / CHARGE_SCALE for q in parm.get("CHARGE"))
    alpha = (electrons + mult - 1) // 2
    return dict(topology=spec["topology"], coordinates=spec["coordinates"],
                natom=natom, atomic_numbers=numbers, qm_indices_0based=qm,
                qm_mask=_amber_mask(parm, qm, "qm_atoms"),
                movable_mm_indices_0based=moving, movable_mm_mask=_amber_mask(parm, moving, "movable_mm_atoms"),
                link_bonds_0based=[list(pair) for pair in cuts], qm_charge=charge, multiplicity=mult,
                expected_qm_atoms_with_links=len(qm) + len(cuts), nalpha=alpha, nbeta=electrons - alpha,
                source_full_partial_charge_sum_e=total_charge, expected_mm_charge_sum_e=total_charge - charge)


def describe_system(model, stream=None):
    """Print the checked system and both AMBER selections using one-based indices."""
    print(f"Topology:    {model['topology']}\nCoordinates: {model['coordinates']}", file=stream)
    print(f"{model['natom']} atoms; {len(model['qm_indices_0based'])} real QM atoms; "
          f"{len(model['link_bonds_0based'])} link H; "
          f"{len(model['movable_mm_indices_0based'])} movable MM atoms.", file=stream)
    print(f"QM atoms:    {model['qm_mask']}\nMovable MM:  {model['movable_mm_mask']}", file=stream)
    print(f"QM charge {model['qm_charge']}, multiplicity {model['multiplicity']}, "
          f"{model['nalpha']} alpha / {model['nbeta']} beta electrons.", file=stream)
