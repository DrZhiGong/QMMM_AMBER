"""Check whether the relaxation gradient belongs to the final single point.
Compare the QM coordinates, MM point charges and ORCA energy before testing convergence.
"""

from bisect import bisect_left, bisect_right
from collections import deque
import math
from pathlib import Path

from common import read, sha, xyz
from energy_precision import energy_matches


# Sander's ASCII restart stores coordinates to 7 decimal places in angstrom.
# Allow for rounding of both geometries and projection of link H atoms.
COORD_TOL = 3e-7
CHARGE_TOL = 1e-10
# A freshly converged SCF at a seven-decimal-place restart can differ slightly.
# 1 microhartree is 0.00063 kcal/mol; a larger change cannot certify the old gradient.
ENERGY_TOL_HARTREE = 1e-6


def _check_event(event, result, expected_spin):
    if (
        Path(result["event"]).resolve() != event
        or int(event.name) != result["callback_index"]
        or read(event / "result.json") != result
    ):
        raise ValueError(f"Inconsistent ORCA callback identity at {event}")
    multiplicity, alpha, beta = expected_spin
    if not math.isfinite(result["energy_hartree"]):
        raise ValueError(f"Nonfinite embedded ORCA energy at {event}")
    if (
        result.get("converged") is not True
        or result.get("coordinate_multiplicity") != multiplicity
        or [result.get("nalpha"), result.get("nbeta")] != [alpha, beta]
    ):
        raise ValueError(f"Expected a converged embedded {alpha}/{beta} checkpoint at {event}")
    for name, key in (
        ("native_qm.xyz", "native_qm_sha256"),
        ("embedding.pc", "pointcharge_sha256"),
    ):
        if sha(event / name) != result[key]:
            raise ValueError(f"ORCA geometry evidence differs from its archive: {event / name}")
    if Path(result["pointcharges_file"]).resolve() != event / "embedding.pc":
        raise ValueError(f"Unexpected MM charge file in {event}")
    native = xyz(event / "native_qm.xyz")
    if _native_distance(native, result["native_qm_atoms"]) > 2e-10:
        raise ValueError(f"Archived native QM geometry disagrees with result.json at {event}")
    return native


def _native_distance(first, second):
    if not first or len(first) != len(second):
        return math.inf
    largest = 0.0
    for (e1, p1), (e2, p2) in zip(first, second):
        if e1 != e2 or len(p1) != 3 or len(p2) != 3:
            return math.inf
        distance = math.dist(p1, p2)
        if not math.isfinite(distance):
            return math.inf
        largest = max(largest, distance)
    return largest


def _points(event, result):
    lines = (
        line.split() for line in (event / "embedding.pc").read_text().splitlines() if line.strip()
    )
    (count,) = next(lines)
    points = [tuple(map(float, row)) for row in lines]
    if (
        not points
        or int(count) != len(points)
        or len(points) != result["pointcharge_count"]
        or any(len(p) != 4 or not all(map(math.isfinite, p)) for p in points)
    ):
        raise ValueError(f"Invalid archived MM point charges at {event}")
    return points


def _compatible(a, b):
    return abs(a[0] - b[0]) <= CHARGE_TOL and math.dist(a[1:], b[1:]) <= COORD_TOL


def _point_assignment(source, target):
    if len(source) != len(target):
        return None
    if all(_compatible(a, b) for a, b in zip(source, target)):
        return range(len(source))
    # Sander can reorder the charges. Use x to find nearby candidates, then compare q and xyz.
    order = sorted(range(len(target)), key=lambda i: target[i][1])
    x_ordered = [target[i][1] for i in order]
    neighbors = []
    for a in source:
        lo = bisect_left(x_ordered, a[1] - COORD_TOL)
        hi = bisect_right(x_ordered, a[1] + COORD_TOL)
        choices = [i for i in order[lo:hi] if _compatible(a, target[i])]
        if not choices:
            return None
        neighbors.append(choices)
    # If two charges share a candidate, look for a path to a free target
    # and reassign the pairs along that path. Each charge gets exactly one partner.
    source_for_target = {}
    target_for_source = {}
    for start in sorted(range(len(source)), key=lambda i: len(neighbors[i])):
        queue = deque([start])
        parent = {}
        visited_sources = {start}
        free_target = None
        while queue and free_target is None:
            source_index = queue.popleft()
            for target_index in neighbors[source_index]:
                if target_index in parent:
                    continue
                parent[target_index] = source_index
                if target_index not in source_for_target:
                    free_target = target_index
                    break
                owner = source_for_target[target_index]
                if owner not in visited_sources:
                    visited_sources.add(owner)
                    queue.append(owner)
        if free_target is None:
            return None
        target_index = free_target
        while target_index is not None:
            source_index = parent[target_index]
            old_target = target_for_source.get(source_index)
            source_for_target[target_index] = source_index
            target_for_source[source_index] = target_index
            target_index = old_target
    return [target_for_source[i] for i in range(len(source))]


def match_final_geometry(
    relax_folder, fresh_bridge, candidate_indices, *, require_match=True, expected_spin
):
    """Match energy-compatible callbacks to native geometry generated from the actual restart."""
    relax = Path(relax_folder).resolve(strict=True)
    fresh = Path(fresh_bridge["event"]).resolve(strict=True)
    fresh_native = _check_event(fresh, fresh_bridge, expected_spin)
    fresh_points = _points(fresh, fresh_bridge)
    matches, observations = [], []
    for index in sorted(set(candidate_indices)):
        event = relax / "evaluations" / f"{index:05d}"
        old = read(event / "result.json")
        native = _check_event(event, old, expected_spin)
        shift = _native_distance(native, fresh_native)
        row = dict(callback_index=index, matched=False)
        if math.isfinite(shift):
            row["max_native_displacement_angstrom"] = shift
        if (
            shift <= COORD_TOL
            and old["native_real_qm_indices"] == fresh_bridge["native_real_qm_indices"]
        ):
            points = _points(event, old)
            assignment = _point_assignment(points, fresh_points)
            if assignment is not None:
                delta = abs(old["energy_hartree"] - fresh_bridge["energy_hartree"])
                row.update(
                    embedded_energy_difference_hartree=delta,
                    max_charge_displacement_angstrom=max(
                        math.dist(p[1:], fresh_points[j][1:]) for p, j in zip(points, assignment)
                    ),
                    max_charge_difference_e=max(
                        abs(p[0] - fresh_points[j][0]) for p, j in zip(points, assignment)
                    ),
                )
                if delta <= ENERGY_TOL_HARTREE:
                    row["matched"] = True
                    matches.append(index)
                else:
                    row["reason"] = "Embedded ORCA energy changed at the same geometry and embedding"
        observations.append(row)
    if not matches and require_match:
        raise ValueError(
            "No source callback matches the final QM/link atoms, MM point charges "
            f"and embedded ORCA energy: {observations}"
        )
    return dict(
        matching_callback_indices=matches,
        candidates=observations,
        native_atoms=len(fresh_native),
        point_charges=len(fresh_points),
        coordinate_tolerance_angstrom=COORD_TOL,
        charge_tolerance_e=CHARGE_TOL,
        embedded_energy_tolerance_hartree=ENERGY_TOL_HARTREE,
        fresh_event=str(fresh),
        scope="Ordered native QM/link geometry and entire native embedding at the Sander input restart",
    )


def match_endpoint(folder, relaxation, final, *, expected_spin):
    final_qm = final["bridge"]
    # At fixed geometry, the first and last QM calls must see the same embedding.
    first = match_final_geometry(Path(final["folder"]), final_qm, [1], require_match=False,
                                 expected_spin=expected_spin)
    if not first["matching_callback_indices"]:
        if "embedded_energy_difference_hartree" in first["candidates"][0]:
            return dict(matching_callback_indices=[], final_single_point_check=first)
        raise ValueError("The final single point changed geometry or embedding between callbacks")
    seen = {}

    def matches(index):
        value = match_final_geometry(folder, final_qm, [index], require_match=False, expected_spin=expected_spin)
        seen[index] = value["candidates"][0]
        return value if value["matching_callback_indices"] else None

    best = matches(relaxation["energy_audit"]["matched_callback_index"])
    if best is not None:
        return best
    content = (folder / "sander.out").read_text(errors="replace")
    for index, _ in energy_matches(content, relaxation["bridge"], all_events=True):
        if index not in seen:
            best = matches(index)
            if best is not None:
                return best
    return dict(matching_callback_indices=[], energy_compatible_callback_indices=sorted(seen),
                candidates=list(seen.values()), embedded_energy_tolerance_hartree=ENERGY_TOL_HARTREE)


def classify_mm_gradient(rms, target):
    if type(target) not in (int, float) or not math.isfinite(target) or target <= 0:
        raise ValueError("Invalid MM gradient target")
    if rms is None:
        return None, "FINAL_SINGLE_POINT_MM_UNVERIFIED"
    if type(rms) not in (int, float) or not math.isfinite(rms) or rms < 0:
        raise ValueError("Invalid MM RMS gradient")
    if rms <= target:
        return True, "RELAXED_UNREVIEWED"
    return False, "MM_NOT_CONVERGED_UNREVIEWED"


def relaxation_status(folder, relaxation, final, target, *, expected_spin):
    """Attach the RMS gradient only if a relaxation callback matches the final single point."""
    audit = match_endpoint(folder, relaxation, final, expected_spin=expected_spin)
    rms = relaxation["sander_rms_gradient"] if audit["matching_callback_indices"] else None
    converged, status = classify_mm_gradient(rms, target)
    return dict(status=status, mm_geometry_converged=converged,
                relaxation_rms_gradient=rms, endpoint_audit=audit)
