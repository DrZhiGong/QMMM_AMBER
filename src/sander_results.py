"""Add the final MM energy terms in sander.out to the embedded QM energy from ORCA.
Check that the sum agrees with the total that sander printed.
"""

from decimal import Decimal, InvalidOperation
import math
import re

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
MM_TERMS = ("BOND", "ANGLE", "DIHED", "VDWAALS", "EEL", "HBOND", "1-4 VDW", "1-4 EEL", "RESTRAINT")
QM_TERM = "EXTERNESCF"


class ExternalEnergyMismatch(ValueError):
    """The reported EXTERNESCF and the ORCA energy differ beyond printing precision."""


def _number(token, label):
    if not re.fullmatch(NUMBER, token):
        raise ValueError("Nonfinite/malformed Sander " + label + ": " + token)
    value = float(token.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError("Nonfinite Sander " + label)
    return value


def _rounding_bound(token):
    """Half the last printed decimal place, including a scientific exponent."""
    try:
        literal = Decimal(token.replace("D", "E").replace("d", "e"))
    except InvalidOperation as exc:
        raise ValueError("Malformed Sander numeric precision") from exc
    if not literal.is_finite():
        raise ValueError("Nonfinite Sander numeric precision")
    return float(Decimal(5).scaleb(literal.as_tuple().exponent - 1))


def _final_table(content):
    sections = re.split(r"^\s*FINAL RESULTS\s*$", content, flags=re.M | re.I)
    if len(sections) != 2:
        raise ValueError("Missing/ambiguous final Sander minimization report")
    section = sections[1]
    # Timing data after the table can also contain "=" signs.
    section = re.split(
        r"^\s*(?:\d+\.\s+TIMINGS|\|?\s*>{3,}PROFILE.*)\s*$", section, maxsplit=1, flags=re.M | re.I
    )[0]
    headers = list(
        re.finditer(r"^\s*NSTEP\s+ENERGY\s+RMS\s+GMAX\s+NAME\s+NUMBER\s*$", section, re.M)
    )
    if len(headers) != 1:
        raise ValueError("Missing/ambiguous final Sander NSTEP table")
    tail = section[headers[0].end() :]
    lines = [line.strip() for line in tail.splitlines() if line.strip()]
    if not lines:
        raise ValueError("Missing final Sander NSTEP values")
    row = re.fullmatch(r"(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)", lines[0])
    if not row:
        raise ValueError("Malformed final Sander NSTEP values")
    return row, lines[1:]


def _energy_components(lines):
    """Keep each printed number as well as its value for the rounding check."""
    terms = {}
    literals = {}
    assignment = re.compile(r"([A-Za-z0-9][A-Za-z0-9 +*/()_-]*?)\s*=\s*(\S+)")
    for line in lines:
        if line.startswith("|"):
            continue
        for match in assignment.finditer(line):
            label = " ".join(match[1].upper().split())
            if label not in MM_TERMS + (QM_TERM, "EAMBER"):
                raise ValueError("Unsupported Sander final energy term: " + label)
            if label in terms:
                raise ValueError("Duplicate Sander final energy term: " + label)
            terms[label] = _number(match[2], label)
            literals[label] = match[2]
    missing = [name for name in MM_TERMS + (QM_TERM,) if name not in terms]
    if missing:
        raise ValueError("Missing Sander final energy terms: " + ", ".join(missing))
    return terms, literals


def parse_sander_results(content, bridge):
    """Read the final energy and MM gradient using the printed conversion factor."""
    return reconstruct_potential(content, bridge)


def reconstruct_potential(content, bridge, *, conversion_literal=None):
    """A factor supplied by the caller uses AMBER's full precision for the energy check."""
    if bridge.get("converged") is not True:
        raise ValueError("Sander energy reconstruction requires a checked ORCA result")
    qm_hartree = float(bridge["energy_hartree"])
    if not math.isfinite(qm_hartree):
        raise ValueError("Nonfinite checked ORCA energy")
    conversions = re.findall(r"^\s*AU_TO_KCAL\s*=\s*(\S+)\s*$", content, re.M | re.I)
    if len(conversions) != 1:
        raise ValueError("Missing/ambiguous native AMBER AU_TO_KCAL conversion")
    if conversion_literal is None:
        conversion_literal = conversions[0]
    conversion = _number(conversion_literal, "AU_TO_KCAL")
    if abs(conversion - 627.50946943) > 1e-4:
        raise ValueError("Unexpected native AMBER Hartree-to-kcal/mol conversion")

    row, component_lines = _final_table(content)
    nstep = int(row[1])
    rounded = _number(row[2], "rounded table energy")
    rms = _number(row[3], "RMS gradient")
    gmax = _number(row[4], "maximum gradient")
    if rms < 0 or gmax < 0:
        raise ValueError("Negative Sander gradient norm")

    terms, literals = _energy_components(component_lines)
    mm_bound = math.fsum(_rounding_bound(literals[name]) for name in MM_TERMS)
    if mm_bound > 0.0004500001:
        raise ValueError("Sander MM component precision is insufficient for relative energies")
    qm_kcal = qm_hartree * conversion
    extern_bound = _rounding_bound(literals[QM_TERM])
    extern_difference = abs(terms[QM_TERM] - qm_kcal)
    if extern_difference > extern_bound + 1e-7:
        raise ExternalEnergyMismatch("Sander EXTERNESCF differs from the checked ORCA total energy")
    mm_sum = math.fsum(terms[name] for name in MM_TERMS)
    total = math.fsum((mm_sum, qm_kcal))
    table_bound = _rounding_bound(row[2])
    if abs(total - rounded) > table_bound + mm_bound + 1e-7:
        raise ValueError("Reconstructed QM/MM energy disagrees with the rounded Sander table")
    if "EAMBER" in terms:
        # EAMBER is a derived diagnostic (epot - restraint), never another
        # additive contribution. Account for it explicitly if it is printed.
        expected = total - terms["RESTRAINT"]
        if abs(terms["EAMBER"] - expected) > mm_bound + _rounding_bound(literals["EAMBER"]) + 1e-7:
            raise ValueError("Sander EAMBER diagnostic differs from the reconstructed total")
    audit = dict(
        status="COMPONENT_RECONSTRUCTED_QMMM_ENERGY_PASS",
        total_source="Nine printed MM components plus checked ORCA total energy",
        mm_components_kcal_mol={name: terms[name] for name in MM_TERMS},
        mm_component_literals={name: literals[name] for name in MM_TERMS},
        mm_component_sum_kcal_mol=mm_sum,
        mm_component_rounding_bound_kcal_mol=mm_bound,
        embedded_qm_energy_hartree=qm_hartree,
        native_au_to_kcal=conversion,
        native_au_to_kcal_literal=conversion_literal,
        embedded_qm_energy_kcal_mol=qm_kcal,
        printed_externescf_kcal_mol=terms[QM_TERM],
        printed_externescf_rounding_bound_kcal_mol=extern_bound,
        externescf_absolute_difference_kcal_mol=extern_difference,
        reconstructed_total_kcal_mol=total,
        rounded_table_energy_kcal_mol=rounded,
        rounded_table_energy_literal=row[2],
        rounded_table_rounding_bound_kcal_mol=table_bound,
        rounded_table_absolute_difference_kcal_mol=abs(total - rounded),
        d3_added_again=False,
        precision_note="MM component rounding bound only; not an estimate of QM/MM method accuracy",
    )
    if "EAMBER" in terms:
        audit["derived_eamber_kcal_mol"] = terms["EAMBER"]
    return dict(
        sander_energy_kcal_mol=total,
        sander_rms_gradient=rms,
        sander_max_gradient=gmax,
        sander_nstep=nstep,
        energy_audit=audit,
    )
