"""Find the ORCA gradient file whose energy matches the final EXTERNESCF in sander.out.
The comparison allows for unit conversion and the rounding in the printed AMBER energy.
"""

import hashlib
import math
import re
from decimal import Decimal
from pathlib import Path

from common import read
from sander_results import NUMBER, ExternalEnergyMismatch, reconstruct_potential


# sander prints a rounded AU_TO_KCAL in the banner but uses this CODATA08 value internally.
REFERENCE = "https://github.com/Amber-MD/AmberClassic/blob/main/src/msander/constants.F90#L158-L176"
K_NATIVE = Decimal("4.3597439422e-18") / Decimal("4.184") / 1000 * Decimal("6.0221417930e23")
ARITHMETIC_SLACK = Decimal("1e-7")  # kcal/mol


class EnergyMismatch(ValueError):
    """This callback's ORCA gradient has a different energy from sander's final step."""


def _decimal(token):
    if not re.fullmatch(NUMBER, token):
        raise ValueError(f"Invalid numeric field: {token}")
    return Decimal(token.upper().replace("D", "E"))


def _half_bin(value):
    """Half a unit in the last printed digit."""
    return Decimal("0.5").scaleb(value.as_tuple().exponent)


def _engrad_energy(bridge):
    event = Path(bridge["event"]).resolve(strict=True)
    path = event / "orca.engrad"
    if Path(bridge["gradient_file"]).resolve(strict=True) != path:
        raise ValueError("Gradient file points outside the selected ORCA event")
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != bridge["gradient_sha256"]:
        raise ValueError("Archived native gradient checksum changed")
    with path.open(encoding="ascii") as stream:
        fields = (v for line in stream if not line.lstrip().startswith("#") for v in line.split())
        natoms = int(next(fields, "0"))
        energy = _decimal(next(fields, ""))
    if natoms != len(bridge["native_qm_atoms"]) or natoms <= 0:
        raise ValueError("Native gradient atom count differs from the verified event")
    return energy


def _printed_energies(content):
    """Read the conversion and final EXTERNESCF exactly as printed, including trailing zeros."""
    factors = re.findall(r"^[ \t]*AU_TO_KCAL[ \t]*=[ \t]*(\S+)[ \t]*$", content, re.M | re.I)
    if len(factors) != 1:
        raise ValueError("Missing or ambiguous Sander AU_TO_KCAL")
    factor = _decimal(factors[0])
    if _half_bin(factor) > Decimal("5e-9") or abs(factor - K_NATIVE) > _half_bin(factor):
        raise ValueError(
            "Printed factor does not identify AMBER's CODATA08_AU_TO_KCAL at the expected precision"
        )
    sections = re.split(r"^\s*FINAL RESULTS\s*$", content, flags=re.M | re.I)
    if len(sections) != 2:
        raise ValueError("Missing or ambiguous final Sander results")
    terms = re.findall(r"\bEXTERNESCF\s*=\s*(\S+)", sections[1], re.I)
    if len(terms) != 1:
        raise ValueError("Missing or ambiguous final Sander EXTERNESCF")
    external = _decimal(terms[0])
    if _half_bin(external) > Decimal("0.00005"):
        raise ValueError("EXTERNESCF has insufficient printed precision")
    return factors[0], terms[0]


def parse_precision_results(content, bridge):
    """Check a callback at AMBER's internal precision, then report using Sander's printed conversion."""
    if bridge.get("converged") is not True:
        raise ValueError("An SCF-converged ORCA callback is required")
    gradient_hartree = _engrad_energy(bridge)
    stdout_hartree = _decimal(str(bridge["energy_hartree"]))
    if abs(gradient_hartree - stdout_hartree) > Decimal("1e-7"):
        raise ValueError("Archived native gradient and ORCA stdout energies disagree in Hartree")

    factor_text, external_text = _printed_energies(content)
    printed_factor = _decimal(factor_text)
    external = _decimal(external_text)
    delta = gradient_hartree * K_NATIVE - external
    bound = _half_bin(external) + ARITHMETIC_SLACK
    try:
        parsed = reconstruct_potential(content, bridge, conversion_literal=str(K_NATIVE))
    except ExternalEnergyMismatch as exc:
        if abs(delta) > bound:
            raise EnergyMismatch(
                f"Event {bridge['callback_index']} differs from EXTERNESCF "
                f"by {delta:+.9f} kcal/mol using its native gradient energy"
            ) from exc
        raise
    if abs(delta) > bound:
        raise ValueError(
            "The native gradient does not reproduce EXTERNESCF while the ORCA stdout energy does"
        )

    # Report the total using the conversion in Sander's output.
    audit = parsed["energy_audit"]
    qm_reported = float(stdout_hartree) * float(printed_factor)
    total = math.fsum((audit["mm_component_sum_kcal_mol"], qm_reported))
    audit["native_precision_energy_check"] = dict(
        status="NATIVE_CODATA08_EXTERNESCF_PASS",
        constants_source=REFERENCE,
        conversion_expression="4.3597439422e-18 / 4.184 / 1000 * 6.0221417930e23",
        conversion_full_precision=str(K_NATIVE),
        conversion_printed_literal=factor_text,
        externescf_printed_literal=external_text,
        rounding_half_bin_kcal_mol=float(_half_bin(external)),
        arithmetic_slack_kcal_mol=float(ARITHMETIC_SLACK),
        engrad_sha256=bridge["gradient_sha256"],
        engrad_energy_hartree=str(gradient_hartree),
        stdout_energy_hartree=str(stdout_hartree),
        engrad_minus_extern_kcal_mol=float(delta),
        stdout_minus_extern_kcal_mol=float(stdout_hartree * K_NATIVE - external),
        stdout_minus_extern_using_printed_factor_kcal_mol=float(stdout_hartree * printed_factor - external),
        reconstructed_total_native_factor_kcal_mol=parsed["sander_energy_kcal_mol"],
    )
    parsed["sander_energy_kcal_mol"] = total
    audit.update(
        native_au_to_kcal=float(printed_factor),
        native_au_to_kcal_literal=factor_text,
        embedded_qm_energy_kcal_mol=qm_reported,
        reconstructed_total_kcal_mol=total,
        externescf_check_au_to_kcal=float(K_NATIVE),
        rounded_table_absolute_difference_kcal_mol=abs(
            audit["rounded_table_energy_kcal_mol"] - total
        ),
        total_source="Nine printed MM components plus checked ORCA stdout total times the printed factor",
    )
    return parsed


def energy_matches(content, last, *, all_events):
    """Try the last ORCA call first. A minimization can finish at an earlier call's geometry."""
    event = Path(last["event"]).resolve(strict=True)
    for index in range(last["callback_index"], 0, -1):
        directory = event.parent / f"{index:05d}"
        callback = read(directory / "result.json")
        if callback["callback_index"] != index or Path(callback["event"]).resolve() != directory:
            raise ValueError(f"Unexpected evaluation history at {directory}")
        if index == last["callback_index"] and callback != last:
            raise ValueError("Final ORCA event differs from its summary")
        try:
            parsed = parse_precision_results(content, callback)
        except EnergyMismatch:
            if not all_events:
                raise
            continue
        parsed["energy_audit"]["matched_callback_index"] = index
        yield index, parsed
        if not all_events:
            return
