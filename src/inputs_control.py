"""Fill in the AMBER and ORCA templates with the atom selections, charge,
spin multiplicity and calculation settings for this run.
"""

from pathlib import Path

TEMPLATES = Path(__file__).resolve().parent.parent / "dat"


def sander_input(*, qm_mask, charge, multiplicity, cycles, drms, moving_mask=None):
    if not qm_mask or len(qm_mask) > 256 or "'" in qm_mask or "\n" in qm_mask:
        raise ValueError("A checked QM mask of at most 256 characters is required")
    belly = "ibelly=0"
    if moving_mask is not None:
        if not moving_mask or len(moving_mask) > 256 or "'" in moving_mask or "\n" in moving_mask:
            raise ValueError("A checked movable-atom mask of at most 256 characters is required")
        belly = f"ibelly=1, bellymask='{moving_mask}'"

    template = (TEMPLATES / "sander.in").read_text()
    return template.format(
        cycles=cycles,
        descent_cycles=min(20, cycles),
        drms=drms,
        belly=belly,
        qm_mask=qm_mask,
        charge=charge,
        multiplicity=multiplicity,
    )


def orca_template(context, *, has_seed):
    """Use the same SCF and basis settings in Sander's template and each ORCA call."""
    template = (TEMPLATES / "orca.inp").read_text()
    rendered = template.format(
        guess_keyword=" MORead" if has_seed else "",
        ranks=context["nprocs"],
        maxcore=context["maxcore_mb"],
        read_seed='%moinp "orbital_seed.gbw"' if has_seed else "",
        spin_guess="" if has_seed else " Guess PModel",
    )
    return "\n".join(line for line in rendered.splitlines() if line) + "\n"
