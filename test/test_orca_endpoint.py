"""Check that the ORCA energies, orbitals, and gradients belong to the step that produced them."""

import os
from contextlib import chdir
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import endpoint
import orca_bridge
from common import sha, write, write_xyz


def orbital_table(spin, occupied):
    rows = [f"{i} {float(i < occupied):.4f} -0.123456 -3.3594" for i in range(7)]
    return f"SPIN {spin} ORBITALS\nNO  OCC  E(Eh)  E(eV)\n" + "\n".join(rows) + "\n"


class OrcaTests(unittest.TestCase):
    def test_final_orbital_tables(self):
        up, down = orbital_table("UP", 6), orbital_table("DOWN", 5)
        self.assertEqual(orca_bridge._spin_electron_counts(up + down), [6, 5])
        for incomplete in (up + up, down + up, down, "SPIN UP ORBITALS\n" + down):
            with self.subTest(incomplete=incomplete), self.assertRaises(ValueError):
                orca_bridge._spin_electron_counts(incomplete)

    def test_gradient_coordinates_in_bohr(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "orc_job.inp"
            path.with_suffix(".pcgrad").write_text("1\n0 0 0\n")
            atoms = [("O", [1.0, 0.0, -0.5])]
            header = "# Number of atoms\n1\n# Energy (Eh)\n-75.0\n# Gradient\n0\n0\n0\n"
            gradient = path.with_suffix(".engrad")
            gradient.write_text(header + "# Atomic numbers and coordinates in Bohr\n8 1.889726 0 -0.944863\n")
            orca_bridge.check_gradient(path, atoms, -75.0, 1)
            gradient.write_text(header + "8 2.0 0 -0.944863\n")
            with self.assertRaisesRegex(ValueError, "geometry"):
                orca_bridge.check_gradient(path, atoms, -75.0, 1)

    def test_callback_requires_fresh_outputs_and_keeps_the_seed(self):
        # An external program may exit normally without writing one of its expected files.
        for missing in (".gbw", ".engrad", ".pcgrad", None):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as tmp, chdir(tmp):
                path = Path("orc_job.inp")
                context = Path("context.json").resolve()
                write(context, dict(charge=1, input_multiplicity=2, nprocs=1,
                                    maxcore_mb=1000, real_qm_atoms=[dict(element="Mg", coord=[1, 0, 0])],
                                    expected_qm_atoms_with_links=1, cap_parent_real_indices=[]))
                Path("charges.pc").write_text("1\n-1 3.0 0 0\n")
                generated = '%pointcharges "charges.pc"\n* xyz 1 2\nMg 1 0 0\n*\n'
                engrad = b"1\n-100.0\n0.0\n0.0\n0.0\n12 1.889726 0 0\n"
                stdout = (
                    "SCF CONVERGED\n" + orbital_table("UP", 6) + orbital_table("DOWN", 5)
                    + "MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS\n0 Mg : 1.0 1.0\n"
                    "Expectation value of <S**2> : 0.75\nNumber of eigenvalues below threshold .... 0\n"
                    "Dispersion correction -0.001\nFINAL SINGLE POINT ENERGY -100.0\nORCA TERMINATED NORMALLY\n"
                )

                def fake_orca(command, output):
                    initial = output.parent.name == "00001"
                    atoms, charge, spin, used_charges, count, charge_sum = orca_bridge.generated_input(path)
                    self.assertTrue(used_charges.samefile(output.parent / "embedding.pc"))
                    self.assertEqual((charge, spin, count, charge_sum), (1, 2, 1, -1.0))
                    if not initial:
                        self.assertEqual(Path("orbital_seed.gbw").read_bytes(), b"seed")
                        self.assertIn('MORead', path.read_text())
                    for suffix, contents in ((".gbw", b"seed" if initial else b"new"),
                                             (".engrad", engrad), (".pcgrad", b"1\n0 0 0\n")):
                        self.assertFalse(path.with_suffix(suffix).exists())
                        if initial or suffix != missing:
                            path.with_suffix(suffix).write_bytes(contents)
                    output.write_text(stdout)

                env = dict(AMBER_ORCA_CONTEXT=str(context), ORCA_BIN="orca", ORCA_MPI_FLAGS="")
                with patch.dict(os.environ, env), patch.object(orca_bridge, "_run_orca", fake_orca):
                    path.write_text(generated)
                    orca_bridge.run(path)
                    path.write_text(generated)
                    if missing:
                        with self.assertRaises((FileNotFoundError, ValueError)):
                            orca_bridge.run(path)
                        self.assertFalse(Path("evaluations/00002/result.json").exists())
                    else:
                        orca_bridge.run(path)
                        self.assertEqual(path.with_suffix(".gbw").read_bytes(), b"new")


def callback(folder, index, energy):
    event = folder / "evaluations" / f"{index:05d}"
    event.mkdir(parents=True)
    native = [["O", [0, 0, 0]], ["H", [1, 0, 0]], ["H", [0, 1, 0]]]
    write_xyz(event / "native_qm.xyz", native, "QM water")
    charges = event / "embedding.pc"
    charges.write_text("2\n0.4 3.0 0.0 0.0\n-0.4 4.0 0.0 0.0\n")
    record = dict(event=str(event), callback_index=index, converged=True, coordinate_multiplicity=1,
                  nalpha=5, nbeta=5, energy_hartree=energy, native_qm_atoms=native,
                  native_qm_sha256=sha(event / "native_qm.xyz"), native_real_qm_indices=[0, 1, 2],
                  pointcharges_file=str(charges), pointcharge_count=2, pointcharge_sha256=sha(charges))
    write(event / "result.json", record)
    return record


class EndpointTests(unittest.TestCase):
    spin = (1, 5, 5)

    def test_restart_energy_must_agree_at_matched_geometry(self):
        for delta, expected in ((0.0000002, True), (0.1, None)):
            with self.subTest(delta=delta), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp).resolve()
                relax = root / "relax"
                old = callback(relax, 1, -76.0)
                final = root / "final"
                fresh = callback(final, 1, -76.0 + delta)
                (relax / "sander.out").write_text("")
                relaxation = dict(bridge=old, energy_audit=dict(matched_callback_index=1), sander_rms_gradient=0.001)
                # Callback 1's match to the source Sander EXTERNESCF was already recorded.
                with patch.object(endpoint, "energy_matches", return_value=[(1, {})]):
                    result = endpoint.relaxation_status(relax, relaxation, dict(folder=str(final), bridge=fresh),
                                                        0.01, expected_spin=self.spin)
                self.assertIs(result["mm_geometry_converged"], expected)
                row = result["endpoint_audit"]["candidates"][0]
                self.assertAlmostEqual(row["embedded_energy_difference_hartree"], delta)
                if expected is None:
                    self.assertIsNone(result["relaxation_rms_gradient"])
                    self.assertIn("energy changed", row["reason"])

    def test_electronic_change_between_single_point_callbacks(self):
        with tempfile.TemporaryDirectory() as tmp:
            final = Path(tmp).resolve() / "final"
            callback(final, 1, -76.0)
            last = callback(final, 2, -75.95)
            result = endpoint.relaxation_status(final.parent / "relax", {}, dict(folder=str(final), bridge=last),
                                                0.01, expected_spin=self.spin)
            self.assertEqual(result["status"], "FINAL_SINGLE_POINT_MM_UNVERIFIED")
            audit = result["endpoint_audit"]["final_single_point_check"]
            self.assertIn("energy changed", audit["candidates"][0]["reason"])


if __name__ == "__main__":
    unittest.main()
