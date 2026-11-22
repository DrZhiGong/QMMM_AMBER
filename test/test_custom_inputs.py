"""Check how the input reader handles AMBER parameters and bonds across the QM/MM boundary."""

import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from common import Prmtop, write_rst
from custom_model import load_system


class CustomInputTests(unittest.TestCase):
    def setUp(self):
        workspace = tempfile.TemporaryDirectory()
        self.addCleanup(workspace.cleanup)
        self.folder = Path(workspace.name) / "system"
        shutil.copytree(ROOT / "examples/water_pair", self.folder)
        self.topology = self.folder / "system.prmtop"
        self.parm = Prmtop(self.topology)
        self.spec = self.folder / "system.json"

    def replace_sections(self, **changes):
        """Write replacements in the same fixed-width format as the complete test topology."""
        lines = self.parm.lines.copy()
        for name in sorted(changes, key=lambda key: self.parm.sections[key]["start"], reverse=True):
            section = self.parm.sections[name]
            start, end = section["start"], section["end"]
            values = changes[name]
            if values is None:
                del lines[start:end]
                continue
            width, repeat = section["width"], section["repeat"]
            if section["kind"].lower() == "e":
                tokens = [f"{v:{width}.{section['precision']}E}" for v in values]
            else:
                tokens = [f"{v:{width}}" for v in values]
            lines[start + 2:end] = ["".join(tokens[i:i + repeat]) + "\n" for i in range(0, len(tokens), repeat)]
        self.topology.write_text("".join(lines))

    def carbon_chain(self, qm, moving, charge):
        # Six carbon sites with a single chain of real C-C boundary candidates.
        pointers = self.parm.get("POINTERS").copy()
        pointers[2], pointers[3], pointers[12] = 0, 5, 5
        self.replace_sections(POINTERS=pointers, ATOMIC_NUMBER=[6]*6, MASS=[12.01]*6,
                              ATOM_NAME=["C1","C2","C3","C4","C5","C6"], CHARGE=[0.0]*6,
                              BONDS_INC_HYDROGEN=[], BONDS_WITHOUT_HYDROGEN=[
                                  x for i in range(5) for x in (3*i, 3*(i+1), 1)])
        write_rst(self.folder / "initial.rst7", [[1.54*i, 0.0, 0.0] for i in range(6)], "Carbon chain")
        spec = json.loads(self.spec.read_text())
        self.spec.write_text(json.dumps({**spec, "qm_atoms": qm, "movable_mm_atoms": moving, "qm_charge": charge}))

    def test_water_pair_has_no_link_atoms(self):
        model = load_system(self.spec)
        self.assertEqual(model["link_bonds_0based"], [])
        self.assertEqual((model["nalpha"], model["nbeta"]), (5, 5))

    def test_valid_two_cuts_two_remaining_mm_atoms(self):
        self.carbon_chain("3-4", "1-2,5-6", 0)
        model = load_system(self.spec)
        self.assertEqual(model["link_bonds_0based"], [[1, 2], [3, 4]])
        self.assertEqual((model["nalpha"], model["nbeta"]), (7, 7))

    def test_shared_mm_link_anchor(self):
        self.carbon_chain("1,3", "2,4-6", 1)
        with self.assertRaisesRegex(ValueError, "MM atom 2 connects to multiple"):
            load_system(self.spec)

    def test_too_few_mm_atoms_outside_link_pairs(self):
        self.carbon_chain("2-5", "1,6", 0)
        with self.assertRaisesRegex(ValueError, "adjust_q=2.*leaves 0"):
            load_system(self.spec)

    def test_missing_and_invalid_mm_parameters(self):
        for name in ("MASS", "ATOM_TYPE_INDEX", "NONBONDED_PARM_INDEX",
                     "LENNARD_JONES_ACOEF", "NUMBER_EXCLUDED_ATOMS", "EXCLUDED_ATOMS_LIST"):
            for replacement in (None, []):
                with self.subTest(name=name, replacement=replacement):
                    self.replace_sections(**{name: replacement})
                    with self.assertRaisesRegex(ValueError, name):
                        load_system(self.spec)

    def test_negative_hbond_index_and_zero_exclusion_entries(self):
        pointers = self.parm.get("POINTERS").copy()
        pointers[19] = 1
        self.replace_sections(POINTERS=pointers, NONBONDED_PARM_INDEX=[1, -1, -1, 3],
                              HBOND_ACOEF=[1.], HBOND_BCOEF=[1.], HBCUT=[0.])
        self.assertEqual(load_system(self.spec)["natom"], 6)

    def test_duplicate_or_malformed_topology_section(self):
        original = "".join(self.parm.lines)
        for extra, message in (("%FLAG CHARGE\n%FORMAT(5E16.8)\n", "Duplicate.*CHARGE"),
                               ("%FLAG UNKNOWN\n", "Missing %FORMAT.*UNKNOWN"),
                               ("%FLAG UNKNOWN\n%FORMAT(10I0)\n", "Invalid.*%FORMAT.*UNKNOWN")):
            with self.subTest(message=message):
                self.topology.write_text(original + extra)
                with self.assertRaisesRegex(ValueError, message):
                    load_system(self.spec)


if __name__ == "__main__":
    unittest.main()
