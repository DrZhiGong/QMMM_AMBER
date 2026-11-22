"""Check what qmmm collect reports for finished, unfinished, and failed runs."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class CollectTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.run_dir = Path(temp.name)
        (self.run_dir / "config.json").write_text(json.dumps({"drms_kcal_mol_angstrom": 0.01}))

    def collect(self):
        return subprocess.run(
            [sys.executable, "-B", str(ROOT / "qmmm.py"), "collect", str(self.run_dir)],
            text=True, capture_output=True,
        )

    def custom_result(self, rms, flag, label, stage="finished"):
        inputs = self.run_dir / "inputs"
        inputs.mkdir(exist_ok=True)
        (inputs / "system.json").write_text("{}")
        status = dict(status=label, mm_geometry_converged=flag)
        (self.run_dir / "setup_status.json").write_text(json.dumps(dict(stage=stage, **status)))
        (self.run_dir / "result.json").write_text(json.dumps(dict(
            full_qmmm_potential_energy_kcal_mol=-100.0, drms_kcal_mol_angstrom=0.01,
            relaxation_rms_gradient=rms, **status,
        )))

    def test_archived_run_uses_its_own_report(self):
        tools = self.run_dir / "tools"
        tools.mkdir()
        (tools / "report.py").write_text("print('Archived report')\nraise SystemExit(2)\n")
        result = self.collect()
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn("Archived report", result.stdout)
        self.assertEqual(result.stderr, "")

    def test_custom_collect_respects_measured_gradient(self):
        for rms, flag, label, code in (
            (0.01, True, "RELAXED_UNREVIEWED", 0),
            (0.1, False, "MM_NOT_CONVERGED_UNREVIEWED", 2),
            (None, None, "FINAL_SINGLE_POINT_MM_UNVERIFIED", 2),
            (0.1, True, "RELAXED_UNREVIEWED", 1),
        ):
            with self.subTest(rms=rms, flag=flag):
                self.custom_result(rms, flag, label)
                result = self.collect()
                self.assertEqual(result.returncode, code, result.stderr)
                if code == 1:
                    self.assertIn("Invalid result", result.stderr)
                    self.assertNotIn("MM converged: True", result.stdout)

    def test_failed_worker_cannot_publish_a_converged_result(self):
        self.custom_result(0.009, True, "RELAXED_UNREVIEWED", stage="FAILED")
        result = self.collect()
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn("FAILED", result.stdout)
        self.assertNotIn("MM converged: True", result.stdout)

    def test_unfinished_custom_run(self):
        inputs = self.run_dir / "inputs"
        inputs.mkdir()
        (inputs / "system.json").write_text("{}")
        result = self.collect()
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn("No result yet", result.stdout)


if __name__ == "__main__":
    unittest.main()
