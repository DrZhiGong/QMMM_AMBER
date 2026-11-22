"""Check which example qmmm opens by default and what it sends to Slurm."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import qmmm


class CommandLineTests(unittest.TestCase):
    def test_default_example(self):
        output = StringIO()
        with patch.object(sys, "argv", ["qmmm", "check"]), redirect_stdout(output):
            qmmm.main()
        self.assertIn("6 atoms; 3 real QM atoms; 0 link H; 3 movable MM atoms.", output.getvalue())

    def test_submit_uses_system_and_configured_task_count(self):
        system = ROOT / "examples/sodium_water/system.json"
        args = ["qmmm", "submit", str(system), "--account", "test-account", "--partition", "test-queue"]
        with (patch.object(sys, "argv", args), patch.object(qmmm.subprocess, "run") as submit,
              patch.object(qmmm, "load_settings", return_value={"orca_ranks": 4}),
              redirect_stdout(StringIO()) as output):
            qmmm.main()
        self.assertIn("QM charge 1, multiplicity 1", output.getvalue())
        submit.assert_called_once()
        self.assertEqual(submit.call_args.args[0], [
            "sbatch", "--ntasks=4", "--account=test-account", "--partition=test-queue",
            str(ROOT / "run.sh"), str(system),
        ])
        self.assertEqual(submit.call_args.kwargs["env"]["QMMM_PROJECT_DIR"], str(ROOT))
        self.assertEqual(submit.call_args.kwargs["env"]["QMMM_PYTHON"], sys.executable)
