# AMBER_ORCA

Start with two water molecules. One stays fixed in the QM region while AMBER relaxes the other. ORCA supplies the QM energy and gradients at each step. Once you've run the example, you can use the same commands with your own AMBER files.

## Install and check the example

Activate a Python 3.11 or newer environment and clone the repository. If you already have a copy, open its folder in a terminal and start with the `pip install` line.

```bash
git clone https://github.com/DrZhiGong/QMMM_AMBER.git AMBER_ORCA
cd AMBER_ORCA
pip install -r requirements.txt
qmmm check examples/water_pair/system.json
```

The input check works on Windows or Linux with Python alone. For the water pair, look for `6 atoms; 3 real QM atoms; 0 link H; 3 movable MM atoms.` To try a charged QM region, change the input path to `examples/sodium_water/system.json`.

## Run the calculation

For the calculation itself, use a Linux Slurm cluster. Install the Python command on shared storage and make sure the compute node can find AMBER `sander`, `cpptraj`, ORCA and its matching MPI installation. If your cluster uses modules, fill in the two module lists at the top of [runtime.sh](src/runtime/runtime.sh). The [settings guide](doc/settings.md) covers software versions and job resources.

```bash
qmmm submit examples/water_pair/system.json
```

Slurm prints a job ID. While the job is queued, use `squeue -u "$USER"`. Once it starts, the `Run:` line in `slurm-<job ID>.out` gives you this job's directory. Use `qmmm status` to follow its progress if it's your most recently prepared run.

## Find your water-pair results

After the job finishes, run this from `AMBER_ORCA` to see a summary of the most recently prepared run. If you've run several jobs, type `qmmm collect` followed by the directory from your water-pair job's `Run:` line.

```bash
qmmm collect
```

1. **Energy.** Read the `Full QM/MM potential` line for the energy of both waters together at the final geometry, in kcal/mol. The same value is saved in `result.json` as `full_qmmm_potential_energy_kcal_mol`. If you also need the embedded QM energy of the first water, use `embedded_qm_energy_hartree` (Hartree).
2. **Convergence.** Check the `MM converged` line. In `result.json`, `mm_geometry_converged` is `true` when the verified RMS gradient met the target, `false` when it was still above target, or `null` when the final gradient could not be verified. The measured gradient is `relaxation_rms_gradient` and the default target is 0.01 kcal/mol/Å.
3. **Final structure.** In a molecular viewer that reads AMBER files, open `inputs/system.prmtop` from this run, then load `relax/relaxed.rst7` as its coordinates. To compare before and after, load `inputs/initial.rst7` on the same topology. Atoms 1–3 are the fixed QM water. Atoms 4–6 are the MM water that was allowed to move.

## Run your own system

Copy an example's `system.json` beside your parameterized AMBER topology and matching coordinate file. Set the two filenames, QM atoms, movable MM atoms, QM charge and spin multiplicity. Atom numbering starts at **1** in the topology. Then repeat `qmmm check` and `qmmm submit` with the path to your JSON.

To change the number of ORCA processes, memory per process, optimization steps or gradient target, edit [settings.json](settings.json). For Slurm wall time and total job memory, edit [run.sh](run.sh). The [settings guide](doc/settings.md) explains each field and where to change the calculation method.

## Where things live

The example inputs stay in `examples/`. Each time you prepare a calculation or a submitted job starts, it gets its own folder under `runs/`. Here's where to look after a water-pair job finishes.

```text
AMBER_ORCA/
├── dat/                            AMBER and ORCA input templates
├── doc/                            Settings and cluster guide
├── examples/
│   ├── water_pair/                 Starting files for the two waters
│   └── sodium_water/               Starting files for sodium and two waters
├── src/                            Calculation code and cluster scripts
├── test/                           Tests for the scripts
├── settings.json                   ORCA resources and MM minimization settings
├── run.sh                          Slurm job script
└── runs/                           Created when a run is prepared
    ├── latest.json                 Points to the most recently prepared run
    └── qm-mm-<date>-<time>-<id>/    One calculation
        ├── inputs/
        │   ├── system.json         QM and movable MM atom selections
        │   ├── system.prmtop       Topology for viewing either structure
        │   └── initial.rst7        Starting coordinates
        ├── config.json             Settings used for this run
        ├── relax/
        │   ├── sander.out          AMBER minimization output
        │   └── relaxed.rst7        Final coordinates of both waters
        ├── final/
        │   ├── sander.out          Final QM/MM single-point output
        │   └── orc_job.gbw         Final ORCA orbitals
        ├── logs/                   Software and atom-selection checks
        └── result.json             Final energies and MM convergence
```

To run the Python tests from this directory, use `python -m unittest discover -s test`. The launcher test runs with `bash test/test_run.sh` and uses a simulated job.
