# Where to change your settings

When you switch molecules, start with that molecule's `system.json`. When you want more memory, more minimization steps, or a tighter convergence target, use the project settings below. Make your edits before submitting the job.

## Choose the QM region and the atoms that can move

Copy the [water pair's system.json](../examples/water_pair/system.json) into the directory with your own AMBER files. All six entries describe the same starting structure.

| Field | What to enter |
|---|---|
| `topology` | Your parameterized, nonperiodic, fixed-charge AMBER topology (`.prmtop`), including the parameters for any ligand or metal site. Filenames are relative to `system.json`. |
| `coordinates` | An ASCII AMBER restart (`.rst7` or `.inpcrd`) containing just the starting coordinates, in topology order. Omit box and velocity records. |
| `qm_atoms` | The atoms ORCA treats quantum mechanically. Count from **1** in the topology. For example, `"1-3,8"` selects atoms 1, 2, 3 and 8. These atoms stay fixed. |
| `movable_mm_atoms` | The MM atoms AMBER can move during minimization. Use the same numbering and range syntax. Leave the QM atoms out of this list. Everything outside the movable set stays fixed. |
| `qm_charge` | The integer charge of the QM fragment, including any link hydrogens. The water pair uses `0`. The [sodium example](../examples/sodium_water/system.json) uses `1` for Na⁺ and its first water. |
| `multiplicity` | The spin multiplicity of that QM fragment, `2S+1`. Both examples use `1` for a singlet. For a metal center, choose the charge and spin state you want to model. ORCA starts from a standard guess at that multiplicity. |

The supplied setup supports QM elements H through Kr and C–C bonds across the QM/MM boundary. AMBER caps each cut with a link hydrogen. Run `qmmm check` with your JSON path and check the reported atom selections, link count and electron counts.

## Adjust the calculation

Open [settings.json](../settings.json) in the project directory.

| Setting | Default | What it controls |
|---|---|---|
| `orca_ranks` | `2` | Number of processes ORCA uses. |
| `maxcore_mb_per_rank` | `1000` | ORCA's memory budget in MB **per process**. |
| `max_cycles` | `200` | Maximum number of MM minimization steps. |
| `drms_kcal_mol_angstrom` | `0.01` | Target RMS gradient in kcal/mol/Å. A smaller value asks for tighter MM convergence. |

## Set up your cluster

Use a Linux Slurm cluster with a shared working directory and a Python 3.11+ environment the compute node can access. The workflow uses AMBER 22 and ORCA 6.1.1 input and output formats. You'll need AMBER's `sander` with its external ORCA interface, `cpptraj`, and ORCA with its matching MPI installation.

- **Time and memory.** The `#SBATCH` lines in [run.sh](../run.sh) request one node, two tasks, 4 GB total and one hour by default. Allow room for AMBER and ORCA overhead above the combined per-process memory budget.
- **Slurm tasks and account.** `qmmm submit` requests the task count in `orca_ranks`. If you use `sbatch run.sh` directly, set `--ntasks` to match. With `qmmm submit`, add `--account YOUR_ACCOUNT` or `--partition YOUR_PARTITION` if your site requires them.
- **Software modules.** Set `amber_modules` and `orca_modules` at the top of [runtime.sh](../src/runtime/runtime.sh). Each nonempty list starts with a clean module environment, so include the compiler and MPI modules your site requires. Leave the lists empty if the programs are already on `PATH`. You can also export `AMBER_MODULES` and `ORCA_MODULES` before submitting.

## Change the method or read the final gradient

The ORCA template uses unrestricted TPSS-D3(0)/def2-SVP. Edit its method or SCF keywords in [orca.inp](../dat/orca.inp). AMBER's minimization and QM/MM options are in [sander.in](../dat/sander.in). Keep the template fields in braces. The result reader uses unrestricted spin populations, the dispersion correction, and gradients for both QM atoms and MM point charges.

After the job finishes, `qmmm collect` reports the final full QM/MM potential and whether the MM gradient met your target. In the run's `result.json`, `mm_geometry_converged` is `true` if it met the target, `false` if it exceeded the target, or `null` if the final gradient could not be verified. The relaxed coordinates are in `relax/relaxed.rst7`, and the final ORCA orbitals are in `final/orc_job.gbw`.
