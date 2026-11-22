# Water pair

Start here to see how a QM/MM calculation works with just two water molecules. The first water stays fixed while AMBER relaxes the second. Their oxygen atoms start 3.22 Å apart.

Open `system.json` to see how the atoms are assigned. Atoms 1–3 form the neutral, singlet QM region. Atoms 4–6 are the movable MM water. The file points to `system.prmtop` for the atom order and TIP3P force-field parameters, and to `initial.rst7` for the starting coordinates.

From the `AMBER_ORCA` directory, run the input check.

```bash
qmmm check examples/water_pair/system.json
```

You should see 3 QM atoms, 3 movable MM atoms and 0 link hydrogens. The QM water has 5 alpha and 5 beta electrons. To submit the job and find its results, follow the [main README](../../README.md).

The topology and coordinates are ready to use. If you want to rebuild them, `leap.in` contains the recipe. Run `tleap -f leap.in` from this example directory with AmberTools loaded. [LEaP input](leap.in) · [AMBER water library](https://github.com/Amber-MD/AmberClassic/blob/main/dat/leap/lib/solvents.lib) · [TIP3P parameters](https://github.com/Amber-MD/AmberClassic/blob/main/dat/leap/parm/frcmod.tip3p)
