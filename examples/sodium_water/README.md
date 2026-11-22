# Sodium ion with two waters

Try this example to see how to set up a charged QM region. Na⁺ and its first water stay fixed. The second water can move, and its oxygen starts 2.8 Å from the ion.

Open `system.json` to see where the charge goes. The QM region includes Na⁺ (atom 1) and the first water (atoms 2–4), with charge +1 and singlet multiplicity. Atoms 5–7 are the movable MM water. The atom order and force-field parameters come from `system.prmtop`. The starting coordinates come from `initial.rst7`.

From the `AMBER_ORCA` directory, run the input check.

```bash
qmmm check examples/sodium_water/system.json
```

You should see 4 QM atoms, 3 movable MM atoms and 0 link hydrogens. With charge +1, the QM region has 10 alpha and 10 beta electrons. To submit the job and find its results, follow the [main README](../../README.md).

The files are ready to use. If you're curious how the model was built, `leap.in` combines TIP3P water with the Joung–Cheatham Na⁺ parameters for TIP3P. To rebuild the files, run `tleap -f leap.in` from this example directory with AmberTools loaded. [LEaP input](leap.in) · [TIP3P parameters](https://github.com/Amber-MD/AmberClassic/blob/main/dat/leap/parm/frcmod.tip3p) · [Ion parameters](https://github.com/Amber-MD/AmberClassic/blob/main/dat/leap/parm/frcmod.ionsjc_tip3p)
