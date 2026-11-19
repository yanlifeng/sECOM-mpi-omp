# Two-phase extrema selection (candidate, not a performance result)

Generate a new two-file overlay against the exact parent source snapshot:

```sh
python3 -B prepare.py --source /snapshot/source --out /fresh/overlay
```

The output directory must not exist. Use both generated files and retain
`extrema2.json`. No changes to the parent, the two dense SIMD loops, either
halo, or the serial/one-rank path. MPI builds assume the parent runtime's
default 32-bit INTEGER/REAL gate and ascending contiguous J ownership.

1. Each owner folds its existing rounded speed and water row summaries in
   ascending J. The four 32-bit words (value bits and J for each) are
   allgathered. The last equal speed wins; the first equal water minimum
   wins. Nonfinite or nonzero subnormal summaries, or unsupported layouts,
   request the original tail via this same message. Classification uses
   integer bits, never an FP comparison on an sNaN. The early water
   summary is private scratch, even when the true global speed aborts.
2. Only the two selected owners re-evaluate the original cells, with the
   original expression/direction: descending I for speed, ascending I
   for water, no water rescan if speed exceeds 100. A three-word integer
   bitwise OR Allreduce publishes their I coordinates and the selected
   water cell's unmodified binary32 bits. Each word has a single writer
   and zero contributions from all other ranks. No MINLOC sign-zero
   assumption, no reduction in floating-point precision.

Normal full-event count: 48,000 iterations, 240,000 -> 96,000 collective
calls per rank (not counting the unchanged halos). The first payload is
16 bytes/rank (512-byte receive at 32 ranks); the second is 12 bytes/rank.
The old pair of row vectors alone received 3,832 bytes/rank. Costs added:
one integer magnitude check per local row minimum and maximum, local
folds, and one rank-record fold. An exceptional case pays for the initial
veto then runs the original tail. Actual speedup must be measured.

Native strict-FP/OMP test in a fresh, explicitly allocated four-core slot:

```sh
taskset -c 20-23 bash /overlay/native_unit.sh /snapshot/source /fixtures/check_extrema_mpi.f
```

The frozen 9 cases and an expanded set (16 normal + 4 exceptional-record
veto tests) passed on p3 with 2x1, 4x1 and 2x2. See
`p3-native-r002.log`; r001 only failed to compile a GCC8-incompatible
BOZ literal in the extra test, corrected using IBSET. Native evidence is
kernel-only: no solver run or timing claim.
