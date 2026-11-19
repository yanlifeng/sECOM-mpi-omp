# Single-collective healthy extrema candidate

Overlay only on an extrema2 parent: `python3 -B prepare.py --source
/snapshot/source --out /fresh/overlay`. Apply it last to the combined
source tree; do not replace a combined runtime with one prepared against
another parent. The public CHECK interface and both dense SIMD loops
remain intact.

When global speed is at most 100 and there is a genuine water minimum,
the original final coordinates are the first minimum-water cell. The
speed winner still rescans its actual global winning row with the old
descending-I SQRT expression; a volatile local integer prevents
dead-code removal even though the later water result supplies both
returned coordinates. No losing speed row is scanned.

Each rank examines at most one additional water row (its own local
minimum); all its speculative loads and comparisons are integer-only.
One native MPI_MAX over three signed INT64 values publishes the global
speed key, water key and veto. No custom MPI datatype or operator.

* Speed prioritizes positive binary32 magnitude, then largest global J,
  then the original sign of zero. Negative speed sentinels are ignored.
* The water score is 2^31 minus the magnitude of a positive float, or
  2^31 plus the magnitude of a negative float. Signed zeros have equal
  score. The low 31 bits encode 2*(IM*JM-global_index)+cell_sign;
  the earlier Fortran-order cell wins, with its actual raw zero sign.
  The range gate IM*JM <= 2^30 prevents signed INT64 overflow.
* Aborts, all-land, no speed/cell match, unsupported layout and
  exceptional/nonzero-subnormal summaries or inspected FSM/ELF values
  leave every public output untouched and request the complete extrema2
  path. The previous fallback from extrema2 to the legacy verdict
  remains intact. Such a case pays an extra small collective; typical
  healthy checks replace the previous two collectives with one.

Validation: original nine solver-level kernel fixtures, expanded 16+4,
and 18 direct fast/veto contract cases passed for p3 2x1, 4x1 and 2x2.
Deliberate sNaN in the true global speed winner raised SIGFPE (136) in
the old and new verdict. A poisoned losing speed row did not trap. See
`p3-native-r001.log`. This is kernel evidence, not full-event timing
or science acceptance.
