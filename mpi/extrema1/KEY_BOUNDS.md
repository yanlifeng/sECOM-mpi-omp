# Independent audit: `SMPI_TRY_HEALTH` in the r014final source

For the qualified binary32/default-32-bit-integer build, let `B=2^31`,
`T=IM*JM`, and `M=0x7f7fffff` (unsigned magnitude of the largest finite
binary32). The production guard admits `1<=T<=2^30`, `1<=J<=JM<2^30`.
Both factors in `T` are promoted before multiplication. Public/no-winner,
abnormal and `>100` cases retain the previous verdict path; zero-initialized
records are never mistaken for valid keys.

- Water: IEEE positive normal magnitudes increase with numerical value;
  negative ones are ordered in reverse. The score `B-mag` or `B+mag` just
  reverses that ordering. The two zero signs have the same score. A row is
  accepted only if strictly below the initial `1.e10` (`0x501502f9`), as
  before. When scores match, `2*(T-loc)+actual_cell_sign` breaks ties in
  favor of smaller row-major `loc`, hence first J, then first I. One cell of
  separation strictly beats either sign; decode retains the winner's sign.
  This suffix is in `[0,B-1]`. The maximum reachable key is
  `(B+M)*B+B-1=0x7fbfffffffffffff`, exactly `2^54` below `INT64_MAX`.
  The decoded magnitude is at most `M`, so conversion back to default
  integer is safe; decoded coordinates are at most `T`.
- Speed: a healthy accepted summary is signed zero or positive normal up
  to 100 (`0x42c80000`). The suffix `2*J+row_sign` fits below `B`; valid
  row scores are positive. Thus larger speed wins, then later J; zero sign
  cannot overturn J. The loose upper bound `0x216400007fffffff` fits signed
  int64. Only the true global winning speed row executes the original
  descending-I `SQRT` search, selecting last I and preserving its floating
  exceptions. This old expression is byte-identical to `SMPI_VERDICT2` and
  the success store remains `VOLATILE`. In a healthy call only, water then
  replaces the returned coordinates, as it did originally.
- The collective is elementwise signed `MPI_MAX` on three `MPI_INTEGER8`
  words. All data keys are non-negative; `0` means absent, the third word
  conveys any rank's veto. If either key is missing or any veto is present,
  the fast path returns before assigning public results.

`key_properties.py SOURCE` extracts and compiles the production encoding
and decoding expressions with integer overflow trapping. It compares output
to an independent Python order on `(numeric_value,J,I)` and checks loaded
cell sign, exact recovered indices and size vetoes. It does not allocate a
`2^30`-cell array or claim to enumerate all `2^32` bit patterns. The proof
above covers all finite binary32 normal values; the executable exercises
representatives from both ends/interior of every exponent bin, both zeros,
threshold ULPs, equal values across cells, and legal/illegal dimensions.
The separate native MPI test (`p3-native-r001.log`) exercised the full code,
including a trapped global-winner `SQRT` and the uninspected losing row.
