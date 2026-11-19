# MPI diagnostic-only timing

```sh
python3 -B mpi/profile/prepare.py --source <frozen-science-source> --out <new-profile-dir> --build
# Use the independent diagnostic runner and fresh rank directories; the generated ELF
# is <new-profile-dir>/s--pom-profile.x. It must never enter an official timing series.
python3 -B mpi/profile/summarize.py <finished-run> --ranks 32 --threads 1 --steps 300 --split 10 --out <new-report.json>
```

Omit `--build` for local source generation. The p3 native mode follows the
campaign's GCC8/Open MPI strict flags, compiles three separate objects and
links the parent's other frozen objects (digests and commands in
`profile-build.json`). The output directory must not already exist.
Run only with normal campaign resource approval; no new collective or
OpenMP team synchronization is added. Non-profile binaries are untouched.

Each rank writes just its own `smpi-profile-RRRRRR.log`, after MPI finalize,
in its existing private working directory. The runner must allow that
additional diagnostic file without relaxing any original-output checks.
A killed/abnormal/partial run must not be used for attribution.

* **life**: five disjoint rank-local intervals: SMPI initialization, subsequent
  model initialization up to the first time-loop statement, the entire time
  loop, post-loop work/output, and SMPI finalization. Time before the first
  main-program SMPI_INIT, MPI launcher/process teardown and profile file I/O
  are excluded. Do not substitute this for launcher wall time.
* **site**: times and dynamic counts of BCOND(7), CBC, BOTTAU, SMAG, BARO,
  ARCHIVE, PUTCDF/TSCDF, PRINTS and the scheduled/final restart windows
  (including any required snapshot); categories are distinguished by lifetime.
  EXTRNL, wet/dry and extrema-check are **inside** BARO, while PUTCDF/TSCDF
  are **inside** ARCHIVE. They are inclusive scope times, not independent
  wall contributions. A writer's PUTCDF includes packing, netCDF and formatted
  text; it is not pure disk duration. The immutable inputs' physically separate
  per-rank initialization remains intact.
* **mpi**: actual time in blocking RAW MPI calls: REAL-plane Allgatherv, halo
  Sendrecv, row-stat and integer-row-stat Allgatherv, coordinate Allreduce,
  row-rooted Bcast, WAD flag Sendrecv, and packed Halo9 Sendrecv if the parent
  has that optional feature. For modern source revisions this includes
  check-record Allgather, fallback picks Allreduce, healthy-key Allreduce,
  one-way generic field Sendrecv (including VAF), WAD paired planes and
  unpacked or packed EXTRNL entry Sendrecv. The manifest records exact
  static callsite counts; the log counts only executed blocking calls,
  not bytes, source lines or high-level helper invocations.
  These **do not overlap each other on a rank**.
  They are already included in the appropriate life/site interval. Multiple
  MPI calls per helper are counted individually; for example an ordinary
  two-direction halo executes two Sendrecvs. MPI setup/finalize collectives
  are included only in their respective life intervals. OpenMP waiting before
  SINGLE, local pack/unpack, local checks and overall load imbalance are not
  reclassified as actual MPI calls.
* **Limitations**: time inside a blocking call can include late peers, protocol
  overhead, transfer, and descheduling; it is not network latency or a pure
  wait/spin fraction. No timer synchronizes the clocks of different ranks.
  In hybrid mode, each thread labels its own MPI context; only team master
  times the EXTRNL/WAD/check span, without a new join. The 32×1 run avoids
  that distinction. Fine-grained timers perturb the profile ELF.

The summary verifies rank ids, normal loop exit, declared OpenMP team size and
observed counts for outer calls and all substeps. It reports rank min/median/max,
sum of rank-seconds, and each rank's MPI fraction of **its own** loop.
Missing optional metrics on a rank count as zero; a maximum from one metric
need not belong to the same rank as a maximum from another. Never add nested
site/mpi timers or sum ranks to estimate wall time.

For an admitted multi-rank production case using unpacked entry+WAD5+keys,
pass `--expect-owned-epochs` during summarization. This verifies **each rank**
executed 9 entry Sendrecvs on first substeps and 5 thereafter, 4 WAD
paired-plane calls, one WAD flag publication and one VAF directional
Sendrecv per substep, and one health-key Allreduce per substep. The
last three WAD categories are distinct: WAD5's five field sends do not
include the separate flag publication. No claim about fallback PICK calls
is inferred from the healthy-key count; use the reported separate counters.

Local smoke: link the generated `secom_mpi.f90` (no `SECOM_USE_MPI`) and
`check_timers.f90` with `-cpp -fopenmp -fcheck=all`; run with one and two
OpenMP lanes in separate `rank000` directories. Pass inputs `--ranks 1
--steps 3 --split 2` and the real lane count to `summarize.py`. This checks
the rank-local accounting only, not native MPI/scientific validation.
