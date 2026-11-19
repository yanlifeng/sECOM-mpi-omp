# sECOM CPU OpenMP and MPI

This local research repository starts from the verified v6 Macau Hato development kit. See [BUILD.md](BUILD.md) for build flags, input and launch requirements. `secom-omp-v6` is the immutable 32-thread baseline; campaign evidence is in [validation/](validation/). These are laboratory changes to the supplied code, not an upstream release.

For new performance reports, use the original same-machine serial executable, the unchanged 4,800-step case and six exact output hashes. Report separate same-revision scaling when available. Include I/O and initialization; do not compare SMT threads to physical cores. No external remote is configured.

The OpenMP research branch is `codex/secom-omp64`; the v6 tag and its
`MANIFEST.json` describe the **original frozen kit**, not a checksum list for
later edits. [The follow-up experiments](validation/OMP64.md) distinguish
accepted source, optional machine tuning and unaccepted candidates.

The latest verified p3 version is tagged **`secom-omp-v46-64`**, with
`secom-omp-v44-64` retained as the rollback baseline. On the unchanged Macau
development case v46 completed in 10.712 seconds at 64 physical cores,
2.95% less time than fresh same-batch v44 (three independent pairwise wins).
The 30.39x figure uses the frozen original serial denominator; a new v46
single-thread performance run was not done. The 40x target is not met.
See [the 2026-11-19 update](validation/OMP64.md) for the six-file byte equality,
transient-field checks, cross-date output metadata, and evidence-package audit.

The MPI development branch is **`codex/secom-mpi-hybrid`**, based on v46.
Its row-owned barotropic path supports MPI plus OpenMP. The latest r014
32 MPI × 1 thread runs took **12.80–13.06 s**, median **12.9677 s** across
five new interleaved comparisons: **25.105×** versus the original frozen
serial executable, exceeding the new **24×** target in all five samples.
The same-day original serial denominator gives **24.124×** at the median
(not in every sample). Against five fresh previous-MPI controls, the median
wall time fell **24.67%**. All six public files and twelve 25-field snapshots
match the original serial byte for byte; snapshots do not mutate solver halos.
This is not a claim of improvement over the existing OpenMP version.
See [MPI results and scope](validationmpi/RESULTS.md) and the
[communication contract](validationmpi/RUNTIME.md). The accepted source
keeps profiling out of official timing. Nine-field entry packing was rejected;
directional/epoch-aware communication and paired WAD fields shorten halos
only on the admitted private path. Multi-rank CHECK independently validates
its exact integer-key healthy reduction, retaining conservative value
fallbacks and public-call input halos. The old `secom-mpi-v1-32` tag is retained.
