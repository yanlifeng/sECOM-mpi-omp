# sECOM CPU / OpenMP

This local research repository starts from the verified v6 Macau Hato development kit. See [BUILD.md](BUILD.md) for build flags, input and launch requirements. `secom-omp-v6` is the immutable 32-thread baseline; campaign evidence is in [validation/](validation/). These are laboratory changes to the supplied code, not an upstream release.

For new performance reports, use the original same-machine serial executable, the unchanged 4,800-step case and six exact output hashes. Report separate same-revision scaling when available. Include I/O and initialization; do not compare SMT threads to physical cores. No external remote is configured.

The active research branch is `codex/secom-omp64`; the v6 tag and its
`MANIFEST.json` describe the **original frozen kit**, not a checksum list for
later edits. [The follow-up experiments](validation/OMP64.md) distinguish
accepted source, optional machine tuning and unaccepted candidates.

The verified p3 baseline is now tagged **`secom-omp-v44-64`**. Its scientific
source matches `87490b3a`; on the unchanged Macau development case it completed
in 10.808 seconds at 64 physical cores (7.37% less time than same-batch v30).
The 30.12x figure uses the frozen original serial denominator, whereas
same-executable 1-to-64-thread scaling is 22.47x; the 40x target is not met.
See [the 2026-11-19 update](validation/OMP64.md) for the six-file byte equality,
transient-field checks, cross-date output metadata, and evidence-package audit.
