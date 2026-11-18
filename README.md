# sECOM CPU / OpenMP

This local research repository starts from the verified v6 Macau Hato development kit. See [BUILD.md](BUILD.md) for build flags, input and launch requirements. `secom-omp-v6` is the immutable 32-thread baseline; campaign evidence is in [validation/](validation/). These are laboratory changes to the supplied code, not an upstream release.

For new performance reports, use the original same-machine serial executable, the unchanged 4,800-step case and six exact output hashes. Report separate same-revision scaling when available. Include I/O and initialization; do not compare SMT threads to physical cores. No external remote is configured.
