# Two isolated halo-epoch experiments

This document preserves the early isolated experiments. The accepted r014
composition and current counts are documented in
[the runtime contract](../../validationmpi/RUNTIME.md#r014-which-ghost-is-current-and-when)
and [the latest results](../../validationmpi/RESULTS.md). In particular,
the table below describes the constraints considered during development;
some items labelled "future" were subsequently implemented and separately
audited by `entry_unpacked.py`, `wad5.py`, `corner.py` and `momentum_south.py`.
The old overlays are not meant to be applied again to the final source.

`prepare.py --source SOURCE --out NEW_DIR --variant checkready|wadcolor|both` writes only changed files into a new directory. Overlay them into a separate build. This does not replace production source. `checkready` independently composes with the extrema-two-collective overlay (the overlap is the `baro_check.f` file; its reduction tail is not modified here). This work does not use timing to establish scientific equivalence.

* `checkready`: the public `SECOM_CHECK_TEAM` remains conservative and refreshes both halos. Only the `BARO_SUBCYCLE` call directly after the WAD return passes the local-boundary flag to the new implementation. At that instant WAD's final `HALO(DVMWAD)` and `HALO(VAF)` completed, including their team exit barriers. No intervening statement, collar write or team nowait store touches either. The public standalone and non-local boundary paths keep the old calls. No extrema expression, ordering, tie recovery, abort decision or index is changed.
* `wadcolor`: only a persistent, multi-rank, local-boundary BARO call with actual arguments `EL,D` enters a new private routine. Every vertical face `(i,j+1)` is tested and mutated by **one** source iteration `(i,j)`; the two possible tests within the iteration retain their original order, including immediate local ghost zeroing. The opposite color never reads or writes that face. Each rank has only one possible outward source row, `j=JHI`, in exactly one color. Do not erase its flags in color two. ET/ETF are written but not read by WAD; both colors write the same target's frozen ETB bits; ETB is read on its owner. Both color loop joins remain, then flags are sent once and the final VA/DVMWAD/VAF halo contract is restored. The non-MPI path, private non-local-boundary path and public WAD procedure are unchanged. This proof would not carry over if a new caller aliased ELEV/TD to ET/ETF/VA/DVMWAD or a third consumer observed the intermediate color.

Per multi-rank substep, from WAD entry through the end of private CHECK, original = 11 two-sided halo operations + 2 north-flags sends; ready only = 9 + 2; both = 7 + 1. Here a halo call uses two `MPI_Sendrecv` and a flags publish uses one, i.e. 24 → 20 → 15 sends/receives. CHECK's value/index collectives are excluded from these counts. These are operation counts, not a speed prediction.

## Additional audited versions (not implemented by these overlays)

| Epoch | Fields and constraints |
|---|---|
| EXTRNL entry | Every ROTATE assigns UAB/VAB/ELB on every owned cell. ADVAVE/stresses need south UAB, north VAB; V pressure/denominator needs south ELB. WAD also mutates DUMWAD, with no final DUMWAD halo; ADVAVE reads its south row. Do not drop these as unchanged. |
| EXTRNL scratch | FLUXVA is reused: free-surface divergence reads the next row of hydraulic flux; U divergence reads the next row of U-advection flux; V divergence reads the previous row of V-advection flux. These are different producers/values with intervening consumers, **not** three duplicate transfers. A future directional transfer, not deletion, may save work. |
| Ordered boundary/corner edits | Before momentum, CURV42D / UA / VA / ELF can be grouped only **after** their producers. LocalBC proves corner RHS and target share a rank, not that target is off a seam: a neighboring momentum stencil must still observe corner updates. If and only if no corner target touches a relevant shared row (NUMECR=0 is sufficient), pre-corner UA/VA ghosts still suffice. |
| WAD entry | ROTATE changed EL,D,VA so previous EXTRNL input versions are insufficient. EXTRNL/BC13 do not write DVMWAD after entry, so a future private input-ready flag could reuse this halo; it must not change public/independently initialized WAD. |
| Next substep | With current WAD’s bilateral input/output contract, D and EL were refreshed on WAD entry, then only read. VA/DVMWAD are current after WAD's final refresh. CHECK is read-only on these. In substeps 2+ of the **same** local-boundary persistent team an independently testable entry9→entry5 `(UA,UAB,VAB,ELB,DUMWAD)` is possible. Do not assume this across outer-step entry: BARO resets both masks, and independent callers may offer stale ghosts. Do not stack this with any alternative that removes the required WAD publications. |

## Focused unit

`test_epochs.py stage NEW_SUBDIR --variant both|checkready` creates a 12×7 kernel from the actual generated sources plus independently checked-in v46 WAD/CHECK. It never starts the solver. Instrumentation counts actual runtime RAW calls, not source statements. `test_epochs.py build EXISTING_DIR --fc mpifort [--mutant public-stale|drop-reset|drop-vaf]` builds it. For a local non-MPI driver smoke, build with `--fc gfortran --smoke` and run `OMP_NUM_THREADS=3 /path/to/qa-none 1 3`. The MPI version poisons all nonowner mutable fields, compares eight fields on owners and final three-field ghosts, and compares signed extrema and both selected indices. It re-poisons only the stand-alone CHECK ghosts to exercise the public/fused distinction. Negative controls lose a north ET/ETF reset, falsely claim stale CHECK halos are ready, or drop WAD's final VAF publication. The 3-rank layout hits a seam at j=3; 7 ranks cover one-row owners. Use `native.py PARENT --cpus <allocated distinct physical cores>` only after explicit allocation. Stage the two children named `checkready`, `both` under PARENT. The launcher binds and validates each rank/thread CPU list; in the p3 campaign it requires 7–32 allocated physical cores within 4..35. No dense gathering is added to production.
