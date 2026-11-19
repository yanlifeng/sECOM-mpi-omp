# Owned-J ARCHIVE candidate (not a general replacement for the old path)

The only scientific source in this candidate is `archive.f`, derived from v46 `1b58f33639538c2d214f14a8f9a575ceffed5267:source/archive.f` (SHA-256 `e38b83224cfb10b38ec1eca494b6fff291966371110e1330ba60f0f1a5d4041b`). The new policy common `/SMPI_OUTER_FLAG/` must be initialized to false in precisely one object (`SMPI_ARC_DEFAULT` here). Set it identically across ranks before the outer evolution and do not change it halfway through an average. It does not control physical output ownership.

## Eligibility and consumers

Use only on a cold-start, internal, CLAMPED, BAROTROPIC, CDF-output case with JTM=0, NCHEMLAY=0, and neglected tracer, sediment, chemistry, particles and waves. The caller owns the JTM/cold-start configuration checks (not in this common deck); retain the narrower forcing, friction, mixing and local-boundary guards from the integrator. Never switch to a legacy consumer of stale rows without publishing the exact fields that consumer reads. In particular, binary history, sediment/chemistry conversion, baroclinic archive, arbitrary reconfiguration/restart are not covered.

* INTX=0: old full-domain public path, including its `IAVGE` and history-window decisions. The main barotropic internal driver also skips INTX=1, as originally. Loop 530 and history resets 570/580 are owned only from INTX>0. On every active window an EL halo supplies the j−1 flux neighbor. H/H1/H2 and other static geometry remain replicated. Cellwise accumulation order and the bed-mask branch are unchanged.
* Per-write publication: the 13 REAL 2D history fields are ARCET, ARCPATM, ARCWSX2D, ARCWSY2D, ARCAIRT, ARCCLD, ARCRHUM, ARCSWOBS, ARCQPREC, ARCQEVAP, ARCWVHT, ARCWVPD, ARCWVDR. Five REAL 3D histories contribute K=1 only: ARCU, ARCV, ARCAM, ARCUX, ARCVX. Ten direct planes are CBC, VARWIF, AREAICE, QPREC2D, QEVAP2D, DUMWAD, DVMWAD, AAM2D, AAM(:,:,2), KM(:,:,2). All 28 are allgathered on a physical history event, not every outer step. This intentionally supports all-rank PUTCDF computation; guard physical I/O separately.
* BAROTROPIC PUTCDF sets KBot=1, PTD/PTS/PTT/PTW/PTKH/PTKM off. The same routine hardwires PTHF/PTQ2 off and PTC off for neglected tracer. Its active source has no ARCTAU, ARCWS1/2, ARCSED1/2, ARCHEM1/2, ARCTHIK or ARCPBED read: these are not full-plane payloads even when zero certification rejects. They are still accumulated/reset on their owners, including the original order and signs. This argument does not apply to the other output formats or physics.
* Certification is not `physics disabled ⇒ +0`. Each owner checks both source and accumulator bits (including negative zero and only the selected bed-mask branch). An all-rank integer MAX of the rejection flag certifies before taking the existing skip path. The old all-grid certificate still applies on the default and cold-start paths.
* At nonempty E/V/F lists publish the current input fields for that list before **every** station accumulation, not just the skill file event. Velocity needs the north-neighbor VA and cross-sections use old stagger and serial REAL addition order. Inactive sediment/chemistry sources may be nonzero, but on this branch are read-only replicated from initialization. At skill samples compute independent per-cell energy contributions on owned rows with a VA halo, gather three REAL fields, then have rank zero perform the original i-fastest serial fold and broadcast the three final bit patterns. No MPI_SUM or regrouping occurs. Preserve the original AREA>0 guard on every eligible step even without a sample.

Frozen `case20m/input/run_data.txt` lines 17–27 specify CDF, ISKILL=999999, EPTS=VPTS=FPTS=0 and JTM=0. Thus this special case avoids station gathers and sampled energy folds altogether; those are implemented because a different skill interval/station list is not equivalent. Window interval is 1200 steps, i.e. four writes per 4800 steps. The historical OpenMP v30 32-thread profile measured roughly 0.999 s/4799 ARCHIVE calls of 12.691 s total. This identifies recurring work, not a prediction of MPI speed.

## Verification

From the outer workspace:

```sh
python3 -B work/secom-omp/mpi/outer/check_archive_local.py \
  --fixture work/secom-omp64-20261118/v39-local/fixture.f \
  --result /tmp/smpi-archive-proof-new.json
```

`archive-local-proof.json` records the local compiler/flags, source fingerprint and per-case output hashes. The 21×19 fixture compares actual v46 and candidate with the policy disabled and with a synthetic single full-domain owner and no MPI transport. It exercises three history writes, positive and negative certification, inactive toxic bed branches, zero AREA STOP 95, sampled energy, nonempty E/V stations and both cross-section directions. It captures raw payloads in a stub PUTCDF/TSCDF; it does **not** validate MPI messages, real NetCDF, or Linux signal behavior.

The separate native test `check_archive_mpi.py` ran on p3, OpenMPI 4.1.1, GNU Fortran 8.5.0, three ranks × one thread on the specifically allocated physical CPUs 20,21,22. `archive-native-proof.json` records the flags, compiler, driver/runtime SHA and all output hashes; `archive-native-proof.log` is stdout. Non-owner values for current dynamics are poisoned **every** step. Baseline and trial write private per-rank files containing only their authoritative history slices at every step, all 28 *global* history fields at the initial INTX=0 write and three later writes; that capture does not publish or modify solver state. The remote row signed-zero rejection, skill-only I/O including original ordered ESUM/TKE/APE, seam velocity station and both transect directions all matched each same-rank v46 raw file exactly. Omitting CBC publication as a negative control produced three differing rank capture files as expected. This verifies raw data delivered to stub PUTCDF/TSCDF, not NetCDF libraries or the whole model. The native archive test changes no solver file. Independent application-scale serial output and trace comparisons are still required. A standalone native rerun from the repo has the following shape (substitute the CPUs assigned to that test, do not reuse an occupied allocation):

```sh
module load mpi/openmpi-x86_64
taskset -c 20-22 python3 -B mpi/outer/check_archive_mpi.py \
  --fixture /path/to/frozen/v39-local/fixture.f \
  --ranks 3 --threads 1 --cpus 20,21,22 --result /tmp/archive-native-new.json
```

On a source-only installation, `--reference /path/to/v46/source/archive.f --source /path/to/r004/source` replaces the MPI harness's Git lookup (original source and fixture SHA are checked). The historical verification fixture `work/secom-omp64-20261118/v39-local/fixture.f` is included here as `fixture-v39.f`, the default fixture if `--fixture` is omitted (SHA-256 `662bd6e05c16448c37dfc36941f40e7aba5ef66e0320d1d0b4e678fb159f363f`). `--smoke --fc gfortran` checks the driver with a single stub rank without claiming transport coverage. The WAD seam unit's independent p3 evidence is in `validationmpi/evidence/wad/`, not part of the archive proof.

Before promoting a change to an event snapshot or the default policy, also run from repo root:

```sh
python3 -B mpi/outer/check_snapshot_shapes.py
python3 -B mpi/outer/check_default_modes.py
```

The first examines GNU Fortran's actual front-end symbol table, including declaration, REAL kind, lower/upper extents, contiguous K slices and admitting only declared COMMON or argument planes; an unrecognized `SMPI_ALL_RAW` call fails closed. It rejects the erroneous implicitly scalar `PRECTEMP`, REAL*8 `ARCTHIK`, and K starting at zero (see `snapshot-shapes-r004b-strict.json`). The second links the production serial stub and archive as separate objects in non-OpenMP and OMP-3 builds. Without setting the policy, it verifies its one default initializer before/after `SMPI_INIT` and from every worker. Actual file bytes in three focused modes match original v46 serial in all four combinations (`default-modes-r004b.json`). This is not a combined MPI+multi-thread full-solver test.
