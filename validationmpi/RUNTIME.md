# Slab runtime contract

In `source`, use `make -f Makefile.portable MODE=mpi MPIFC=mpifort`
with the same strict floating-point, architecture and netCDF settings as
the serial control. MPI mode links OpenMP as well. The serial/omp modes
compile the same runtime file without MPI headers, MPI libraries or
additional OpenMP barriers. Build each mode in a different OBJDIR.
The include file `mpideck` is compatible with fixed and free form;
do not include it twice if the scientific routine already includes comdeck.

Call `SMPI_INIT(im,jm)` and `SMPI_FINALIZE()` outside parallel regions,
in the initial thread on all ranks. A private duplicated communicator
isolates the halo tags. If the embedding application already initialized
MPI, it retains ownership of finalization; it must have requested at least
MPI_THREAD_SERIALIZED. Do not call from independently concurrent/nested
teams, on a subset of the ranks or from differently ordered call paths.
All ranks must use the same positive dimensions with IM×JM in range of
default INTEGER, a homogeneous default REAL32/INTEGER32 build and
no more ranks than rows. Contract violations or MPI errors are fatal,
not recoverable fallbacks to a local solve.

All planes still have full global (IM,JM) shape, I contiguous. **Only owned
rows are current** unless explicitly exchanged. Each rank owns an ascending
interval of global J, inclusive SMPI_JLO..SMPI_JHI. SMPI_RANK is zero-based,
SMPI_SIZE is the rank count. Before INIT and after FINALIZE defaults are
rank 0, size 1, bounds 1..HUGE; normal pre-init legacy loop bounds work.
The read-only functions SMPI_FIRST(lo), SMPI_LAST(hi), SMPI_OWNS(j) and
SMPI_COLOR(lo,hi,positive_stride) intersect the original loop with ownership.
COLOR preserves the original global residue class or returns >hi if empty.
SMPI_OWNER_RAW(j) is a read-only noncollective integer query (zero-based
owner rank) valid only after INIT for global 1..JM. SMPI_LOCALBC is an
independent COMMON flag, false by default and reset at INIT/FINALIZE.

| Interface | Effect |
| --- | --- |
| SMPI_HALO(A,IM,JM) | Adjacent rank on each side gets one boundary row. Only existing immediate ghost rows are received; not a full gather. |
| SMPI_ALL(A,IM,JM) | In-place Allgatherv of owned plane rows; replicated full result. |
| SMPI_ROWS(vmax,emin,JM) / SMPI_IROWS(imax,imin,JM) | In-place Allgatherv of two **real** or **integer** per-row arrays. Does not recompute or reduce any values. |
| SMPI_INDEX(i) | Integer max: the single owner supplies the positive index, all peers supply 0. It is the caller's job to identify the unique owner. |
| SMPI_RBCAST_ROW(x,j) | Broadcast the owner rank's default REAL scalar as four bytes, preserving e.g. negative zero. |
| SMPI_WAD_PUBLISH(flags,et,etf,va,vaf,dvmwad,etb,im,jm) | Send IM integer event flags north; a receiver applies them only at its first owned row. Bit 0: force VA/VAF/DVMWAD to +0; bit 1: ET=ETF=local ETB. Unmarked entries are untouched. Not a halo refresh. |

All seven have a matching `_RAW` suffix: call raw **once per rank** from a
pre-existing SINGLE/serial context, only when producer data are already
published and the caller provides the subsequent team join. The ordinary
MPI wrappers require **every thread of a team** to enter once with the same
shared buffers/scalars. They perform an entry barrier then one MPI-calling
thread, and an exit join. Do not call a wrapper from within SINGLE/MASTER,
worksharing loops, tasks or a thread-conditional branch. Outside a parallel
region the orphaned directives have one participant. `MPI_THREAD_SERIALIZED`
permits whichever worker executes SINGLE to call MPI, not simultaneous MPI
calls. MPI builds retain those joins even at communicator size 1. Ordinary
serial/OpenMP stubs have no collective semantics and no extra joins.

## p3 microtest, not a solver benchmark

Use the OS-supplied Open MPI 4.1.1 with GNU Fortran 8.5:
`module load mpi/openmpi-x86_64`. In a unique scratch directory place
`secom_mpi.f90`, `mpideck`, `check_runtime.f90`, `runtime_place.c` and
`runtime_rank.sh`. Run all commands under the same allowed CPU mask
(the fixture uses **24–31**, within the campaign's 4–35 budget):

```sh
mpicc -O2 -c runtime_place.c
mpifort -cpp -DSECOM_USE_MPI -fopenmp -fcheck=all -O2 -g -c secom_mpi.f90
mpifort -fopenmp -fcheck=all -O2 check_runtime.f90 secom_mpi.o runtime_place.o -o check-runtime
for args in '5 1 1' '5 2 1' '7 4 1' '4 4 1' '7 2 2' '5 1 2'; do
  set -- $args  # JM, ranks, OpenMP lanes per rank
  timeout 25s mpirun --host localhost:4 -np "$2" --nooversubscribe --bind-to none \
    --mca plm isolated --mca pml ob1 --mca btl self,vader \
    sh runtime_rank.sh "$@"
done
```

The wrapper pins before exec and the compiled test checks each application
thread's actual singleton affinity and execution CPU. It tests full-size and
uneven/one-row slabs, exterior boundaries, repeated halo and owner gathers,
an MPI-calling OpenMP SINGLE, signed-zero broadcast, integer index reduction,
and selected/unselected seam updates. The runtime reports the provided thread
level. The fixture is not evidence that real solver ranks make equal progress
or that their internal MPI runtime has no auxiliary threads.

Remote original evidence: p3
`/home/user_home/ylf/secom-mpi-20261119/runtime-unit-r002/test.log`,
SHA-256 `6433278ba2a6b814b6880671179a19746afd8df2bd019a8ce94da4af0313252f`.
All six MPI cases passed. The tested source hash was
`14f5f37b412121f3e2efaab7b28834c200df1fad94bf3c8e04edb8887bbc325b`;
the subsequent non-MPI-only guard removes stub OpenMP synchronization
without changing MPI-mode statements.

The r002 solver source revision, including the LOCALBC flag and row-owner
query, was then tested in `runtime-unit-r003` on the same host and mask.
All six cases passed; each rank additionally checks every row owner, and
tests the default, initialization and finalization flag. The complete local
log copies are `runtime-unit-r002.log` and `runtime-unit-r003.log`.
The latter SHA-256 is
`0ebc0757f534d87df616f2ea618c3d0dedda5daaaddf3b6f8ea7f85023b8252b`;
it records the exact inputs and binary hashes and application-thread CPUs.
