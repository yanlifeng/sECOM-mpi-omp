# Optional nine-field halo (not a measured solver improvement)

Generate the two-file overlay against the **chosen parent science snapshot**:

```sh
python3 -B mpi/packed/prepare.py --source source --out packed-new
```

Use both generated files, `extrnl.f` and `secom_mpi.f90`, in a fresh build.
No mpideck, outer-kernel, I/O, flags or legacy boundary fallback edits.
The source generator requires the nine uninterrupted HALOs in the
SMPI_LOCALBC arm and rejects repeated application or an unexpected layout.
Preserve later parent edits by regenerating; do not transplant the r002 QA
runtime binary or an earlier generated source over a new runtime.

Interface: `SMPI_HALO9(A,B,C,D,E,F,G,H,V,IM,JM)` and identical
`SMPI_HALO9_RAW`. The caller supplies nine distinct contiguous default REAL
planes. The public wrapper requires all threads, with an entry barrier and a
SINGLE exit barrier; RAW is serialized by its caller as for other SMPI RAW
operations. MPI p=1 retains the public joins; non-MPI builds have neither
communication nor new joins. The other eight intermediate entry/exit joins
are unnecessary at this exact call site: the arrays are distinct and no
work or consumer occurs between the original exchanges.

At INIT, MPI size>1 creates two reusable IM×9 REAL buffers (26,784 bytes
total for IM=372); at FINALIZE it frees them. Low rows travel south, high
rows north; their recipient updates only the high/low adjacent ghost when
that neighbour exists. Each direction uses one blocking MPI_Sendrecv,
private communicator, tags 19 and 20, and no overlap of receive/send buffers.
Single-row slabs pack the same authoritative row for both neighbours.
Dummy PROC_NULL receives never index outside the plane. If 9×IM cannot fit
the MPI INTEGER count all ranks use nine original RAW calls instead. An
allocation failure aborts; independent per-rank fallbacks could deadlock.

## Verification

On p3 in a fresh scratch directory put the two generated files, the parent
`mpideck`, `prepare.py`, `halo9.f90.inc`, `check_packed.f90`,
`native_unit.sh`, and `validationmpi/{runtime_place.c,runtime_rank.sh}`.
With an approved CPU24–31 window, use `taskset -c 24-31 bash native_unit.sh`.
The ordinary launcher uses explicit singleton OpenMP places, and a C probe
checks the actual application lane masks/CPUs.

The [p3 test log](p3-packed-r002.log) records the tested source hashes and
seven passes (1/2/4 ranks, OMP 1/2, uneven and one-row slabs). The fixture
checks 4,901,472 field words against **both** nine ordinary HALOs and an
independent expected-bit oracle, including signed zeros, signaling/quiet
NaN payloads, subnormals and poisoned far rows, across twelve replacements.
It exercises both public and RAW/single entries under strict FP traps and
bounds checks. The earlier [compile-only log](p3-packed-r001-compile.log)
records GCC8 rejecting an out-of-range hex constructor in the fixture;
the portable IBSET fix did not change the candidate runtime or call site.

The per-direction payload is 13,392 bytes instead of nine 1,488-byte
messages at the production width. This trades latency for local packing
and can cross MPI transport protocol thresholds. No bandwidth, overlap or
whole-solver speedup is established. For acceptance, run the established
scientific output/field gates and independently interleave fresh real
32-rank solver comparisons against the same parent. A separate communication
microbenchmark, if required, should synchronize outside each timed batch,
alternate methods, report each batch's slowest rank and never count that
isolated result as model throughput.
