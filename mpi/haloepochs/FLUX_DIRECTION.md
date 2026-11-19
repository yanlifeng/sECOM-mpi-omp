# Directed FLUXVA exchanges: three generations, not one

Use `python3 -B flux_direction.py --source PARENT --out FRESH`. Only
`extrnl.f`, `advave.f`, and `secom_mpi.f90` are emitted, with their
parent/candidate SHA-256 on stdout. Apply sequentially to the already
combined parent when composing overlapping generated source files.

| Producer -> next consumer | Read beyond owned rows | Send edge / receive ghost |
|---|---|---|
| EXTRNL 405 -> ELF 410 | FLUXVA(i,j+1) | first row south / receive from north |
| ADV34 first generation -> U 480 | FLUXVA(i,j+1) | first row south / receive from north |
| ADV34 second generation -> V 880 | FLUXVA(i,j-1) | last row north / receive from south |

None of the three generations can be reused across a produce/reset. In
the admitted path the neighboring-row reads above are the only active
off-slab references to FLUXVA. Producer-side reads within a row remain
local; face velocity and geometry stencils are distinct input planes and
keep their existing halos. Later diagnostic/output publication gathers
owned rows; it does not treat a discarded ghost as authoritative.
Public EXTRNL (ROTATE=false), public ADVAVE (EDGE_LATE=false), and all
failed SMPI_LOCALBC plans retain the original bidirectional publication.

One MPI_Sendrecv per generation replaces two; OMP entry/join semantics,
wire datatype, vector size, and model expressions remain unchanged. The
opposite ghost is neither read nor updated. Across 48,000 substeps this
removes 144,000 Sendrecv calls per rank; effect on wall time is unmeasured.

Native microtest only, in a fresh directory with a specifically assigned
contiguous four-core subset: `taskset -c LO-HI bash /overlay/flux_native.sh
/snapshot/source`. It verifies all three generations, every owner and
untouched opposite/far ghost byte (including case-specific NaN payloads
and signed zero), reconstructed finite differences, uneven and one-row
slabs, 1/2 OpenMP threads, and a deliberately reversed second generation.
