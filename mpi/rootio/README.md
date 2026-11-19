# r003 physical output ownership candidate

Do not apply this to an existing run/build. Produce two files in a fresh
directory, overlay them on a new source snapshot, keep all compiler flags,
objects, and solver input unchanged. No runtime options or MPI calls added.

```sh
python3 -B mpi/rootio/prepare.py --source source --out rootio-new
python3 -B mpi/rootio/audit.py --source source --candidate rootio-new --compiler gfortran
```

The output rank is always rank 0, including serial and MPI size 1. Other
ranks still execute PUTCDF's first-call setup and its original control flow:
it sets ten COMMON print flags, calls RDCORLOC (shared grid coordinates and
gcmprt), computes dates, packs arrays and advances its saved ICDF. They
skip only netCDF output functions/procedures, fort.1002/1006 writes and
their formatted-output helper. Private netCDF IDs that would be returned
to a nonwriter are zeroed; they are not used for I/O. No netCDF handle is
shared between processes. The unchanged wind-text helper only reads its
caller-supplied arrays; formatting-only arithmetic is no longer performed
on peers. Restart time/counter changes, internal filenames and branches are
unchanged: only IUWRS OPEN (4 source locations), WRITE (24) and CLOSE (2)
are guarded. Archive accumulation, diagnostics, input handling and temporary
file cleanup are not modified.

The audit uses the GNU Fortran front end, not textual grep: for both files
the root-specialized executable IR must match the parent. No targeted
physical sinks may remain on the peer path, and the ten direct COMMON writes,
RDCORLOC, ICDF and date calls must remain. The source digest and code-count
guards reject unexpected PUTCDF revisions. This is a structural check, not
a replacement for the real early/late/final scientific comparisons.

For the 24-input cold/barotropic campaign, in each clean **separate rank
directory**, rank 0 must produce all six historical files with the same-day
serial SHA-256: `fort.1002`, `fort.1006`, `gcmplt.cdf`, `startup`, `gcmprt`,
`heat.dat`. In an uninstrumented run, each nonroot may contain exactly the
24 unchanged inputs, launch metadata `launch.json` and `stdout.log`, and
only two scientific outputs: `gcmprt` and `heat.dat` (zero length for the
frozen case). Compare these two to the serial proof as well, unless a later
separately reviewed change expressly alters their scope. Reject any
`startup*`, `fort.*`, `gcmplt.cdf` on nonroot, and unexpected leftover
scratch. Trace runs may additionally contain only their declared frame files
on their designated rank(s); do not relax the six-file gate on root.
Old prototype runners expecting six outputs *on every rank* will correctly
mark this layout as different; use a new gate, keep historical runs intact.

Do not share run directories or symlink writable scratch across ranks:
legacy BCDATA/FIRST create, rewind and reopen rank-local temporary inputs,
and ECOM removes `gcm_temp*` in its current directory. Root-only output is
not root-only initialization and is not an MPI scaling claim.
