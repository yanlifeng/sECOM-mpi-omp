# Post-boundary CURV42D and ELF: receive the south ghost only

`python3 -B momentum_south.py --source /frozen/r013/source --out /fresh/overlay`
emits only `extrnl.f`. It reuses the existing one-message
`SMPI_FLUX_HALO_S`; operator, wire datatype, OMP joins and public
bidirectional path are unchanged. Only SMPI_LOCALBC + ROTATE +
SMPI_OUTER enters the new branch. This removes two redundant north-ghost
Sendrecv calls per substep (96,000 over the frozen 48,000-substep event).

Audit: the only declaration of /BLK2D/ is in comdeck. CURV42D is the
full plane after COR, and ELF is the distinct full plane between VAB
and EL. Neither occurs in an EQUIVALENCE. At this epoch, the 310
CURV42D copies (including JE+1) and ordered hard corners already
finished; the corner routine modifies only UA/VA. Momentum U 425 reads
CURV42D and ELF at J (including I-1); V 435 reads J and J-1. No
J+1 consumer follows. BCOND(13)'s sole ELF reference is in the zero-trip
NUMQBC loop; the admitted LINEAR boundary loop only zeros UAF/VAF based
on static FSM. BCOND(13) never reads CURV42D. ROTATE copies/filters
ELF at owned J only. Subsequent WAD receives EL and D, not ELF, and has
no direct reference to either audited plane; CHECK and both verdict
paths read ELF only on owned rows. Full snapshots and captured frames
publish owners before use. The next substep overwrites the interior
owner values and refreshes its south ghosts at the same boundary.

Physics microtest: 372x19, two outer steps x four substeps, changing
advection, real EXTRNL (including its two actual BCOND bodies), WAD and
CHECK; 24 captured owner planes + scalar verdict. Before each update,
poison the south ghosts with distinct wrong finite values and the north
ghosts with distinct signaling-NaN payloads. The selected path must
leave the north payloads untouched. Changing just CURV or ELF to the
opposite direction must produce a scientific owner mismatch. For the
normal ROTATE path the oracle is one rank. An artificial, out-of-scope
public ROTATE=false -> WAD sequence has a pre-existing r013 2-rank vs
1-rank VAF discrepancy; it is measured separately, comparing this
candidate to the unchanged 2-rank parent. This is not a claim of
generic public MPI equivalence or whole-event performance.
