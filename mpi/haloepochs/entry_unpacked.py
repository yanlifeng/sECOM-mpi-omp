"""Experimental same directional/epoch contract, separate plane messages."""
from pathlib import Path
import argparse
import hashlib


def build(src):
    path = Path(src) / 'secom_mpi.f90'
    text = path.read_bytes().decode('latin1')
    anchor = '  if(im>huge(im)/6)then\n'
    if text.count(anchor) != 1:
        raise ValueError('expected exactly one directional entry fallback')
    # INIT/dims have already checked IM>0. Use the exact existing
    # per-plane directional path, never the old bidirectional HALO.
    text = text.replace(anchor, '  ! Deliberate unpacked experiment: positive grid always takes this path.\n  if(im>0)then\n')
    return {'secom_mpi.f90': text}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    files = build(a.source)
    a.out.mkdir()
    for name, text in files.items():
        b = text.encode('latin1')
        (a.out / name).write_bytes(b)
        print(hashlib.sha256(b).hexdigest(), name)
