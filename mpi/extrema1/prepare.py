"""Overlay only on an extrema2 source tree. The previous verdict remains the fallback."""
import argparse
import hashlib
import json
from pathlib import Path


def build(src):
    here = Path(__file__).resolve().parent
    b = (src/'baro_check.f').read_text(encoding='latin1')
    r = (src/'secom_mpi.f90').read_text(encoding='latin1')
    old = '        CALL SMPI_VERDICT2(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN,OLDTAIL2)'
    if b.count(old) != 1 or 'SMPI_TRY_HEALTH' in b or 'subroutine smpi_checkkeys_raw' in r:
        raise ValueError('requires a fresh extrema2 parent')
    b = b.replace(old, '''        CALL SMPI_TRY_HEALTH(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN,OLDTAIL2)
        IF(OLDTAIL2) THEN
''' + old + '''
        END IF''')
    return {'baro_check.f': b + (here/'health.f.inc').read_text(),
            'secom_mpi.f90': r + (here/'keys.f90.inc').read_text()}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    files = build(a.source)
    a.out.mkdir()
    meta = {}
    for n, s in files.items():
        data = s.encode('latin1')
        with (a.out/n).open('xb') as f:
            f.write(data)
        meta[n] = {'parent': hashlib.sha256((a.source/n).read_bytes()).hexdigest(),
                   'candidate': hashlib.sha256(data).hexdigest()}
    with (a.out/'extrema1.json').open('x') as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta, sort_keys=True))
