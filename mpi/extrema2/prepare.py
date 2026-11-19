"""Overlay for the two-message extrema verdict; do not edit the parent source."""
import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def once(s, old, new):
    if s.count(old) != 1:
        raise ValueError('missing/ambiguous parent anchor: ' + old[:72])
    return s.replace(old, new)


def build(src):
    b = (src / 'baro_check.f').read_text(encoding='latin1')
    r = (src / 'secom_mpi.f90').read_text(encoding='latin1')
    if 'SMPI_VERDICT2' in b.upper() or 'SMPI_CHECKREC_RAW' in r.upper():
        raise ValueError('already applied')
    b = once(b, '! Wet/dry owns and commits', '      LOGICAL OLDTAIL2\n! Wet/dry owns and commits')
    anchor = '      IF(SMPI_SIZE.GT.1) CALL SMPI_ROWS(ROWVMAX,ROWEMIN,JM)'
    b = once(b, anchor, '''      IF(SMPI_SIZE.GT.1) THEN
!$OMP SINGLE
        CALL SMPI_VERDICT2(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN,OLDTAIL2)
!$OMP END SINGLE COPYPRIVATE(OLDTAIL2)
        IF(.NOT.OLDTAIL2) RETURN
      END IF
''' + anchor)
    return {'baro_check.f': b + (HERE / 'verdict.f.inc').read_text(),
            'secom_mpi.f90': r + (HERE / 'collect.f90.inc').read_text()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    files = build(args.source)
    args.out.mkdir()
    manifest = {}
    for name, data in files.items():
        data = data.encode('latin1')
        with (args.out / name).open('xb') as f:
            f.write(data)
        manifest[name] = {'parent_sha256': hashlib.sha256((args.source/name).read_bytes()).hexdigest(),
                          'candidate_sha256': hashlib.sha256(data).hexdigest()}
    with (args.out / 'extrema2.json').open('x') as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps(manifest, sort_keys=True))


if __name__ == '__main__':
    main()
