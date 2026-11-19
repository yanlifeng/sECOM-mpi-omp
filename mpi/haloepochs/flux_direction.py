"""Generate only three directed FLUXVA handoffs; preserve public/bidirectional fallback."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def _live(text):
    return '\n'.join(s.split('!', 1)[0].upper() for s in text.splitlines()
                     if s and s[0] not in 'cC*!')


def _verify(text, label, offsets, calls):
    # Fail closed if the consumer starts reading either the opposite ghost
    # or via a new alias. Header declarations live in the untouched comdeck.
    s = _live(text)
    halo = r'CALL\s+SMPI_HALO\(FLUXVA,IM,JM\)'
    if len(re.findall(halo, s)) != calls:
        raise ValueError(label + ': unexpected FLUXVA communication')
    s = re.sub(halo, '', s)
    expr = r'\bFLUXVA\(I,J([+-]\d+)?\)'
    seen = [v for v in re.findall(expr, s) if v]
    if seen != offsets or re.search(r'\bFLUXVA\b', re.sub(expr, '', s)):
        raise ValueError(label + ': FLUXVA aliases or nonlocal row accesses changed')


def build(src):
    src = Path(src)
    names = ('extrnl.f', 'advave.f', 'secom_mpi.f90')
    ext, adv, rt = ((src / n).read_bytes().decode('latin1') for n in names)
    if 'smpi_flux_one_raw' in rt.lower():
        raise ValueError('directed flux already installed')
    _verify(ext, 'EXTRNL_TEAM', ['+1'], 1)
    _verify(adv, 'SECOM_ADV34_TEAM', ['+1', '-1'], 2)
    if ('CALL SECOM_ADV34_TEAM(ADVUA,ADVVA,.FALSE.)' not in adv or
            "EDGE_LATE=ROTATE.AND.TOR.EQ.'BAROTROPIC'.AND.BCTYPE.EQ.'CLAMPED'" not in ext):
        raise ValueError('public/private ADV entry contract changed')
    old = '      CALL SMPI_HALO(FLUXVA,IM,JM)\n'
    if ext.count(old) != 1:
        raise ValueError('cannot identify EXTRNL flux call')
    st = '''      IF(SMPI_LOCALBC.AND.ROTATE) THEN
        CALL SMPI_FLUX_HALO_N(FLUXVA,IM,JM)
      ELSE
        CALL SMPI_HALO(FLUXVA,IM,JM)
      END IF
'''
    ext = ext.replace(old, st)
    old = '      IF(SMPI_SIZE.GT.1) CALL SMPI_HALO(FLUXVA,IM,JM)\n'
    parts = adv.split(old)
    if len(parts) != 3 or any('DO %d J=' % n not in segment
                               for n, segment in zip((480, 880), parts[1:])):
        raise ValueError('cannot identify two distinct ADV flux generations')
    st = '''      IF(SMPI_LOCALBC.AND.EDGE_LATE) THEN
        CALL SMPI_FLUX_HALO_%s(FLUXVA,IM,JM)
      ELSE IF(SMPI_SIZE.GT.1) THEN
        CALL SMPI_HALO(FLUXVA,IM,JM)
      END IF
'''
    adv = parts[0] + st % 'N' + parts[1] + st % 'S' + parts[2]
    rt += Path(__file__).with_name('flux_direction.f90.inc').read_text()
    return dict(zip(names, (ext, adv, rt)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    files = build(a.source)
    a.out.mkdir()
    record = {}
    for n, s in files.items():
        data = s.encode('latin1')
        with (a.out / n).open('xb') as f:
            f.write(data)
        record[n] = {'parent': hashlib.sha256((a.source / n).read_bytes()).hexdigest(),
                     'candidate': hashlib.sha256(data).hexdigest()}
    print(json.dumps(record, sort_keys=True))


if __name__ == '__main__':
    main()
