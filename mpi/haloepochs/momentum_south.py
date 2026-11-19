"""Publish only the consumed south CURV42D/ELF ghosts in private rotating external mode."""
import argparse
import hashlib
from pathlib import Path
import re


def live(text):
    return '\n'.join(s.split('!', 1)[0].upper() for s in text.splitlines()
                     if s and s[0] not in '*!cC')


def build(src):
    src = Path(src)
    e = (src/'extrnl.f').read_text(encoding='latin1')
    rt = (src/'secom_mpi.f90').read_text(encoding='latin1').lower()
    if 'subroutine smpi_flux_halo_s(' not in rt:
        raise ValueError('requires the already qualified south-ghost helper')
    if 'MOMENTUM_SOUTH:' in e:
        raise ValueError('already applied')
    c = live((src/'comdeck').read_text()).replace(' ', '')
    if c.count('CURV42D(IM,JM)') != 1 or c.count('ELF(IM,JM)') != 1 or not all(
            x in c for x in ('COR(IM,JM),CURV42D(IM,JM)', 'VAB(IM,JM),ELF(IM,JM),EL(IM,JM)')):
        raise ValueError('unexpected BLK2D plane layout')
    tail = live(e.split(' 5000 CONTINUE', 1)[1].split('      RETURN', 1)[0])
    tail = re.sub(r'CALL\s+SMPI_ALL_RAW\(ELF,IM,JM\)', '', tail)
    for field, allowed, south_n in (('CURV42D', {'I,J','I-1,J','I,J-1'}, 1),
                                    ('ELF', {'I,J','I-1,J','I,J-1'}, 2)):
        ids = re.findall(r'\b'+field+r'\(([^()]*)\)', tail)
        if set(ids) != allowed or ids.count('I,J-1') != south_n or tail.count(field) != len(ids):
            raise ValueError(field + ' momentum/rotation consumers changed')
    b = live((src/'bcond.f').read_text(encoding='latin1'))
    b13 = b.split('  20  CONTINUE',1)[1].split('IF(IDX.EQ.13) RETURN',1)[0]
    if b13.count('ELF(') != 1 or 'ELF(IC,JC)' not in b13 or 'DO 200 N=1,NUMQBC' not in b13:
        raise ValueError('BCOND(13) ELF access/zero-discharge gate changed')
    if 'CURV42D' in b13 or "OPTEBC(15:20).EQ.'LINEAR'" not in b13:
        raise ValueError('BCOND(13) curvature access or LINEAR branch changed')
    if re.search(r'\b(ELF|CURV42D)\b', live((src/'wandd_kky.f').read_text(encoding='latin1'))):
        raise ValueError('WAD has a new direct read of an audited field')
    for field in ('CURV42D','ELF'):
        pattern = r'(?m)^([ ]*)CALL SMPI_HALO\('+field+r',IM,JM\)$'
        def sub(m):
            sp = m.group(1)
            return ('! MOMENTUM_SOUTH: only private V momentum reads this field across J.\n'
                    + sp + 'IF(SMPI_LOCALBC.AND.ROTATE.AND.SMPI_OUTER) THEN\n'
                    + sp + '  CALL SMPI_FLUX_HALO_S('+field+',IM,JM)\n'
                    + sp + 'ELSE\n' + sp + '  CALL SMPI_HALO('+field+',IM,JM)\n'
                    + sp + 'END IF')
        e, count = re.subn(pattern, sub, e)
        if count != 1:
            raise ValueError('missing/ambiguous '+field+' midstep publication')
    return e


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    out = build(a.source).encode('latin1')
    a.out.mkdir()
    with (a.out/'extrnl.f').open('xb') as f:
        f.write(out)
    print(hashlib.sha256(out).hexdigest(), 'extrnl.f')
