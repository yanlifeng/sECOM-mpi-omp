"""Remove the unused packed branch after selecting directional per-plane entry."""
from pathlib import Path
import argparse
import hashlib


def build(src):
    text = (Path(src) / 'secom_mpi.f90').read_bytes().decode('latin1')
    a = text.index('module smpi_entry_workspace\n')
    b = text.index('subroutine smpi_entry_epoch(d,', a)
    text = text[:a] + text[b:]
    if text.count('  call smpi_entry_epoch_release()\n') != 1:
        raise ValueError('unexpected entry workspace release')
    text = text.replace('  call smpi_entry_epoch_release()\n', '')
    a = text.index('subroutine smpi_entry_epoch_raw(')
    b = text.index('\ncontains\n', a)
    part = text[a:b]
    if part.count('  if(im>0)then\n') != 1:
        raise ValueError('requires already selected unpacked entry variant')
    part = part.replace('  use smpi_entry_workspace\n', '')
    part = part.replace('  integer::lo,hi,ji,ns,nr,rc,nf\n', '  integer::lo,hi\n')
    head = part.split('  ! Preserve exact *directional*')[0]
    body = part.split('  if(im>0)then\n', 1)[1].split('    return\n  endif', 1)[0]
    body = '\n'.join(line[2:] if line.startswith('  ') else line for line in body.split('\n'))
    part = head + '  ! Only needed directions and live field versions are sent.\n' + body
    return {'secom_mpi.f90': text[:a] + part + text[b:]}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    files = build(a.source)
    a.out.mkdir()
    for n, s in files.items():
        b = s.encode('latin1')
        (a.out / n).write_bytes(b)
        print(hashlib.sha256(b).hexdigest(), n)
