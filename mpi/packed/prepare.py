"""Generate the optional nine-plane MPI halo overlay; never write to the parent."""
import argparse
import hashlib
from pathlib import Path


def replace_once(txt, needle, new):
    if txt.count(needle) != 1:
        raise ValueError('expected exactly one parent anchor: ' + repr(needle[:100]))
    return txt.replace(needle, new)


def build(src):
    rt = (src / 'secom_mpi.f90').read_bytes().decode('latin1')
    if 'smpi_halo9' in rt.lower():
        raise ValueError('parent already has packed halo')
    rt = replace_once(rt, '  integer, save :: south,north\n',
                      '  integer, save :: south,north\n  real, allocatable, save :: send9(:,:),recv9(:,:)\n')
    rt = replace_once(rt, '  do k=0,smpi_size-1\n',
                      '''  ! A uniform oversized count uses the existing separate path, not a partial message.
  if(smpi_size>1.and.im<=huge(im)/9)then
    allocate(send9(im,9),recv9(im,9),stat=ierr)
    if(ierr/=0)call bad('cannot allocate nine-plane halo buffers')
  endif
  do k=0,smpi_size-1
''')
    rt = replace_once(rt, '  deallocate(rows,offs,cells,starts,wad_recv)\n',
                      '  deallocate(rows,offs,cells,starts,wad_recv)\n  if(allocated(send9))deallocate(send9,recv9)\n')
    rt += Path(__file__).with_name('halo9.f90.inc').read_text()
    ext = (src / 'extrnl.f').read_bytes().decode('latin1')
    planes = 'D UA VA UAB VAB EL ELB DUMWAD DVMWAD'.split()
    old = '      IF(SMPI_LOCALBC) THEN\n' + ''.join('      CALL SMPI_HALO(%s,IM,JM)\n' % p for p in planes) + '      ELSE\n'
    new = '''      IF(SMPI_LOCALBC) THEN
! All nine input snapshots share the same owner/ghost boundary. No
! intermediate consumer observes any of the individual transfers.
      CALL SMPI_HALO9(D,UA,VA,UAB,VAB,EL,ELB,DUMWAD,DVMWAD,IM,JM)
      ELSE
'''
    return {'secom_mpi.f90': rt, 'extrnl.f': replace_once(ext, old, new)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    files = build(a.source)
    a.out.mkdir()
    for name, text in files.items():
        data = text.encode('latin1')
        with (a.out / name).open('xb') as f:
            f.write(data)
        print(hashlib.sha256(data).hexdigest(), name)


if __name__ == '__main__':
    main()
