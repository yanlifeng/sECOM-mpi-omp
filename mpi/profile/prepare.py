"""Generate, or independently relink, the diagnostic-only profile snapshot."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess


def edit(txt, pat, sub, n=1):
    txt, count = re.subn(pat, sub, txt, flags=re.M)
    if count != n:
        raise ValueError(f'profile anchor expected {n}, found {count}: {pat}')
    return txt


def surround(t, pat, id, n=1):
    return edit(t, pat, lambda m: f'      CALL SD_SITE({id},1)\n' + m[0] + f'\n      CALL SD_SITE({id},0)', n)


def mainf(t):
    t = edit(t, r'^      CALL SMPI_INIT\(IM,JM\)$',
             '      CALL SD_STAGE(1)\n\\g<0>\n      CALL SD_ID(SMPI_RANK,SMPI_SIZE)\n      CALL SD_STAGE(2)')
    t = edit(t, r'^ +DO 9000 INTX=ISTART,IEND$', '      CALL SD_STAGE(3)\n\\g<0>\n      CALL SD_STEP()')
    t = edit(t, r'^ 9000 +CONTINUE$', '\\g<0>\n      CALL SD_STAGE(4)')
    t = edit(t, r'^      CALL SMPI_FINALIZE\(\)\n      STOP\b',
             '      CALL SD_STAGE(5)\n      CALL SMPI_FINALIZE()\n      CALL SD_DUMP()\n      STOP')
    t = surround(t, r'^      CALL BCOND\(7,DTI2,0\)$', 1)
    t = surround(t, r'^ +CALL UPDATE_CBCM?\b.*$', 2, 4)  # startup and loop separated by stage
    t = surround(t, r'^      CALL BOTTAU\b.*$', 3)
    t = surround(t, r'^ +CALL (?:ARCHIVE\b|SECOM_ARCHIVE_BODY\().*$', 4, 4)
    t = surround(t, r'^ +CALL BARO_SUBCYCLE\([^\n]*\n     &[^\n]*\)$', 5)
    for pat in [r'^      IF \(INTX.EQ.IRESTART\) THEN$', r'^      IF \(IEND.LE.IRESTART.OR.NRESEND.EQ.0\) THEN$']:
        t = edit(t, pat, '\\g<0>\n      CALL SD_SITE(6,1)')
    t = edit(t, r"^(      WRITE\(\*,\*\) 'RESTART FILE FINISHED,[^\n]*\n)(      endif)$",
             '\\1      CALL SD_SITE(6,0)\n\\2')
    t = edit(t, r'^(      ENDIF\n!NG\n! *\n 9100 CONTINUE)$', '      CALL SD_SITE(6,0)\n\\1')
    # Conditional statement must stay conditional; startup/abnormal print paths
    # remain inside their lifetime or are stopped before a diagnostic is accepted.
    t = edit(t, r"^      IF\(RESTAR.EQ.'COLD START'\) CALL (PRINTS\([^\n]+\))$",
             "      IF(RESTAR.EQ.'COLD START') THEN\n      CALL \\1\n      ENDIF")
    t = surround(t, r'^ +CALL PRINTS\(DRHOX,DRHOY,TRNU,TRNV\)$', 8, 4)
    t = edit(t, r"^      IF\(HORZMIX.EQ.'CLOSURE   '\) CALL SMAG[^\n]*$",
             "      IF(HORZMIX.EQ.'CLOSURE   ') THEN\n      CALL SD_SITE(14,1)\n      CALL SMAG\n"
             "      CALL SD_SITE(14,0)\n      ENDIF")
    t = surround(t, r'^        CALL EXTRNL_TEAM\([^\n]+$', 10)
    t = surround(t, r'^        CALL WANDD_TEAM\(EL,D\)$', 11)
    t = surround(t, r'^        CALL SECOM_CHECK_TEAM\([^\n]+$', 12)
    return t


def archive(t):
    for pat,id in [(r'^ +CALL PUTCDF\(TMIDDLE\)$',7),(r'^ +CALL TSCDF\(TMIDDLE,icdfcheck\)$',9)]:
        t=surround(t,pat,id)
    return t


def runtime(t):
    kinds = {'halo':1,'all':2,'rows':3,'irows':4,'index':5,'rbcast_row':6,'wad_publish':7,'halo9':8}
    got, out, raw = {}, [], None
    for line in t.splitlines(keepends=True):
        if line.startswith('subroutine ') or line.startswith('end subroutine'):
            m = re.match(r'subroutine smpi_(\w+)_raw\(', line)
            raw = m[1] if m else None
        if raw and re.match(r'\s*call MPI_\w+\(', line):
            if raw not in kinds:
                raise ValueError('unclassified RAW communication: ' + raw)
            id = kinds[raw]
            got[raw] = got.get(raw,0) + 1
            # checked(ierr) belongs outside the actual MPI call timer.
            call, sep, checked = line.rstrip('\n').partition(';call checked(')
            out.extend([f'  call sd_mp({id},1)\n', call + '\n', f'  call sd_mp({id},0)\n'])
            if sep:
                out.append('  call checked(' + checked + '\n')
        else:
            out.append(line)
    exp = {'halo':4,'all':1,'rows':2,'irows':2,'index':1,'rbcast_row':1,'wad_publish':1}
    if 'subroutine smpi_halo9_raw(' in t:
        exp['halo9']=1
    if got != exp:
        raise ValueError('unexpected runtime callsites: ' + repr(got))
    return ''.join(out) + Path(__file__).with_name('timers.inc').read_text()


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--build',action='store_true',help='p3 native relink with the existing parent science objects')
    ap.add_argument('--fc',default='/usr/lib64/openmpi/bin/mpifort')
    ap.add_argument('--nf-config',default='/home/user_home/ylf/kt5-licom-secom-ylf/deps/bin/nf-config')
    a=ap.parse_args()
    s,d=a.source.resolve(),a.out.resolve()
    methods={'ecom3d.f':mainf, 'archive.f':archive, 'secom_mpi.f90':runtime}
    rendered={n:fn((s/n).read_bytes().decode('latin1')).encode('latin1') for n,fn in methods.items()}
    d.mkdir()
    for n,b in rendered.items():
        with (d/n).open('xb') as f: f.write(b)
    proof={'kind':'diagnostic_not_eligible_for_timing','parent':{n:sha(s/n) for n in rendered},
           'generated':{n:sha(d/n) for n in rendered},'commands':[]}
    if a.build:
        flags=shlex.split('-O3 -fno-fast-math -ffp-contract=off -std=legacy -ffixed-line-length-132 '
                          '-g -fbacktrace -ffpe-trap=invalid,zero,overflow -march=native -mprefer-vector-width=256 -fopenmp')
        conf=lambda k:shlex.split(subprocess.check_output([a.nf_config,k],universal_newlines=True))
        inc,libs=conf('--fflags'),conf('--flibs')
        objs=sorted((s/'build/mpi').glob('*.o'))
        swaps={Path(n).with_suffix('.o').name for n in rendered}
        if len(objs)<90 or not swaps<={p.name for p in objs}:raise ValueError('incomplete parent object set')
        proof['reused_objects']={p.name:sha(p) for p in objs if p.name not in swaps}
        proof['parent_elf_sha']=sha(s/'build/mpi/s--pom.x')
        for n in ['secom_mpi.f90','ecom3d.f','archive.f']:
            extra=['-cpp','-DSECOM_USE_MPI','-J'+str(d)] if n.endswith('.f90') else []
            proof['commands'].append([a.fc,*flags,*inc,'-I'+str(s),*extra,'-c',str(d/n),'-o',str(d/Path(n).with_suffix('.o'))])
        proof['commands'].append([a.fc,*flags,'-o',str(d/'s--pom-profile.x'),
                                  *[str(d/p.name if p.name in swaps else p) for p in objs],*libs])
        with (d/'build.log').open('x') as log:
            for cmd in proof['commands']:
                print(' '.join(shlex.quote(w) for w in cmd),file=log,flush=True)
                subprocess.run(cmd,cwd=d,stdout=log,stderr=subprocess.STDOUT,check=True)
        proof['binary_sha']=sha(d/'s--pom-profile.x')
    with (d/'profile-build.json').open('x') as f: json.dump(proof,f,indent=2)


if __name__=='__main__':
    main()
