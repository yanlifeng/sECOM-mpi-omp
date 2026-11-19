"""Verify r014 trace-source sidecar against the immutable v2 bundle; do not compile/run a solver."""
import argparse
import importlib.util
import json
from pathlib import Path
import shlex
import sys
import tempfile
from unittest.mock import patch

import audit_evidence as ref


BUILDER='2e9ca8821ae6589e348fdd44414d73c29d462e17b484bd9f68218a548a333273'
TREE='46f5da3de7d5904ed20e902e1c1ba501f2ca40fa090e11c07a8d7754a6b9d1a4'
ELFS={'early':'58e632e2ff9c001f14572e0dc63ff6f57c8f8ede39cc6a178a3da7983fcac554',
      'late':'8f20fc11e7f6ea0f76f40f73a1cce947cfbc40a5f76dda9af27e9399cbe1089c'}
REMOTE=Path(ref.REMOTE)/'r014final'


def audit(bundle,side):
    src=bundle/'source'; build=bundle/'build'
    generator=bundle/'scripts/build_trace.py'
    ref.need(ref.sha(generator)==BUILDER,'refusing to execute an unpinned generator')
    mapping={f.name:ref.sha(f) for f in src.iterdir() if f.is_file() and not f.is_symlink()}
    ref.need(len(mapping)==107 and ref.canon(mapping)==TREE,'frozen source tree mismatch')
    receipt=ref.read(side/'MANIFEST.json'); remote=receipt['remote_sha256']
    files={m+'/'+n for m in ELFS for n in ('main.f','capture.f','build.log')}
    ref.need(set(remote)==files and receipt['remote_directory']==str(REMOTE) and
             all((side/n).is_file() and not (side/n).is_symlink() and ref.sha(side/n)==remote[n] for n in files),
             'downloaded source/log differs from explicit remote hashes')
    spec=importlib.util.spec_from_file_location('r014_pinned_trace_builder',generator)
    builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
    result={}
    class Rendered(Exception): pass
    with tempfile.TemporaryDirectory() as tmp:
        for mode,expected_elf in ELFS.items():
            dest=Path(tmp)/mode; staged=side/mode; remote_dir=REMOTE/('trace-'+mode)
            with patch.object(sys,'argv',['trace',str(src),str(dest)]+(['--late'] if mode=='late' else [])), \
                 patch.object(builder.subprocess,'check_output',side_effect=Rendered):
                try:builder.main()
                except Rendered:pass  # both sources are now on disk; no compiler nor network process invoked
                else:raise ValueError('trace rendering did not stop before the environment-dependent compile')
            for n in ('main.f','capture.f'):
                ref.need((dest/n).read_bytes()==(staged/n).read_bytes(),mode+' regenerated '+n+' differs')
            main=(staged/'main.f').read_bytes()
            hook=(b'!$OMP SINGLE\n        CALL CAPTURE_MPI_SUBSTEP(ADVUA,ADVVA,TRNU,TRNV,VAMAX,ELMINX,IAMAX,JAMAX)\n'
                  b'!$OMP END SINGLE\n')
            ref.need(main.count(hook)==1 and main.replace(hook,b'')==(src/'ecom3d.f').read_bytes(),
                     mode+' main not exactly one capture hook')
            source=(staged/'capture.f').read_text()
            gathers=[x.strip() for x in source.splitlines() if 'CALL SMPI_' in x]
            ref.need(gathers==['CALL SMPI_ALL_RAW(BUFFER,IM,JM)'] and
                     'DO J=SMPI_JLO,SMPI_JHI\n        BUFFER(:,J)=CAPFIELD(:,J)\n      ENDDO' in source,
                     mode+' did not gather a dedicated copy of owned rows')
            commands=[shlex.split(x) for x in (staged/'build.log').read_text().splitlines()]
            b=ref.read(build/('trace-'+mode+'.json')); actual=ref.sha(build/('trace-'+mode+'-executable'))
            ref.need(len(commands)==3 and commands==b['commands'] and
                     b['snapshot_mode']=='owned_copy_nondestructive' and b['binary_sha']==actual==expected_elf,
                     mode+' actual ELF/build commands mismatch')
            for unit,c in zip(('main','capture'),commands[:2]):
                ref.need(str(remote_dir/(unit+'.f'))==c[c.index('-c')+1] and
                         str(remote_dir/(unit+'.o'))==c[c.index('-o')+1],mode+' wrong compiled sources')
            ref.need(all(c[0]=='/usr/lib64/openmpi/bin/mpifort' and
                         all(flag in c for flag in ('-fno-fast-math','-ffp-contract=off',
                         '-ffpe-trap=invalid,zero,overflow','-fopenmp')) for c in commands),
                     mode+' compiler/safety flags differ')
            link=commands[-1];objects=[x for x in link if x.endswith('.o')]
            science={str(REMOTE/'source/build/mpi'/(f.stem+'.o')) for f in src.iterdir() if
                     f.suffix.lower() in ('.f','.f90','.c') and f.stem!='ecom3d'}
            ref.need(link[link.index('-o')+1]==str(remote_dir/'s--pom.x') and len(objects)==len(set(objects))
                     and set(objects)==science|{str(remote_dir/(n+'.o')) for n in ('main','capture')},
                     mode+' contains an unrelated link object')
            result[mode]=dict(regenerated_sources_byte_identical=True,unchanged_science_objects_in_link=len(science),
                              actual_elf_sha256=actual,commands_match_build_receipt=True,only_owned_scratch_allgather=True)
    return dict(passed=True,scope='source reconstruction, build command and pre-pinned ELF evidence; not a native recompile',
                manifest_sha256=ref.sha(side/'MANIFEST.json'),bundle_manifest_sha256=ref.sha(bundle/'MANIFEST.json'),
                science_tree_sha256=TREE,builder_sha256=BUILDER,remote_source_and_log_sha256=remote,modes=result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('evidence',type=Path);p.add_argument('sidecar',type=Path)
    p.add_argument('--output',type=Path);a=p.parse_args()
    if a.output and a.evidence.resolve() in a.output.resolve().parents:p.error('cannot write to the immutable main bundle')
    s=json.dumps(audit(a.evidence.resolve(),a.sidecar.resolve()),indent=2)+'\n'
    if a.output:
        with a.output.open('x') as f:f.write(s)
    print(s,end='')
