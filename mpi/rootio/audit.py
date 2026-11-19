"""Verify compiler-parsed rank-zero paths against the unedited parent (GNU Fortran)."""
import argparse
from pathlib import Path
import re
import subprocess


def parse(compiler, includes, filename, omp):
    args = [compiler, '-w', '-ffixed-line-length-132', '-std=legacy',
            '-fsyntax-only', '-fdump-fortran-original', '-I' + str(includes)]
    if omp:
        args.append('-fopenmp')
    p = subprocess.run(args + [str(filename)], stdout=subprocess.PIPE, check=True, text=True)
    return '\n'.join(re.findall(r'^  code:\n(.*?)^---+', p.stdout, re.M | re.S)).splitlines()


def project(lines, proc, writer):
    guards, stack, ans = 0, [], []
    start = re.compile(r'^( *)IF \(= ' + proc + r':smpi_rank 0\)$')
    for line in lines:
        m = start.fullmatch(line)
        if m:
            stack.append([len(m[1]), writer])
            guards += 1
        elif stack and line == ' ' * stack[-1][0] + 'ELSE':
            stack[-1][1] = not writer
        elif stack and line == ' ' * stack[-1][0] + 'ENDIF':
            stack.pop()
        elif all(arm for _, arm in stack):
            ans.append(line.strip())
    if stack:
        raise ValueError('unmatched rank selection in compiler IR')
    return ans, guards


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--candidate', type=Path, required=True)
    ap.add_argument('--compiler', default='gfortran')
    args = ap.parse_args()
    for name, proc, delta in [('putcdf.f', 'putcdf', 445), ('ecom3d.f', 'ecom_3d', 30)]:
        old = parse(args.compiler, args.source, args.source / name, proc != 'putcdf')
        new = parse(args.compiler, args.source, args.candidate / name, proc != 'putcdf')
        before, n0 = project(old, proc, True)
        after, n1 = project(new, proc, True)
        if not before or before != after or n1 - n0 != delta:
            raise ValueError(name + ': rank-zero executable IR differs, or unreviewed guard count')
        peer, _ = project(new, proc, False)
        sinks = (r'CALL (?:nc\w+|secom_windtext)\b|^ASSIGN putcdf:\w+ nc\w+\[|'
                 r'WRITE UNIT=100[26]\b') if proc == 'putcdf' else r'(?:WRITE|OPEN|CLOSE) UNIT=ecom_3d:iuwrs\b'
        if any(re.search(sinks, line) for line in peer):
            raise ValueError(name + ': physical output call remains in nonroot path')
        if proc == 'putcdf':
            # Direct COMMON writes found from GCC8's IN-COMMON symbol attribute;
            # RDCORLOC also modifies shared geometry indirectly.
            for field in 'pts ptt ptw ptkh ptkm ptuf ptvf ptq2 pthf ptcd'.split():
                pattern = 'ASSIGN putcdf:' + field + ' '
                if [l for l in before if l.startswith(pattern)] != [l for l in peer if l.startswith(pattern)]:
                    raise ValueError('missing common mutation: ' + field)
            for prefix in ['CALL rdcorloc ', 'CALL datehr ', 'CALL getdate ', 'ASSIGN putcdf:icdf ']:
                actions = [l for l in before if l.startswith(prefix)]
                if not actions or actions != [l for l in peer if l.startswith(prefix)]:
                    raise ValueError('missing original side effects: ' + prefix)
        print(f'{name}: root executable IR unchanged; {delta} physical-output guards; nonroot sinks absent')


if __name__ == '__main__':
    main()
