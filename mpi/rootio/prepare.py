"""Generate two candidate source files in a fresh directory; never edit science in place.

Physical output only. The untouched parent build must have mpideck in comdeck.
"""
import argparse
from collections import Counter
import hashlib
from pathlib import Path
import re

PUT_SHA = 'af03fcec1b66d3c9a994049ed1b01944bb0d8bc7fd0221b364cbfa6d80b22e8b'
FUN = {'NCCRE': 1, 'NCSFIL': 2, 'NCDDEF': 6, 'NCVDEF': 47, 'NCOPN': 1, 'NCVID': 29}
SUB = {'NCAPTC': 180, 'NCAPT': 127, 'NCVPT': 44, 'NCVPT1': 2, 'NCCLOS': 1, 'NCENDF': 1,
       'SECOM_WINDTEXT': 1}
NAME = r'[A-Z][A-Z0-9_]*'
ROOT = '      IF (SMPI_RANK.EQ.0)'


def fixed_groups(src):
    """Original statement starts + continuations; keep every source byte."""
    current = []
    for line in src.splitlines(keepends=True):
        if len(line) > 6 and line[:1] == ' ' and line[5] in ('0', ' ') and line[6:].strip():
            if current:
                yield current
            current = []
        current.append(line)
    if current:
        yield current


def rewrite(src, nc):
    count, out = Counter(), []
    for lines in fixed_groups(src):
        start = lines[0]
        code = start[6:].split('!', 1)[0].strip().upper()
        symbol = None
        if nc:
            fn = re.match('(' + NAME + r')\s*=\s*(' + '|'.join(FUN) + r')\s*\(', code)
            sub = re.match(r'CALL\s+(' + '|'.join(SUB) + r')\s*\(', code)
            text = re.match(r'WRITE\s*\(\s*(1002|1006)\s*,', code)
            if fn:
                symbol, kind = fn.groups()
            elif sub or text:
                kind = (sub or text)[1]
            else:
                out.extend(lines)
                continue
        else:
            f = re.match(r'(OPEN|WRITE|CLOSE)\s*\(IUWRS[,)]', code)
            if not f:
                out.extend(lines)
                continue
            kind = f[1]
        if start[:6] != '      ':
            raise ValueError('label, continuation, or unexpected form: ' + start)
        # The function result is a private file handle/id. A later unchanged
        # DIMS=... assignment may copy it, so never leave it undefined.
        if symbol:
            out.extend([ROOT + ' THEN\n', *lines,
                        '      ELSE\n', '        ' + symbol + '=0\n', '      ENDIF\n'])
        else:
            out.extend([ROOT + '\n', '     &' + start[6:], *lines[1:]])
        count[kind] += 1
    want = {**FUN, **SUB, '1002': 1, '1006': 2} if nc else {'OPEN': 4, 'WRITE': 24, 'CLOSE': 2}
    if count != Counter(want):
        raise ValueError(f'unreviewed statement inventory: {count}')
    return ''.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True, help='new, nonexistent output directory')
    a = ap.parse_args()
    put = (a.source / 'putcdf.f').read_bytes()
    if hashlib.sha256(put).hexdigest() != PUT_SHA:
        raise ValueError('PUTCDF changed; repeat common-symbol and external-call audit before adapting')
    deck = (a.source / 'comdeck').read_text().lower()
    if not re.search(r'''(?m)^\s*include\s*['"]mpideck['"]''', deck):
        raise ValueError('missing rank declaration')
    mainf = (a.source / 'ecom3d.f').read_bytes()
    edits = {'putcdf.f': rewrite(put.decode('latin1'), True),
             'ecom3d.f': rewrite(mainf.decode('latin1'), False)}
    a.out.mkdir()
    for name, data in edits.items():
        with (a.out / name).open('xb') as f:
            b = data.encode('latin1')
            f.write(b)
            print(hashlib.sha256(b).hexdigest(), name)


if __name__ == '__main__':
    main()
