"""Independent byte/field audit, read-only to both simulation directories."""
import argparse
import array
import hashlib
import json
from pathlib import Path
import struct

FIELDS = ('DUMWAD DVMWAD EL ELB ELF ET ETF UA UAB UAF VA VAB VAF D '
          'CURV42D CBC TXBOT TYBOT FLUXUA FLUXVA TPS TRNU TRNV ADVUA ADVVA').split()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('reference', type=Path)
    ap.add_argument('candidate', type=Path)
    ap.add_argument('report', type=Path)
    ap.add_argument('--count', type=int, required=True)
    args = ap.parse_args()
    ref = sorted(args.reference.glob('transient-*.bin'))
    cand = sorted(args.candidate.glob('transient-*.bin'))
    assert len(ref) == len(cand) == args.count
    assert [p.name for p in ref] == [p.name for p in cand]
    report = dict(reference=str(args.reference), candidate=str(args.candidate),
                  field_order=FIELDS, frames={}, all_bitwise_equal=True)
    for a, b in zip(ref, cand):
        x, y = a.read_bytes(), b.read_bytes()
        h = struct.unpack('<4i2f2i', x[:32])
        frame = dict(header=list(h), header_equal=x[:32] == y[:32], fields={},
                     reference_sha=hashlib.sha256(x).hexdigest(), candidate_sha=hashlib.sha256(y).hexdigest())
        n = h[2] * h[3] * 4
        assert len(x) == len(y) == 32 + len(FIELDS) * n
        for k, name in enumerate(FIELDS):
            lo = 32 + k * n
            u, v = x[lo:lo + n], y[lo:lo + n]
            s = dict(bitwise_equal=u == v)
            if u != v:
                ui, vi = array.array('I'), array.array('I')
                uf, vf = array.array('f'), array.array('f')
                ui.frombytes(u); vi.frombytes(v); uf.frombytes(u); vf.frombytes(v)
                mismatch = [i for i, (r, z) in enumerate(zip(ui, vi)) if r != z]
                s.update(bit_mismatches=len(mismatch), first_i=mismatch[0] % h[2] + 1,
                         first_j=mismatch[0] // h[2] + 1,
                         max_abs_error=max(abs(float(uf[i]) - float(vf[i])) for i in mismatch))
            frame['fields'][name] = s
        frame['bitwise_equal'] = x == y
        report['all_bitwise_equal'] &= frame['bitwise_equal']
        report['frames'][a.name] = frame
    with args.report.open('x') as f:
        json.dump(report, f, indent=2)
    print(json.dumps(dict(all_bitwise_equal=report['all_bitwise_equal'], frames=len(ref),
                         fields_per_frame=len(FIELDS), report=str(args.report))))
    if not report['all_bitwise_equal']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
