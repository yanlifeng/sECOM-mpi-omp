"""Read-only local check of selected 2026-11-19 receipts; does not launch solvers.

The package contains SHA receipts, not raw full fields. Physical byte/field
re-reads were performed on p3; this independently checks their reported values.
Usage: python3 -B campaign/audit_v44_package.py PATH_TO_TAR_GZ
"""
import hashlib
import json
import math
from pathlib import Path
from statistics import median
import sys
import tarfile

PIN = 'ef3c94e1e7eb6d7a773775733c205fb4a71a80d3dbfd056448a49aed011f6c82'
CPUS = list(range(4, 36)) + [i for i in range(39, 72) if i != 64]


def audit(path):
    raw = Path(path).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == PIN, 'changed package'
    with tarfile.open(path, 'r:gz') as tar:
        items = tar.getmembers()
        assert len(items) == len({m.name for m in items}) == 44
        assert all(m.isfile() and '/' not in m.name and m.size < 1048576 for m in items)
        def read(name):
            return json.load(tar.extractfile(name))
        plan = read('day19-v44-plan.json')
        ref = read('full64recover19serial-00-s1.json')
        result = read('day19-v44-result.json')
        tags = [m.name for m in items if m.name.startswith('full64day19') and m.name.endswith('.json')]
        assert len(tags) == 20
        for name in tags:
            d = read(name)
            n = d['threads']
            cpus = CPUS if n == 64 else [39]
            assert n in (1, 64) and d['max_threads'] == n
            assert d['steps'] == 4800 and d['trace_frames'] == {}
            assert d['input_sha'] == ref['input_sha'] and d['outputs'] == ref['outputs']
            assert len(d['input_sha']) == 24 and len(d['outputs']) == 6
            assert d['binary_sha'] == plan['ids'][d['variant']]
            assert d['single_cpu_bindings_seen'] == cpus
            assert 0 <= d['estimated_other_busy_core_equivalents_including_smt'] <= (1.5 if n == 64 else .05)
            assert d['matches_p3_serial'] and d['layout'] == 'normal' and d['memory_policy'] == 'default'
            shots = d['complete_affinity_snapshots']
            assert len(shots) == 2 and shots[0]['tid_to_cpu'] == shots[1]['tid_to_cpu']
            for s in shots:
                assert sorted(s['tid_to_cpu'].values()) == cpus
                assert str(d['solver_pid']) in s['tid_to_cpu']
            assert shots[1]['at_s'] - shots[0]['at_s'] - shots[0]['scan_s'] >= 1.0
        confirm = result['confirmation']
        a, b = [median(confirm['times'][k]) for k in ('EdgeHist', 'Cj')]
        assert confirm['candidate_median'] == a and confirm['control_median'] == b
        assert sum(x < y for x, y in zip(confirm['times']['EdgeHist'], confirm['times']['Cj'])) == 3
        assert a <= .99*b and result['median64'] == a
        for mode, count, variants in [('early', 5, 2), ('late', 7, 1)]:
            d = read('day19-v44-'+mode+'-fields.json')
            assert d['all_matches'] and d['mismatches'] == []
            assert len(d['frames']) == count and d['field_count'] == 25 and len(d['candidates']) == variants
        return dict(samples=20, v44_seconds=a, v30_seconds=b, reduction=1-a/b,
                    reference_scope='same-day raw SHA; no field files in local package')


if __name__ == '__main__':
    print(json.dumps(audit(sys.argv[1]), indent=2))
