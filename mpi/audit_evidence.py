"""Offline r004b 32 x 1 acceptance. Read-only bundle; optional fresh JSON outside it.

Usage: python3 -B audit_evidence.py DIR [--output NEW.json]
Exit 0 = accepted; 1 = failed or incomplete evidence. No solver, network or dependencies.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import struct


FILES = ('fort.1002', 'fort.1006', 'startup', 'gcmplt.cdf', 'gcmprt', 'heat.dat')
FIELDS = ('DUMWAD DVMWAD EL ELB ELF ET ETF UA UAB UAF VA VAB VAF D '
          'CURV42D CBC TXBOT TYBOT FLUXUA FLUXVA TPS TRNU TRNV ADVUA ADVVA').split()
RANKS = 32
CPU = list(range(4, 36))
PREFIX = 'r004b-p32-t1-'
REMOTE = '/home/user_home/ylf/secom-mpi-20261119/'
TAGS = [PREFIX + ('full%02d' % n) for n in range(1, 6)]
REP = TAGS[2]
ELF = {'science': ('source/build/mpi', 'b2206db02406aa2dd2ea08e670b213a4aa989d6845cec5dcb775909dff18c155'),
       'early': ('trace-early', '8275712b8a89ed0dc910b5a9eb2ae8d8859bdb05a8017ec13f820aa426e731c0'),
       'late': ('trace-late', '6e84ca4b72df5aa4563b5318fb90ec1d052f61967b4cb0bb980b8710e37092f4')}
# Frozen, independently obtained from day19-selected-evidence.tar.gz; timing denominators
# from day19-human-readback.json SHA 1ee14fc2afe1fa5b9d70d620d1618480f880ae091dd01d9f5ba1e0dc1bd98d68.
PROOF_SHA = 'bc05ced6408653dd95cc7af58672fa4226421cb81a9668f8cd8803b70b91c7a8'
FROZEN, SAME_DAY = 325.5552825294435, 312.834679977037
# Canonical mapping SHA-256 from ORIGINAL serial traceref-00-11.json / deepref-00-81.json.
INPUT_SHA = {'short': '9dca8005e6123886712306c3e26131dbd7e44881d62a3d5274c3fc28eabf5d7b',
             'long': '31db4ae7a8531cab1f2b3baee0e050b7b995f14881a04afe3285f79391e2056b'}
FRAME_SHA = {'early': '63bb782533be9654eaeded43651da613bb71f8ef121b31518cc7c24baf7acfab',
             'late': '133b673fca85d499681fa622bef192279e4b5c0ef15f3478b8b4b86238a623dd'}
POSITIONS = {'early': [(1, 1), (1, 10), (17, 3), (153, 5), (300, 10)],
             'late': [(601, 2), (1199, 10), (1200, 10), (1201, 1), (2400, 7), (3600, 4), (4800, 10)]}
PICK = {'early': 'transient-153-05.bin', 'late': 'transient-3600-04.bin'}
WIDTH, HEIGHT = 372, 479
PLANE = WIDTH * HEIGHT * 4


def need(condition, what):
    if not condition:
        raise ValueError(what)


def read(path):
    def unique(pairs):
        out = {}
        for k, v in pairs:
            need(k not in out, 'duplicate JSON key: ' + k)
            out[k] = v
        return out
    return json.loads(path.read_text(), object_pairs_hook=unique)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for b in iter(lambda: stream.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def canon(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def inventory(root):
    m = read(root / 'MANIFEST.json')
    need(m['representative'] == REP, 'representative is not full03')
    files = {}
    for p in root.rglob('*'):
        need(not p.is_symlink(), 'symlink: ' + str(p))
        if not p.is_dir():
            need(p.is_file(), 'not regular: ' + str(p))
            if p != root / 'MANIFEST.json':
                files[p.relative_to(root).as_posix()] = p
    need(set(files) == set(m['members']), 'manifest/file inventory mismatch')
    for n, p in files.items():
        v = m['members'][n]
        need(v == {'bytes': p.stat().st_size, 'sha256': sha(p)}, 'manifest bytes or digest: ' + n)
    return len(files)


def run(root, tag, kind, proof):
    d = root / 'runs' / tag
    a, b = read(d / 'started.json'), read(d / 'result.json')
    check = lambda ok, why: need(ok, tag + ': ' + why)
    check(all(k in b and b[k] == v for k, v in a.items()), 'started/result disagree')
    size = 'short' if kind == 'early' else 'long'
    check(a['ranks'] == RANKS and a['threads'] == 1 and a['steps'] == (300 if size == 'short' else 4800),
          'wrong ranks, threads or steps')
    check(a['root_only_large_outputs'] is True and a['cpus'] == CPU, 'wrong I/O or CPU set')
    check(len(a['input_sha']) == 24 and canon(a['input_sha']) == INPUT_SHA[size],
          '24 input hashes do not match frozen serial')
    binsub, bindigest = ELF[kind]
    binary = REMOTE + 'r004b/' + binsub + '/s--pom.x'
    cmd = a['command']
    def opt(k):
        check(cmd.count(k) == 1, 'missing/duplicate launch option ' + k)
        return cmd[cmd.index(k) + 1]
    check(a['binary_sha'] == bindigest and opt('--binary') == binary and
          opt('--output') == REMOTE + 'runs/' + tag and opt('-np') == '32' and
          opt('--ranks') == '32' and opt('--threads') == '1' and '--nooversubscribe' in cmd,
          'wrong binary, tag or launch command')
    check(b['passed'] is True and b['returncode'] == 0 and b['error'] is None, 'solver failure')
    elapsed = b['elapsed']
    check(type(elapsed) in (int, float) and math.isfinite(elapsed) and elapsed > 0, 'elapsed')
    topo = a['topology']
    check(set(topo) == set(map(str, CPU)) and len({(x['physical_package_id'], x['core_id'])
          for x in topo.values()}) == RANKS, 'not 32 distinct physical cores')
    ranks = set(map(str, range(RANKS)))
    check(set(b['checks']) == set(b['rank_binding_samples']) == ranks, 'incomplete rank receipts')
    pids = set()
    for i in range(RANKS):
        rk, cpu = str(i), [CPU[i]]
        launch = read(d / ('rank%03d' % i) / 'launch.json')
        check((launch['rank'], launch['size'], launch['threads'], launch['cpus'],
               launch['binary'], launch['binary_sha']) == (i, RANKS, 1, cpu, binary, bindigest),
              'launch receipt rank ' + rk)
        obs = b['rank_binding_samples'][rk]
        masks = obs['masks']
        check(obs['pid'] == launch['pid'] and obs['pid'] not in pids and str(obs['pid']) in masks and
              obs['application_omp_threads'] == 1 and obs['os_threads_including_mpi'] == len(masks) and
              all(v == cpu for v in masks.values()), 'observed CPU masks rank ' + rk)
        pids.add(obs['pid'])
        q = b['checks'][rk]
        keys = FILES if i == 0 else ('gcmprt', 'heat.dat')
        check(q['outputs'] == {k: proof['output_sha'][size][k] for k in keys} and
              sorted(q['expected_output_files']) == sorted(keys) and q['unchanged_inputs'] is True and
              q['matches_same_day_serial'] is True, 'output/input mismatch rank ' + rk)
    check((d / 'trace-audit.json').is_file() == (kind != 'science'), 'untimed/trace separation')
    return elapsed, a['started']


def trace(root, label):
    d = root / 'runs' / (PREFIX + label)
    t = read(d / 'trace-audit.json')
    need(t['all_bitwise_equal'] is True and t['field_order'] == FIELDS, label + ': trace schema')
    fmt = 'transient-%03d-%02d.bin' if label == 'early' else 'transient-%04d-%02d.bin'
    names = {fmt % pt: pt for pt in POSITIONS[label]}
    need(set(t['frames']) == set(names), label + ': expected frame schedule')
    orig = {}
    for n, (step, sub) in names.items():
        f = t['frames'][n]
        h = f['header']
        need(len(h) == 8 and h[:4] == [step, sub, WIDTH, HEIGHT] and
             f['header_equal'] is True and f['bitwise_equal'] is True and
             set(f['fields']) == set(FIELDS) and all(x == {'bitwise_equal': True} for x in f['fields'].values()),
             'wrong frame/header/field result: ' + n)
        need(f['reference_sha'] == f['candidate_sha'], 'whole frame SHA mismatch: ' + n)
        orig[n] = {'bytes': 32 + PLANE * len(FIELDS), 'sha256': f['reference_sha']}
    need(canon(orig) == FRAME_SHA[label], label + ': trace hashes not original serial')
    n = PICK[label]
    f = t['frames'][n]
    ps = [root / 'frames' / label / (k + '.bin') for k in ('reference', 'candidate')]
    need(all(p.stat().st_size == orig[n]['bytes'] and sha(p) == orig[n]['sha256'] for p in ps),
         label + ': selected frame SHA/size')
    with ps[0].open('rb') as a, ps[1].open('rb') as b:
        head = a.read(32)
        need(head == b.read(32) and list(struct.unpack('<4i2f2i', head)) == f['header'], n + ': local header')
        individual = {}
        for field in FIELDS:
            x, y = a.read(PLANE), b.read(PLANE)
            need(len(x) == PLANE and x == y, n + ': local field ' + field)
            individual[field] = hashlib.sha256(x).hexdigest()
        need(not a.read(1) and not b.read(1), n + ': trailing bytes')
    return dict(receipt_frames=len(names), fields_per_frame=len(FIELDS), local_frame=n,
                local_frame_sha256=orig[n]['sha256'], local_field_sha256=individual)


def audit(root):
    members = inventory(root)
    proofpath = root / 'reference/day19-serial-proof.json'
    need(sha(proofpath) == PROOF_SHA, 'serial proof differs from frozen independent copy')
    proof = read(proofpath)
    need(proof['date'] == '20261119' and proof['both_verified'] is True and
         all(set(proof['output_sha'][n]) == set(FILES) for n in ('short', 'long')), 'invalid serial proof')
    times = {}
    prev = None
    for t in TAGS:
        elapsed, epoch = run(root, t, 'science', proof)
        need(type(epoch) in (int, float) and math.isfinite(epoch) and
             (prev is None or epoch >= prev), 'non-independent timing intervals: ' + t)
        prev = epoch + elapsed
        times[t] = elapsed
    outputs = {n: sha(root / 'outputs' / n) for n in FILES}
    need(outputs == proof['output_sha']['long'], 'actual local six outputs != day19 serial')
    traces = {}
    for label in POSITIONS:
        run(root, PREFIX + label, label, proof)
        traces[label] = trace(root, label)
    vals = list(times.values())
    m = statistics.median(vals)
    return dict(accepted=FROZEN / m >= 16, manifest_members=members,
                manifest_sha256=sha(root / 'MANIFEST.json'), serial_proof_sha256=PROOF_SHA,
                scope='r004b; five 4800-step 32MPI x 1OMP uninstrumented samples only; unrelated historical failures retained',
                input_mapping_sha256=INPUT_SHA, binaries_sha256={k: v[1] for k, v in ELF.items()},
                observed_affinity='every sampled OS task on recorded rank-specific CPU; 32 distinct physical cores (4..35)',
                science_seconds=times, median_seconds=m, range_seconds=[min(vals), max(vals)],
                frozen_original_seconds=FROZEN, same_day_original_seconds=SAME_DAY,
                frozen_speedup=FROZEN/m, same_day_speedup=SAME_DAY/m, frozen_16x_seconds=FROZEN/16,
                frozen_16x_met=FROZEN/m >= 16, same_day_16x_met=SAME_DAY/m >= 16,
                outputs_sha256=outputs, traces=traces)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('evidence', type=Path)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    root = args.evidence.resolve()
    if args.output and (root == args.output.resolve() or root in args.output.resolve().parents):
        p.error('--output must be a fresh file outside the immutable evidence directory')
    try:
        report = audit(root)
    except (OSError, KeyError, TypeError, IndexError, ValueError, struct.error) as e:
        report = dict(accepted=False, error=type(e).__name__ + ': ' + str(e))
    result = json.dumps(report, indent=2, allow_nan=False) + '\n'
    if args.output:
        with args.output.open('x') as f:
            f.write(result)
    print(result, end='')
    return 0 if report['accepted'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
