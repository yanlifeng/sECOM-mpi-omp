"""Offline compatibility checks, not performance samples.

The caller first verifies the bundle inventory and frozen serial proof, then
passes that proof to audit(). No solver, subprocess, or remote access.
"""
import math
import re

import audit_evidence as ref


_R014_SHA = {'omp': '1fcea9cdc48c8a589154a2943adfde3e1c7f1d937fb7dd0eaaffb95e42335191',
             'serial': 'bdacf51d0349955e84bc6dedc5a285203c424245cfc98d424aea9d64e39f0307'}


def check_run(root, revision, tag, build, ranks, threads, steps, proof, digest):
    d = root / 'runs' / tag
    st = ref.read(d / 'started.json')
    res = ref.read(d / 'result.json')
    def need(ok, msg):
        ref.need(ok, tag + ': ' + msg)
    size = 'short' if steps == 300 else 'long'
    cpus = list(range(4, 4 + ranks * threads))
    base = ref.REMOTE + revision + '/source/build/' + build + '/s--pom.x'
    need(all(k in res and res[k] == v for k, v in st.items()), 'started/result disagree')
    # r014 registered unpartitioned one-rank serial/OMP launches without
    # --root-only-io, unlike its rank-partitioned hybrid launch. Both one-rank
    # modes must still produce (and below match) the same six actual hashes.
    io_expected = (ranks > 1) if revision == 'r014final' else True
    need((st['ranks'], st['threads'], st['steps'], st['cpus']) ==
         (ranks, threads, steps, cpus) and st['root_only_large_outputs'] is io_expected, 'shape/CPU bounds')
    need(len(st['input_sha']) == 24 and ref.canon(st['input_sha']) == ref.INPUT_SHA[size], '24 original inputs')
    need(res['passed'] is True and res['error'] is None and res['returncode'] == 0, 'solver return')
    need(all(type(x) in (float, int) and math.isfinite(x) and x > 0 for x in
             (st['started'], res['elapsed'])), 'non-positive or non-finite clock')
    need(len(st['topology']) == len(cpus) and set(st['topology']) == set(map(str, cpus)) and
         len({(v['physical_package_id'], v['core_id']) for v in st['topology'].values()}) == len(cpus),
         'duplicate or missing physical core')
    cmd = st['command']
    need(st['binary_sha'] == digest and cmd[0] == '/usr/lib64/openmpi/bin/mpirun' and
         cmd.count('--rank-exec') == 1 and cmd.count('--nooversubscribe') == 1 and not
         (d / 'trace-audit.json').exists(), 'wrong binary/launch mode')
    for flag, want in [('-np', str(ranks)), ('--host', 'localhost:' + str(ranks)),
                       ('--ranks', str(ranks)), ('--threads', str(threads)), ('--binary', base),
                       ('--output', ref.REMOTE + 'runs/' + tag)]:
        need(cmd.count(flag) == 1 and cmd.index(flag) + 1 < len(cmd) and
             cmd[cmd.index(flag) + 1] == want, 'launch flag ' + flag)
    rs = set(map(str, range(ranks)))
    need(set(res['checks']) == set(res['rank_binding_samples']) == rs and
         {p.name for p in d.glob('rank[0-9][0-9][0-9]')} ==
         {'rank%03d' % i for i in range(ranks)}, 'missing/extra ranks')
    all_tids = set()
    singles_seen = 0
    for rank in range(ranks):
        rd = d / ('rank%03d' % rank)
        obs = res['rank_binding_samples'][str(rank)]
        launch = ref.read(rd / 'launch.json')
        assigned = cpus[rank*threads:(rank+1)*threads]
        need((launch['rank'], launch['size'], launch['threads'], launch['cpus'],
              launch['binary'], launch['binary_sha']) ==
             (rank, ranks, threads, assigned, base, digest), 'launch receipt rank ' + str(rank))
        masks = obs['masks']
        need(launch['pid'] == obs['pid'] and str(launch['pid']) in masks and
             not all_tids.intersection(masks) and obs['os_threads_including_mpi'] == len(masks)
             and obs['application_omp_threads'] == threads, 'thread/PID receipts rank ' + str(rank))
        need(all(k.isdigit() and int(k) > 0 and v and set(v) <= set(assigned) and len(set(v)) == len(v)
                 for k, v in masks.items()), 'observed OS affinity outside rank assignment')
        all_tids.update(masks)
        singles_seen += len({v[0] for v in masks.values() if len(v) == 1})
        z = res['checks'][str(rank)]
        names = ref.FILES if rank == 0 else ('heat.dat', 'gcmprt')
        need(z['unchanged_inputs'] is True and z['matches_same_day_serial'] is True and
             sorted(z['expected_output_files']) == sorted(names) and
             z['outputs'] == {n: proof['output_sha'][size][n] for n in names},
             'wrong output/input receipts rank ' + str(rank))
        if ranks > 1:
            # We require the rank-specific domain and active scientific fast
            # path, not particular geometry-dependent corner-veto decisions.
            logs = (rd / 'stdout.log').read_text().splitlines()
            ds = [s for s in logs if 'MPI DOMAIN' in s]
            q, r = divmod(ref.HEIGHT, ranks)
            lo = 1 + rank*q + min(rank, r)
            hi = lo + q + int(rank < r) - 1
            m = re.fullmatch(r'\s*MPI DOMAIN\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+LOCALBC\s+T\s+OWNED_OUTER\s+T\s*',
                             ds[0]) if len(ds) == 1 else None
            need(m and tuple(map(int, m.groups())) == (rank, ranks, lo, hi), 'MPI domain ' + str(rank))
    return dict(tag=tag, build=build, ranks=ranks, threads=threads, steps=steps, started=st['started'],
                elapsed_seconds=res['elapsed'], actual_elf_sha256=digest, passed=True, timed_for_acceptance=False,
                physical_cores=len(cpus), observed_singleton_cores=singles_seen, observed_os_tasks=len(all_tids),
                all_os_masks_within_rank_assignment=True, all_input_and_output_receipts_match=True)


def audit(root, revision, prefix, proof, science_sha):
    """Verify OMP32, one 300-step serial run, and 8MPI x4 hybrid."""
    results = {}
    for suffix, build, ranks, threads, steps in [('omp32', 'omp', 1, 32, 4800),
                                                 ('serial300', 'serial', 1, 1, 300),
                                                 ('hybrid8x4', 'mpi', 8, 4, 4800)]:
        elf = root / 'build' / (build + '-executable')
        ref.need(elf.is_file() and not elf.is_symlink(), 'missing actual ' + build + ' ELF')
        with elf.open('rb') as f:
            ref.need(f.read(4) == b'\x7fELF', 'missing ELF signature: ' + build)
        sha = ref.sha(elf)
        if build == 'mpi':
            ref.need(sha == science_sha, 'hybrid is not the frozen science binary')
        elif revision == 'r014final':
            ref.need(sha == _R014_SHA[build], build + ' ELF differs from the independently supplied build SHA')
        results[suffix] = check_run(root, revision, prefix + suffix, build, ranks, threads, steps, proof, sha)
    ref.need(len({r['actual_elf_sha256'] for r in results.values()}) == 3, 'regression build alias')
    return results
