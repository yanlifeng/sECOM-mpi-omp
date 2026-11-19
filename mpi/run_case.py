"""p3-only correctness development runner. No scheduling, old-run mutation, or retry.

MPI ranks currently have separate legacy I/O directories. This explicitly
measures a transitional bulk-publication prototype, NOT a scaling result.
CPU 4..35 is the entire allowed physical-core budget for this campaign.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time

OLD = Path('/home/user_home/ylf/secom-omp64-20261118.0SVeEm')
INPUT = Path('/home/user_home/ylf/secom-parallel-check-20261118.jNqsOL/code-delivery/case20m/input')
PRELOAD = '/home/user_home/ylf/secom-parallel-check-20261118.jNqsOL/h5_preinit.so'
DEPS = '/home/user_home/ylf/kt5-licom-secom-ylf/deps/lib'
MPIRUN = '/usr/lib64/openmpi/bin/mpirun'
OUTPUTS = ('fort.1002', 'fort.1006', 'startup', 'gcmplt.cdf', 'gcmprt', 'heat.dat')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2)


def rank_exec(args):
    rank = int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK'])
    size = int(os.environ['OMPI_COMM_WORLD_SIZE'])
    assert size == args.ranks and 0 <= rank < size
    cpus = list(range(4 + rank * args.threads, 4 + (rank + 1) * args.threads))
    assert set(cpus) <= set(range(4, 36))
    os.sched_setaffinity(0, cpus)
    run = args.output / ('rank%03d' % rank)
    os.chdir(run)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC='FALSE',
               OMP_PROC_BIND='TRUE', OMP_PLACES=','.join('{%d}' % c for c in cpus),
               OPENBLAS_NUM_THREADS='1', LD_PRELOAD=PRELOAD,
               LD_LIBRARY_PATH=DEPS + ':' + os.environ.get('LD_LIBRARY_PATH', ''))
    write(run / 'launch.json', dict(pid=os.getpid(), rank=rank, size=size,
          threads=args.threads, cpus=cpus, binary=str(args.binary), binary_sha=sha(args.binary)))
    fd = os.open(str(run / 'stdout.log'), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    os.dup2(fd, 1)
    os.dup2(fd, 2)
    os.close(fd)
    os.execve(str(args.binary), [str(args.binary)], env)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--binary', type=lambda p: Path(p).resolve(), required=True)
    ap.add_argument('--output', type=lambda p: Path(p).resolve(), required=True)
    ap.add_argument('--ranks', type=int, required=True)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--steps', type=int, choices=(300, 4800), default=300)
    ap.add_argument('--timeout', type=int, default=600)
    ap.add_argument('--rank-exec', action='store_true')
    ap.add_argument('--root-only-io', action='store_true', help='expect only rank0 large public outputs')
    args = ap.parse_args()
    assert args.ranks > 0 and args.threads > 0 and args.ranks * args.threads <= 32
    if args.rank_exec:
        rank_exec(args)
        return
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    root = args.output.parent
    root.mkdir(parents=True, exist_ok=True)
    with (root / '.solver.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        args.output.mkdir()  # Never overwrite/retry a previous sample.
        proof = json.loads((OLD / 'day19-serial-proof.json').read_text())
        assert proof['both_verified'] and proof['date'] == time.strftime('%Y%m%d')
        expected = proof['output_sha']['short' if args.steps == 300 else 'long']
        input_hashes = {}
        for rank in range(args.ranks):
            d = args.output / ('rank%03d' % rank)
            d.mkdir()
            for p in sorted(INPUT.iterdir()):
                if not p.is_file():
                    raise ValueError('unexpected input entry ' + str(p))
                if p.name == 'run_data.txt':
                    b = p.read_bytes()
                    if args.steps == 300:
                        for a, z in [(b'      4800     21600', b'       300     21600'),
                                     (b'1200.0000', b' 100.0000')]:
                            assert b.count(a) == 1
                            b = b.replace(a, z)
                    with (d / p.name).open('xb') as f:
                        f.write(b)
                else:
                    shutil.copyfile(p, d / p.name)
            hs = {p.name: sha(d / p.name) for p in INPUT.iterdir() if p.is_file()}
            if rank:
                assert hs == input_hashes
            input_hashes = hs
        assert len(input_hashes) == 24
        cpus = list(range(4, 4 + args.ranks * args.threads))
        topology = {c: {k: Path('/sys/devices/system/cpu/cpu%d/topology/%s' % (c, k)).read_text().strip()
                       for k in ('physical_package_id', 'core_id', 'thread_siblings_list')} for c in cpus}
        assert len({(v['physical_package_id'], v['core_id']) for v in topology.values()}) == len(cpus)
        argsv = [MPIRUN, '-np', str(args.ranks), '--host', 'localhost:%d' % args.ranks,
                 '--nooversubscribe', '--bind-to', 'none', '--mca', 'plm', 'isolated',
                 '--mca', 'pml', 'ob1', '--mca', 'btl', 'self,vader',
                 sys.executable, str(Path(__file__).resolve()), '--rank-exec',
                 '--binary', str(args.binary), '--output', str(args.output),
                 '--ranks', str(args.ranks), '--threads', str(args.threads)]
        info = dict(kind='correctness_prototype_not_scaling', ranks=args.ranks, threads=args.threads,
                    root_only_large_outputs=args.root_only_io,
                    steps=args.steps, cpus=cpus, topology=topology, input_sha=input_hashes,
                    binary_sha=sha(args.binary), command=argsv, load0=os.getloadavg(),
                    shared_with_other_jobs_allowed=True, started=time.time())
        write(args.output / 'started.json', info)
        seen = {}
        start = time.monotonic()
        env = dict(os.environ, LD_LIBRARY_PATH=DEPS + ':' + os.environ.get('LD_LIBRARY_PATH', ''))
        error = None
        with (args.output / 'mpirun.log').open('x') as log:
            proc = subprocess.Popen(argsv, stdout=log, stderr=subprocess.STDOUT, env=env,
                    start_new_session=True, preexec_fn=lambda: os.sched_setaffinity(0, cpus))
            try:
                while proc.poll() is None:
                    for rank in range(args.ranks):
                        p = args.output / ('rank%03d/launch.json' % rank)
                        try:
                            launch = json.loads(p.read_text())
                            pid = launch['pid']
                            if Path('/proc/%d/exe' % pid).resolve() != args.binary:
                                continue
                            masks = {int(t.name): sorted(os.sched_getaffinity(int(t.name)))
                                     for t in Path('/proc/%d/task' % pid).iterdir()}
                            assert masks and all(set(m) <= set(launch['cpus']) for m in masks.values())
                            if args.threads == 1:
                                assert all(m == launch['cpus'] for m in masks.values())
                            seen[rank] = dict(pid=pid, masks=masks, application_omp_threads=args.threads,
                                              os_threads_including_mpi=len(masks))
                        except (FileNotFoundError, ProcessLookupError, json.JSONDecodeError):
                            pass
                    if time.monotonic() - start > args.timeout:
                        raise TimeoutError('bounded prototype timeout')
                    time.sleep(0.1)
            except BaseException as ex:
                error = repr(ex)
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=8)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
        info.update(elapsed=time.monotonic() - start, returncode=proc.returncode,
                    load1=os.getloadavg(), rank_binding_samples=seen, error=error)
        checks = {}
        for rank in range(args.ranks):
            d = args.output / ('rank%03d' % rank)
            outputs = {n: sha(d / n) for n in OUTPUTS if (d / n).exists()}
            subset = ('gcmprt', 'heat.dat') if args.root_only_io and rank > 0 else OUTPUTS
            checks[rank] = dict(outputs=outputs, expected_output_files=list(subset),
                       matches_same_day_serial=outputs == {n: expected[n] for n in subset},
                       unchanged_inputs=all(sha(d / n) == s for n, s in input_hashes.items()))
        info['checks'] = checks
        info['passed'] = (error is None and proc.returncode == 0 and len(seen) == args.ranks and
                          all(x['matches_same_day_serial'] and x['unchanged_inputs'] for x in checks.values()))
        write(args.output / 'result.json', info)
        print(json.dumps(dict(passed=info['passed'], elapsed=info['elapsed'], ranks=args.ranks,
                             output=str(args.output),
                             rank_matches={r: x['matches_same_day_serial'] for r, x in checks.items()}), indent=2), flush=True)
        if not info['passed']:
            raise SystemExit(1)


if __name__ == '__main__':
    main()
