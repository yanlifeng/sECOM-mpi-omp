"""Collect selected immutable receipts and outputs, only after all solvers exit.

No deletion, full-directory mirroring, or mutation of historical references.
The archive manifest names every regular member and records SHA-256 and size.
"""
import argparse
import fcntl
import hashlib
import io
import json
from pathlib import Path
import tarfile


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--archive', type=Path, required=True)
    ap.add_argument('--representative', default='r004b-p32-t1-full03')
    ap.add_argument('--trace-prefix', default='r004b-p32-t1-',
                    help='Exact tested candidate run prefix, including trailing dash')
    args = ap.parse_args()
    root, source = args.root.resolve(), args.source.resolve()
    assert root == Path('/home/user_home/ylf/secom-mpi-20261119')
    assert source.is_relative_to(root) if hasattr(source, 'is_relative_to') else str(source).startswith(str(root) + '/')
    assert args.archive.resolve().parent == root and args.archive.name.endswith('.tar.gz')
    old = Path('/home/user_home/ylf/secom-omp64-20261118.0SVeEm')
    files = {}
    def add(name, p):
        p = Path(p)
        assert name not in files and p.is_file() and not p.is_symlink(), (name, p)
        assert not name.startswith('/') and '..' not in Path(name).parts
        files[name] = p
    with (root / 'runs/.solver.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for p in sorted((root / 'runs').iterdir()):
            if not p.is_dir():
                continue
            assert (p / 'result.json').is_file(), 'unfinished run: ' + str(p)
            for n in ('started.json', 'result.json', 'trace-audit.json', 'mpirun.log'):
                if (p / n).is_file():
                    add('runs/' + p.name + '/' + n, p / n)
            for rank in sorted(p.glob('rank[0-9][0-9][0-9]')):
                for n in ('launch.json', 'stdout.log'):
                    if (rank / n).is_file():
                        add('runs/' + p.name + '/' + rank.name + '/' + n, rank / n)
                for log in sorted(rank.glob('smpi-profile-*.log')):
                    add('runs/' + p.name + '/' + rank.name + '/' + log.name, log)
        add('reference/day19-serial-proof.json', old / 'day19-serial-proof.json')
        add('collector/collect_evidence.py', Path(__file__).resolve())
        for n in ('build.log', 'build-omp.log', 'build-serial.log'):
            if (source.parent / n).is_file():
                add('build/' + n, source.parent / n)
        add('build/r004-initial-build.log', root / 'r004outer/build.log')
        add('build/mpi-executable', source / 'build/mpi/s--pom.x')
        for mode in ('omp', 'serial'):
            p = source / 'build' / mode / 's--pom.x'
            if p.is_file():
                add('build/' + mode + '-executable', p)
        plan = source.parent / 'formal-series-plan.json'
        if plan.is_file():
            add('build/formal-series-plan.json', plan)
            add('build/control-executable', root / 'r004b/source/build/mpi/s--pom.x')
        # Diagnostic evidence is separate from the timed scientific ELF.
        profile = root / 'profile-r014final-v1'
        if source.parent.name == 'r014final' and profile.is_dir():
            for p in sorted(profile.glob('*.json')):
                add('profile/' + p.name, p)
            for n in ('s--pom-profile.x', 'build.log'):
                if (profile / n).is_file():
                    add('profile/' + n, profile / n)
        for label in ('early', 'late'):
            add('build/trace-' + label + '.json', source.parent / ('trace-' + label) / 'build.json')
            add('build/trace-' + label + '-executable', source.parent / ('trace-' + label) / 's--pom.x')
        for p in sorted(source.iterdir()):
            if p.is_file():
                add('source/' + p.name, p)
        for p in sorted((source.parent / 'mpi').rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:
                add('scripts/' + str(p.relative_to(source.parent / 'mpi')), p)
        for n in ('fort.1002', 'fort.1006', 'startup', 'gcmplt.cdf', 'gcmprt', 'heat.dat'):
            add('outputs/' + n, root / 'runs' / args.representative / 'rank000' / n)
        # Two full 25-field frames independently compared on the Mac.
        # All twelve frames remain fully audited on p3 and have receipt hashes.
        for label, refdir, run, name in [
            ('early', 'traceref-00-11', args.trace_prefix + 'early', 'transient-153-05.bin'),
            ('late', 'deepref-00-81', args.trace_prefix + 'late', 'transient-3600-04.bin')]:
            add('frames/' + label + '/reference.bin', old / refdir / name)
            add('frames/' + label + '/candidate.bin', root / 'runs' / run / 'rank000' / name)
        inventory = {n: dict(bytes=p.stat().st_size, sha256=digest(p)) for n, p in files.items()}
        manifest = json.dumps(dict(representative=args.representative, members=inventory), indent=2).encode()
        with args.archive.open('xb') as out, tarfile.open(fileobj=out, mode='w:gz') as tar:
            for n, p in files.items():
                tar.add(str(p), arcname=n, recursive=False)
            info = tarfile.TarInfo('MANIFEST.json')
            info.size = len(manifest)
            info.mode = 0o600
            tar.addfile(info, io.BytesIO(manifest))
        archive_sha = digest(args.archive)
        with Path(str(args.archive) + '.sha256').open('x') as f:
            f.write(archive_sha + '  ' + args.archive.name + '\n')
        print(json.dumps(dict(archive=str(args.archive), bytes=args.archive.stat().st_size,
                              sha256=archive_sha, members=len(files))))


if __name__ == '__main__':
    main()
