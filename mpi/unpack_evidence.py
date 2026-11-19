"""Verify archive SHA, strict regular-member manifest, then extract to a new directory."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile


def sha(stream):
    h = hashlib.sha256()
    for b in iter(lambda: stream.read(1 << 20), b''):
        h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('archive', type=Path)
    ap.add_argument('sidecar', type=Path)
    ap.add_argument('destination', type=Path)
    ap.add_argument('--receipt', type=Path, required=True)
    a = ap.parse_args()
    assert not a.destination.exists() and not a.receipt.exists()
    expected = a.sidecar.read_text().split()[0]
    assert len(expected) == 64
    with a.archive.open('rb') as f:
        actual = sha(f)
    assert actual == expected, 'archive SHA mismatch'
    with tarfile.open(a.archive, 'r:gz') as tar:
        members = tar.getmembers()
        names = [m.name for m in members]
        assert len(names) == len(set(names))
        assert sum(m.size for m in members) <= 1 << 30, 'oversized evidence'
        for m in members:
            p = PurePosixPath(m.name)
            assert m.isfile() and not m.issym() and not m.islnk()
            assert not p.is_absolute() and '..' not in p.parts and str(p) == m.name
            assert '\\' not in m.name and 0 <= m.size <= 200 << 20
        manifest = json.load(tar.extractfile('MANIFEST.json'))
        assert set(manifest['members']) == set(names) - {'MANIFEST.json'}
        for m in members:
            if m.name == 'MANIFEST.json':
                continue
            record = manifest['members'][m.name]
            assert record == {'bytes': m.size, 'sha256': sha(tar.extractfile(m))}, m.name
        # Only validated regular members can be written, with no tar ownership,
        # mode, link, path traversal or special-file metadata being applied.
        a.destination.mkdir()
        for m in members:
            p = a.destination.joinpath(*PurePosixPath(m.name).parts)
            p.parent.mkdir(parents=True, exist_ok=True)
            with p.open('xb') as out:
                shutil.copyfileobj(tar.extractfile(m), out, 1 << 20)
    receipt = dict(archive_sha256=actual, members=len(members),
                   uncompressed_bytes=sum(m.size for m in members),
                   destination=str(a.destination.resolve()), all_members_verified=True)
    with a.receipt.open('x') as f:
        json.dump(receipt, f, indent=2)
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
