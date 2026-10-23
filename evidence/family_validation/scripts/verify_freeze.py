#!/usr/bin/env python3
"""Verify all frozen files and reject unexpected or missing files."""
import hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
expected=set()
for line in (root/'SHA256SUMS.txt').read_text().splitlines():
    sha,rel=line.split('  ',1)
    path=root/rel
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=sha:
        raise SystemExit('CHECKSUM FAILED: '+rel)
    expected.add(rel)
actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p!=root/'SHA256SUMS.txt'}
if actual!=expected:raise SystemExit('File inventory differs from frozen manifest')
print('PASS: all '+str(len(expected))+' frozen files match their SHA256 checksums')
