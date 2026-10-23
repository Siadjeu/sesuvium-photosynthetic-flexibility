#!/usr/bin/env python3
"""Verify the curated publication's SHA256 manifest; no external packages."""
import hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
count=0
for line in (root/'SHA256SUMS.txt').read_text().splitlines():
 expected,name=line.split(maxsplit=1);p=(root/name).resolve();p.relative_to(root)
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise SystemExit('FAIL: '+name)
 count+=1
print(f'PASS: {count} publication files match their SHA256 checksums')
