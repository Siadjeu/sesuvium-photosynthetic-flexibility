#!/usr/bin/env python3
"""Run frozen SS analysis using published input copies in a temporary workspace."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--component',required=True,choices=['metabolite_comparisons','transcript_metabolite_comparison']);p.add_argument('--out',required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[1];src=root/'evidence/ss_metabolite/evidence'/a.component;out=Path(a.out).resolve()
if out.exists():raise SystemExit('Choose a new output directory')
with tempfile.TemporaryDirectory(prefix='ss_reproduction_') as d:
 work=Path(d);shutil.copytree(src/'inputs',work/'data');shutil.copyfile(src/'run.py',work/'run.py')
 hp=work/'data/input_hashes.json'
 if hp.exists():
  for name,h in json.loads(hp.read_text()).items():
   f=(work/'data'/name).resolve();f.relative_to(work/'data')
   if hashlib.sha256(f.read_bytes()).hexdigest()!=h:raise SystemExit('Input checksum mismatch: '+name)
 subprocess.run([sys.executable,str(work/'run.py'),'--out',str(out)],check=True)
