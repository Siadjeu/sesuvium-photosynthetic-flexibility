#!/usr/bin/env python3
"""Copy verified completed releases and documentation; standard library only."""
import argparse,datetime,hashlib,json,shutil,tarfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--project-root',required=True);a=p.parse_args()
root=Path(a.project_root).resolve();base=root/'08_focused_families/downstream_supported/metabolomics';package=Path(__file__).resolve().parent
sources={'metabolite_comparisons':base/'release_20681419_20261021T134025Z','transcript_metabolite_comparison':base/'concordance_20732199_20261021T155959Z'}
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def verify(folder):
 manifest=folder/'SHA256SUMS.txt'
 if not manifest.is_file():raise RuntimeError('Missing source manifest: '+folder.name)
 records=[]
 for line in manifest.read_text().splitlines():
  h,n=line.split('  ',1); f=folder/n
  if Path(n).is_absolute() or '..' in Path(n).parts or f.is_symlink() or not f.resolve().is_relative_to(folder.resolve()):raise RuntimeError('Unsafe manifest path')
  if not f.is_file() or digest(f)!=h:raise RuntimeError('Checksum mismatch: '+folder.name+'/'+n)
  records.append((n,h))
 if any(f.is_symlink() for f in folder.rglob('*')):raise RuntimeError('Source has symlinks')
 return records
checks={k:verify(v) for k,v in sources.items()}
qc=json.loads((sources['transcript_metabolite_comparison']/'QC.json').read_text())
assert qc['genes']==99 and qc['assessable']==83 and qc['families']==18
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
out=base/'frozen_releases'/('SS_transcript_metabolite_evidence_v1_'+stamp);out.mkdir(parents=True,exist_ok=False)
for k,v in sources.items():
 shutil.copytree(v,out/'evidence'/k)
 verify(out/'evidence'/k)
shutil.copytree(package/'docs',out/'docs');shutil.copy2(__file__,out/'freeze.py');shutil.copy2(package/'freeze.slurm',out/'freeze.slurm');shutil.copy2(package/'README.txt',out/'README.txt')
provenance={k:{'source_relative_path':str(v.relative_to(root)),'source_manifest_sha256':digest(v/'SHA256SUMS.txt'),'verified_files':len(checks[k])} for k,v in sources.items()}
(out/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2))
(out/'RELEASE_STATUS.json').write_text(json.dumps({'status':'FROZEN_COMPUTATIONAL_EVIDENCE','version':'1.0.0','laboratory_metadata_status':'OPEN','github_publication_performed':False,'new_analyses':False},indent=2))
# Nested manifests are included; only the new top-level manifest excludes itself.
manifest=out/'SHA256SUMS.txt'
with manifest.open('w') as f:
 for item in sorted(out.rglob('*')):
  if item.is_file() and item!=manifest:f.write(digest(item)+'  '+str(item.relative_to(out))+'\n')
verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname=out.name)
bundle.with_name(bundle.name+'.sha256').write_text(digest(bundle)+'  '+bundle.name+'\n')
print('CONSOLIDATION COMPLETE');print('Release:',out);print('Upload this bundle:',bundle)
