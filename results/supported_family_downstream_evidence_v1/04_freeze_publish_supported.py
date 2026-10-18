#!/usr/bin/env python3
"""Freeze verified supported-family evidence. Run on the HPC login node.
Default: prepare and verify local branch/tag. --push: also publish both refs.
Requires Python >=3.9 and git. No existing release or working tree is edited.
"""
import argparse,hashlib,json,os,re,subprocess,tempfile
from pathlib import Path
BASE_TAG='photosynthesis-family-evidence-v1.0.0'
BASE_COMMIT='2e6e7a0f93b813cc5d4c58f60e635e1a82f215b1'
BRANCH='freeze/supported-family-downstream-evidence-v1'
TAG='supported-family-downstream-evidence-v1.0.0'
DEST='results/supported_family_downstream_evidence_v1'
MANIFEST_SHA='fc8457b0c0303c0a804935b27aaaa5e0fa07c28f0c01dc90f31b3c0892cedf2b'
ALLOWED_REMOTES={'git@github.com:Siadjeu/sesuvium-photosynthetic-flexibility.git','https://github.com/Siadjeu/sesuvium-photosynthetic-flexibility.git','https://github.com/Siadjeu/sesuvium-photosynthetic-flexibility'}
def require(ok,message):
 if not ok:raise RuntimeError(message)
def digest(data):return hashlib.sha256(data).hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--push',action='store_true',help='Publish verified branch and annotated tag, without force');args=ap.parse_args()
 project=Path(os.environ['PROJECT_ROOT']).resolve()
 repo=Path(os.environ.get('GITHUB_REPO_DIR',str(project/'github_release'))).resolve()
 release=Path(os.environ.get('SUPPORTED_RELEASE',str(project/'08_focused_families/downstream_supported/integration/release_20559738_20261018T192606Z'))).resolve()
 def git(*argv,data=None,env=None,check=True):
  p=subprocess.run(['git','-C',str(repo),*argv],input=data,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
  if check and p.returncode:raise RuntimeError('git '+argv[0]+' failed: '+p.stderr.decode(errors='replace'))
  return p.stdout.decode().strip()
 def ref(name):
  return git('rev-parse','--verify',name,check=False)
 require(repo.is_dir(),'Repository directory missing')
 require(git('rev-parse','--is-inside-work-tree')=='true','Use the existing GitHub checkout')
 require(ref('refs/tags/'+BASE_TAG+'^{commit}')==BASE_COMMIT,'Family freeze tag does not match the verified base commit')
 require(not git('status','--porcelain'),'Repository has uncommitted/untracked files; keep them safe and resolve before running')
 for opts in [('remote','get-url','--all','origin'),('remote','get-url','--push','--all','origin')]:
  urls=git(*opts).splitlines();require(len(urls)==1 and urls[0] in ALLOWED_REMOTES,'Origin must be the expected GitHub repository with one URL')
 original_head=ref('HEAD');original_symbolic=git('symbolic-ref','-q','HEAD',check=False)
 protected={x:ref(x) for x in ['refs/heads/main','refs/tags/'+BASE_TAG,'refs/heads/freeze/photosynthesis-family-evidence-v1']}
 manifest=(release/'SHA256SUMS.txt').read_bytes()
 require(digest(manifest)==MANIFEST_SHA,'Source manifest differs from the verified job 20559738 release; stop for review')
 payload={};seen=set()
 for line in manifest.decode().splitlines():
  h,name=line.split(maxsplit=1);name=name.lstrip('*');p=Path(name)
  require(not p.is_absolute() and '..' not in p.parts and len(p.parts)==1,'Unsafe source manifest path')
  require(name not in seen,'Duplicate manifest entry');seen.add(name)
  f=release/name;require(f.is_file() and not f.is_symlink(),'Missing or symlinked input: '+name)
  b=f.read_bytes();require(digest(b)==h,'Source checksum mismatch: '+name);payload[name]=b
 require(len(payload)==25,'Unexpected verified payload size')
 payload['SHA256SUMS.txt']=manifest
 qc=json.loads(payload['QC.json']);require(qc['supported_records']==657 and qc['mapped_loci']==657 and qc['new_analyses'] is False,'Unexpected integration QC')
 require(not git('ls-tree',BASE_COMMIT,'--',DEST),'Destination already exists in the family freeze')
 payload['04_freeze_publish_supported.py']=Path(__file__).read_bytes()
 payload['README.md']=('''# Supported-family downstream evidence v1

This separate release integrates existing WGD/segmental duplication, CAFE5,
aBSREL positive-selection and PO/SS differential-expression results for 657
supported family loci. Read METHODS_RESULTS.txt before using the tables.

- 02_family_species_summary.tsv: family/species counts and denominators.
- 09_family_species_gene_lists.tsv: genes grouped by evidence category.
- matrix_*.tsv: Wang-style species-by-family tables.
- 04_DE_locus_contrast_evidence.tsv and 05_DE_family_contrast_summary.tsv:
  individual authorized contrasts; pooled DE totals mix biological comparisons.
- 01_supported_locus_evidence.tsv: identifiers and assessment coverage.

NA is not zero. Nine PO PEPC loci belong to a primary-model expanded
orthogroup; sensitivity-model branch tests were unavailable for this subset.
The 13 orthogroup-coverage exceptions are not failed gene-ID mappings.
Supported family membership does not establish a C4/CAM-specific role.

SHA256SUMS.txt verifies the original integration payload.
FREEZE_SHA256SUMS.txt additionally covers this README and freeze script.
SOURCE_BUNDLE_SHA256SUMS.txt applies to the earlier source collection.
Source data remain in the existing HPC collection; source hashes and relative
paths are retained in SOURCE_MANIFEST.tsv. No new biological analysis was run.

## Reproduce the integration

Use Python >=3.9 and the collected source release:

```bash
python3 integrate_supported.py --source-release "$SOURCE_RELEASE" --output "$NEW_OUTPUT"
```

For SLURM, set PROJECT_ROOT, SOURCE_RELEASE and INTEGRATION_SCRIPT_DIR to the
directory containing the integration script, then submit 03_integrate_supported.slurm.
The separate freeze script uses PROJECT_ROOT, optionally GITHUB_REPO_DIR and
SUPPORTED_RELEASE. Run without arguments to prepare; add --push to publish.
The expected verified integration manifest and base commit are pinned.

The parent is the unchanged photosynthesis-family-evidence-v1.0.0 freeze.
This branch/tag does not update main. Version tags are immutable by convention;
future corrections require a new version, not replacement of this tag.
''').encode()
 metadata=dict(version=TAG,parent_family_freeze=BASE_TAG,parent_commit=BASE_COMMIT,source_job='20559738',source_manifest_sha256=MANIFEST_SHA,supported_loci=657,privacy_scope='New release files scanned; inherited base is the previously audited exact family-freeze commit. Historical Git commits are not rescanned.',new_analyses=False)
 payload['FREEZE_METADATA.json']=(json.dumps(metadata,indent=2)+'\n').encode()
 # Detect private filesystem roots and cluster hostnames. Portable OS paths are allowed.
 private=re.compile(r'/(?:fs/)?(?:dss|home|scratch|workspace|lustre|gpfs|mnt|Users)/[^\s\"\x27<>]+|\b(?:hpcl|mpcb|mpcs)\d{3}\b')
 findings=[]
 for name,b in payload.items():
  try:text=b.decode('utf-8')
  except UnicodeDecodeError:raise RuntimeError('Unexpected binary payload: '+name)
  if private.search(text):findings.append(name)
 require(not findings,'Private-path/hostname audit failed in: '+', '.join(findings))
 payload['FREEZE_SHA256SUMS.txt']=''.join(digest(b)+'  '+name+'\n' for name,b in sorted(payload.items())).encode()
 # A private temporary Git index preserves the user's checkout and normal index.
 with tempfile.TemporaryDirectory(prefix='supported-freeze-') as temp:
  env=os.environ.copy();env['GIT_INDEX_FILE']=str(Path(temp)/'index')
  git('read-tree',BASE_COMMIT,env=env)
  for name,b in sorted(payload.items()):
   blob=git('hash-object','-w','--stdin',data=b)
   git('update-index','--add','--cacheinfo','100644',blob,DEST+'/'+name,env=env)
  tree=git('write-tree',env=env)
  branch_ref='refs/heads/'+BRANCH;tag_ref='refs/tags/'+TAG
  existing=ref(branch_ref);existing_tag=ref(tag_ref)
  if existing or existing_tag:
   require(bool(existing) and bool(existing_tag),'Only one target ref exists; stop to inspect instead of overwriting')
   require(ref(branch_ref+'^{tree}')==tree,'Existing freeze content differs; use a new version')
   require(git('rev-list','--parents','-n','1',existing).split()==[existing,BASE_COMMIT],'Unexpected freeze parent')
   require(git('cat-file','-t',existing_tag)=='tag' and ref(tag_ref+'^{commit}')==existing,'Tag does not identify the prepared freeze')
   commit=existing;tag_object=existing_tag
  else:
   commit=git('commit-tree',tree,'-p',BASE_COMMIT,data=b'Freeze supported-family downstream evidence v1\n')
   ident=git('var','GIT_COMMITTER_IDENT')
   tag_text=f'object {commit}\ntype commit\ntag {TAG}\ntagger {ident}\n\nVerified supported-family downstream evidence v1.\n'
   tag_object=git('mktag',data=tag_text.encode())
   transaction=f'start\ncreate {branch_ref} {commit}\ncreate {tag_ref} {tag_object}\nprepare\ncommit\n'
   git('update-ref','--stdin',data=transaction.encode())
 require(ref('HEAD')==original_head and git('symbolic-ref','-q','HEAD',check=False)==original_symbolic,'Checkout changed unexpectedly')
 require(all(ref(k)==v for k,v in protected.items()),'Protected local ref changed unexpectedly')
 changed=git('diff-tree','--no-commit-id','--name-status','-r',BASE_COMMIT,commit).splitlines()
 require(len(changed)==len(payload) and all(x.startswith('A\t'+DEST+'/') for x in changed),'Unexpected change outside new release')
 print(f'PASS: 25 source checksums; {len(payload)} new freeze files; private-path audit.')
 print('Prepared commit: '+commit+'\nBranch: '+BRANCH+'\nTag: '+TAG)
 if args.push:
  remote={line.split('\t')[1]:line.split('\t')[0] for line in git('ls-remote','origin',branch_ref,tag_ref).splitlines()}
  for name,value in [(branch_ref,commit),(tag_ref,tag_object)]:require(name not in remote or remote[name]==value,'Remote ref differs; refusing to replace '+name)
  git('push','--atomic','origin',branch_ref+':'+branch_ref,tag_ref+':'+tag_ref)
  remote={line.split('\t')[1]:line.split('\t')[0] for line in git('ls-remote','origin',branch_ref,tag_ref).splitlines()}
  require(remote.get(branch_ref)==commit and remote.get(tag_ref)==tag_object,'Remote verification failed; inspect before retrying')
  print('PUBLISHED and remote refs verified. Main and family freeze were not changed.')
  print('https://github.com/Siadjeu/sesuvium-photosynthetic-flexibility/tree/'+TAG+'/'+DEST)
 else:print('PREPARED locally. Re-run the same script with --push to publish.')
if __name__=='__main__':
 try:main()
 except (RuntimeError,OSError,KeyError,ValueError) as e:raise SystemExit('ERROR: '+str(e))
