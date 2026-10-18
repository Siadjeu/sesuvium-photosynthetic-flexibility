#!/usr/bin/env bash
#SBATCH --job-name=beta_pepc_audit
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=beta_pepc_audit_%j.out
#SBATCH --error=beta_pepc_audit_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
for tool in python3; do command -v "$tool"; done
python3 - <<'PY'
import gzip,csv,hashlib,json,os,re,subprocess,tarfile,shutil
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def release(step,env):
 paths=[root/os.environ[env]] if os.environ.get(env) else [p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
 if len(paths)!=1:raise SystemExit('Expected one release; set '+env+' to the intended release directory')
 p=paths[0].resolve(strict=True)
 for line in (p/'SHA256SUMS.txt').read_text().splitlines():
  h,n=line.split(maxsplit=1);f=(p/n).resolve()
  if not f.is_relative_to(p) or sha(f)!=h:raise SystemExit('Checksum mismatch: '+n)
 return p
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta(p):
 d={};key=None
 for line in p.read_text().splitlines():
  if line.startswith('>'):
   key=line[1:].split()[0]
   if key in d:raise SystemExit('Duplicate FASTA ID: '+key)
   d[key]=''
  elif line.strip():
   if key is None:raise SystemExit('Invalid FASTA')
   d[key]+=line.strip()
 return d
def writefa(p,d):p.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
src=release('step25kz17.38','RECONCILIATION_RELEASE');prior=release('step25kz17.37','INTEGRATED_RELEASE')
assert json.loads((src/'PROVENANCE.json').read_text())['integrated_manifest_sha256']==sha(prior/'SHA256SUMS.txt')
beta=read(src/'01_BETA_CA_protein_evidence.tsv');pepc=read(src/'02_PEPC_protein_evidence.tsv');ledger=read(prior/'01_integrated_locus_evidence.tsv')
assert len(beta)==110 and len(pepc)==55 and len(ledger)==710
assert all(r['mapping_status']=='UNAMBIGUOUS' for r in beta+pepc)
expected={'BETA_SUPPORTED_SCREEN_EVIDENCE':45,'BETA_SUPPORTED_WITH_REVIEW_CAVEAT':2,'BETA_PARTIAL_UNRESOLVED':2,'OUTSIDE_SUPPORTED_BETA_SET':61}
assert dict(Counter(r['evidence_category'] for r in beta))==expected
out=base/'step25kz17.39'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'sources').mkdir()
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
tier={'BETA_SUPPORTED_SCREEN_EVIDENCE':'SUPPORTED','BETA_SUPPORTED_WITH_REVIEW_CAVEAT':'MODEL_CAVEAT','BETA_PARTIAL_UNRESOLVED':'UNRESOLVED','OUTSIDE_SUPPORTED_BETA_SET':'OUTSIDE_FAMILY_SET'}
for r in beta:
 ledger.append(dict(wang_family='beta-CA',species=r['species'],gene_id=r['gene_id'],protein_id=r['protein_id'],subunit='',evidence_category=r['evidence_category'],count_tier=tier[r['evidence_category']],source_file='step25kz17.38/01_BETA_CA_protein_evidence.tsv',model_or_scope_limitation='Original class='+r['original_classification']+'; residue='+r['residue_check']+'; selection='+r['sequence_selection_reason']+'; '+r['limitation'],exhaustive_count='NOT_ESTABLISHED',C4_CAM_specific_function='NOT_ESTABLISHED'))
assert len(ledger)==820 and len({(r['wang_family'],r['species'],r['gene_id']) for r in ledger})==820
table('01_integrated_locus_evidence.tsv',ledger,list(ledger[0]))
shutil.copy2(src/'01_BETA_CA_protein_evidence.tsv',out/'02_BETA_CA_source_decisions.tsv')
# Trace existing PEPC validation within bounded, explicitly recorded project scopes.
# Read metadata/scripts only, never execute them or follow arbitrary shell commands.
scopes=[base/'step25kz16.8b',base/'step25kz16.8',base/'step25kz16.7',base/'step25kn/results/PEPC',base/'step25kz17.2/release/PEPC']
if os.environ.get('PEPC_AUDIT_DIR'):scopes.append(root/os.environ['PEPC_AUDIT_DIR'])
ids={r['protein_id'] for r in pepc};alias_to_ids=defaultdict(set)
for r in pepc:
 alias_to_ids[r['protein_id']].add(r['protein_id'])
 # Local IDs are used only when a row explicitly provides the matching species.
 local=r['protein_id'].split('__',1)[1];alias_to_ids[r['species']+'__'+local].add(r['protein_id'])
files=set()
for scope in scopes:
 if not scope.is_dir():continue
 for folder,dirs,names in os.walk(scope):
  dirs[:]=[d for d in dirs if d not in ['work','tmp','logs','source_checkpoints'] and not d.startswith('.')]
  for n in names:
   p=Path(folder)/n;rel=str(p.relative_to(scope))
   if len(p.relative_to(scope).parts)>5:continue
   if p.suffix.lower() not in ['.tsv','.csv','.json','.txt','.py','.sh','.slurm']:continue
   if 'PEPC_KINASE' in rel.upper():continue
   if p.suffix in ['.py','.sh','.slurm'] or re.search('pepc|registry|validation|provenance|manifest',rel,re.I):files.add(p.resolve())
manifest=[];matches=[];contexts=[];issues=[]
for i,p in enumerate(sorted(files),1):
 rel=str(p.relative_to(root)) if p.is_relative_to(root) else p.name
 if i>200 or p.stat().st_size>2_000_000:
  issues.append(dict(source=rel,issue='SKIPPED_BOUNDED_FILE_COUNT_OR_SIZE_LIMIT'));continue
 try:text=p.read_text()
 except (UnicodeError,OSError):issues.append(dict(source=rel,issue='UNREADABLE_TEXT'));continue
 digest=sha(p);label='S'+str(i).zfill(3);dest=out/'sources'/(label+'_'+p.name)
 # Anonymize project/home paths in snapshots; hash records the original bytes.
 clean=text.replace(str(root),'<PROJECT_ROOT>')
 clean=re.sub(r'/(?:fs/dss/home|home)/[^/\s]+','<HOME>',clean)
 dest.write_text(clean)
 manifest.append(dict(source_id=label,source=rel,original_sha256=digest,snapshot_sha256=sha(dest),bytes=p.stat().st_size,snapshot=str(dest.relative_to(out))))
 for line_no,line in enumerate(text.splitlines(),1):
  if re.search(r'tree_validated_registry|registry|PRIMARY_FULL_LENGTH_PEPC|candidate_validation|read_csv|read_table',line,re.I):
   if len([x for x in contexts if x['source_id']==label])<150:contexts.append(dict(source_id=label,line_number=line_no,text=line[:2000].replace(str(root),'<PROJECT_ROOT>')))
 if p.suffix.lower() not in ['.tsv','.csv']:continue
 try:
  import io
  rr=csv.DictReader(io.StringIO(text),delimiter='\t' if p.suffix.lower()=='.tsv' else ',')
  for rownum,r in enumerate(rr,2):
   vals=[str(v) for v in r.values() if v is not None];found=set()
   for v in vals:
    if v in ids:found.add(v)
   sp=r.get('species') or r.get('species_code')
   if sp:
    for v in vals:
     candidate=sp+'__'+v
     if candidate in ids:found.add(candidate)
   for pid in found:matches.append(dict(protein_id=pid,source_id=label,row_number=rownum,record_json=json.dumps(r,sort_keys=True)))
 except (csv.Error,TypeError):issues.append(dict(source=rel,issue='TABLE_PARSE_ERROR'))
# Presence in a registry is provenance evidence, not automatic proof of its validation criteria.
audit=[]
source_lookup={r['source_id']:r for r in manifest}
for r in pepc:
 mm=[x for x in matches if x['protein_id']==r['protein_id']]
 sourceids=sorted({x['source_id'] for x in mm})
 non_tree=[x for x in sourceids if not re.search(r'tree_tip_metadata|tree_tip_map',source_lookup[x]['source'],re.I)]
 audit.append(dict(protein_id=r['protein_id'],species=r['species'],gene_id=r['gene_id'],historical_classification=r['original_classification'],matched_source_files=len(sourceids),non_tree_source_files=len(non_tree),source_ids=';'.join(sourceids),provenance_status='NON_TREE_RECORDS_LOCATED_REVIEW_CRITERIA' if non_tree else 'TREE_ONLY_OR_NO_ADDITIONAL_RECORDS',membership_decision='PENDING_PROVENANCE_REVIEW',completeness='NOT_ESTABLISHED'))
table('03_PEPC_provenance_audit.tsv',audit,list(audit[0]));table('04_PEPC_exact_matching_records.tsv',matches,['protein_id','source_id','row_number','record_json']);table('05_PEPC_producer_and_registry_context.tsv',contexts,['source_id','line_number','text']);table('06_PEPC_source_manifest.tsv',manifest,['source_id','source','original_sha256','snapshot_sha256','bytes','snapshot']);table('07_bounded_collection_issues.tsv',issues,['source','issue'])
# Refresh six-family evidence matrices; preserve final-count NA for all 23 categories.
counts=read(prior/'02_species_family_evidence_counts.tsv');tiers=['SUPPORTED','LIMITED_PANEL_CAVEAT','MODEL_CAVEAT','UNRESOLVED','OUTSIDE_FAMILY_SET']
for r in counts:
 if r['wang_family']=='beta-CA':
  rr=[x for x in ledger if x['wang_family']=='beta-CA' and x['species']==r['species']];c=Counter(x['count_tier'] for x in rr)
  for t in tiers:r[t]=c[t]
  r['records_integrated']=len(rr);r['coverage_status']='RECENT_LEDGER_INTEGRATED'
table('08_species_family_evidence_counts.tsv',counts,list(counts[0]))
coverage=read(prior/'05_refreshed_family_coverage.tsv')
for r in coverage:
 if r['wang_family']=='beta-CA':r.update(current_status='EVIDENCE_CATEGORIES_INTEGRATED_WITH_LIMITATIONS',integrated_records=110,next_action='RETAIN_SUPPORTED_CAVEAT_PARTIAL_AND_OUTSIDE_CATEGORIES')
 elif r['wang_family']=='PEPC':r.update(current_status='BOUNDED_PROVENANCE_REVIEW_REQUIRED',next_action='REVIEW_REGISTRY_ROWS_AND_PRODUCER_CRITERIA')
table('09_refreshed_family_coverage.tsv',coverage,list(coverage[0]))
shutil.copy2(prior/'04_final_Wang_counts_PENDING.tsv',out/'10_final_Wang_counts_PENDING.tsv')
(out/'README.txt').write_text('BETA-CA integration: 45 supported screen records, two supported with caveats, two unresolved partial records, 61 outside the supported beta set. Alpha-family evidence is not accepted as beta-CA. The AT alignment isoform is not an extra locus. Prior five-family ledger unchanged. Final counts remain NA; supported does not imply confirmed biochemical function or exhaustive membership.\nPEPC provenance audit searches only recorded bounded scopes, at most 200 text files <=2 MB each and depth <=5, excluding work/tmp/logs/source_checkpoints and PEPC_KINASE. File presence and exact row matches do not establish quality, independence or acceptance. Multiple releases remain visible; no latest-file assumption. Sources have original hashes and anonymized snapshots. No scripts executed, no analyses rerun, no PEPC promotion.\nUnmatched IDs can reflect different identifier formats; only exact full IDs or local IDs with explicit row species are matched. Missing/limited collection is recorded, not treated as evidence of biological absence. Use source context to trace upstream producer/registry criteria before deciding membership.\n')
(out/'QC.json').write_text(json.dumps(dict(status='BETA_INTEGRATED_PEPC_PROVENANCE_REVIEW_REQUIRED',integrated_records=len(ledger),beta_categories=expected,PEPC_loci=len(audit),PEPC_sources=len(manifest),PEPC_matching_rows=len(matches),PEPC_provenance_status=dict(Counter(r['provenance_status'] for r in audit)),collection_issues=len(issues),analysis_rerun=False,final_counts=False),indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(reconciliation_manifest_sha256=sha(src/'SHA256SUMS.txt'),integration_manifest_sha256=sha(prior/'SHA256SUMS.txt'),scopes=[str(p.relative_to(root)) if p.is_relative_to(root) else p.name for p in scopes]),indent=2)+'\n')
allfiles=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in allfiles))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='BETA_integrated_PEPC_provenance')
print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
