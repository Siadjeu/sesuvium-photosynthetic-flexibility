#!/usr/bin/env bash
#SBATCH --job-name=dit_close_almt
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=dit_close_almt_%j.out
#SBATCH --error=dit_close_almt_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
python3 - <<'PY'
import urllib.request,time,csv,gzip,hashlib,json,os,re,shutil,subprocess,tarfile
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True);base=root/'08_focused_families/validation'
def need(v,msg):
 if not v:raise SystemExit(msg)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def verify(p):
 for line in (p/'SHA256SUMS.txt').read_text().splitlines():
  h,n=line.split(maxsplit=1);f=(p/n.lstrip('*')).resolve()
  need(f.is_relative_to(p) and f.is_file() and sha(f)==h,'Checksum mismatch: '+n)
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta(p):
 d={};headers={};key=None
 with p.open() as f:
  for line in f:
   if line.startswith('>'):
    key=line[1:].split()[0];need(key not in d,'Duplicate FASTA ID '+key);d[key]='';headers[key]=line[1:].strip()
   elif line.strip():
    need(key is not None,'Invalid FASTA');d[key]+=''.join(line.split()).upper()
 return d,headers
def norm(s):return s[:-1] if s.endswith('*') else s
def valid(s):return bool(s) and set(s)<=set('ACDEFGHIKLMNPQRSTVWYBXZJUO') and not set(s)<=set('ACGTN')
def release(step,env):
 pp=[root/os.environ[env]] if os.environ.get(env) else [p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
 need(len(pp)==1,'Set '+env+' to the intended release');p=pp[0].resolve(strict=True);verify(p);return p




src=release('step25kz17.58','DIT_ANNOTATION_RELEASE')
old=read(src/'source_prior_01_integrated_locus_evidence.tsv');need(len(old)==946 and not any(r['wang_family']=='DIT2' for r in old),'Unexpected prior ledger')
rows=read(src/'01_all_candidate_domain_review.tsv');exceptions=read(src/'10_combined_exception_summary.tsv');checks=read(src/'08_protein_CDS_consistency.tsv');tree=read(src/'source_02_DIT_candidate_reconciliation.tsv');tm={r['protein_id']:r for r in tree};need(len(rows)==38 and len(exceptions)==9,'Unexpected review scope')
for r in exceptions:
 k=r['protein_id'];need(r['category']=='ANNOTATION_CONSISTENT_MODEL_CAVEAT_REMAINS','Exception status changed '+k);need(any(c['candidate_id']==k and c['protein_matches_candidate']=='True' and c['translation_matches_candidate']=='True' and c['frame_remainder']=='0' and c['internal_stops']=='0' and c['ambiguous_amino_acids']=='0' for c in checks),'Consistency evidence missing '+k)
new=[];comparators=[];ff=list(old[0])
for r in rows:
 k=r['protein_id'];fam=r['tree_group'];need(tm[k]['stable_across_four_trees']=='True' and tm[k]['screen_family']==fam,'Unstable assignment '+k);flag=r['flags']!='NO_HEURISTIC_EXCEPTION';tier=('PARTIAL' if fam=='DIT2' else 'MODEL_CAVEAT') if flag else 'SUPPORTED'
 record=dict(wang_family=fam,species=r['species'],gene_id=r['gene_id'],protein_id=k,subunit='',evidence_category='PHYLOGENETIC_DOMAIN_SUPPORTED' if not flag else ('DIT2_ASSOCIATED_PARTIAL' if fam=='DIT2' else 'DIT1_MODEL_CAVEAT'),count_tier=tier,source_file='step25kz17.58/01_all_candidate_domain_review.tsv;step25kz17.58/10_combined_exception_summary.tsv;step25kz17.58/source_02_DIT_candidate_reconciliation.tsv',model_or_scope_limitation='BOUNDED_REFERENCE_PANEL;SEARCH_COMPLETENESS_AND_TRANSPORT_FUNCTION_NOT_ESTABLISHED'+(';' +r['flags']+';CDS_CONSISTENT_NOT_BIOLOGICAL_VALIDATION' if flag else ''),exhaustive_count='NOT_ESTABLISHED',C4_CAM_specific_function='NOT_ESTABLISHED')
 (new if fam=='DIT2' else comparators).append(record)
need(Counter(r['count_tier'] for r in new)=={'SUPPORTED':17,'PARTIAL':4} and Counter(r['count_tier'] for r in comparators)=={'SUPPORTED':12,'MODEL_CAVEAT':5},'Unexpected evidence categories')
ledger=old+new;need(len({(r['wang_family'],r['species'],r['gene_id']) for r in ledger})==len(ledger),'Duplicate family/locus')
out=base/'step25kz17.59'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'ALMT_evidence').mkdir()
def table(n,rr,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
table('01_integrated_locus_evidence.tsv',ledger,ff);table('02_DIT2_closed_evidence_ledger.tsv',new,ff);table('03_DIT1_comparators_NOT_WANG_COUNTS.tsv',comparators,ff)
summary=[]
for fam,rr in [('DIT2',new),('DIT1_COMPARATOR',comparators)]:
 for sp in sorted({r['species'] for r in old}):
  subset=[r for r in rr if r['species']==sp];summary.append(dict(family=fam,species=sp,supported_loci=sum(r['count_tier']=='SUPPORTED' for r in subset),partial_loci=sum(r['count_tier']=='PARTIAL' for r in subset),model_caveat_loci=sum(r['count_tier']=='MODEL_CAVEAT' for r in subset),final_exhaustive_count='NA'))
table('04_DIT_species_evidence_summary.tsv',summary,list(summary[0]))
# Read-only ALMT inventory within focused-family tree. Never follow directory symlinks.
searchroot=root/'08_focused_families';inventory=[];budget=25*1024*1024;copied=0
for folder,dirs,files in os.walk(searchroot,followlinks=False):
 dirs[:]=[d for d in dirs if not d.startswith('.') and d not in {'work','tmp','temp'} and (Path(folder)/d).resolve()!=out.parent.resolve()]
 for name in sorted(files):
  p=Path(folder)/name;rel=p.relative_to(root)
  if 'ALMT' not in str(rel).upper():continue
  if p.is_symlink() or not p.is_file():continue
  size=p.stat().st_size;eligible=p.suffix.lower() in {'.tsv','.csv','.txt','.json','.faa','.fa','.fasta','.treefile','.nwk','.log','.domtbl','.tbl'};digest=sha(p) if size<=10*1024*1024 else 'NOT_HASHED_SIZE_LIMIT';dest='';status='METADATA_ONLY'
  if eligible and size<=5*1024*1024 and copied+size<=budget:
   target=out/'ALMT_evidence'/('E%04d_'%(len(inventory)+1)+name);shutil.copyfile(p,target);dest=str(target.relative_to(out));copied+=size;status='COPIED_COMPLETE_FILE'
  inventory.append(dict(source_file=str(rel),bytes=size,sha256=digest,collection_status=status,bundle_file=dest))
table('05_ALMT_existing_file_inventory.tsv',inventory,['source_file','bytes','sha256','collection_status','bundle_file'])
coverage=read(src/'source_07_refreshed_family_coverage.tsv');need(len(coverage)==23,'Expected 23 families')
for r in coverage:
 fam=r['wang_family'];r['integrated_records']=str(sum(x['wang_family']==fam for x in ledger))
 if fam=='DIT2':r.update(current_status='BOUNDED_REVIEW_CLOSED_LEDGER_INTEGRATED',search_completeness='NOT_ESTABLISHED',priority='4',next_action='RETAIN_17_SUPPORTED_AND_4_PARTIAL_RECORDS;EXHAUSTIVE_COUNTS_PENDING')
 if fam=='ALMT':r.update(current_status='EXISTING_FILES_INVENTORIED_NOT_VALIDATED' if inventory else 'NO_ALMT_NAMED_FILES_LOCATED_IN_BOUNDED_SCAN',priority='1',next_action='RECONCILE_ALMT_REFERENCES_AND_CANDIDATES_AGAINST_CORRECTED_INPUT;DOMAIN_AND_SPECIFICITY_VALIDATION')
table('06_refreshed_family_coverage.tsv',coverage,list(coverage[0]));shutil.copyfile(src/'source_prior_07_final_Wang_counts_PENDING.tsv',out/'07_final_Wang_counts_PENDING.tsv');unresolved=[r for r in coverage if int(r['integrated_records'])==0];table('08_families_without_integrated_ledger.tsv',unresolved,list(coverage[0]))
bylocus=defaultdict(list)
for r in ledger:bylocus[r['species'],r['gene_id']].append(r)
conflicts=[dict(species=sp,gene_id=g,families=';'.join(sorted({r['wang_family'] for r in rr}))) for (sp,g),rr in bylocus.items() if len({r['wang_family'] for r in rr})>1];table('09_cross_family_overlap_review.tsv',conflicts,['species','gene_id','families'])
for n in ['01_all_candidate_domain_review.tsv','10_combined_exception_summary.tsv','08_protein_CDS_consistency.tsv','11_existing_domain_hits.tsv','12_existing_model_coverage.tsv','source_01_existing_tree_DIT_group_audit.tsv','source_02_DIT_candidate_reconciliation.tsv']:shutil.copyfile(src/n,out/('evidence_'+n))
(out/'README.txt').write_text('DIT bounded review closed:21 DIT2 evidence loci (17 supported,4 partial).17 DIT1 comparators (12 supported,5 model caveats) kept separate from Wang counts. Partial models have consistent source CDS but lack completeness; no pseudogene or functional claims. Prior946 records unchanged plus21 DIT2=967 records across11 families. Final exhaustive counts remain pending. ALMT filename/path inventory restricted to08_focused_families, excluding work/tmp/temp and this step; no symlink traversal. Copies complete eligible files<=5MiB within25MiB budget; large/uncollected files remain in inventory. Absence in filename scan is not biological absence. Historical coverage fields retained; new inventory is separate. No searches, trees or repairs rerun.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step58_manifest_sha256=sha(src/'SHA256SUMS.txt'),prior_ledger_sha256=sha(src/'source_prior_01_integrated_locus_evidence.tsv')),indent=2)+'\n');qc=dict(status='DIT_CLOSED_ALMT_INVENTORY_READY',integrated_records=len(ledger),integrated_families=len({r['wang_family'] for r in ledger}),DIT2_tiers=dict(Counter(r['count_tier'] for r in new)),DIT1_comparator_tiers=dict(Counter(r['count_tier'] for r in comparators)),ALMT_files_located=len(inventory),ALMT_files_copied=sum(r['collection_status']=='COPIED_COMPLETE_FILE' for r in inventory),families_without_integrated_ledger=len(unresolved),cross_family_overlap_flags=len(conflicts),final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out);bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='DIT_closed_ALMT_inventory/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
