#!/usr/bin/env bash
#SBATCH --job-name=almt_close_ppt
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=almt_close_ppt_%j.out
#SBATCH --error=almt_close_ppt_%j.err
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




src=release('step25kz17.63','ALMT_GROUP_REVIEW_RELEASE');old=read(src/'source_prior_01_integrated_locus_evidence.tsv');need(len(old)==967 and not any(r['wang_family']=='ALMT' for r in old),'Unexpected prior ledger')
rows=read(src/'09_candidate_evidence_categories.tsv');iso=read(src/'10_SS_alternative_review_outcomes.tsv');need(len(rows)==142,'Expected142 candidate records');need({r['alternative_protein_id'] for r in iso}=={'SS__g16509.t1','SS__g17774.t1'} and all(r['review_outcome']=='NO_ALMT_SUPPORT_IN_ALTERNATIVE' for r in iso),'Alternative isoform conclusions changed')
new=[];ff=list(old[0]);outside={'SLP_Hap1__SP4G00453','SS__g16509.t2','SS__g17774.t2'};short='SP_Hex__gene_Spor0056610.1'
for r in rows:
 k=r['protein_id'];cat=r['evidence_category'];flags=r['retained_completeness_flags'];lim='BOUNDED_REFERENCE_PANEL;FUNCTION_AND_EXHAUSTIVE_COUNTS_NOT_ESTABLISHED;'+flags+';'+r['architecture_category']
 if cat=='MAIN_TREE_DOMAIN_EVIDENCE':
  need(r['smallest_stable_reference_group'],'Main candidate lacks supported stable group');tier='SUPPORTED' if flags=='NO_RECALIBRATED_HEURISTIC_EXCEPTION' and r['architecture_category']=='NO_TARGETED_TLD_ARCHITECTURE_CAVEAT' else 'MODEL_CAVEAT';category='PHYLOGENETIC_DOMAIN_SUPPORTED' if tier=='SUPPORTED' else 'ALMT_ASSOCIATED_COMPLETENESS_OR_ARCHITECTURE_CAVEAT'
 elif cat=='PARTIAL_SENSITIVITY_ONLY':tier='PARTIAL';category='ALMT_DOMAIN_PARTIAL_SENSITIVITY_ONLY';lim+=';JOINT_FRAGMENT_PLACEMENT_UNCERTAIN'
 elif k in outside:tier='EXCLUDED';category='OUTSIDE_SUPPORTED_ALMT_SET';lim+=';NO_PF11744_OR_POSITIVE_REFERENCE_SIMILARITY;NO_GENOMIC_ABSENCE_CLAIM'
 elif k==short:tier='UNRESOLVED';category='SHORT_SIMILARITY_ONLY_NO_PF11744';lim+=';131AA;NO_TREE_PLACEMENT'
 else:need(False,'Unexpected outside-tree record '+k)
 new.append(dict(wang_family='ALMT',species=r['species'],gene_id=r['gene_id'],protein_id=k,subunit='',evidence_category=category,count_tier=tier,source_file='step25kz17.63/09_candidate_evidence_categories.tsv;step25kz17.63/10_SS_alternative_review_outcomes.tsv',model_or_scope_limitation=lim,exhaustive_count='NOT_ESTABLISHED',C4_CAM_specific_function='NOT_ESTABLISHED'))
counts=Counter(r['count_tier'] for r in new);need(counts=={'SUPPORTED':95,'MODEL_CAVEAT':32,'PARTIAL':11,'EXCLUDED':3,'UNRESOLVED':1},'Unexpected category totals');ledger=old+new;need(len({(r['wang_family'],r['species'],r['gene_id']) for r in ledger})==len(ledger),'Duplicate family/locus')
out=base/'step25kz17.64'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'PPT_evidence').mkdir()
def table(n,rr,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
table('01_integrated_locus_evidence.tsv',ledger,ff);table('02_ALMT_closed_evidence_ledger.tsv',new,ff)
summary=[]
for sp in sorted({r['species'] for r in old}):
 rr=[r for r in new if r['species']==sp];summary.append(dict(species=sp,**{tier.lower()+'_records':sum(r['count_tier']==tier for r in rr) for tier in ['SUPPORTED','MODEL_CAVEAT','PARTIAL','EXCLUDED','UNRESOLVED']},total_evidence_records=len(rr),final_exhaustive_count='NA'))
table('03_ALMT_species_evidence_summary.tsv',summary,list(summary[0]))
# Read-only PPT inventory within focused-family tree. Never follow directory symlinks.
searchroot=root/'08_focused_families';inventory=[];budget=25*1024*1024;copied=0
for folder,dirs,files in os.walk(searchroot,followlinks=False):
 dirs[:]=[d for d in dirs if not d.startswith('.') and d not in {'work','tmp','temp'} and (Path(folder)/d).resolve()!=out.parent.resolve()]
 for name in sorted(files):
  p=Path(folder)/name;rel=p.relative_to(root)
  if not re.search(r'(?<![A-Z0-9])PPT(?:[12])?(?![A-Z0-9])',str(rel).upper()):continue
  if p.is_symlink() or not p.is_file():continue
  size=p.stat().st_size;eligible=p.suffix.lower() in {'.tsv','.csv','.txt','.json','.faa','.fa','.fasta','.treefile','.nwk','.log','.domtbl','.tbl'};digest=sha(p) if size<=10*1024*1024 else 'NOT_HASHED_SIZE_LIMIT';dest='';status='METADATA_ONLY'
  if eligible and size<=5*1024*1024 and copied+size<=budget:
   target=out/'PPT_evidence'/('E%04d_'%(len(inventory)+1)+name);shutil.copyfile(p,target);dest=str(target.relative_to(out));copied+=size;status='COPIED_COMPLETE_FILE'
  inventory.append(dict(source_file=str(rel),bytes=size,sha256=digest,collection_status=status,bundle_file=dest))
table('05_PPT_existing_file_inventory.tsv',inventory,['source_file','bytes','sha256','collection_status','bundle_file'])
coverage=read(src/'source_prior_06_refreshed_family_coverage.tsv');need(len(coverage)==23,'Expected23 families')
for r in coverage:
 fam=r['wang_family'];r['integrated_records']=str(sum(x['wang_family']==fam for x in ledger))
 if fam=='ALMT':r.update(current_status='BOUNDED_REVIEW_CLOSED_LEDGER_INTEGRATED',search_completeness='NOT_ESTABLISHED',priority='4',next_action='RETAIN_SUPPORTED_CAVEAT_PARTIAL_EXCLUDED_UNRESOLVED_CATEGORIES;EXHAUSTIVE_COUNTS_PENDING')
 if fam=='PPT':r.update(current_status='EXISTING_FILES_INVENTORIED_NOT_VALIDATED' if inventory else 'NO_PPT_NAMED_FILES_LOCATED_IN_BOUNDED_SCAN',priority='1',next_action='REVIEW_PPT_REFERENCE_DEFINITION_AND_EXISTING_CANDIDATE_PROVENANCE')
table('06_refreshed_family_coverage.tsv',coverage,list(coverage[0]));shutil.copyfile(src/'source_prior_07_final_Wang_counts_PENDING.tsv',out/'07_final_Wang_counts_PENDING.tsv');unresolved=[r for r in coverage if int(r['integrated_records'])==0];table('08_families_without_integrated_ledger.tsv',unresolved,list(coverage[0]))
bylocus=defaultdict(list)
for r in ledger:bylocus[r['species'],r['gene_id']].append(r)
conflicts=[dict(species=sp,gene_id=g,families=';'.join(sorted({r['wang_family'] for r in rr})),tiers=';'.join(sorted({r['count_tier'] for r in rr}))) for (sp,g),rr in bylocus.items() if len({r['wang_family'] for r in rr})>1];table('09_cross_family_overlap_review.tsv',conflicts,['species','gene_id','families','tiers'])
for n in ['07_backbone_split_stability.tsv','08_stable_reference_groups.tsv','09_candidate_evidence_categories.tsv','10_SS_alternative_review_outcomes.tsv','22_protein_CDS_consistency.tsv']:shutil.copyfile(src/n,out/('evidence_'+n))
(out/'README.txt').write_text('Bounded ALMT review closed:95 supported,32 main-tree model/completeness caveats,11 partial sensitivity-only,3 excluded from supported set,1 unresolved short similarity-only record. All142 historical loci retained once; categories are not summed into a final gene count. No functional subtype or exhaustive-count claims. Alternative SS isoforms do not restore ALMT evidence. Exclusion is bounded evidence, not genomic absence or pseudogenization. Previous967 records unchanged plus142=1109 evidence records across12 families. Final Wang matrix remains pending. PPT inventory uses standalone PPT/PPT1/PPT2 path tokens within08_focused_families, skips work/tmp/temp and symlinks; absence in this bounded scan does not establish absence of evidence elsewhere. Complete small files collected under25MiB cap; larger files retained as metadata. No analyses rerun.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step63_manifest_sha256=sha(src/'SHA256SUMS.txt'),prior_ledger_sha256=sha(src/'source_prior_01_integrated_locus_evidence.tsv')),indent=2)+'\n');qc=dict(status='ALMT_CLOSED_PPT_INVENTORY_READY',integrated_records=len(ledger),integrated_families=len({r['wang_family'] for r in ledger}),ALMT_tiers=dict(counts),PPT_files_located=len(inventory),PPT_files_copied=sum(r['collection_status']=='COPIED_COMPLETE_FILE' for r in inventory),families_without_integrated_ledger=len(unresolved),cross_family_overlap_flags=len(conflicts),final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out);bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='ALMT_closed_PPT_inventory/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
