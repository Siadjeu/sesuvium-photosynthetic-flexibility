#!/usr/bin/env bash
#SBATCH --job-name=beta_pepc_map
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=beta_pepc_map_%j.out
#SBATCH --error=beta_pepc_map_%j.err
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
integrated=release('step25kz17.37','INTEGRATED_RELEASE')
def choose_file(step,name,env):
 paths=[root/os.environ[env]] if os.environ.get(env) else sorted((base/step).rglob(name))
 groups=defaultdict(list)
 for p in paths:
  if p.is_file():groups[sha(p)].append(p)
 if len(groups)!=1:raise SystemExit('Expected one distinct '+name+'; set '+env+' explicitly. Distinct files='+str(len(groups)))
 return next(iter(groups.values()))[0]
beta=choose_file('step25kz17.7','02_CA_family_specific_review.tsv','BETA_SOURCE')
residue=choose_file('step25kz17.8','01_BETA_CA_residue_completeness_review.tsv','BETA_RESIDUE_SOURCE')
selection=choose_file('step25kz17.12','01_selection_manifest.tsv','BETA_SELECTION_SOURCE')
pep=root/os.environ.get('PEPC_METADATA','08_focused_families/validation/step25kz16.8b/release/01_PEPC_tree_tip_metadata.tsv')
assert pep.is_file(),'PEPC metadata missing; set PEPC_METADATA'
prior=integrated/'sources/PRIOR_WANG__04_gene_evidence.tsv'
assert prior.is_file()
species=['AT','BV','MC','PO','SLP_Hap1','SP_Hex','SP_Tet','SS']
bs=read(beta);res={r['protein_id']:r for r in read(residue)};sel={r['locus_id']:r for r in read(selection)}
assert len(bs)==110 and len({r['protein_id'] for r in bs})==110
out=base/'step25kz17.38'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'sources').mkdir()
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
provenance=[]
for label,p in [('BETA_SCREEN',beta),('BETA_RESIDUE',residue),('BETA_SELECTION',selection),('PEPC_METADATA',pep),('PRIOR_EVIDENCE',prior)]:
 dest=out/'sources'/(label+'__'+p.name);shutil.copy2(p,dest)
 provenance.append(dict(label=label,source=str(p.relative_to(root)) if p.is_relative_to(root) else p.name,sha256=sha(p),snapshot=str(dest.relative_to(out))))
# Exact aliases only: no suffix stripping, no sequence-based locus merging.
index={}
for sp in species:
 p=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv')
 assert p.is_file(),'Missing isoform mapping '+sp
 rr=read(p);ix=defaultdict(set)
 for r in rr:
  for k in ['gene_id','transcript_id','protein_id','cds_id']:
   if r.get(k) and r.get('gene_id'):ix[r[k]].add(r['gene_id'])
 index[sp]=ix;provenance.append(dict(label=sp+'_MAPPING',source=str(p.relative_to(root)),sha256=sha(p),snapshot='NOT_COPIED_FULL_MAPPING'))
rows=[];issues=[];evidence=[]
def mapid(pid,sp):
 if sp not in index or not pid.startswith(sp+'__'):return set()
 return index[sp].get(pid.split('__',1)[1],set())
for r in bs:
 pid=r['protein_id'];sp=r['species'];s=sel.get(pid,{});a=res.get(pid,{});genes=mapid(pid,sp);sc=r['screen_class']
 if sc=='BETA_SUPPORTED_REVIEW':
  cat='BETA_SUPPORTED_WITH_REVIEW_CAVEAT' if s.get('reason')=='HOLD_UNRESOLVED_RESIDUES' else 'BETA_SUPPORTED_SCREEN_EVIDENCE'
 elif 'BETA' in sc and 'PARTIAL' in sc:cat='BETA_PARTIAL_UNRESOLVED'
 else:cat='OUTSIDE_SUPPORTED_BETA_SET'
 rows.append(dict(wang_family='beta-CA',species=sp,protein_id=pid,gene_id=';'.join(sorted(genes)),mapping_status='UNAMBIGUOUS' if len(genes)==1 else 'UNRESOLVED',evidence_category=cat,original_classification=sc,tree_or_alignment_inclusion=s.get('alignment_selection','NOT_RECORDED'),tree_tip=s.get('tip_id',''),sequence_selection_reason=s.get('reason',''),residue_check=a.get('reference_residue_check','NOT_ASSESSED'),prior_reviewer_decision=r.get('prior_reviewer_decision',''),latest_recorded_reviewer_decision=r.get('new_reviewer_decision',''),count_decision='PENDING_LOCUS_EVIDENCE_REVIEW',limitation='Specificity screen is not a final functional count; tree exclusion is not automatic family exclusion.'))
 evidence.append(dict(wang_family='beta-CA',protein_id=pid,source='BETA_SCREEN',record_json=json.dumps(r,sort_keys=True)))
# Union evidence records, not inferred accepted memberships.
pe=defaultdict(list)
for r in read(prior):
 if r.get('wang_category')=='PEPC':pe[r['protein_id']].append(('PRIOR_EVIDENCE',r))
for r in read(pep):
 if r.get('species')=='REF' or r.get('tip_id','').startswith('REF'):continue
 pid=r.get('tip_id','')
 if not any(pid.startswith(s+'__') for s in species):
  issues.append(dict(family='PEPC',protein_id=pid,issue='NO_EXACT_PROJECT_ID_IN_TREE_METADATA'));continue
 pe[pid].append(('PEPC_TREE_METADATA',r))
for pid,rr in sorted(pe.items()):
 sp=pid.split('__',1)[0];genes=mapid(pid,sp);tree=[r for source,r in rr if source=='PEPC_TREE_METADATA'];cl=';'.join(sorted({r.get('classification','') for r in tree}))
 rows.append(dict(wang_family='PEPC',species=sp,protein_id=pid,gene_id=';'.join(sorted(genes)),mapping_status='UNAMBIGUOUS' if len(genes)==1 else 'UNRESOLVED',evidence_category='PEPC_TREE_EVIDENCE_REQUIRES_MEMBERSHIP_REVIEW' if tree else 'PEPC_HISTORICAL_CANDIDATE_REQUIRES_REVIEW',original_classification=cl,tree_or_alignment_inclusion='TREE_METADATA_MEMBER' if tree else 'NOT_RECORDED',tree_tip=pid if tree else '',sequence_selection_reason='',residue_check='NOT_ASSESSED_HERE',prior_reviewer_decision='',latest_recorded_reviewer_decision='',count_decision='PENDING_LOCUS_EVIDENCE_REVIEW',limitation='Tree inclusion and legacy candidate status do not establish complete family inventory or C4/CAM-specific function.'))
 for source,r in rr:evidence.append(dict(wang_family='PEPC',protein_id=pid,source=source,record_json=json.dumps(r,sort_keys=True)))
for r in rows:
 if r['mapping_status']!='UNAMBIGUOUS':issues.append(dict(family=r['wang_family'],protein_id=r['protein_id'],issue='UNRESOLVED_OR_AMBIGUOUS_LOCUS_MAPPING'))
table('01_BETA_CA_protein_evidence.tsv',[r for r in rows if r['wang_family']=='beta-CA'],list(rows[0]));table('02_PEPC_protein_evidence.tsv',[r for r in rows if r['wang_family']=='PEPC'],list(rows[0]))
# Collapse only exact mapped loci. Preserve mixed evidence and all protein IDs.
groups=defaultdict(list)
for r in rows:
 key=(r['wang_family'],r['species'],r['gene_id'] if r['mapping_status']=='UNAMBIGUOUS' else 'UNRESOLVED:'+r['protein_id']);groups[key].append(r)
loci=[]
for (fam,sp,gene),rr in sorted(groups.items()):
 cats=sorted({r['evidence_category'] for r in rr});mapped=not gene.startswith('UNRESOLVED:')
 loci.append(dict(wang_family=fam,species=sp,gene_id=gene if mapped else '',protein_ids=';'.join(sorted({r['protein_id'] for r in rr})),mapping_status='UNAMBIGUOUS' if mapped else 'UNRESOLVED',evidence_categories=';'.join(cats),protein_records=len(rr),category_conflict='YES' if len(cats)>1 else 'NO',tree_or_alignment_status=';'.join(sorted({r['tree_or_alignment_inclusion'] for r in rr})),count_decision='PENDING_REVIEW',confirmed_function='NOT_ESTABLISHED'))
table('03_locus_reconciliation_ledger.tsv',loci,list(loci[0]));table('04_evidence_records.tsv',evidence,['wang_family','protein_id','source','record_json']);table('05_reconciliation_issues.tsv',issues,['family','protein_id','issue'])
summary=[]
for fam in ['beta-CA','PEPC']:
 for sp in species:
  for cat in sorted({r['evidence_category'] for r in rows if r['wang_family']==fam}):
   rr=[r for r in rows if r['wang_family']==fam and r['species']==sp and r['evidence_category']==cat]
   summary.append(dict(wang_family=fam,species=sp,evidence_category=cat,protein_records=len(rr),mapped_loci=len({r['gene_id'] for r in rr if r['mapping_status']=='UNAMBIGUOUS'}),unresolved_records=sum(r['mapping_status']!='UNAMBIGUOUS' for r in rr),count_type='EVIDENCE_ONLY_NOT_FINAL'))
table('06_species_evidence_summary.tsv',summary,list(summary[0]))
# Keep five-family completed integration unchanged; add provisional ledgers separately.
shutil.copy2(integrated/'01_integrated_locus_evidence.tsv',out/'07_prior_completed_family_ledger_UNCHANGED.tsv')
coverage=read(integrated/'05_refreshed_family_coverage.tsv')
for r in coverage:
 if r['wang_family'] in ['beta-CA','PEPC']:
  r['current_status']='LOCUS_RECONCILIATION_COMPLETE_DECISIONS_PENDING';r['next_action']='REVIEW_LOCUS_LEDGER_AND_CATEGORY_CONFLICTS_BEFORE_COUNT_INTEGRATION'
table('08_refreshed_family_coverage.tsv',coverage,list(coverage[0]))
table('09_source_provenance.tsv',provenance,['label','source','sha256','snapshot'])
(out/'README.txt').write_text('BETA-CA and PEPC locus reconciliation using existing evidence only. All source classifications and prior reviewer decisions are preserved; earlier acceptance of alpha-CA as beta-CA is not propagated over the later specificity screen. BETA main-alignment selection is not final family acceptance. The short AT alignment isoform is represented through the existing selection record, not counted as another locus.\nPEPC prior evidence and current tree metadata are retained together. Presence in a tree does not establish family completeness or subtype function. Exact mapping aliases are used without guessed suffix removal. Only unambiguous identical species/family/gene keys are collapsed; identical sequences at distinct loci remain distinct.\nSummary categories can overlap when multiple protein records at a locus have different evidence; inspect category_conflict before summing. Unresolved mappings are not counted as resolved loci. No sequence analysis rerun, no final count promotion, no source changes. The earlier five-family integrated ledger remains unchanged.\n')
(out/'QC.json').write_text(json.dumps(dict(status='BETA_PEPC_LOCUS_RECONCILIATION_COMPLETE_REVIEW_REQUIRED',protein_records=len(rows),locus_ledger_rows=len(loci),unresolved_mappings=sum(r['mapping_status']!='UNAMBIGUOUS' for r in rows),category_conflicts=sum(r['category_conflict']=='YES' for r in loci),evidence_categories=dict(Counter(r['evidence_category'] for r in rows)),analysis_rerun=False,final_counts=False),indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(integrated_manifest_sha256=sha(integrated/'SHA256SUMS.txt')),indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='BETA_PEPC_locus_reconciliation')
print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
