#!/usr/bin/env bash
#SBATCH --job-name=dit_reconcile
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=dit_reconcile_%j.out
#SBATCH --error=dit_reconcile_%j.err
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





closed=release('step25kz17.55','COUNTING_LEDGER_RELEASE');screen=release('step25kz17.46','TRANSPORTER_RECOVERY_RELEASE');trees=release('step25kz17.48','TDT_TREE_RELEASE');sensitivity=release('step25kz17.49','TD029_REVIEW_RELEASE')
need(json.loads((trees/'PROVENANCE.json').read_text())['step46_manifest_sha256']==sha(screen/'SHA256SUMS.txt'),'Tree/screen provenance mismatch');need(json.loads((sensitivity/'PROVENANCE.json').read_text())['source_manifest_sha256']==sha(trees/'SHA256SUMS.txt'),'Sensitivity provenance mismatch')
rows=[r for r in read(screen/'07_all_candidate_evidence_REVIEW_ONLY.tsv') if r['best_panel_family'] in {'DIT1','DIT2'}];need(Counter(r['best_panel_family'] for r in rows)=={'DIT1':18,'DIT2':21},'Unexpected DIT screen set');need(all(r['mapping_status']=='UNAMBIGUOUS' for r in rows),'Ambiguous mapping')
seqs,_=fasta(screen/'05_all_discovered_and_prior_candidates.faa');refs,_=fasta(screen/'02_reference_panel.faa');meta=read(trees/'02_tree_id_map.tsv');byprotein={r['protein_id']:r for r in meta};refids={r['tree_id'] for r in meta if r['role']=='CURATED_REFERENCE'};need(refids=={'REF__Q8LG88','REF__Q9LXV3','REF__Q9FMF7'},'Reference panel changed')
def edges(path):
 t=re.findall(r'[^\s(),:;]+|[(),:;]',path.read_text());pos=0;ee=[]
 def node():
  nonlocal pos
  label=''
  if t[pos]=='(':
   pos+=1;z=set()
   while True:
    z|=node()
    if t[pos]==',':pos+=1;continue
    need(t[pos]==')','Malformed tree');pos+=1;break
   if t[pos] not in {':',',',')',';'}:label=t[pos];pos+=1
  else:z={t[pos]};pos+=1
  if t[pos]==':':pos+=2
  ee.append((z,label));return z
 universe=node();need(t[pos]==';','Invalid tree ending');return universe,[(z,l) for a,l in ee for z in (a,universe-a) if 1<len(z)<len(universe)-1]
out=base/'step25kz17.56'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
assignments={};splitrows=[];presence={};sources=[]
for source,label in [(trees,'original'),(sensitivity,'without_TD029')]:
 for mode in ['main_without_short','sensitivity_with_short']:
  name=label+'__'+mode;path=source/mode/'tree.treefile';universe,ee=edges(path);presence[name]=universe;sources.append(dict(analysis=name,tree_sha256=sha(path)))
  for fam,anchor in [('DIT1','REF__Q9LXV3'),('DIT2','REF__Q9FMF7')]:
   options=[(z,l) for z,l in ee if anchor in z and z&refids=={anchor}]
   if options:z,support=max(options,key=lambda a:len(a[0]));z=z-{anchor}
   else:z=set();support='NO_REFERENCE_EXCLUSIVE_EDGE'
   assignments[name,fam]=z;splitrows.append(dict(analysis=name,family=fam,candidates_on_reference_exclusive_side=len(z),support_SHaLRT_UFBoot=support,tree_ids=';'.join(sorted(z)),interpretation='TOPOLOGICAL_EVIDENCE_NOT_AUTOMATIC_ACCEPTANCE'))
table('01_existing_tree_DIT_group_audit.tsv',splitrows,list(splitrows[0]));decisions=[];excluded='SP_Hex__gene_Spor0030220.1'
for r in rows:
 k=r['protein_id'];need(k in byprotein,'Candidate missing tree map');tid=byprotein[k]['tree_id'];a={name:next((fam for fam in ['DIT1','DIT2'] if tid in assignments[name,fam]),'OUTSIDE_DIT_REFERENCE_GROUPS') if tid in ids else 'NOT_IN_TREE' for name,ids in presence.items()};available=[v for v in a.values() if v!='NOT_IN_TREE'];stable=len(set(available))==1 and len(available)==4
 status='REVIEW_REQUIRED'
 if k==excluded:status='PREVIOUSLY_EXCLUDED_SUPPORTED_TDT_DIT_SET'
 elif stable and available[0] in {'DIT1','DIT2'}:status='CONSISTENT_TREE_GROUP_REQUIRES_DOMAIN_COMPLETENESS_REVIEW'
 decisions.append(dict(protein_id=k,species=r['species'],gene_id=r['gene_id'],tree_id=tid,screen_family=r['best_panel_family'],screen_flags=r['flags'],**a,stable_across_four_trees=stable,review_status=status,domain_completeness='NOT_VALIDATED_FOR_FULL_CURRENT_DIT_SET',final_count='NA'))
table('02_DIT_candidate_reconciliation.tsv',decisions,list(decisions[0]));active=[r for r in rows if r['protein_id']!=excluded];need(len(active)==38,'Expected 38 retained candidates');writefa('03_DIT_candidates_REVIEW_ONLY.faa',{r['protein_id']:seqs[r['protein_id']] for r in active});writefa('04_DIT_TDT_reference_controls.faa',{k:refs[k] for k in sorted(refids)});writefa('05_excluded_TD029.faa',{excluded:seqs[excluded]})
shutil.copyfile(screen/'01_reference_definitions.tsv',out/'04_source_reference_definitions.tsv');table('06_source_DIT_screen_evidence.tsv',rows,list(rows[0]))
coverage=read(closed/'06_refreshed_family_coverage.tsv')
for r in coverage:
 if r['wang_family']=='DIT2':r.update(current_status='EXISTING_SCREEN_AND_COMPARATIVE_PHYLOGENY_AVAILABLE_NOT_INTEGRATED',integrated_records='0',search_completeness='NOT_ESTABLISHED',priority='1',next_action='DOMAIN_COMPLETENESS_VALIDATION_OF_38_DIT_CANDIDATES_WITH_DIT1_DIT2_TDT_CONTROLS;RETAIN_TREE_SUPPORT_LIMITATIONS')
table('07_refreshed_family_coverage.tsv',coverage,list(coverage[0]))
for n in ['01_integrated_locus_evidence.tsv','07_final_Wang_counts_PENDING.tsv','04_all_carrier_candidates_disposition.tsv','05_remaining_carrier_candidates.tsv']:shutil.copyfile(closed/n,out/('prior_'+n))
# Explicit outstanding work, no implied absence when older evidence is incomplete.
requirements=[dict(item='DOMAIN_COMPLETENESS',scope='38 retained DIT1/DIT2 candidates plus three controls',action='Reuse or obtain full Pfam evidence for exact current sequences; inspect overlaps, terminal loss and additional domains'),dict(item='PHYLOGENETIC_SUPPORT',scope='Four existing TDT-versus-DIT trees',action='Inspect DIT1/DIT2 edge supports and membership; topology consistency alone is not strong support'),dict(item='HISTORICAL_EXCEPTION',scope='PO__evm.TU.LG08.1714',action='Retain historical long-model/additional-domain caveat and reassess against current sequence'),dict(item='COUNTING_SCOPE',scope='Wang DIT2; DIT1 comparators retained separately',action='Do not pool DIT1 with DIT2 or count the excluded TD029')];table('08_next_validation_requirements.tsv',requirements,list(requirements[0]))
(out/'README.txt').write_text('Existing-evidence reconciliation, not a new discovery or validation run. Original 39 best-panel DIT candidates:18 DIT1 and21 DIT2. TD029 retained as previously excluded;38 candidates exported for domain/completeness review. Four existing unrooted tree edge groups examined using reference-exclusive sides; support values retained without forced acceptance. Stability alone is not strong support. DIT2 stale coverage corrected; historical path counts are inherited audit fields, not a current inventory. DIT1 remains a comparator, not pooled into Wang DIT2. Full current-set domain completeness unresolved; prior PO LG08.1714 architecture caveat retained as review priority. Integrated946-record ledger and final-count matrix unchanged. No searches, alignments, trees or annotation repairs rerun.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step55_manifest_sha256=sha(closed/'SHA256SUMS.txt'),step46_manifest_sha256=sha(screen/'SHA256SUMS.txt'),step48_manifest_sha256=sha(trees/'SHA256SUMS.txt'),step49_manifest_sha256=sha(sensitivity/'SHA256SUMS.txt'),trees=sources),indent=2)+'\n');qc=dict(status='DIT_EXISTING_EVIDENCE_RECONCILED',screen_candidates=39,previously_excluded=1,retained_for_validation=38,screen_groups=dict(Counter(r['best_panel_family'] for r in active)),review_status_counts=dict(Counter(r['review_status'] for r in decisions)),coverage_corrected='DIT2',ledger_changed=False,final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out);bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='DIT_evidence_reconciliation/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
