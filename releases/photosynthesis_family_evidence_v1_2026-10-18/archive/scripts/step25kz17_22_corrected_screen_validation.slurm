#!/usr/bin/env bash
#SBATCH --job-name=corrected_screen
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=corrected_screen_%j.out
#SBATCH --error=corrected_screen_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
command -v python3
command -v diamond
python3 - <<'PY'
import csv,hashlib,json,os,subprocess,tarfile
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
def choose(env,step):
 choices=[root/os.environ[env]] if os.environ.get(env) else sorted(p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file())
 if len(choices)!=1:raise SystemExit('Set '+env+' explicitly; found '+str(len(choices)))
 return choices[0]
cor=choose('CORRECTED_RELEASE','step25kz17.21');old=choose('SCREEN_RELEASE','step25kz17.18')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def verify(d):
 for l in (d/'SHA256SUMS.txt').read_text().splitlines():
  h,n=l.split(maxsplit=1)
  if sha(d/n)!=h:raise SystemExit('Checksum mismatch: '+n)
verify(cor);verify(old)
proteome=cor/'full_outputs/corrected_protein_candidates.faa'
expect=[r['sha256'] for r in read(cor/'08_full_output_checksums.tsv') if r['file']=='full_outputs/corrected_protein_candidates.faa']
if len(expect)!=1 or sha(proteome)!=expect[0]:raise SystemExit('Corrected FASTA checksum mismatch')
out=base/'step25kz17.22'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir();(out/'references').mkdir()
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rr)
def fasta(p):
 h=None;s=[]
 with p.open() as f:
  for l in f:
   if l.startswith('>'):
    if h is not None:yield h,''.join(s)
    h=l[1:].strip();s=[]
   elif h is not None:s.append(l.strip())
 if h is not None:yield h,''.join(s)
def writefa(p,rr):
 with p.open('w') as f:
  for h,s in rr:f.write('>'+h+'\n'+s+'\n')
seqs={}
for h,s in fasta(proteome):
 pid=h.split()[0]
 if pid in seqs:raise SystemExit('Duplicate protein ID')
 seqs[pid]=s
ledger={r['record_id']:r for r in read(cor/'02_retained_record_ledger.tsv') if r['action']=='RETAIN'}
if set(ledger)!=set(seqs):raise SystemExit('Ledger/FASTA mismatch')
for p,s in seqs.items():
 if hashlib.sha256(s.encode()).hexdigest()!=ledger[p]['output_sequence_sha256']:raise SystemExit('Sequence hash mismatch: '+p)
removed={r['tree_id'] for r in read(cor/'01_removed_CDS_provenance.tsv')}
oldrows=read(old/'08_family_subunit_candidates.tsv')
oldseq={h.split()[0]:s for h,s in fasta(old/'11_candidate_sequences.faa')}
reconcile=[]
for pid in sorted({r['protein_id'] for r in oldrows}):
 status='REMOVED_CONFIRMED_CDS' if pid in removed else ('RETAINED_IDENTICAL' if pid in seqs and oldseq.get(pid)==seqs[pid] else 'MISSING_OR_CHANGED')
 reconcile.append(dict(protein_id=pid,status=status))
if any(r['status']=='MISSING_OR_CHANGED' for r in reconcile):raise SystemExit('Unexpected candidate sequence change or missing record')
table('01_previous_candidate_reconciliation.tsv',reconcile,['protein_id','status'])
meta={}
for r in read(old/'01_reference_curation.tsv'):meta[r['accession']]=r
targets={h.split()[0]:s for h,s in fasta(old/'references/curated_targets.faa')}
refs=dict(targets)
for h,s in fasta(old/'references/BACKGROUND_REFERENCE.faa'):refs[h.split()[0]]=s
writefa(out/'references/targets.faa',sorted(targets.items()))
writefa(out/'references/target_and_background.faa',sorted(refs.items()))
table('02_reference_metadata.tsv',list(meta.values()),['accession','assigned_category','subunit','header'])
cpu=os.environ.get('SLURM_CPUS_PER_TASK','4');commands=[]
def run(args,name):
 commands.append([a.replace(str(root),'<PROJECT_ROOT>') for a in args])
 with (out/name).open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
run(['diamond','version'],'diamond_version.txt')
db=out/'work/corrected'
run(['diamond','makedb','--in',str(proteome),'--db',str(db),'--threads',cpu],'03_makedb.log')
fmt=['qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen']
def search(query,db,output,log):
 run(['diamond','blastp','--query',str(query),'--db',str(db),'--out',str(output),'--threads',cpu,'--sensitive','--evalue','1e-5','--max-target-seqs','0','--outfmt','6']+fmt,log)
hits=out/'04_corrected_target_hits.tsv'
search(out/'references/targets.faa',db,hits,'04_search.log')
def hitrows(p):
 with p.open() as f:return list(csv.DictReader(f,fieldnames=fmt,delimiter='\t'))
best={}
for r in hitrows(hits):
 ref=meta[r['qseqid']];pid=r['sseqid']
 if pid not in seqs:raise SystemExit('Unexpected hit ID')
 key=(ref['assigned_category'],ref['subunit'],pid)
 if key not in best or float(r['bitscore'])>float(best[key]['bitscore']):best[key]=r
candidateids={k[2] for k in best}
writefa(out/'05_candidate_sequences.faa',[(p,seqs[p]) for p in sorted(candidateids)])
oldkeys={(r['category'],r['subunit'],r['protein_id']) for r in oldrows};newkeys=set(best)
delta=[]
for key in sorted(oldkeys|newkeys):
 c,s,p=key
 status='SHARED' if key in oldkeys and key in newkeys else ('NEW_HIT' if key in newkeys else ('REMOVED_CDS' if p in removed else 'OLD_HIT_NOT_RECOVERED'))
 delta.append(dict(category=c,subunit=s,protein_id=p,status=status))
table('06_old_new_hit_comparison.tsv',delta,['category','subunit','protein_id','status'])
# Reverse search is limited-panel evidence, not reciprocal orthology or complete specificity.
rdb=out/'work/reference_panel'
run(['diamond','makedb','--in',str(out/'references/target_and_background.faa'),'--db',str(rdb),'--threads',cpu],'07_reference_db.log')
reverse=out/'08_candidate_reference_hits.tsv'
if candidateids:search(out/'05_candidate_sequences.faa',rdb,reverse,'08_reverse_search.log')
else:reverse.write_text('')
rev=defaultdict(list)
for r in hitrows(reverse):rev[r['qseqid']].append(r)
evidence=[]
for (category,subunit,pid),r in sorted(best.items()):
 rr=rev[pid]
 # Competing group includes other target families/subunits as well as background.
 same=[h for h in rr if meta[h['sseqid']]['assigned_category']==category and meta[h['sseqid']]['subunit']==subunit]
 other=[h for h in rr if not(meta[h['sseqid']]['assigned_category']==category and meta[h['sseqid']]['subunit']==subunit)]
 same=sorted(same,key=lambda x:-float(x['bitscore']));other=sorted(other,key=lambda x:-float(x['bitscore']))
 top=same[0] if same else None;comp=other[0] if other else None
 delta_score=float(top['bitscore'])-float(comp['bitscore']) if top and comp else ''
 flag='NO_REVERSE_TARGET_HIT' if not top else ('NO_COMPETITOR_HIT_IN_PANEL' if not comp else ('TARGET_HIGHER_IN_PANEL' if delta_score>0 else 'COMPETITOR_TIED_OR_HIGHER'))
 evidence.append(dict(category=category,subunit=subunit,protein_id=pid,species=pid.split('__')[0],
  length_aa=len(seqs[pid]),reference=r['qseqid'],evalue=r['evalue'],bitscore=r['bitscore'],
  reference_span=round((abs(int(r['qend'])-int(r['qstart']))+1)/int(r['qlen']),4),
  protein_span=round((abs(int(r['send'])-int(r['sstart']))+1)/int(r['slen']),4),
  reverse_target=top['sseqid'] if top else '',competitor=comp['sseqid'] if comp else '',
  competitor_category=meta[comp['sseqid']]['assigned_category'] if comp else '',
  competitor_subunit=meta[comp['sseqid']]['subunit'] if comp else '',
  reverse_bitscore_margin=delta_score,limited_panel_status=flag,
  composition_review=ledger[pid]['composition_review'],final_decision='PENDING'))
fields=['category','subunit','protein_id','species','length_aa','reference','evalue','bitscore','reference_span','protein_span','reverse_target','competitor','competitor_category','competitor_subunit','reverse_bitscore_margin','limited_panel_status','composition_review','final_decision']
table('09_specificity_evidence_REVIEW_ONLY.tsv',evidence,fields)
table('10_family_screen_counts_NOT_FINAL.tsv',[dict(category=c,species=s,candidate_ids=len({r['protein_id'] for r in evidence if r['category']==c and r['species']==s})) for c in sorted({meta[a]['assigned_category'] for a in targets}) for s in ['AT','BV','MC','PO','SLP_Hap1','SP_Hex','SP_Tet','SS']],['category','species','candidate_ids'])
table('11_priority_PPDK_RP_review.tsv',[r for r in evidence if r['category'] in ['PPDK','PPDK_RP']],fields)
table('12_input_provenance.tsv',[dict(role=k,path=str(p.relative_to(root)),sha256=sha(p)) for k,p in [('corrected_proteome',proteome),('previous_hits',old/'08_family_subunit_candidates.tsv'),('reference_curation',old/'01_reference_curation.tsv')]],['role','path','sha256'])
(out/'COMMANDS.json').write_text(json.dumps(commands,indent=2)+'\n')
(out/'README.txt').write_text("""Corrected-database reference screening and limited-panel specificity evidence.
Verifies corrected FASTA bytes and each retained sequence against correction ledger.
Reconciles previous step17.18 candidates by exact ID and sequence. No suffix aliases.
Recomputes target search statistics against the corrected database with original settings.
06 records shared, added, lost and removed-CDS group hits.
09 compares candidate-query hits to the same group and competing reference groups.
A score advantage in this 114-reference panel is NOT validated family membership.
The panel has sparse kinase competitors: PPCK specificity remains unresolved.
V-ATPase related subunit hits require domain and subunit-specific review.
Alignment span is not domain coverage. No arbitrary margin/coverage threshold accepts genes.
No-hit does not prove absence. All counts are candidate IDs, not unique validated loci.
Protein sequences and IDs are unchanged; prior alignments can be reused only after
checking their sequence mapping. Prior trees, HMMs and other families are not rerun.
BETA-CA/PEPC validation remains separate; this refresh covers the five step17.18 categories.
AT composition flags remain visible. The orphan SS CDS is not added.
Next use 11 for focused PPDK/RP completion and 09 to select domain/specificity checks.
No domains are asserted or automatically assigned by this script.
""")
qc=dict(status='CORRECTED_SCREEN_COMPLETE_REVIEW_ONLY',protein_records=len(seqs),candidates=len(candidateids),
 previous_reconciliation=dict(Counter(r['status'] for r in reconcile)),
 hit_changes=dict(Counter(r['status'] for r in delta)),final_gene_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Corrected_screen_validation/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
