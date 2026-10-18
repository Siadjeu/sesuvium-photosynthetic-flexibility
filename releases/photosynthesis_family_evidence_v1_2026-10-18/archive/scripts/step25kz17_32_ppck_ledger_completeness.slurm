#!/usr/bin/env bash
#SBATCH --job-name=ppck_ledger
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=ppck_ledger_%j.out
#SBATCH --error=ppck_ledger_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
for tool in python3 diamond hmmbuild hmmsearch; do command -v "$tool"; done
python3 - <<'PY'
import csv,hashlib,json,os,re,subprocess,tarfile,shutil
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
tree=release('step25kz17.31','TREE_RELEASE');screen=release('step25kz17.22','SCREEN_RELEASE')
prov=json.loads((tree/'PROVENANCE.json').read_text())
if sha(screen/'SHA256SUMS.txt')!=prov['screen_manifest_sha256']:raise SystemExit('Tree and screen provenance mismatch')
tips=read(tree/'02_tip_map.tsv');high={r['tip'] for r in tips if r['role']=='HIGH_PRIORITY'}
refs={r['tip'] for r in tips if r['role']=='PPCK_REFERENCE'}
assert len(high)==27 and refs=={'REF_Q93VK0','REF_Q9SPK4'}
# Verify the exact supported bipartition; never assume that HIGH_PRIORITY means accepted.
tokens=re.findall(r'[^(),:;\s]+|[(),:;]',(tree/'10_PPCK_kinase.treefile').read_text());i=0;nodes=[]
def parse():
 global i
 leaves=set();label=''
 if tokens[i]=='(':
  i+=1
  while True:
   leaves|=parse()
   if tokens[i]==',':i+=1;continue
   assert tokens[i]==')';i+=1;break
  if tokens[i] not in ':,);':label=tokens[i];i+=1
 else:leaves={tokens[i]};i+=1
 if tokens[i]==':':i+=2
 nodes.append((leaves,label));return leaves
alltips=parse();assert alltips=={r['tip'] for r in tips}
clade=high|refs
support=next((s for leaves,s in nodes if leaves==clade or alltips-leaves==clade),None)
assert support is not None,'Supported PPCK split absent'
a,b=map(float,support.split('/'));assert a>=99 and b>=99
loc={r['tip']:r for r in read(tree/'08_locus_reconciliation.tsv')}
assert all(loc[t]['mapping_status']=='UNAMBIGUOUS' for t in high)
assert len({(loc[t]['species'],loc[t]['gene_id']) for t in high})==27
rows=[r for r in read(screen/'09_specificity_evidence_REVIEW_ONLY.tsv') if r['category']=='PPCK']
assert len(rows)==5121 and len({r['protein_id'] for r in rows})==5121
seq=fasta(screen/'05_candidate_sequences.faa');queries={r['protein_id']:seq[r['protein_id']] for r in rows}
full=fasta(tree/'01_full_proteins.faa');aln=fasta(tree/'07_reference_core.faa')
for r in tips:
 if r['tip'] not in refs and not r['tip'].startswith('REF_'):
  assert queries[r['protein_id']]==full[r['tip']],'Changed sequence'
accepted={r['protein_id'] for r in tips if r['tip'] in high}
out=base/'step25kz17.32'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False)
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
commands=[];cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
def run(args,log):
 commands.append(args)
 with (out/log).open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
ledger=[]
for r in tips:
 if r['tip'] in high:
  ledger.append(dict(protein_id=r['protein_id'],species=r['species'],gene_id=loc[r['tip']]['gene_id'],tree_tip=r['tip'],decision='PHYLOGENETICALLY_SUPPORTED_PPCK_LOCUS',split_support=support,subtype='UNRESOLVED',functional_validation='NOT_ESTABLISHED',exhaustive_count='NOT_ESTABLISHED'))
table('01_PPCK_evidence_ledger.tsv',ledger,list(ledger[0]))
counts=Counter(r['species'] for r in ledger)
table('02_PPCK_supported_counts.tsv',[dict(species=s,supported_loci=counts[s],exhaustive_count='NOT_ESTABLISHED') for s in sorted(counts)],['species','supported_loci','exhaustive_count'])
writefa(out/'03_PPCK_training_alignment.faa',{t:aln[t] for t in sorted(clade)})
writefa(out/'04_all_screened_PPCK_candidates.faa',queries)
# Use the same aligned kinase core for every comparison panel sequence.
writefa(out/'05_comparison_panel.faa',{t:s.replace('-','').replace('.','') for t,s in aln.items()})
shutil.copy2(tree/'02_tip_map.tsv',out/'05_comparison_panel_map.tsv')
run(['hmmbuild','--amino','--cpu',cpu,str(out/'PPCK.hmm'),str(out/'03_PPCK_training_alignment.faa')],'06_hmmbuild.log')
# --max removes acceleration filters for this small, bounded 5,121-protein search.
run(['hmmsearch','--cpu',cpu,'--max','-E','10','--domE','10','--tblout',str(out/'07_profile_hits.tbl'),'--domtblout',str(out/'07_profile_domains.tbl'),str(out/'PPCK.hmm'),str(out/'04_all_screened_PPCK_candidates.faa')],'07_hmmsearch.log')
run(['diamond','makedb','--threads',cpu,'--in',str(out/'05_comparison_panel.faa'),'--db',str(out/'panel')],'08_makedb.log')
run(['diamond','blastp','--threads',cpu,'--very-sensitive','--query',str(out/'04_all_screened_PPCK_candidates.faa'),'--db',str(out/'panel'),'--evalue','10','--max-target-seqs','0','--outfmt','6','qseqid','sseqid','pident','length','evalue','bitscore','qlen','slen','--out',str(out/'08_panel_hits.tsv')],'08_diamond.log')
profile={}
for l in (out/'07_profile_hits.tbl').read_text().splitlines():
 if l and not l.startswith('#'):
  x=l.split();profile[x[0]]=(float(x[5]),float(x[4]))
assert accepted<=profile.keys(),'Supported locus missing from profile output'
# Empirical in-sample reference band: discovery heuristic, NOT a validated cutoff.
floor=min(profile[p][0] for p in accepted)-20.0
best=defaultdict(dict);idmap={r['tip']:r['protein_id'] for r in tips}
for l in (out/'08_panel_hits.tsv').read_text().splitlines():
 x=l.split('\t');q,t=x[:2]
 if idmap[t]==q:continue # exclude the exact self match during score calibration
 group='PPCK' if t in clade else 'COMPARISON'
 hit=(float(x[5]),t,float(x[4]),float(x[2]),int(x[3]))
 if group not in best[q] or hit[0]>best[q][group][0]:best[q][group]=hit
results=[];queue=[]
for r in rows:
 p=r['protein_id'];ph=profile.get(p);pos=best[p].get('PPCK');neg=best[p].get('COMPARISON');reasons=[]
 if ph and ph[0]>=floor:reasons.append('PROFILE_WITHIN_SUPPORTED_BAND_MINUS_20_BITS')
 if pos and pos[2]<=1e-5 and (not neg or pos[0]>=neg[0]-20):reasons.append('PPCK_PANEL_HIT_WITHIN_20_BITS_OF_COMPARISON')
 if ph and ph[1]<=1e-5 and pos and pos[2]<=1e-5 and not neg:reasons.append('PROFILE_AND_PPCK_HIT_NO_COMPARISON_HIT')
 status='ALREADY_SUPPORTED' if p in accepted else ('REVIEW_TRIGGERED' if reasons else 'NO_TRIGGER_IN_BOUNDED_CHECK')
 row=dict(protein_id=p,species=r['species'],length_aa=r['length_aa'],status=status,profile_bitscore=ph[0] if ph else '',profile_evalue=ph[1] if ph else '',best_PPCK_tip=pos[1] if pos else '',PPCK_bitscore=pos[0] if pos else '',best_comparison_tip=neg[1] if neg else '',comparison_bitscore=neg[0] if neg else '',score_margin=pos[0]-neg[0] if pos and neg else '',triggers=';'.join(reasons),previous_screen_status=r['limited_panel_status'])
 results.append(row)
 if status=='REVIEW_TRIGGERED':queue.append(row)
queue.sort(key=lambda r:(-float(r['profile_bitscore']) if r['profile_bitscore']!='' else float('inf'),r['protein_id']))
table('09_all_candidates_completeness_evidence.tsv',results,list(results[0]));table('10_remaining_review_queue.tsv',queue,list(results[0]))
writefa(out/'11_remaining_review_queue.faa',{r['protein_id']:queries[r['protein_id']] for r in queue})
# Show empirical calibration on accepted loci and all 98 comparison candidates.
roles={r['protein_id']:r['role'] for r in tips if not r['tip'].startswith('REF_')}
cal=[dict(**r,tree_role=roles[r['protein_id']]) for r in results if r['protein_id'] in roles]
table('12_panel_calibration.tsv',cal,list(cal[0]))
(out/'QC.json').write_text(json.dumps(dict(status='BOUNDED_CHECK_COMPLETE_REVIEW_REQUIRED',supported_loci=27,screened_candidates=len(rows),remaining_candidates=len(rows)-27,review_queue=len(queue),profile_discovery_floor=floor,score_margin_allowance_bits=20,final_counts=False,whole_proteome_completeness_tested=False),indent=2)+'\n')
(out/'LIMITATIONS.txt').write_text('27 distinct loci supported by the reviewed PPCK split; no PPCK1/2 subtype or biochemical function assigned.\nTree has short-alignment parameter warning; internal reference attachment is weak.\nSearch universe is the 5121 previously screened candidates, not the whole proteome or genome.\nProfile is trained on accepted candidates and references; score calibration is in-sample and not independent validation.\nDiscovery rules use a 20-bit allowance below the weakest supported profile score OR a PPCK panel hit within 20 bits of the best comparison hit (E<=1e-5). These are permissive prioritization heuristics, not family acceptance/rejection thresholds.\nComparison panel consists of 98 non-PPCK-group tree candidates and four kinase controls, not an exhaustive kinase classification.\nNo-trigger candidates are not excluded; short/divergent/missing gene models may be missed. Queue entries require review and cannot increase counts automatically.\nDistinct loci are preserved even if protein sequences are identical. No original release is modified.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(tree_manifest_sha256=sha(tree/'SHA256SUMS.txt'),screen_manifest_sha256=sha(screen/'SHA256SUMS.txt'),verified_split_support=support),indent=2)+'\n')
for tool,args in [('diamond',['version']),('hmmbuild',['-h']),('hmmsearch',['-h'])]:run([tool]+args,tool+'_version.txt')
(out/'COMMANDS.json').write_text(json.dumps([[v.replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for v in c] for c in commands],indent=2)+'\n')
(out/'panel.dmnd').unlink(missing_ok=True)
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.iterdir()):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+p.name+'\n')
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='PPCK_ledger_completeness')
print('RELEASE='+str(out));print('BUNDLE='+str(bundle));print('SUPPORTED_LOCI=27; REMAINING_REVIEW_QUEUE='+str(len(queue)))
PY
