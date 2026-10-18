#!/usr/bin/env bash
#SBATCH --job-name=ppdk_review
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=ppdk_review_%j.out
#SBATCH --error=ppdk_review_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from project directory}"
command -v python3
command -v diamond
command -v mafft
python3 - <<'PY'
import csv,hashlib,json,os,re,subprocess,tarfile
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
choices=[root/os.environ['SCREEN_RELEASE']] if os.environ.get('SCREEN_RELEASE') else sorted(p for p in (base/'step25kz17.18').glob('release_*') if p.is_dir() and (p/'QC.json').is_file())
if len(choices)!=1:raise SystemExit('Set SCREEN_RELEASE to the intended step17.18 release; found '+str(len(choices)))
src=choices[0]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for line in (src/'SHA256SUMS.txt').read_text().splitlines():
 h,n=line.split(maxsplit=1);p=src/n
 if not p.is_file() or sha(p)!=h:raise SystemExit('Source checksum mismatch: '+n)
out=base/'step25kz17.19'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
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
def writefa(p,records):
 with p.open('w') as f:
  for h,s in records:f.write('>'+h+'\n'+s+'\n')
proteome=root/os.environ.get('PROTEOME','08_focused_families/validation/step25kz16.7/release/work/all_species_primary_proteins.faa')
seqs={}
for h,s in fasta(proteome):
 pid=h.split()[0]
 if pid in seqs:raise SystemExit('Duplicate primary ID: '+pid)
 seqs[pid]=s
# Alias candidates are explicitly reported, never accepted as proven locus mappings.
aliases=defaultdict(list)
for pid in seqs:aliases[re.sub(r'\.t\d+$','',pid)].append(pid)
current=[r for r in read(src/'08_family_subunit_candidates.tsv') if r['category'] in ['PPDK','PPDK_RP']]
selected={r['protein_id'] for r in current};screenfam={r['protein_id']:r['category'] for r in current}
if not selected<=seqs.keys():raise SystemExit('Current screen IDs absent from selected proteome')
oldfile=base/'step25kn/results/PPDK/candidate_validation.tsv'
old=read(oldfile);mapping=[]
for r in old:
 raw=r['sequence_id'];pid=raw.split('|')[-1]
 matches=[pid] if pid in seqs else aliases.get(pid,[])
 status='EXACT_ID' if pid in seqs else ('UNIQUE_TRANSCRIPT_LABEL_CANDIDATE' if len(matches)==1 else 'UNRESOLVED')
 if len(matches)==1:selected.add(matches[0])
 mapping.append(dict(old_id=raw,old_status=r['status'],old_length=r['length_or_alignment'],
  current_id=matches[0] if len(matches)==1 else '',mapping_status=status,
  sequence_matches_old_accepted='NOT_TESTED'))
oldaccepted=base/'step25kn/results/PPDK/PPDK.accepted.proteins.fa'
oldseq={h.split()[0].split('|')[-1]:s for h,s in fasta(oldaccepted)}
for r in mapping:
 pid=r['old_id'].split('|')[-1]
 if pid in oldseq and r['current_id']:
  r['sequence_matches_old_accepted']='YES' if oldseq[pid]==seqs[r['current_id']] else 'NO'
table('01_old_to_current_ID_audit.tsv',mapping,['old_id','old_status','old_length','current_id','mapping_status','sequence_matches_old_accepted'])
# Preserve missing legacy sequence as a separate query when the current primary sequence is absent/changed.
querymap=[];qseq={}
for pid in sorted(selected):
 q='Q'+str(len(querymap)+1).zfill(4);qseq[q]=seqs[pid]
 querymap.append(dict(query_id=q,protein_id=pid,sequence_source='CURRENT_PRIMARY',length_aa=len(seqs[pid]),screen_family=screenfam.get(pid,'NOT_RECOVERED_AS_PPDK_OR_RP')))
for pid,s in sorted(oldseq.items()):
 rr=[r for r in mapping if r['old_id'].split('|')[-1]==pid]
 if not rr or not rr[0]['current_id'] or rr[0]['sequence_matches_old_accepted']=='NO':
  q='Q'+str(len(querymap)+1).zfill(4);qseq[q]=s
  querymap.append(dict(query_id=q,protein_id=pid,sequence_source='LEGACY_ACCEPTED',length_aa=len(s),screen_family='LEGACY_PPDK'))
table('02_query_ID_map.tsv',querymap,['query_id','protein_id','sequence_source','length_aa','screen_family'])
writefa(out/'03_review_queries.faa',sorted(qseq.items()))
refmeta={}
for r in read(src/'01_reference_curation.tsv'):refmeta[r['accession']]=r
refs={}
for name in ['curated_targets.faa','BACKGROUND_REFERENCE.faa']:
 for h,s in fasta(src/'references'/name):
  acc=h.split()[0]
  if acc in refs and refs[acc]!=s:raise SystemExit('Reference sequence conflict')
  refs[acc]=s
writefa(out/'04_competition_references.faa',sorted(refs.items()))
table('04_reference_metadata.tsv',list(refmeta.values()),['accession','assigned_category','subunit','header'])
cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
def run(args,log,stdout=None):
 with (out/log).open('w') as f:
  if stdout:
   with stdout.open('w') as g:subprocess.run(args,stdout=g,stderr=f,check=True)
  else:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
db=out/'work/refs'
run(['diamond','makedb','--in',str(out/'04_competition_references.faa'),'--db',str(db),'--threads',cpu],'05_makedb.log')
hits=out/'06_competitive_hits.tsv'
run(['diamond','blastp','--query',str(out/'03_review_queries.faa'),'--db',str(db),'--out',str(hits),'--threads',cpu,'--sensitive','--evalue','1e-5','--max-target-seqs','0','--outfmt','6','qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen'],'06_diamond.log')
fields=['query','ref','identity','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen']
byq=defaultdict(list)
with hits.open() as f:
 for r in csv.DictReader(f,fieldnames=fields,delimiter='\t'):byq[r['query']].append(r)
summary=[];align={'PPDK':[],'PPDK_RP':[]}
for m in querymap:
 q=m['query_id'];hh=sorted(byq[q],key=lambda r:-float(r['bitscore']))
 target=[r for r in hh if refmeta[r['ref']]['assigned_category'] in ['PPDK','PPDK_RP']]
 other=[r for r in hh if refmeta[r['ref']]['assigned_category'] not in ['PPDK','PPDK_RP']]
 best=target[0] if target else None
 row=dict(m,best_target_reference=best['ref'] if best else '',best_target_family=refmeta[best['ref']]['assigned_category'] if best else '',
  target_bitscore=best['bitscore'] if best else '',best_other_reference=other[0]['ref'] if other else '',other_bitscore=other[0]['bitscore'] if other else '',
  target_query_span='',target_reference_span='',review_status='NO_PPDK_RP_HIT' if not best else 'REVIEW_REQUIRED')
 if best:
  row['target_query_span']=round((abs(int(best['qend'])-int(best['qstart']))+1)/int(best['qlen']),4)
  row['target_reference_span']=round((abs(int(best['send'])-int(best['sstart']))+1)/int(best['slen']),4)
  row['review_status']='OTHER_REFERENCE_SCORES_AT_LEAST_AS_HIGH' if other and float(other[0]['bitscore'])>=float(best['bitscore']) else 'TARGET_FAVORED_IN_LIMITED_PANEL'
  align[row['best_target_family']].append((q,qseq[q]))
 summary.append(row)
table('07_specificity_review.tsv',summary,list(summary[0]))
for fam,rr in align.items():
 rr=[('REF_'+a,s) for a,s in refs.items() if refmeta[a]['assigned_category']==fam]+rr
 writefa(out/(fam+'_alignment_input.faa'),rr)
 run(['mafft','--auto','--thread',cpu,str(out/(fam+'_alignment_input.faa'))],fam+'_mafft.log',out/(fam+'_alignment.faa'))
table('08_source_provenance.tsv',[dict(source=str(p.relative_to(root)),sha256=sha(p)) for p in [oldfile,oldaccepted,src/'08_family_subunit_candidates.tsv',src/'01_reference_curation.tsv']],['source','sha256'])
(out/'README.txt').write_text("""PPDK/PPDK-RP focused specificity and identifier review.
All legacy PPDK validation rows are audited, including failed records.
Exact identifiers and unique transcript-label candidates are distinguished.
A unique suffix alias alone is not proof of locus identity; old accepted protein
sequence identity is recorded where available. Original IDs remain in the map.
Legacy accepted sequences missing/changed in primary data are separately queried.
The 114-reference panel contains current target and background references.
Target-favored means only within this limited panel, not definitive orthology.
No domain architecture validation or final gene counting is performed here.
Separate MAFFT alignments use provisional best-target assignment for review only;
they are not phylogenetic membership proof. Q IDs map through 02_query_ID_map.tsv.
Inspect the missing SP_Hex candidate and partial RP candidates first, then resolve
locus mapping and domain completeness. PPCK/ATPase specificity is a separate task.
""")
(out/'QC.json').write_text(json.dumps(dict(status='REVIEW_COMPLETE_NOT_FINAL',queries=len(querymap),reference_count=len(refs),old_rows=len(mapping),unresolved_old_ids=sum(r['mapping_status']=='UNRESOLVED' for r in mapping)),indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='PPDK_specificity_review/'+str(p.relative_to(out)))
print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
