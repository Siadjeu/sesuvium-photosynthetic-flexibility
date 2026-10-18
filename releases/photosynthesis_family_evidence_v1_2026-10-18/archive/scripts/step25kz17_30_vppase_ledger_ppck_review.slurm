#!/usr/bin/env bash
#SBATCH --job-name=ppck_review
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=ppck_review_%j.out
#SBATCH --error=ppck_review_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
command -v python3
command -v hmmsearch
command -v mafft
python3 - <<'PY'
import csv,gzip,hashlib,json,os,re,shutil,subprocess,tarfile
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
choices=[root/os.environ['SCREEN_RELEASE']] if os.environ.get('SCREEN_RELEASE') else sorted(p for p in (base/'step25kz17.22').glob('release_*') if p.is_dir() and (p/'QC.json').is_file())
if len(choices)!=1:raise SystemExit('Set SCREEN_RELEASE to the intended step17.22 release')
src=choices[0]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta(p):
 d={};h=None;s=[]
 with p.open() as f:
  for l in f:
   if l.startswith('>'):
    if h is not None:
     if h in d:raise ValueError('Duplicate FASTA ID '+h)
     d[h]=''.join(s)
    h=l[1:].split()[0];s=[]
   elif h is not None:s.append(l.strip())
 if h is not None:
  if h in d:raise ValueError('Duplicate FASTA ID '+h)
  d[h]=''.join(s)
 return d
for l in (src/'SHA256SUMS.txt').read_text().splitlines():
 h,n=l.split(maxsplit=1)
 if sha(src/n)!=h:raise SystemExit('Screen checksum mismatch: '+n)
# Use installed Pfam profiles; never silently select between database versions.
pfam=os.environ.get('PFAM_HMM')
if pfam:
 db=root/pfam
else:
 paths=set()
 homes=[root/'databases',root/'00_databases',root.parent/'my_interproscan']
 if os.environ.get('INTERPROSCAN_HOME'):homes.append(Path(os.environ['INTERPROSCAN_HOME']))
 exe=shutil.which('interproscan.sh')
 if exe:homes.append(Path(exe).resolve().parent)
 names={'pfam-a.hmm','pfam_a.hmm','pfam-a.hmm.gz','pfam_a.hmm.gz'}
 for home in homes:
  if home.is_dir():
   for folder,dirs,files in os.walk(home):
    dirs[:]=[d for d in dirs if not d.startswith('.') and d not in ['temp','tmp','work']]
    for name in files:
     if name.lower() in names:paths.add((Path(folder)/name).resolve())
 if len(paths)!=1:raise SystemExit('Set PFAM_HMM to an installed Pfam-A ASCII .hmm or .hmm.gz file. Found '+str(len(paths))+' candidates: '+'; '.join(map(str,sorted(paths))))
 db=next(iter(paths))
if not db.is_file():raise SystemExit('PFAM_HMM does not exist')
out=base/'step25kz17.30'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rr)
def writefa(p,d):
 with p.open('w') as f:
  for h,s in d.items():f.write('>'+h+'\n'+s+'\n')
# Record the agreed bounded V-PPase review without declaring unresolved models repaired.
def choose(step,env):
 pp=[root/os.environ[env]] if os.environ.get(env) else [p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
 if len(pp)!=1:raise SystemExit('Set '+env+' to the intended '+step+' directory')
 p=pp[0].resolve(strict=True)
 for line in (p/'SHA256SUMS.txt').read_text().splitlines():
  h,n=line.split(maxsplit=1);f=(p/n).resolve()
  if not f.is_relative_to(p) or sha(f)!=h:raise SystemExit('Checksum failure '+n)
 return p
bounded=choose('step25kz17.29','BOUNDED_RELEASE')
if json.loads((bounded/'QC.json').read_text()).get('consistency_issues'):raise SystemExit('Unexpected consistency issues in bounded review')
vrows=read(bounded/'08_candidate_ledger_unchanged.tsv')
if len(vrows)!=48:raise SystemExit('Expected 48 V-PPase ledger rows')
for r in vrows:
 flag=r['review_flag'];pid=r['protein_id']
 r['decision']='SUPPORTED_NO_PRIORITY_MODEL_FLAG' if flag=='NO_PRIORITY_EXCEPTION' else ('UNRESOLVED_FRAGMENT' if flag in ['SHORT_FRAGMENT','PARTIAL_N_TERMINAL','PARTIAL_C_TERMINAL'] else 'FAMILY_SUPPORTED_MODEL_CAVEAT')
 r['suspected_region']='PO_LG07_1805_1806' if pid in ['PO__evm.TU.LG07.1805','PO__evm.TU.LG07.1806'] else ''
 r['count_status']='EVIDENCE_CATEGORY_NOT_FINAL';r['read_validation']='NOT_AVAILABLE_FOR_MODEL_RESOLUTION'
if Counter(r['decision'] for r in vrows)!=Counter({'SUPPORTED_NO_PRIORITY_MODEL_FLAG':40,'FAMILY_SUPPORTED_MODEL_CAVEAT':4,'UNRESOLVED_FRAGMENT':4}):raise SystemExit('Unexpected V-PPase category totals')
table('11_VPPASE_decision_ledger.tsv',vrows,list(vrows[0]))
vc=[]
for sp in sorted({r['species'] for r in vrows}):
 rr=[r for r in vrows if r['species']==sp]
 if any(r['mapping_status']!='UNAMBIGUOUS' for r in rr):raise SystemExit('Unresolved V-PPase locus mapping')
 if len({r['gene_id'] for r in rr})!=len(rr):raise SystemExit('Duplicate mapped V-PPase genes')
 vc.append(dict(species=sp,annotated_candidate_loci=len(rr),supported_no_priority_flag=sum(r['decision']=='SUPPORTED_NO_PRIORITY_MODEL_FLAG' for r in rr),family_supported_model_caveat=sum(r['decision']=='FAMILY_SUPPORTED_MODEL_CAVEAT' for r in rr),unresolved_fragment_records=sum(r['decision']=='UNRESOLVED_FRAGMENT' for r in rr),suspected_split_regions=len({r['suspected_region'] for r in rr if r['suspected_region']}),status='NOT_FINAL_GENE_COUNTS'))
table('12_VPPASE_count_categories.tsv',vc,list(vc[0]))
(out/'13_VPPASE_limitations.txt').write_text('V-PPase screening identified 48 annotated candidate loci. Forty had no priority model flag; four additional family-supported loci retained model or domain-architecture caveats. Four fragment records remained unresolved, including two adjacent Portulaca annotations bridged by protein-to-genome alignments and treated as a suspected split-gene region without automatic merging. Suitable independent read evidence was unavailable to resolve these models. No confirmed pseudogene, corrected gene model, or final functional gene count was inferred.\n')
(out/'LEDGER_PROVENANCE.json').write_text(json.dumps({'bounded_release':os.path.relpath(bounded,root),'bounded_manifest_sha256':sha(bounded/'SHA256SUMS.txt'),'screen_manifest_sha256':sha(src/'SHA256SUMS.txt')},indent=2)+'\n')

if db.suffix=='.gz':
 expanded=out/'work/Pfam-A.hmm'
 with gzip.open(db,'rb') as a,expanded.open('wb') as b:shutil.copyfileobj(a,b)
 dbuse=expanded
else:dbuse=db
# Validate database identity and threshold annotations before spending compute.
models=ga=pf=0
with dbuse.open() as f:
 for line in f:
  if line.startswith('HMMER3/'):models+=1
  elif line.startswith('GA '):ga+=1
  elif line.startswith('ACC ') and 'PF' in line:pf+=1
if models<1000 or ga!=models or pf!=models:raise SystemExit('Expected a complete Pfam-A HMMER3 ASCII database with GA cutoffs; inspect PFAM_HMM')
cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
allrows=[r for r in read(src/'09_specificity_evidence_REVIEW_ONLY.tsv') if r['category']=='PPCK']
if len(allrows)!=5121 or len({r['protein_id'] for r in allrows})!=5121:raise SystemExit('Expected 5121 unique PPCK screen hits')
# Bounded prioritization, not an acceptance rule. Retain the full screen ledger.
# Two complementary rankings prevent compact models being buried among long kinases.
selected=set();ranks={}
for sp in sorted({r['species'] for r in allrows}):
 rr=[r for r in allrows if r['species']==sp]
 byscore=sorted(rr,key=lambda r:(-float(r['bitscore']),r['protein_id']))
 compact=sorted(rr,key=lambda r:(abs(int(r['length_aa'])-281),-float(r['bitscore']),r['protein_id']))
 for i,r in enumerate(byscore,1):ranks.setdefault(r['protein_id'],{})['score_rank']=i
 for i,r in enumerate(compact,1):ranks.setdefault(r['protein_id'],{})['length_rank']=i
 selected.update(r['protein_id'] for r in byscore[:20]);selected.update(r['protein_id'] for r in compact[:20])
rows=[r for r in allrows if r['protein_id'] in selected]
queue=[dict(r,**ranks[r['protein_id']],domain_batch='SELECTED' if r['protein_id'] in selected else 'DEFERRED_NOT_EXCLUDED') for r in allrows]
table('09_PPCK_all_screen_hits_REVIEW_ONLY.tsv',queue,list(queue[0]))
table('10_PPCK_priority_batch.tsv',rows,list(rows[0]))
seqs=fasta(src/'05_candidate_sequences.faa');refs=fasta(src/'references/targets.faa')
refgroups={'PPCK':['Q93VK0','Q9SPK4']}
query={};mapping=[];groups=defaultdict(dict)
for i,r in enumerate(sorted(rows,key=lambda x:(x['category'],x['protein_id'])),1):
 q='C'+str(i).zfill(3);s=seqs[r['protein_id']]
 query[q]=s;groups[r['category']][q]=s
 mapping.append(dict(query_id=q,protein_id=r['protein_id'],family=r['category'],role='CANDIDATE',length_aa=len(s)))
for fam,ids in refgroups.items():
 for acc in ids:
  q='REF_'+acc;query[q]=refs[acc];groups[fam][q]=refs[acc]
  mapping.append(dict(query_id=q,protein_id=acc,family=fam,role='REFERENCE',length_aa=len(refs[acc])))
table('01_sequence_map.tsv',mapping,['query_id','protein_id','family','role','length_aa'])
writefa(out/'02_candidates_and_references.faa',query)
def run(args,log,output=None):
 with (out/log).open('w') as f:
  if output:
   with output.open('w') as g:subprocess.run(args,stdout=g,stderr=f,check=True)
  else:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
run(['hmmsearch','-h'],'hmmsearch_version.txt')
run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'03_pfam.domtbl'),str(dbuse),str(out/'02_candidates_and_references.faa')],'03_hmmsearch.log')
domains=[];byq=defaultdict(list)
for l in (out/'03_pfam.domtbl').read_text().splitlines():
 if not l or l.startswith('#'):continue
 x=l.split(maxsplit=22)
 r=dict(query_id=x[0],model=x[3],accession=x[4],model_length=int(x[5]),full_evalue=x[6],
  independent_evalue=x[12],domain_bitscore=x[13],hmm_start=int(x[15]),hmm_end=int(x[16]),
  protein_start=int(x[17]),protein_end=int(x[18]),envelope_start=int(x[19]),envelope_end=int(x[20]),
  model_span=round((int(x[16])-int(x[15])+1)/int(x[5]),4))
 domains.append(r);byq[x[0]].append(r)
fields=['query_id','model','accession','model_length','full_evalue','independent_evalue','domain_bitscore','hmm_start','hmm_end','protein_start','protein_end','envelope_start','envelope_end','model_span']
table('04_domain_evidence.tsv',domains,fields)
table('05_reference_domain_baselines.tsv',[r for r in domains if r['query_id'].startswith('REF_')],fields)
# Preserve every overlapping hit: this is raw Pfam evidence, not clan-resolved annotation.
comparisons=[];completeness=[]
for fam,records in groups.items():
 inp=out/(fam+'_input.faa');writefa(inp,records)
 aln=out/(fam+'_alignment.faa')
 run(['mafft','--auto','--thread',cpu,str(inp)],fam+'_mafft.log',aln)
 aligned=fasta(aln)
 if set(aligned)!=set(records) or len(set(map(len,aligned.values())))!=1:raise SystemExit('Invalid MAFFT alignment')
 for q,s in aligned.items():
  if s.replace('-','').upper()!=records[q].upper():raise SystemExit('Alignment sequence changed')
 for acc in refgroups[fam]:
  ref='REF_'+acc;reference=aligned[ref];pos=[i for i,a in enumerate(reference) if a!='-']
  rd=sorted(byq[ref],key=lambda x:x['protein_start'])
  signature=[r['accession'] for r in rd]
  for q in records:
   if q.startswith('REF_'):continue
   cd=sorted(byq[q],key=lambda x:x['protein_start']);candidate_signature=[r['accession'] for r in cd]
   comparisons.append(dict(query_id=q,reference=acc,family=fam,reference_domain_order=';'.join(signature),
    candidate_domain_order=';'.join(candidate_signature),
    missing_reference_models=';'.join(sorted(set(signature)-set(candidate_signature))),
    additional_models=';'.join(sorted(set(candidate_signature)-set(signature))),
    exact_raw_domain_order_match='YES' if signature and candidate_signature==signature else 'NO',
    status='REFERENCE_NO_GA_DOMAINS' if not signature else 'REVIEW_NOT_FINAL'))
   occupied=[j for j,i in enumerate(pos) if aligned[q][i]!='-']
   longest=gap=0
   for i in pos:
    gap=gap+1 if aligned[q][i]=='-' else 0;longest=max(longest,gap)
   completeness.append(dict(query_id=q,family=fam,reference=acc,
    reference_length=len(pos),candidate_length=len(records[q]),
    reference_positions_occupied=len(occupied),reference_position_coverage=round(len(occupied)/len(pos),4),
    missing_N_reference_positions=occupied[0] if occupied else len(pos),
    missing_C_reference_positions=len(pos)-1-occupied[-1] if occupied else len(pos),
    longest_reference_gap=longest,review_status='PENDING'))
table('06_domain_architecture_comparison.tsv',comparisons,['query_id','reference','family','reference_domain_order','candidate_domain_order','missing_reference_models','additional_models','exact_raw_domain_order_match','status'])
table('07_alignment_completeness.tsv',completeness,['query_id','family','reference','reference_length','candidate_length','reference_positions_occupied','reference_position_coverage','missing_N_reference_positions','missing_C_reference_positions','longest_reference_gap','review_status'])
table('08_manual_decision_template.tsv',[dict(**r,domain_decision='PENDING',completeness_decision='PENDING',locus_mapping='PENDING',final_decision='PENDING',reason='') for r in mapping if r['role']=='CANDIDATE'],['query_id','protein_id','family','length_aa','domain_decision','completeness_decision','locus_mapping','final_decision','reason'])
(out/'DATABASE_PROVENANCE.json').write_text(json.dumps(dict(database_filename=db.name,compressed_or_source_sha256=sha(db),uncompressed_sha256=sha(dbuse),profiles=models,method='hmmsearch --cut_ga',database_release='Confirm from local installation metadata; hash pins exact bytes'),indent=2)+'\n')
(out/'README.txt').write_text("""V-PPase decision consolidation and bounded PPCK specificity review.
V-PPase: 40 unflagged loci, four family-supported loci with model caveats,
and four unresolved fragment records. The PO fragment pair shares a suspected
split-region identifier; original loci are preserved and no merging is performed.
Counts are evidence categories, not final functional or Wang-comparable totals.
PPCK: all 5121 screen hits retained. Per species, select the union of the top
20 DIAMOND bit-score hits and 20 proteins closest in length to 281 aa (midpoint
of the 278/284-aa reference lengths). This is a practical batch-selection rule,
not a biological threshold. Deferred sequences are not excluded.
Pfam domains and reference alignments help identify incompatible architectures
and incompleteness, but shared kinase domains do not establish PPCK specificity.
The existing limited competitor-panel evidence is retained in the queue.
No PPCK final membership/count is assigned. Expanded kinase comparison and
phylogenetic review may still be required after this batch.
No sequences or existing releases are changed. No additional V-PPase rescue.
""")
(out/'QC.json').write_text(json.dumps(dict(status='DOMAIN_REVIEW_COMPLETE_NOT_FINAL',candidates=len(rows),references=2,
 pfam_models=models,domain_hits=len(domains),reference_controls_without_hits=[q for q in query if q.startswith('REF_') and not byq[q]],final_counts=False),indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='VPPASE_ledger_PPCK_review/'+str(p.relative_to(out)))
print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
