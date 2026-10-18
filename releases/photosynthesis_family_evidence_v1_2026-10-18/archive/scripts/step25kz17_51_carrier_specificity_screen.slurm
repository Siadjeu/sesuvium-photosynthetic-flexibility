#!/usr/bin/env bash
#SBATCH --job-name=carrier_spec
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=carrier_spec_%j.out
#SBATCH --error=carrier_spec_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v hmmsearch
command -v diamond
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

src=release('step25kz17.50','CARRIER_PANEL_RELEASE');need(json.loads((src/'QC.json').read_text())['carrier_panel_status']=='BROAD_PANEL_READY_FOR_DEFINITION_REVIEW','Expanded panel is not ready')
queue=read(src/'08_mitochondrial_carrier_queue_REVIEW_ONLY.tsv');proteins,_=fasta(src/'09_mitochondrial_carrier_candidates.faa');refs,_=fasta(src/'references/04_carrier_panel.faa');meta=read(src/'references/03_broad_panel_annotations_REVIEW_REQUIRED.tsv');refmeta={r['reference_id']:r for r in meta}
need(len(queue)==len(proteins)==724 and len(refs)==43,'Unexpected screen size');need(set(proteins)=={r['protein_id'] for r in queue} and len({(r['species'],r['gene_id']) for r in queue})==724,'Candidate/locus mismatch');need(set(refs)==set(refmeta),'Panel metadata mismatch')
for k,v in refs.items():need(hashlib.sha256(v.encode()).hexdigest()==refmeta[k]['sequence_sha256'],'Reference hash mismatch')
targets={'REF__Q9C5M0':'DTC','REF__Q9SJY5':'DIC','REF__Q9SB52':'DIC','REF__Q9FY68':'DIC'}
for k,v in targets.items():need(refmeta[k]['assigned_family']==v,'Target definition mismatch')
need(all(valid(v) for v in list(proteins.values())+list(refs.values())),'Invalid sequence')
out=base/'step25kz17.51'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
writefa('01_candidates.faa',proteins);writefa('02_reference_panel.faa',refs);writefa('03_domain_queries.faa',{**proteins,**refs});table('02_reference_annotations.tsv',meta,list(meta[0]))
if os.environ.get('PFAM_HMM'):db=(root/os.environ['PFAM_HMM']).resolve(strict=True)
else:
 paths=set();homes=[root/'databases',root/'00_databases',root.parent/'my_interproscan']
 if os.environ.get('INTERPROSCAN_HOME'):homes.append(Path(os.environ['INTERPROSCAN_HOME']))
 exe=shutil.which('interproscan.sh')
 if exe:homes.append(Path(exe).resolve().parent)
 for home in homes:
  if home.is_dir():
   for folder,dirs,files in os.walk(home):
    dirs[:]=[d for d in dirs if not d.startswith('.') and d not in ['tmp','temp','work']]
    for n in files:
     if n.lower() in {'pfam-a.hmm','pfam_a.hmm','pfam-a.hmm.gz','pfam_a.hmm.gz'}:paths.add((Path(folder)/n).resolve())
 need(len(paths)==1,'Set PFAM_HMM to the intended installed Pfam-A .hmm or .hmm.gz');db=next(iter(paths))
dbuse=db
if db.suffix=='.gz':
 dbuse=out/'work/Pfam-A.hmm'
 with gzip.open(db,'rb') as a,dbuse.open('wb') as b:shutil.copyfileobj(a,b)
models=ga=0
with dbuse.open() as f:
 for line in f:
  if line.startswith('HMMER3/'):models+=1
  elif line.startswith('GA '):ga+=1
need(models>1000 and models==ga,'Expected complete Pfam with GA thresholds')
cpu=os.environ.get('SLURM_CPUS_PER_TASK','4');commands=[]
def run(args,log,stdout=None):
 commands.append([str(x).replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for x in args])
 with (out/log).open('w') as f:
  if stdout:
   with (out/stdout).open('w') as a:subprocess.run(args,stdout=a,stderr=f,check=True)
  else:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)

run(['diamond','makedb','--in',str(out/'02_reference_panel.faa'),'--db',str(out/'work/panel'),'--threads',cpu],'04_makedb.log')
fmt=['qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen']
run(['diamond','blastp','--query',str(out/'01_candidates.faa'),'--db',str(out/'work/panel'),'--out',str(out/'04_all_panel_hits.tsv'),'--threads',cpu,'--very-sensitive','--evalue','1e-3','--max-target-seqs','0','--max-hsps','1','--outfmt','6']+fmt,'04_diamond.log')
run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'05_pfam.domtbl'),str(dbuse),str(out/'03_domain_queries.faa')],'05_hmmsearch.log')
domains=defaultdict(list);dh=[]
for line in (out/'05_pfam.domtbl').read_text().splitlines():
 if not line or line.startswith('#'):continue
 x=line.split(maxsplit=22);need(x[0] in proteins or x[0] in refs,'Unknown domain sequence');r=dict(sequence_id=x[0],model=x[3],accession=x[4],model_length=int(x[5]),i_evalue=float(x[12]),score=float(x[13]),hmm_from=int(x[15]),hmm_to=int(x[16]),ali_from=int(x[17]),ali_to=int(x[18]));dh.append(r);domains[x[0]].append(r)
table('06_all_domain_hits.tsv',dh,['sequence_id','model','accession','model_length','i_evalue','score','hmm_from','hmm_to','ali_from','ali_to'])
dsummary=[];ds={}
for k,seq in {**proteins,**refs}.items():
 hh=domains[k];car=[h for h in hh if h['accession'].split('.')[0]=='PF00153'];covered=set()
 for h in car:covered.update(range(h['ali_from'],h['ali_to']+1))
 r=dict(sequence_id=k,length_aa=len(seq),PF00153_hit_count=len(car),PF00153_union_protein_coverage=round(len(covered)/len(seq),6),PF00153_intervals=';'.join(str(h['ali_from'])+'-'+str(h['ali_to']) for h in car),other_Pfam_models=';'.join(sorted({h['accession'] for h in hh if h not in car})),all_Pfam_models=';'.join(sorted({h['accession'] for h in hh})));dsummary.append(r);ds[k]=r
need(all(ds[k]['PF00153_hit_count']>0 for k in targets),'Target reference domain control failed')
table('07_domain_architecture_summary.tsv',dsummary,list(dsummary[0]))
hits=defaultdict(list);annot=[]
for line in (out/'04_all_panel_hits.tsv').read_text().splitlines():
 x=line.split('\t');need(len(x)==12,'Unexpected DIAMOND format');q,k=x[:2];need(q in proteins and k in refs,'Unknown similarity ID');need(int(x[10])==len(proteins[q]) and int(x[11])==len(refs[k]),'Hit sequence length mismatch')
 r=dict(protein_id=q,reference_id=k,reference_group=targets.get(k,'OTHER_CARRIER'),reference_name=refmeta[k]['protein_name'],pident=float(x[2]),alignment_length=int(x[3]),query_span=round((int(x[5])-int(x[4])+1)/int(x[10]),6),reference_span=round((int(x[7])-int(x[6])+1)/int(x[11]),6),evalue=float(x[8]),bitscore=float(x[9]));hits[q].append(r);annot.append(r)
table('08_annotated_panel_hits.tsv',annot,['protein_id','reference_id','reference_group','reference_name','pident','alignment_length','query_span','reference_span','evalue','bitscore'])
review=[];familyrows=[]
for original in queue:
 q=original['protein_id'];hh=sorted(hits[q],key=lambda h:(-h['bitscore'],h['evalue'],h['reference_id']));best=hh[0] if hh else None;flags=['BOUNDED_PANEL_NOT_FUNCTIONAL_ASSIGNMENT'];priority=[]
 for fam in ['DTC','DIC']:
  own=[h for h in hh if h['reference_group']==fam];other=[h for h in hh if h['reference_group']!=fam];a=own[0] if own else None;b=other[0] if other else None;margin=a['bitscore']-b['bitscore'] if a and b else None
  # Inclusive triage only: top or within 20 bits of any competing reference.
  selected=bool(a and (not b or margin>=-20))
  if selected:priority.append(fam)
  familyrows.append(dict(protein_id=q,species=original['species'],gene_id=original['gene_id'],target_family=fam,best_target_reference=a['reference_id'] if a else 'NONE',target_bitscore=a['bitscore'] if a else 'NA',target_evalue=a['evalue'] if a else 'NA',target_query_span=a['query_span'] if a else 'NA',target_reference_span=a['reference_span'] if a else 'NA',best_competing_reference=b['reference_id'] if b else 'NONE',competing_group=b['reference_group'] if b else 'NONE',competing_bitscore=b['bitscore'] if b else 'NA',target_minus_competitor_bits=round(margin,3) if margin is not None else 'NA',phylogeny_triage=selected,decision='REVIEW_ONLY'))
 if not best:flags.append('NO_PANEL_HIT')
 elif min(best['query_span'],best['reference_span'])<.7:flags.append('LIMITED_BEST_HIT_SPAN')
 if ds[q]['PF00153_hit_count']==0:flags.append('NO_PF00153_AT_GA')
 elif ds[q]['PF00153_hit_count']!=3:flags.append('PF00153_HIT_COUNT_NOT_THREE_REVIEW_ONLY')
 if ds[q]['other_Pfam_models']:flags.append('ADDITIONAL_DOMAIN_MODELS')
 if best and len({h['reference_group'] for h in hh if best['bitscore']-h['bitscore']<=20})>1:flags.append('CLOSE_COMPETING_GROUPS')
 review.append(dict(protein_id=q,species=original['species'],gene_id=original['gene_id'],prior_best_panel_family=original['best_panel_family'],best_reference=best['reference_id'] if best else 'NONE',best_group=best['reference_group'] if best else 'NONE',best_bitscore=best['bitscore'] if best else 'NA',best_query_span=best['query_span'] if best else 'NA',best_reference_span=best['reference_span'] if best else 'NA',PF00153_hit_count=ds[q]['PF00153_hit_count'],other_Pfam_models=ds[q]['other_Pfam_models'],phylogeny_priority_families=';'.join(priority) or 'NONE',flags=';'.join(flags),decision='REVIEW_ONLY',final_count='NA'))
table('09_all_candidates_specificity_REVIEW_ONLY.tsv',review,list(review[0]));table('10_target_competitor_comparisons.tsv',familyrows,list(familyrows[0]));prioritized=[r for r in review if r['phylogeny_priority_families']!='NONE'];table('11_phylogeny_priority_queue.tsv',prioritized,list(review[0]));writefa('11_phylogeny_priority_candidates.faa',{r['protein_id']:proteins[r['protein_id']] for r in prioritized})
for n in ['01_integrated_locus_evidence.tsv','06_refreshed_family_coverage.tsv','07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/('prior_'+n))
run(['diamond','version'],'diamond_version.txt');run(['hmmsearch','-h'],'hmmer_version.txt')
(out/'README.txt').write_text('Bounded specificity screen: 724 existing candidates versus 43 references; not a new whole-proteome discovery. Accession-defined targets DTC Q9C5M0; DIC Q9SJY5/Q9SB52/Q9FY68, regardless of PUMP naming. Other references are comparison sequences, not automatically validated subtype labels. DIAMOND very-sensitive E<=1e-3, all reported reference targets, one HSP per pair; spans are single-HSP spans, not aggregated coverage. All reported hits retained. Pfam GA search covers candidates and all controls. PF00153 repeat hit count and overlaps retained; three hits is a review heuristic, not proof of completeness. Family-specific triage includes a target hit with score within 20 bits of best other-group hit or no competitor; no domain/coverage filter removes triage candidates. This heuristic prioritizes review, does not accept/exclude genes. All other candidates retained in full table/FASTA. No trees or count updates. Additional domains, fragments and weak/competing hits need bounded review before phylogeny. No internet needed.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),pfam_sha256=sha(db),pfam_models=models,target_accession_definitions=targets,commands=commands),indent=2)+'\n')
qc=dict(status='CARRIER_SPECIFICITY_REVIEW_READY',candidates=len(proteins),references=len(refs),reported_similarity_hits=len(annot),candidates_with_hits=len(hits),candidate_best_groups=dict(Counter(r['best_group'] for r in review)),candidates_with_PF00153=sum(ds[k]['PF00153_hit_count']>0 for k in proteins),priority_unique_candidates=len(prioritized),priority_family_rows=sum(r['phylogeny_triage'] for r in familyrows),final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
for p in out.glob('*'):
 if p.suffix in {'.log','.domtbl'}:p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(dbuse),'<PFAM_HMM>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts);(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Carrier_specificity_screen/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
