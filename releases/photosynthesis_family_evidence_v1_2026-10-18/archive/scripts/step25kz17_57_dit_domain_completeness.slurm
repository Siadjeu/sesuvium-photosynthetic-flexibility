#!/usr/bin/env bash
#SBATCH --job-name=dit_domains
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=dit_domains_%j.out
#SBATCH --error=dit_domains_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v hmmsearch
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
src=release('step25kz17.56','DIT_RECONCILIATION_RELEASE');trees=release('step25kz17.48','TDT_TREE_RELEASE')
need(json.loads((src/'PROVENANCE.json').read_text())['step48_manifest_sha256']==sha(trees/'SHA256SUMS.txt'),'Alignment provenance mismatch')
proteins,_=fasta(src/'03_DIT_candidates_REVIEW_ONLY.faa');controls,_=fasta(src/'04_DIT_TDT_reference_controls.faa');queries={**proteins,**controls}
need(len(proteins)==38 and set(controls)=={'REF__Q9LXV3','REF__Q9FMF7','REF__Q8LG88'},'Unexpected candidate/control set')
need(all(valid(s) for s in queries.values()),'Invalid amino-acid sequences')
selected=[r for r in read(src/'02_DIT_candidate_reconciliation.tsv') if r['protein_id'] in proteins];need(len(selected)==38,'Candidate table mismatch')
oldaln,_=fasta(trees/'05_all_alignment.faa');meta=read(trees/'02_tree_id_map.tsv');lookup={r['protein_id']:r['tree_id'] for r in meta};aln={k:oldaln[lookup[k]] for k in queries}
need(len({len(s) for s in aln.values()})==1 and all(aln[k].replace('-','')==s for k,s in queries.items()),'Alignment/sequence provenance mismatch')
out=base/'step25kz17.57'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
writefa('01_domain_queries.faa',queries);writefa('02_existing_reference_alignment.faa',aln)
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
run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'04_pfam.domtbl'),str(dbuse),str(out/'01_domain_queries.faa')],'04_hmmsearch.log')
dh=[];domains=defaultdict(list)
for line in (out/'04_pfam.domtbl').read_text().splitlines():
 if not line or line.startswith('#'):continue
 x=line.split(maxsplit=22);r=dict(sequence_id=x[0],model=x[3],accession=x[4],model_length=int(x[5]),i_evalue=float(x[12]),score=float(x[13]),hmm_from=int(x[15]),hmm_to=int(x[16]),ali_from=int(x[17]),ali_to=int(x[18]));dh.append(r);domains[x[0]].append(r)
table('06_all_Pfam_domain_hits.tsv',dh,['sequence_id','model','accession','model_length','i_evalue','score','hmm_from','hmm_to','ali_from','ali_to'])
expected={'PF00939','PF03600'};modelrows=[];review=[];gaps=[]
for k in queries:
 for model in sorted(expected|{h['accession'].split('.')[0] for h in domains[k]}):
  hh=[h for h in domains[k] if h['accession'].split('.')[0]==model];covered=set()
  for h in hh:covered.update(range(h['hmm_from'],h['hmm_to']+1))
  modelrows.append(dict(sequence_id=k,model=model,hit_count=len(hh),union_model_coverage=round(len(covered)/hh[0]['model_length'],6) if hh else 0,protein_intervals=';'.join(str(h['ali_from'])+'-'+str(h['ali_to']) for h in hh)))
table('07_domain_model_coverage.tsv',modelrows,list(modelrows[0]))
for r in selected:
 k=r['protein_id'];ref={'DIT1':'REF__Q9LXV3','DIT2':'REF__Q9FMF7'}[r['screen_family']];cols=[i for i,c in enumerate(aln[ref]) if c!='-'];present=[j for j,i in enumerate(cols) if aln[k][i]!='-'];n=len(present);mk={h['accession'].split('.')[0] for h in domains[k]};flags=[]
 if n/len(cols)<.9:flags.append('REFERENCE_OCCUPANCY_LT_90_PERCENT')
 if expected-mk:flags.append('MISSING_REFERENCE_PFAM_MODEL')
 if mk-expected:flags.append('ADDITIONAL_PFAM_MODELS')
 if any(x['union_model_coverage']<.8 for x in modelrows if x['sequence_id']==k and x['model'] in expected):flags.append('MODEL_COVERAGE_LT_80_PERCENT')
 if 'X' in proteins[k]:flags.append('AMBIGUOUS_RESIDUES')
 if k=='PO__evm.TU.LG08.1714':flags.append('PRIOR_LONG_MODEL_ARCHITECTURE_CAVEAT')
 start=None
 for j in range(len(cols)+1):
  gap=j<len(cols) and aln[k][cols[j]]=='-'
  if gap and start is None:start=j
  if not gap and start is not None:gaps.append(dict(protein_id=k,reference_id=ref,reference_start=start+1,reference_end=j,gap_length=j-start));start=None
 review.append(dict(protein_id=k,species=r['species'],gene_id=r['gene_id'],tree_group=r['screen_family'],reference_id=ref,length_aa=len(proteins[k]),reference_occupancy=round(n/len(cols),6),N_terminal_reference_gaps=present[0] if present else len(cols),C_terminal_reference_gaps=len(cols)-present[-1]-1 if present else len(cols),internal_reference_gaps=present[-1]-present[0]+1-n if present else 0,inserted_residues=sum(a!='-' and b=='-' for a,b in zip(aln[k],aln[ref])),Pfam_models=';'.join(sorted(mk)),flags=';'.join(flags) or 'NO_HEURISTIC_EXCEPTION',decision='REVIEW_ONLY',final_count='NA'))
table('08_DIT_completeness_review.tsv',review,list(review[0]));exceptions=[r for r in review if r['flags']!='NO_HEURISTIC_EXCEPTION'];table('09_exception_queue.tsv',exceptions,list(review[0]));writefa('09_exception_sequences.faa',{r['protein_id']:proteins[r['protein_id']] for r in exceptions});table('10_reference_gap_intervals.tsv',gaps,['protein_id','reference_id','reference_start','reference_end','gap_length'])
for n in ['01_existing_tree_DIT_group_audit.tsv','02_DIT_candidate_reconciliation.tsv','07_refreshed_family_coverage.tsv','prior_01_integrated_locus_evidence.tsv','prior_07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/('source_'+n))
run(['hmmsearch','-h'],'hmmer_version.txt')
qc=dict(status='DIT_DOMAIN_COMPLETENESS_REVIEW_READY',candidates=38,controls=3,exception_candidates=len(exceptions),candidate_groups=dict(Counter(r['tree_group'] for r in review)),control_domains={k:sorted({h['accession'].split('.')[0] for h in domains[k]}) for k in controls},final_counts=False)
need(all(expected<={h['accession'].split('.')[0] for h in domains[k]} for k in controls),'Expected control domains missing; inspect Pfam database/results')
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),alignment_manifest_sha256=sha(trees/'SHA256SUMS.txt'),pfam_sha256=sha(db),pfam_models=models,commands=commands),indent=2)+'\n')
(out/'README.txt').write_text('Bounded full-Pfam GA search of38 candidates and3 controls. Existing exact-sequence-verified alignment reused; no new phylogeny or discovery search. DIT1 and DIT2 assessed against respective anchors. Shared PF00939/PF03600 domains do not distinguish transporter specificity. Overlapping models and split hits retained; union HMM coverage is descriptive, not proof of intact architecture. Occupancy90% and model coverage80% are review heuristics, not biological acceptance criteria. Terminal alignment gaps can reflect targeting-region divergence. Additional domains and prior PO LG08.1714 caveat remain explicit. All candidates REVIEW_ONLY; no count promotion, gene repair or functional claims.\n')
for p in out.glob('*'):
 if p.suffix in {'.log','.domtbl'}:p.write_text(p.read_text().replace(str(dbuse),'<PFAM_HMM>').replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts);(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='DIT_domain_completeness/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
