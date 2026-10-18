#!/usr/bin/env bash
#SBATCH --job-name=td029_review
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=td029_review_%j.out
#SBATCH --error=td029_review_%j.err
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

src=release('step25kz17.48','TDT_TREE_RELEASE')
iq=next((shutil.which(x) for x in ['iqtree3','iqtree2','iqtree'] if shutil.which(x)),None);need(iq,'IQ-TREE not found')
meta=read(src/'02_tree_id_map.tsv');target=next((r for r in meta if r['tree_id']=='TD029'),None);need(target and target['protein_id']=='SP_Hex__gene_Spor0030220.1','TD029 identity mismatch')
seqs,_=fasta(src/'01_all_sequences.faa');refids={r['tree_id'] for r in meta if r['role']=='CURATED_REFERENCE'};need(len(refids)==3,'Expected three controls')
out=base/'step25kz17.49'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
writefa('01_domain_queries.faa',{k:seqs[k] for k in ['TD029']+sorted(refids)})
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

run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'02_pfam.domtbl'),str(dbuse),str(out/'01_domain_queries.faa')],'02_hmmsearch.log')
hits=[]
for line in (out/'02_pfam.domtbl').read_text().splitlines():
 if not line or line.startswith('#'):continue
 x=line.split(maxsplit=22);hits.append(dict(sequence_id=x[0],model=x[3],accession=x[4],model_length=int(x[5]),i_evalue=x[12],score=x[13],hmm_from=int(x[15]),hmm_to=int(x[16]),ali_from=int(x[17]),ali_to=int(x[18])))
table('03_domain_hits.tsv',hits,['sequence_id','model','accession','model_length','i_evalue','score','hmm_from','hmm_to','ali_from','ali_to'])
summary=[]
for k in ['TD029']+sorted(refids):
 hh=[h for h in hits if h['sequence_id']==k];covered=set()
 for h in hh:covered.update(range(h['ali_from'],h['ali_to']+1))
 summary.append(dict(sequence_id=k,length_aa=len(seqs[k]),domain_hits=len(hh),models=';'.join(sorted({h['accession'] for h in hh})) or 'NONE',union_protein_coverage=round(len(covered)/len(seqs[k]),6),interpretation='REVIEW_ONLY; shared domains do not establish specificity; no hit does not prove absence'))
table('04_domain_summary.tsv',summary,list(summary[0]))
# Parse unrooted edge bipartitions, avoiding root-dependent clade assignments.
def edges(path):
 tokens=re.findall(r'[^\s(),:;]+|[(),:;]',path.read_text());pos=0;ee=[]
 def node():
  nonlocal pos
  label=''
  if tokens[pos]=='(':
   pos+=1;leaves=set()
   while True:
    leaves|=node()
    if tokens[pos]==',':pos+=1;continue
    need(tokens[pos]==')','Malformed Newick');pos+=1;break
   if tokens[pos] not in {':',',',')',';'}:label=tokens[pos];pos+=1
  else:leaves={tokens[pos]};pos+=1
  if tokens[pos]==':':pos+=2
  ee.append((set(leaves),label));return leaves
 allids=node();need(tokens[pos]==';','Invalid Newick ending');return allids,ee
comparisons=[];runs=[]
for mode,seed in [('main_without_short',174901),('sensitivity_with_short',174902)]:
 alignment,_=fasta(src/mode/'alignment.faa');need('TD029' in alignment,'TD029 missing from source tree');need(len({len(v) for v in alignment.values()})==1,'Invalid source alignment')
 report=(src/mode/'tree.iqtree').read_text();mm=re.search(r'Best-fit model according to BIC:\s*(\S+)',report);need(mm,'Source model missing');model=mm.group(1)
 dd={k:v for k,v in alignment.items() if k!='TD029'};folder=out/mode;folder.mkdir();writefa(mode+'/alignment.faa',dd)
 prefix=folder/'tree';run([iq,'-s',str(folder/'alignment.faa'),'-st','AA','-m',model,'-B','1000','-bnni','-alrt','1000','-T',cpu,'-seed',str(seed),'-keep-ident','--prefix',str(prefix)],mode+'/iqtree.stdout.log')
 newtree=Path(str(prefix)+'.treefile');need(newtree.is_file(),'Tree missing');ids,ee=edges(newtree);need(ids==set(dd),'Tree tips differ from input')
 expected={r['tree_id'] for r in meta if r['screen_group']=='TDT'} & ids
 for stage,path in [('source_with_TD029',src/mode/'tree.treefile'),('new_without_TD029',newtree)]:
  universe,es=edges(path);support=[label for side,label in es if side==expected or universe-side==expected]
  comparisons.append(dict(analysis=mode,stage=stage,TDT_group_size=len(expected),exact_TDT_bipartition=len(support)>0,support_SHaLRT_UFBoot=';'.join(support) if support else 'NOT_FOUND',TDT_members=';'.join(sorted(expected))))
 labels=['LABELS','SEPARATOR TAB','DATA']+[r['tree_id']+'\t'+r['protein_id']+' | '+r['screen_group']+' | '+r['role'] for r in meta if r['tree_id'] in dd];(folder/'itol_labels.txt').write_text('\n'.join(labels)+'\n')
 runs.append(dict(analysis=mode,sequences=len(dd),columns=len(next(iter(dd.values()))),fixed_source_model=model,seed=seed,rooting='UNROOTED'))
table('05_TDT_split_comparison.tsv',comparisons,list(comparisons[0]));table('06_tree_run_summary.tsv',runs,list(runs[0]));table('07_tree_id_map.tsv',meta,list(meta[0]))
for n in ['03_candidate_screen_evidence.tsv','08_short_candidate_CDS_check.tsv','source_08_TDT_completeness_review.tsv','source_prior_01_integrated_locus_evidence.tsv','source_prior_06_refreshed_family_coverage.tsv','source_prior_07_final_Wang_counts_PENDING.tsv','source_10_mitochondrial_carrier_queue_UNCHANGED.tsv']:shutil.copyfile(src/n,out/n)
run([iq,'--version'],'iqtree_version.txt');run(['hmmsearch','-h'],'hmmer_version.txt')
(out/'README.txt').write_text('Bounded TD029 domain and tree sensitivity review. Pfam GA search: TD029 plus three curated TDT/DIT controls. Both source tree alignments retained exactly, removing only TD029; no realignment or remasking. Source model fixed, parameters re-estimated. 1000 SH-aLRT and 1000 UFBoot with bnni; unrooted. Exact TDT edge bipartitions compared with source trees, including/excluding short TD037 as previously. Branch support order SH-aLRT/UFBoot. No automated functional assignment, candidate exclusion, annotation repair or final counting. Original ledger and mitochondrial queue carried unchanged. Shared Pfam models are not specificity proof. Upload tree.treefile and accompanying itol_labels.txt to iTOL.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),pfam_sha256=sha(db),pfam_models=models,commands=commands),indent=2)+'\n')
qc=dict(status='TD029_DOMAIN_SENSITIVITY_REVIEW_READY',domain_queries=4,TD029_domain_hits=sum(h['sequence_id']=='TD029' for h in hits),trees=runs,TDT_split_comparisons=comparisons,final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
for p in out.rglob('*'):
 if p.is_file() and 'work' not in p.relative_to(out).parts and p.suffix in {'.log','.iqtree','.txt','.domtbl'}:p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(dbuse),'<PFAM_HMM>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts);(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='TD029_domain_tree_sensitivity/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
