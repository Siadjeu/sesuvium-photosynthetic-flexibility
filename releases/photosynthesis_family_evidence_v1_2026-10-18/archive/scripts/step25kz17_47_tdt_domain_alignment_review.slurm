#!/usr/bin/env bash
#SBATCH --job-name=tdt_validate
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=tdt_validate_%j.out
#SBATCH --error=tdt_validate_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v hmmsearch
command -v mafft
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
src=release('step25kz17.46','TRANSPORTER_RECOVERY_RELEASE')
allrows=read(src/'07_all_candidate_evidence_REVIEW_ONLY.tsv');selected=[r for r in allrows if r['best_panel_family']=='TDT']
need(len(selected)==11 and all(r['mapping_status']=='UNAMBIGUOUS' for r in selected),'Expected 11 mapped TDT candidates')
seqs,_=fasta(src/'05_all_discovered_and_prior_candidates.faa');proteins={r['protein_id']:seqs[r['protein_id']] for r in selected}
refs,_=fasta(src/'02_reference_panel.faa');refmeta=read(src/'01_reference_definitions.tsv');refids={r['family']:r['reference_id'] for r in refmeta if r['family'] in {'TDT','DIT1','DIT2'}}
need(set(refids)=={'TDT','DIT1','DIT2'},'Missing TDT/DIT controls');tdt=refids['TDT'];controls={k:refs[k] for k in refids.values()}
for r in refmeta:
 if r['reference_id'] in controls:need(hashlib.sha256(controls[r['reference_id']].encode()).hexdigest()==r['sequence_sha256'],'Reference hash mismatch')
out=base/'step25kz17.47'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir();(out/'annotation').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+s+'\n' for k,s in d.items()))
writefa('01_TDT_candidates.faa',proteins);writefa('02_domain_queries.faa',{**proteins,**controls});writefa('03_TDT_alignment_input.faa',{**proteins,tdt:refs[tdt]})
table('01_candidate_screen_evidence.tsv',selected,list(selected[0]))
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
run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'04_pfam.domtbl'),str(dbuse),str(out/'02_domain_queries.faa')],'04_hmmsearch.log')
run(['mafft','--amino','--thread',cpu,'--localpair','--maxiterate','1000',str(out/'03_TDT_alignment_input.faa')],'05_mafft.log','05_TDT_reference_alignment.faa')
aln,_=fasta(out/'05_TDT_reference_alignment.faa');need(set(aln)==set(proteins)|{tdt} and len({len(s) for s in aln.values()})==1,'Invalid TDT alignment')
need(all(aln[k].replace('-','')==s for k,s in {**proteins,tdt:refs[tdt]}.items()),'Alignment changed sequences')
dh=[];domains=defaultdict(list)
for line in (out/'04_pfam.domtbl').read_text().splitlines():
 if not line or line.startswith('#'):continue
 x=line.split(maxsplit=22);r=dict(sequence_id=x[0],model=x[3],accession=x[4],model_length=int(x[5]),i_evalue=float(x[12]),score=float(x[13]),hmm_from=int(x[15]),hmm_to=int(x[16]),ali_from=int(x[17]),ali_to=int(x[18]));dh.append(r);domains[x[0]].append(r)
table('06_all_Pfam_domain_hits.tsv',dh,['sequence_id','model','accession','model_length','i_evalue','score','hmm_from','hmm_to','ali_from','ali_to'])
reference_models={r['accession'].split('.')[0] for r in domains[tdt]};need(reference_models,'TDT reference has no Pfam support')
coverage=[]
for k in list(proteins)+list(controls):
 for model in sorted({h['accession'].split('.')[0] for h in domains[k]}|reference_models):
  hh=[h for h in domains[k] if h['accession'].split('.')[0]==model];covered=set()
  for h in hh:covered.update(range(h['hmm_from'],h['hmm_to']+1))
  coverage.append(dict(sequence_id=k,model=model,hit_count=len(hh),union_model_coverage=round(len(covered)/hh[0]['model_length'],6) if hh else 0,protein_intervals=';'.join(str(h['ali_from'])+'-'+str(h['ali_to']) for h in hh),note='Overlapping Pfam models and repeated/split hits retained; not independent domain counts'))
table('07_domain_model_coverage.tsv',coverage,list(coverage[0]))
cols=[i for i,c in enumerate(aln[tdt]) if c!='-'];review=[];gaps=[]
for r in selected:
 k=r['protein_id'];s=aln[k];present=[j for j,i in enumerate(cols) if s[i]!='-'];n=len(present);matches=sum(s[i]==aln[tdt][i] for i in cols if s[i]!='-')
 models_k={h['accession'].split('.')[0] for h in domains[k]};flags=['LIMITED_REFERENCE_PANEL_NO_FUNCTIONAL_VALIDATION']
 if n/len(cols)<.9:flags.append('REFERENCE_OCCUPANCY_LT_90_PERCENT')
 if reference_models-models_k:flags.append('MISSING_REFERENCE_PFAM_MODEL')
 if models_k-reference_models:flags.append('ADDITIONAL_PFAM_MODELS')
 start=None
 for j in range(len(cols)+1):
  gap=j<len(cols) and s[cols[j]]=='-'
  if gap and start is None:start=j
  if not gap and start is not None:
   gaps.append(dict(protein_id=k,reference_start=start+1,reference_end=j,gap_length=j-start,alignment_start=cols[start]+1,alignment_end=cols[j-1]+1));start=None
 review.append(dict(protein_id=k,species=r['species'],gene_id=r['gene_id'],length_aa=len(proteins[k]),reference_occupancy=round(n/len(cols),6),identity_over_paired=round(matches/n,6) if n else 0,N_terminal_reference_gaps=present[0] if present else len(cols),C_terminal_reference_gaps=len(cols)-present[-1]-1 if present else len(cols),internal_reference_gaps=sum(s[cols[j]]=='-' for j in range(present[0],present[-1]+1)) if present else 0,candidate_inserted_residues=sum(a!='-' and b=='-' for a,b in zip(s,aln[tdt])),Pfam_models=';'.join(sorted(models_k)),flags=';'.join(flags),decision='REVIEW_ONLY',final_count='NA'))
table('08_TDT_completeness_review.tsv',review,list(review[0]));table('09_reference_gap_intervals.tsv',gaps,['protein_id','reference_start','reference_end','gap_length','alignment_start','alignment_end'])
# Bounded annotation collection for the pre-identified shorter SP_Hex protein.
priority='SP_Hex__gene_Spor0444110.1';need(priority in proteins,'Expected shorter SP_Hex candidate missing')
rr=next(r for r in selected if r['protein_id']==priority);sp=rr['species'];gene=rr['gene_id'];mp=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv')
need(sha(mp)==json.loads((src/'PROVENANCE.json').read_text())['mapping_hashes'][sp],'Mapping changed since discovery')
mm=[r for r in read(mp) if r['gene_id']==gene];need(mm and len({(r['seqid'],r['strand']) for r in mm})==1,'Invalid exception mapping')
table('annotation/01_all_isoform_mapping.tsv',mm,list(mm[0]));aliases=set()
for m in mm:
 for f in ['gene_id','transcript_id','protein_id','cds_id']:
  if m.get(f) and m[f] not in {'.','NA'}:aliases|={m[f],sp+'__'+m[f]}
need(priority in aliases,'Candidate not confirmed by exact mapping')
chrom=mm[0]['seqid'];start=max(1,min(int(m['start']) for m in mm)-5000);end=max(int(m['end']) for m in mm)+5000
manifest=root/'00_admin/manifests/input_sets.tsv';config=[r for r in read(manifest) if r['input_set']==sp];need(len(config)==1,'Ambiguous input manifest');status=[];provenance=[]
def op(p):return gzip.open(p,'rt') if p.suffix=='.gz' else p.open()
for kind,col in [('gff','gff3'),('protein','protein_all'),('cds','cds_all')]:
 p=root/config[0][col];provenance.append(dict(kind=kind,source=str(p.relative_to(root)) if p.is_relative_to(root) else '<EXTERNAL>/'+p.name,sha256=sha(p) if p.is_file() else 'MISSING'))
 if not p.is_file():status.append(dict(kind=kind,status='MISSING',records=0));continue
 n=0
 if kind=='gff':
  with op(p) as f,(out/'annotation/02_context.gff3').open('w') as w:
   w.write('##gff-version 3\n')
   for line in f:
    if line.startswith('##FASTA'):break
    if line.startswith('#'):continue
    c=line.rstrip('\n').split('\t')
    if len(c)==9 and c[0]==chrom and int(c[3])<=end and int(c[4])>=start:w.write(line.rstrip('\n')+'\n');n+=1
 else:
  found={};key=None;parts=[]
  def commit():
   if key in aliases:
    need(key not in found,'Duplicate source sequence ID');found[key]=''.join(parts).upper()
  with op(p) as f:
   for line in f:
    if line.startswith('>'):commit();key=line[1:].split()[0];parts=[]
    elif key in aliases:parts.append(''.join(line.split()))
   commit()
  writefa('annotation/03_'+kind+'_isoforms.fasta',found);n=len(found)
  if kind=='protein':
   status.append(dict(kind='priority_protein_sequence',status='EXACT_MATCH' if any(norm(s)==proteins[priority] for s in found.values()) else 'NO_EXACT_SEQUENCE_MATCH',records=sum(norm(s)==proteins[priority] for s in found.values())))
 status.append(dict(kind=kind,status='COLLECTED' if n else 'NO_EXACT_ID_MATCH',records=n))
table('annotation/04_collection_status.tsv',status,['kind','status','records']);table('annotation/05_input_provenance.tsv',provenance,['kind','source','sha256'])
# Carry forward all mitochondrial candidates; no broad-carrier subtype assignment here.
mito=[r for r in allrows if r['best_panel_family'] in {'DTC','DIC','UCP1'}]
table('10_mitochondrial_carrier_queue_UNCHANGED.tsv',mito,list(allrows[0]))
for n in ['prior_01_integrated_locus_evidence.tsv','prior_06_refreshed_family_coverage.tsv','prior_07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/n)
run(['hmmsearch','-h'],'hmmer_version.txt');run(['mafft','--version'],'mafft_version.txt')
(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),pfam_sha256=sha(db),pfam_models=models,mapping_sha256=sha(mp),input_manifest_sha256=sha(manifest),commands=commands),indent=2)+'\n')
(out/'README.txt').write_text('Focused TDT domain/completeness review: 11 best-panel TDT candidates plus TDT control aligned de novo. DIT1/DIT2 controls included in Pfam search only; shared domains do not distinguish TDT from DIT. All Pfam overlaps and split/repeated hits retained; union model coverage is descriptive, not proof of an intact domain. Alignment occupancy and 90% flag are review heuristics, not functional acceptance. Exact sequences/IDs verified; distinct loci retained. Shorter SP_Hex candidate annotation collected using exact mapped aliases with +/-5kb GFF context; neighbors are not assigned to the candidate. Missing files/IDs remain explicit. No genome/read validation, repair, phylogeny, new carrier panel retrieval or final counting. DTC/DIC/UCP1 candidate queue carried unchanged for subsequent broader mitochondrial-carrier comparator curation. No internet required.\n')
qc=dict(status='TDT_DOMAIN_ALIGNMENT_REVIEW_READY',TDT_candidates=len(proteins),alignment_sequences=len(aln),alignment_columns=len(next(iter(aln.values()))),reference_Pfam_models=sorted(reference_models),candidates_with_additional_flags=sum(r['flags']!='LIMITED_REFERENCE_PANEL_NO_FUNCTIONAL_VALIDATION' for r in review),mitochondrial_queue_retained=len(mito),annotation_status=status,final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
for p in out.glob('*'):
 if p.suffix in {'.log','.domtbl'}:p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(dbuse),'<PFAM_HMM>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts);(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='TDT_domain_alignment_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
