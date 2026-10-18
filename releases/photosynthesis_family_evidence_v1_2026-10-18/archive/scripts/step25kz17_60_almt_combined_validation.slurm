#!/usr/bin/env bash
#SBATCH --job-name=almt_combined
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=almt_combined_%j.out
#SBATCH --error=almt_combined_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v diamond
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
src=release('step25kz17.59','DIT_ALMT_RELEASE');cor=release('step25kz17.21','CORRECTED_INPUT_RELEASE')
proteome=cor/'full_outputs/corrected_protein_candidates.faa';expected=[r['sha256'] for r in read(cor/'08_full_output_checksums.tsv') if r['file']=='full_outputs/corrected_protein_candidates.faa'];need(len(expected)==1 and proteome.is_file() and sha(proteome)==expected[0],'Corrected proteome provenance mismatch')
inv=read(src/'05_ALMT_existing_file_inventory.tsv')
def source(path):
 rr=[r for r in inv if r['source_file']==path];need(len(rr)==1 and rr[0]['bundle_file'],'Required inventoried file unavailable '+path);p=src/rr[0]['bundle_file'];need(sha(p)==rr[0]['sha256'],'Inventory provenance mismatch');return p
reg=read(source('08_focused_families/validation/step25kn/results/ALMT/candidate_validation.tsv'));need(len(reg)==142,'Expected142 historical records')
refs,headers=fasta(source('08_focused_families/validation/step25kn/references/ALMT.faa'));historical,_=fasta(source('08_focused_families/validation/step25kz17.2/release/ALMT/primary_candidates.faa'))
ids={}
for r in reg:
 sp,sep,k=r['sequence_id'].partition('|');need(sep and k.startswith(sp+'__') and k not in ids,'Unexpected historical ID');ids[k]=(sp,r)
# Stream selected corrected sequences without loading the whole proteome.
current={};key=None;parts=[]
def commit():
 if key in ids:need(key not in current,'Duplicate current ID');current[key]=''.join(parts).upper()
with proteome.open() as f:
 for l in f:
  if l.startswith('>'):commit();key=l[1:].split()[0];parts=[]
  elif key in ids:parts.append(''.join(l.split()))
 commit()
out=base/'step25kz17.60'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
positive={};other={};refrows=[]
for k,v in refs.items():
 acc=k.split('|')[1];h=headers[k];role='ALMT_LABELLED_POSITIVE' if re.search(r'GN=ALMT\d+\b',h) else 'OTHER_REFERENCE_CONTROL';rid='REF__'+acc;(positive if role=='ALMT_LABELLED_POSITIVE' else other)[rid]=v;refrows.append(dict(reference_id=rid,original_id=k,role=role,source_header=h,sequence_sha256=hashlib.sha256(v.encode()).hexdigest(),curation_status='HISTORICAL_ANNOTATION_NOT_INDEPENDENT_FUNCTIONAL_VERIFICATION'))
need(len(positive)==14 and len(other)==5,'Reference composition changed');table('01_reference_role_audit.tsv',refrows,list(refrows[0]));writefa('02_ALMT_positive_references.faa',positive);writefa('02_other_reference_controls.faa',other)
mapping={};mhash={}
for sp in sorted({v[0] for v in ids.values()}):
 p=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv');need(p.is_file(),'Missing mapping '+sp);mhash[sp]=sha(p);mapping[sp]=read(p)
audit=[];proteins={}
for k,(sp,r) in ids.items():
 local=k[len(sp)+2:];mm=[m for m in mapping[sp] if local in {m.get(f,'') for f in ['gene_id','transcript_id','protein_id','cds_id']} or k in {m.get(f,'') for f in ['gene_id','transcript_id','protein_id','cds_id']}];genes={m['gene_id'] for m in mm};seq=norm(current.get(k,''));old=historical.get('ALMT__'+r['sequence_id'],'');status='READY_FOR_SCREEN' if len(genes)==1 and valid(seq) else 'UNRESOLVED_MAPPING_OR_SEQUENCE'
 if status=='READY_FOR_SCREEN':proteins[k]=seq
 audit.append(dict(protein_id=k,species=sp,gene_id=next(iter(genes)) if len(genes)==1 else '',mapping_gene_count=len(genes),historical_status=r['status'],historical_sequence_available=bool(old),historical_equals_corrected=norm(old)==seq if old and seq else 'NA',corrected_sequence_found=k in current,length_aa=len(seq),input_status=status,final_count='NA'))
need(proteins,'No reconciled candidates');table('03_all_142_candidate_reconciliation.tsv',audit,list(audit[0]));writefa('03_reconciled_candidates.faa',proteins);writefa('01_domain_queries.faa',{**proteins,**positive,**other})
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
run(['diamond','makedb','--threads',cpu,'--in',str(out/'02_ALMT_positive_references.faa'),'--db',str(out/'work/ALMT')],'05_diamond_makedb.log')
run(['diamond','blastp','--threads',cpu,'--query',str(out/'03_reconciled_candidates.faa'),'--db',str(out/'work/ALMT'),'--out',str(out/'05_positive_reference_hits.tsv'),'--very-sensitive','--evalue','1e-3','--max-target-seqs','0','--max-hsps','1','--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore'],'05_diamond.log')
best={}
for l in (out/'05_positive_reference_hits.tsv').read_text().splitlines():
 x=l.split('\t');need(len(x)==12,'Unexpected similarity format');k=x[0]
 if k not in best or float(x[11])>float(best[k][11]):best[k]=x
review=[];architecture=[]
for k in list(proteins)+list(positive)+list(other):
 hh=domains[k];models_k=sorted({h['accession'].split('.')[0] for h in hh})
 for m in models_k:
  h=[h for h in hh if h['accession'].split('.')[0]==m];covered=set()
  for d in h:covered.update(range(d['hmm_from'],d['hmm_to']+1))
  architecture.append(dict(sequence_id=k,model=m,hit_count=len(h),union_model_coverage=round(len(covered)/h[0]['model_length'],6),protein_intervals=';'.join(str(d['ali_from'])+'-'+str(d['ali_to']) for d in h)))
# Derive shared positive-reference models from measured controls, avoiding guessed identifiers.
freq=Counter(m for k in positive for m in {h['accession'].split('.')[0] for h in domains[k]});shared={m for m,n in freq.items() if n>=len(positive)*.8};need(shared,'No Pfam model shared by>=80% positive references; inspect controls')
for r in audit:
 k=r['protein_id'];flags=[];x=best.get(k);mm={h['accession'].split('.')[0] for h in domains[k]}
 if r['input_status']!='READY_FOR_SCREEN':flags.append(r['input_status'])
 else:
  if shared-mm:flags.append('MISSING_COMMON_POSITIVE_REFERENCE_MODEL')
  if mm-shared:flags.append('ADDITIONAL_MODELS_REVIEW')
  if not x:flags.append('NO_POSITIVE_REFERENCE_SIMILARITY')
  elif (int(x[9])-int(x[8])+1)/int(x[5])<.8:flags.append('REFERENCE_SPAN_LT_80_PERCENT')
  if any(a['union_model_coverage']<.8 for a in architecture if a['sequence_id']==k and a['model'] in shared):flags.append('MODEL_COVERAGE_LT_80_PERCENT')
  if 'X' in proteins[k]:flags.append('AMBIGUOUS_RESIDUES')
 if r['historical_status']!='PASS':flags.append('HISTORICAL_FAILURE_RETAINED')
 review.append(dict(r,models=';'.join(sorted(mm)),best_positive_reference=x[1] if x else '',best_bitscore=x[11] if x else '',best_evalue=x[10] if x else '',query_span=round((int(x[7])-int(x[6])+1)/int(x[4]),6) if x else '',reference_span=round((int(x[9])-int(x[8])+1)/int(x[5]),6) if x else '',flags=';'.join(flags) or 'NO_HEURISTIC_EXCEPTION',decision='REVIEW_ONLY'))
table('07_domain_architecture.tsv',architecture,['sequence_id','model','hit_count','union_model_coverage','protein_intervals']);table('08_combined_ALMT_review.tsv',review,list(review[0]));table('09_exception_queue.tsv',[r for r in review if r['flags']!='NO_HEURISTIC_EXCEPTION'],list(review[0]))
for n in ['01_integrated_locus_evidence.tsv','06_refreshed_family_coverage.tsv','07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/('prior_'+n))
run(['hmmsearch','-h'],'hmmer_version.txt');run(['diamond','version'],'diamond_version.txt')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step59_manifest_sha256=sha(src/'SHA256SUMS.txt'),corrected_manifest_sha256=sha(cor/'SHA256SUMS.txt'),corrected_proteome_sha256=sha(proteome),mapping_hashes=mhash,pfam_sha256=sha(db),commands=commands),indent=2)+'\n')
qc=dict(status='ALMT_RECONCILIATION_DOMAIN_SCREEN_READY',historical_records=len(audit),screened_candidates=len(proteins),unresolved_inputs=len(audit)-len(proteins),positive_references=len(positive),other_controls=len(other),shared_positive_reference_models=sorted(shared),exception_records=sum(r['flags']!='NO_HEURISTIC_EXCEPTION' for r in review),final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'README.txt').write_text('All142 historical records retained including3 failures. Exact species-scoped mapping aliases only; no transcript suffix stripping or nearest-ID recovery. Current corrected proteins screened when mapping resolves to one gene and sequence passes amino-acid validation. Historical139 sequences compared where available.14 historically ALMT-labelled references separate from STOP1/DTX controls; labels are not independent functional verification. Full Pfam GA search includes candidates and all19 controls; DIAMOND compares candidates only against14 positive references. Shared Pfam models derived from>=80% positive controls.80% span/model thresholds are review heuristics; HSP span is not alignment occupancy or proof of completeness. Additional models retained for review. No new phylogeny, genome repair, count promotion or whole-proteome discovery; search completeness remains unknown.\n')
for p in out.glob('*'):
 if p.suffix in {'.log','.domtbl'}:p.write_text(p.read_text().replace(str(dbuse),'<PFAM_HMM>').replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts);(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='ALMT_combined_validation/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
