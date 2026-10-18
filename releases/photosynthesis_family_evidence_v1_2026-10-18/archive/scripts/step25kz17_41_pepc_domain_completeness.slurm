#!/usr/bin/env bash
#SBATCH --job-name=pepc_domains
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=pepc_domains_%j.out
#SBATCH --error=pepc_domains_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v hmmsearch
python3 - <<'PY'
import csv,gzip,hashlib,json,os,re,shutil,subprocess,tarfile
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
pp=[root/os.environ['PEPC_CORRECTION_RELEASE']] if os.environ.get('PEPC_CORRECTION_RELEASE') else [p for p in (base/'step25kz17.40').glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
need(len(pp)==1,'Set PEPC_CORRECTION_RELEASE to the intended step25kz17.40 release')
src=pp[0].resolve(strict=True);verify(src)
audit=read(src/'01_PEPC_corrected_provenance.tsv');comp=read(src/'03_sequence_comparison.tsv');seq,_=fasta(src/'05_collected_sequences.fasta')
need(len(audit)==55 and len({r['protein_id'] for r in audit})==55,'Expected 55 unique loci')
byid={r['protein_id']:r for r in audit}
good=[r for r in comp if r['sequence_relation']=='EXACT_NORMALIZED_MATCH'];bad=[r for r in comp if r['sequence_relation']=='HISTORICAL_SEQUENCE_EQUALS_RECOVERY_CDS']
need(len(good)==53 and len(bad)==2 and len(comp)==55,'Unexpected provenance categories')
need({r['protein_id'] for r in bad}=={'SS__g24599','SS__g8156'},'Unexpected CDS exceptions')
proteins={}
for r in good:
 k=r['protein_id'];s=seq['CORRECTED_PROTEIN|'+k]
 need(hashlib.sha256(s.encode()).hexdigest()==r['corrected_sha256'],'Protein hash mismatch: '+k)
 need(norm(s)==norm(seq['HISTORICAL|'+k]),'Sequence comparison mismatch')
 need(valid(norm(s)),'Invalid protein sequence '+k);proteins[k]=norm(s)
# Select only explicitly annotated PPC1-4 controls from the original mixed reference file.
refpath=(root/os.environ.get('PEPC_REFERENCE_FASTA','08_focused_families/validation/step25kn/references/PEPC.faa')).resolve(strict=True)
refhash=sha(refpath);refseq,headers=fasta(refpath);refs={};refmeta=[]
for k,s in refseq.items():
 m=re.search(r'\bGN=(PPC[1-4])(?:\s|$)',headers[k],re.I)
 if m:
  s=norm(s);need(valid(s),'Invalid selected reference '+k)
  rid='REF__'+k;refs[rid]=s;refmeta.append(dict(reference_id=rid,source_id=k,gene=m.group(1).upper(),length_aa=len(s),header=headers[k]))
need(len(refs)==4 and {r['gene'] for r in refmeta}=={'PPC1','PPC2','PPC3','PPC4'},'Expected exactly one explicitly GN=PPC1-4 annotated control each; inspect PEPC_REFERENCE_FASTA')
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
 need(len(paths)==1,'Set PFAM_HMM to the intended installed Pfam-A ASCII .hmm or .hmm.gz file');db=next(iter(paths))
out=base/'step25kz17.41'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+s+'\n' for k,s in d.items()))
writefa('01_PEPC_corrected_53_candidates.faa',proteins);writefa('02_PPC_reference_controls.faa',refs)
table('03_reference_metadata.tsv',refmeta,['reference_id','source_id','gene','length_aa','header'])
# Translate the two documented CDS only as separately labelled diagnostic fragments.
code=dict(zip((a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'),'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
fragments={};exceptions=[]
for r in bad:
 k=r['protein_id'];s=seq['HISTORICAL|'+k];cds=seq['RECOVERY_CDS|'+k]
 need(s==cds and len(cds)%3==0 and set(cds)<=set('ACGT'),'Unexpected recovery CDS')
 aa=norm(''.join(code[cds[i:i+3]] for i in range(0,len(cds),3)))
 need(valid(aa),'Diagnostic translation invalid');fragments['DIAGNOSTIC_ONLY__'+k]=aa
 exceptions.append(dict(protein_id=k,species=byid[k]['species'],gene_id=byid[k]['gene_id'],historical_input_type='CDS_MISUSED_AS_PROTEIN',cds_length=len(cds),translation_length=len(aa),decision='EXCLUDE_HISTORICAL_TREE_RECORD_LOCUS_UNRESOLVED',formal_count='NA'))
writefa('04_diagnostic_fragment_translations_NOT_TREE_INPUT.faa',fragments)
table('05_CDS_exception_ledger.tsv',exceptions,list(exceptions[0]))
queries={**proteins,**refs,**fragments};writefa('06_domain_search_queries.faa',queries)
dbhash=sha(db);dbuse=db
if db.suffix=='.gz':
 dbuse=out/'work/Pfam-A.hmm'
 with gzip.open(db,'rb') as a,dbuse.open('wb') as b:shutil.copyfileobj(a,b)
models=ga=acc=0;target={};meta={}
with dbuse.open() as f:
 for line in f:
  if line.startswith('HMMER3/'):models+=1;meta={}
  elif line.startswith('GA '):ga+=1
  elif line.startswith('ACC '):acc+=1;meta['accession']=line.split()[1]
  elif line.startswith('NAME '):meta['name']=line.split()[1]
  elif line.startswith('LENG '):meta['length']=int(line.split()[1])
  elif line.strip()=='//' and meta.get('accession','').split('.')[0]=='PF00311':target=meta.copy()
need(models>1000 and models==ga==acc and target.get('name')=='PEPcase','Expected complete Pfam-A with GA cutoffs and PF00311 PEPcase')
def run(args,name):
 with (out/name).open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
run(['hmmsearch','-h'],'hmmsearch_version.txt')
run(['hmmsearch','--cpu',os.environ.get('SLURM_CPUS_PER_TASK','4'),'--cut_ga','--noali','--domtblout',str(out/'07_pfam.domtbl'),str(dbuse),str(out/'06_domain_search_queries.faa')],'07_hmmsearch.log')
hits=[];byseq=defaultdict(list)
for line in (out/'07_pfam.domtbl').read_text().splitlines():
 if not line or line.startswith('#'):continue
 x=line.split(maxsplit=22);need(len(x)>=22,'Malformed domtbl row');need(x[0] in queries,'Unexpected domain target')
 r=dict(sequence_id=x[0],protein_length=int(x[2]),model=x[3],accession=x[4],model_length=int(x[5]),full_evalue=float(x[6]),domain_i_evalue=float(x[12]),domain_score=float(x[13]),hmm_from=int(x[15]),hmm_to=int(x[16]),ali_from=int(x[17]),ali_to=int(x[18]),env_from=int(x[19]),env_to=int(x[20]))
 hits.append(r);byseq[x[0]].append(r)
fields=['sequence_id','protein_length','model','accession','model_length','full_evalue','domain_i_evalue','domain_score','hmm_from','hmm_to','ali_from','ali_to','env_from','env_to']
table('08_all_Pfam_domain_hits.tsv',hits,fields)
def union(intervals):
 covered=set()
 for a,b in intervals:covered.update(range(a,b+1))
 return len(covered)
review=[]
for k,s in queries.items():
 hh=[h for h in byseq[k] if h['accession'].split('.')[0]=='PF00311']
 best=max(hh,key=lambda h:h['domain_score']) if hh else None
 cov=(best['hmm_to']-best['hmm_from']+1)/best['model_length'] if best else 0
 ucov=union([(h['hmm_from'],h['hmm_to']) for h in hh])/target['length']
 flags=[]
 if not hh:flags.append('NO_PF00311_GA_HIT')
 elif cov<0.9:flags.append('BEST_DOMAIN_MODEL_COVERAGE_LT_90_PERCENT')
 if len(hh)>1:flags.append('MULTIPLE_PROFILE_HITS_REVIEW_SPLIT_OR_REPEAT')
 if not 800<=len(s)<=1300:flags.append('OUTSIDE_HISTORICAL_LENGTH_RANGE')
 if 'X' in s:flags.append('AMBIGUOUS_X_RESIDUES')
 other=sorted({h['accession']+':'+h['model'] for h in byseq[k] if h['accession'].split('.')[0]!='PF00311'})
 if other:flags.append('ADDITIONAL_PFAM_MODELS_REVIEW')
 role='REFERENCE_CONTROL' if k in refs else 'DIAGNOSTIC_FRAGMENT_NOT_COUNTABLE' if k in fragments else 'CANDIDATE'
 review.append(dict(sequence_id=k,role=role,length_aa=len(s),PF00311_hit_count=len(hh),best_domain_model_coverage=round(cov,5),union_model_coverage=round(ucov,5),best_domain_i_evalue=best['domain_i_evalue'] if best else '',best_domain_score=best['domain_score'] if best else '',best_hmm_from=best['hmm_from'] if best else '',best_hmm_to=best['hmm_to'] if best else '',additional_models=';'.join(other),review_flags=';'.join(flags),evidence_category='DOMAIN_SUPPORTED_NO_PRIORITY_FLAG' if not flags else 'REVIEW_REQUIRED',final_membership='PENDING',C4_CAM_function='NOT_ESTABLISHED'))
table('09_PEPC_domain_completeness_review.tsv',review,list(review[0]))
controls_ok=all(any(h['accession'].split('.')[0]=='PF00311' for h in byseq[k]) for k in refs)
shutil.copyfile(src/'01_PEPC_corrected_provenance.tsv',out/'10_prior_provenance.tsv')
shutil.copyfile(src/'06_integrated_locus_evidence_UNCHANGED.tsv',out/'11_integrated_locus_evidence_UNCHANGED.tsv')
shutil.copyfile(src/'07_final_Wang_counts_PENDING.tsv',out/'12_final_Wang_counts_PENDING.tsv')
(out/'HISTORICAL_TREE_STATUS.txt').write_text('Historical 55-candidate PEPC tree superseded for this counting workflow: two inputs are CDS, not proteins. Earlier files are preserved. Rebuild alignment and tree from validated proteins after reviewing this release; pruning the old tree alone does not correct its input alignment. Diagnostic translations are not main-tree inputs or accepted loci.\n')
need(sha(db)==dbhash and sha(refpath)==refhash,'Input database changed during run')
qc=dict(status='DOMAIN_REVIEW_READY' if controls_ok else 'REFERENCE_CONTROL_FAILURE_REVIEW_REQUIRED',candidate_proteins=len(proteins),reference_controls=len(refs),diagnostic_fragments=len(fragments),all_reference_controls_have_PF00311=controls_ok,candidate_categories=dict(Counter(r['evidence_category'] for r in review if r['role']=='CANDIDATE')),candidate_PF00311_hits=sum(r['PF00311_hit_count']>0 for r in review if r['role']=='CANDIDATE'),whole_genome_completeness='NOT_TESTED',final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(prior_manifest_sha256=sha(src/'SHA256SUMS.txt'),reference_file_sha256=sha(refpath),pfam_source_sha256=dbhash,pfam_uncompressed_sha256=sha(dbuse),pfam_models=models,PEPC_model=target,method='hmmsearch --cut_ga against full Pfam-A',review_thresholds=dict(best_domain_model_coverage=0.9,historical_length_range=[800,1300])),indent=2)+'\n')
(out/'README.txt').write_text('53 provenance-verified proteins, four explicitly GN=PPC1-4 annotated controls, and two diagnostic CDS translations searched against complete installed Pfam-A using GA cutoffs. PF00311 PEPcase is the target family model. Domain coverage measures model coordinates, not proof of full-length gene models or function. 90% best-domain coverage and historical 800-1300 aa range are review heuristics, not acceptance thresholds. Union coverage is supplementary and can conceal split/repeated hits; retain all raw hits and flags. Other Pfam hits are retained without clan-overlap filtering and are not automatic exclusions. PPC controls do not assign plant/bacterial subtype or C4/CAM function. Diagnostic fragments remain separate, non-countable, with locus status unresolved even if a domain matches. No tree inference or whole-proteome completeness search; no formal count promotion. Integrated ledger unchanged.\n')
# Scrub paths from external-tool logs before bundling; record tool versions.
for p in [out/'07_hmmsearch.log',out/'07_pfam.domtbl',out/'hmmsearch_version.txt']:
 s=p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(dbuse),'<PFAM_HMM>').replace(str(Path.home()),'<HOME>');p.write_text(s)
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='PEPC_domain_completeness/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
