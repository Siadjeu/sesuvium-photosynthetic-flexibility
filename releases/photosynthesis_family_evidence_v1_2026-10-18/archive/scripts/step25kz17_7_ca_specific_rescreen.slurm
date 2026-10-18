#!/usr/bin/env bash
#SBATCH --job-name=ca_specific_screen
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=ca_specific_screen_%j.out
#SBATCH --error=ca_specific_screen_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
: "${BETA_REVIEW:?Export BETA_REVIEW as the completed Step25kz17.6 release directory}"
command -v python3 >/dev/null
command -v mafft >/dev/null
command -v hmmbuild >/dev/null
command -v hmmsearch >/dev/null
export CA_PROTEOME="${CA_PROTEOME:-$PROJECT_ROOT/08_focused_families/validation/step25kz16.7/release/work/all_species_primary_proteins.faa}"
export CA_REFERENCES="${CA_REFERENCES:-$PROJECT_ROOT/08_focused_families/validation/step25kn/references/BETA_CA.faa}"
python3 - <<'PY'
import os,csv,json,hashlib,subprocess,tarfile,shutil,math
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
review_dir=Path(os.environ['BETA_REVIEW']).resolve(strict=True)
proteome=Path(os.environ['CA_PROTEOME']).resolve(strict=True)
refs_path=Path(os.environ['CA_REFERENCES']).resolve(strict=True)
review_path=review_dir/'07_BETA_CA_manual_review_corrected.tsv'
for p in [proteome,refs_path,review_path]:
 if not p.is_file() or p.stat().st_size==0:raise SystemExit('Missing or empty required input')
# Explicit accession allowlists. Unknown/missing/duplicate references stop the job.
accessions={
 'BETA_CA':{'P27140','P42737','Q94CE3','Q94CE4','Q9C6F5','Q9ZUC2'},
 'ALPHA_CA':{'F4HUC4','F4IHR4','F4JIK2','O04846','Q8L817','Q9FM99','Q9FYE3','Q9SUB4'}}
threads=int(os.environ.get('SLURM_CPUS_PER_TASK','4'))
E_CUTOFF=1e-5
COVERAGE=0.50

def fasta(path):
 seqs={};key=None;parts=[]
 def save():
  if key is None:return
  seq=''.join(parts).upper()
  if not seq:raise SystemExit('Empty sequence: '+key)
  if key in seqs:raise SystemExit('Duplicate FASTA identifier: '+key)
  seqs[key]=seq
 with path.open() as f:
  for line in f:
   if line.startswith('>'):
    save();key=line[1:].split()[0];parts=[]
   elif line.strip():
    if key is None:raise SystemExit('Sequence before FASTA header')
    parts.append(''.join(line.split()))
 save();return seqs

def write_fasta(path,records):
 with path.open('w') as f:
  for k,s in records:
   f.write('>'+k+'\n')
   for start in range(0,len(s),80):f.write(s[start:start+80]+'\n')

def tsv(path,fields,rows):
 with path.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def run(args,stdout,stderr):
 with stdout.open('w') as o,stderr.open('w') as e:subprocess.run(args,stdout=o,stderr=e,check=True)

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()

refs=fasta(refs_path);groups={g:{} for g in accessions};ref_records=[]
for pid,seq in refs.items():
 tokens=pid.split('|')
 if len(tokens)<3:raise SystemExit('Expected accession-bearing reference ID: '+pid)
 acc=tokens[1];matches=[g for g,allowed in accessions.items() if acc in allowed]
 if len(matches)!=1:raise SystemExit('Unrecognized reference accession: '+acc)
 g=matches[0]
 if acc in groups[g]:raise SystemExit('Duplicate reference accession: '+acc)
 if any(c not in 'ACDEFGHIKLMNPQRSTVWYBXZJUO' for c in seq):raise SystemExit('Non-protein symbol in reference: '+acc)
 groups[g][acc]=(pid,seq)
 ref_records.append({'accession':acc,'reference_id':pid,'family':g,'length_aa':len(seq)})
for g in groups:
 if set(groups[g])!=accessions[g]:raise SystemExit('Missing expected reference accession(s) for '+g)
with review_path.open(encoding='utf-8-sig',newline='') as f:
 reader=csv.DictReader(f,delimiter='\t');review_fields=reader.fieldnames or [];review_rows=list(reader)
if 'protein_id' not in review_fields:raise SystemExit('Missing review protein_id column')
old={r['protein_id']:r for r in review_rows}
if len(old)!=len(review_rows):raise SystemExit('Duplicate review identifier')
seqs=fasta(proteome)
if not set(old).issubset(seqs):raise SystemExit('Some reviewed proteins are absent from the complete proteome')
stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
out=root/'08_focused_families/validation/step25kz17.7'/('release_'+os.environ.get('SLURM_JOB_ID','manual')+'_'+stamp)
out.mkdir(parents=True,exist_ok=False)
tsv(out/'01_reference_family_audit.tsv',['accession','reference_id','family','length_aa'],ref_records)
version={}
for name,args in [('mafft',['mafft','--version']),('hmmbuild',['hmmbuild','-h']),('hmmsearch',['hmmsearch','-h'])]:
 p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=True)
 version[name]=p.stdout.splitlines()[:5]
all_hits={};all_domains={};coverage_rows=[]
for g in ['BETA_CA','ALPHA_CA']:
 d=out/g;d.mkdir()
 write_fasta(d/'references.faa',[(pid,s) for acc,(pid,s) in sorted(groups[g].items())])
 run(['mafft','--thread',str(threads),'--auto',str(d/'references.faa')],d/'reference_alignment.faa',d/'mafft.log')
 run(['hmmbuild','--cpu',str(threads),'--amino','-n',g,str(d/'profile.hmm'),str(d/'reference_alignment.faa')],d/'hmmbuild.log',d/'hmmbuild.err')
 # Broad reporting retains weak hits for inspection; significance is evaluated below.
 run(['hmmsearch','--cpu',str(threads),'--noali','-E','10','--domE','10',
      '--incE',str(E_CUTOFF),'--incdomE',str(E_CUTOFF),
      '--tblout',str(d/'sequence_hits.tbl'),'--domtblout',str(d/'domain_hits.domtbl'),
      str(d/'profile.hmm'),str(proteome)],d/'search.txt',d/'search.err')
 hits={}
 for line in (d/'sequence_hits.tbl').read_text().splitlines():
  if not line.strip() or line.startswith('#'):continue
  t=line.split()
  if len(t)<18 or t[2]!=g:raise SystemExit('Unexpected sequence table format')
  if t[0] not in seqs or t[0] in hits:raise SystemExit('Unknown or repeated HMMER target')
  ev,score=float(t[4]),float(t[5])
  if not math.isfinite(ev) or ev<0 or not math.isfinite(score):raise SystemExit('Invalid HMMER value')
  hits[t[0]]={'evalue':ev,'score':score}
 domains={}
 for line in (d/'domain_hits.domtbl').read_text().splitlines():
  if not line.strip() or line.startswith('#'):continue
  t=line.split()
  if len(t)<22 or t[3]!=g:raise SystemExit('Unexpected domain table format')
  pid=t[0]
  if pid not in seqs:raise SystemExit('Unknown domain target')
  mlen=int(t[5]);ie=float(t[12]);start,end=int(t[15]),int(t[16]);a,b=int(t[17]),int(t[18])
  if int(t[2])!=len(seqs[pid]) or not (1<=start<=end<=mlen) or not (1<=a<=b<=len(seqs[pid])):raise SystemExit('Invalid domain coordinates')
  row={'protein_id':pid,'family':g,'independent_evalue':ie,'domain_score':float(t[13]),'model_length':mlen,
       'hmm_from':start,'hmm_to':end,'alignment_from':a,'alignment_to':b,
       'model_coverage':(end-start+1)/mlen}
  domains.setdefault(pid,[]).append(row);coverage_rows.append(row)
 all_hits[g]=hits;all_domains[g]=domains

def evidence(g,pid):
 h=all_hits[g].get(pid);ds=all_domains[g].get(pid,[])
 significant=[d for d in ds if d['independent_evalue']<=E_CUTOFF]
 maxcov=max((d['model_coverage'] for d in significant),default=0)
 # Coverage is per single domain; separated fragments/repeats cannot sum into a pass.
 sig=h is not None and h['evalue']<=E_CUTOFF and bool(significant)
 return h,ds,maxcov,sig,sig and maxcov>=COVERAGE

# Union includes every old candidate plus newly significant alpha/beta hits.
ids=set(old)
for g in all_hits:
 ids.update(pid for pid in all_hits[g] if evidence(g,pid)[3])
rows=[];beta_provisional=[];counts=Counter();species=Counter()
for pid in sorted(ids):
 bh,bd,bc,bs,bfull=evidence('BETA_CA',pid)
 ah,ad,ac,ass,afull=evidence('ALPHA_CA',pid)
 if bs and ass:classification='AMBIGUOUS_BOTH_PROFILES'
 elif bs:classification='BETA_SUPPORTED_REVIEW' if bfull else 'BETA_PARTIAL_REVIEW'
 elif ass:classification='ALPHA_SUPPORTED_REVIEW' if afull else 'ALPHA_PARTIAL_REVIEW'
 else:classification='NO_SIGNIFICANT_FAMILY_SUPPORT'
 flags=[]
 seq=seqs[pid]
 if '*' in seq:flags.append('STOP_SYMBOL')
 if '-' in seq or '.' in seq:flags.append('GAP_SYMBOL')
 if 'X' in seq:flags.append('AMBIGUOUS_X')
 if any(c not in 'ACDEFGHIKLMNPQRSTVWYBXZJUO*-.' for c in seq):flags.append('UNEXPECTED_SYMBOL')
 prior=old.get(pid,{})
 row={'protein_id':pid,'species':pid.split('__')[0],'length_aa':len(seq),
      'previous_candidate':'YES' if pid in old else 'NO','screen_class':classification,
      'sequence_flags':';'.join(flags) or 'NONE'}
 for g,h,ds,cov in [('beta',bh,bd,bc),('alpha',ah,ad,ac)]:
  row[g+'_full_evalue']=format(h['evalue'],'.8g') if h else 'NA'
  row[g+'_full_score']=format(h['score'],'.8g') if h else 'NA'
  row[g+'_best_significant_domain_model_coverage']=format(cov,'.6f')
  row[g+'_best_domain_i_evalue']=format(min(d['independent_evalue'] for d in ds),'.8g') if ds else 'NA'
 for k in ['reviewer_decision','catalytic_region_check','truncation_check','reviewer_reason']:
  row['prior_'+k]=prior.get(k,'NOT_PREVIOUSLY_REVIEWED')
 # Family-specific evidence is new; previous approvals are history, not fresh validation.
 row['new_reviewer_decision']='PENDING';row['new_catalytic_region_check']='PENDING';row['new_truncation_check']='PENDING';row['new_reviewer_reason']=''
 row['formal_status']='REVIEW_ONLY_NOT_EXPORTED'
 rows.append(row);counts[classification]+=1;species[(row['species'],classification)]+=1
 if classification=='BETA_SUPPORTED_REVIEW' and not flags:beta_provisional.append(pid)
if not rows:raise SystemExit('No records recovered')
tsv(out/'02_CA_family_specific_review.tsv',list(rows[0]),rows)
tsv(out/'03_domain_evidence.tsv',['protein_id','family','independent_evalue','domain_score','model_length','hmm_from','hmm_to','alignment_from','alignment_to','model_coverage'],coverage_rows)
tsv(out/'04_species_screen_counts.tsv',['species','screen_class','candidate_records'],
 [{'species':s,'screen_class':c,'candidate_records':n} for (s,c),n in sorted(species.items())])
write_fasta(out/'05_all_screened_candidates.faa',[(pid,seqs[pid]) for pid in sorted(ids)])
write_fasta(out/'06_beta_provisional_candidates.faa',[(pid,seqs[pid]) for pid in beta_provisional])
write_fasta(out/'07_beta_provisional_with_references.faa',
 [(pid,seqs[pid]) for pid in beta_provisional]+ [('REF__'+pid,s) for acc,(pid,s) in sorted(groups['BETA_CA'].items())])
run(['mafft','--thread',str(threads),'--auto',str(out/'07_beta_provisional_with_references.faa')],out/'08_beta_only_review_alignment.faa',out/'08_mafft.log')
aligned=fasta(out/'08_beta_only_review_alignment.faa')
expected=set(beta_provisional)|{'REF__'+pid for pid,s in groups['BETA_CA'].values()}
if set(aligned)!=expected or len(set(map(len,aligned.values())))!=1:raise SystemExit('Beta-only alignment QC failed')
for pid,seq in aligned.items():
 original=seqs[pid] if pid in seqs else refs[pid.removeprefix('REF__')]
 if seq.replace('-','')!=original:raise SystemExit('Alignment altered a sequence: '+pid)
qc={'status':'COMPLETE_FAMILY_SPECIFIC_SCREEN_REVIEW_ONLY','proteome_sequences':len(seqs),
 'previous_candidates':len(old),'screened_review_records':len(rows),'new_significant_candidates':len(ids-set(old)),
 'reference_counts':{g:len(v) for g,v in groups.items()},'classification_counts':dict(counts),
 'beta_provisional_alignment_candidates':len(beta_provisional),'alignment_records':len(aligned),
 'alignment_columns':len(next(iter(aligned.values()))),
 'screen_thresholds':{'sequence_evalue':E_CUTOFF,'independent_domain_evalue':E_CUTOFF,'single_domain_model_coverage':COVERAGE},
 'policies':[
 'Thresholds are provisional screening rules, not calibrated family-specific acceptance cutoffs.',
 'No alpha/beta score subtraction: independent profile support is evaluated separately.',
 'No reported hit means not reported under this search configuration, not proven biological absence.',
 'Partial and dual-profile hits remain unresolved. No automatic exclusion or formal acceptance.',
 'Previous review selections are preserved separately. New family-specific checks start PENDING.',
 'Beta provisional FASTA/alignment are review inputs, not formal gene counts or tree inputs.',
 'No catalytic residue pass is inferred from an E-value, length, or domain coverage.',
 'No other project module or previous release was changed.'],
 'input_sha256':{'proteome':sha(proteome),'references':sha(refs_path),'prior_review':sha(review_path)},
 'software':version}
(out/'SCREEN_QC.json').write_text(json.dumps(qc,indent=2)+'\n')
# A small review bundle includes every prior candidate, including alpha/weak rows.
bundle=out/'review_bundle';bundle.mkdir()
for name in ['01_reference_family_audit.tsv','02_CA_family_specific_review.tsv','03_domain_evidence.tsv','04_species_screen_counts.tsv','08_beta_only_review_alignment.faa','SCREEN_QC.json']:
 shutil.copyfile(out/name,bundle/name)
(bundle/'COMPLETE.txt').write_text('Family-specific screen complete; formal acceptance pending.\n')
with (bundle/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(bundle.iterdir()):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+p.name+'\n')
with tarfile.open(out/'CA_specific_review_bundle.tar.gz','w:gz') as t:t.add(bundle,arcname='CA_specific_review')
(out/'COMPLETE.txt').write_text('Family-specific screen complete. No formal proteins exported.\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p!=out/'SHA256SUMS.txt':f.write(sha(p)+'  '+p.relative_to(out).as_posix()+'\n')
print('Status:',qc['status']);print('Screen classes:',dict(counts))
print('Output: <PROJECT_ROOT>/'+out.relative_to(root).as_posix())
print('Upload CA_specific_review_bundle.tar.gz for review.')
PY
