#!/usr/bin/env bash
#SBATCH --job-name=vppase_domains
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=vppase_domains_%j.out
#SBATCH --error=vppase_domains_%j.err
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
out=base/'step25kz17.27'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rr)
def writefa(p,d):
 with p.open('w') as f:
  for h,s in d.items():f.write('>'+h+'\n'+s+'\n')
# Consolidate prior PPDK decisions, preserving unresolved loci and mapping uncertainty.
def choose(step,env):
 candidates=[root/os.environ[env]] if os.environ.get(env) else [p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
 if len(candidates)!=1:raise SystemExit('Set '+env+' to the intended '+step+' release')
 p=candidates[0].resolve(strict=True)
 for line in (p/'SHA256SUMS.txt').read_text().splitlines():
  digest,name=line.split(maxsplit=1);f=(p/name).resolve()
  if not f.is_relative_to(p) or sha(f)!=digest:raise SystemExit('Checkpoint checksum failure: '+name)
 return p
domain=choose('step25kz17.23','PPDK_DOMAIN_RELEASE');audit=choose('step25kz17.26','AUDIT_RELEASE')
if json.loads((audit/'QC.json').read_text()).get('parser_discrepancies'):raise SystemExit('Resolve audit parser discrepancies first')
prior=[r for r in read(domain/'01_sequence_map.tsv') if r['role']=='CANDIDATE']
if Counter(r['family'] for r in prior)!=Counter({'PPDK':15,'PPDK_RP':17}):raise SystemExit('Unexpected PPDK checkpoint candidate set')
flagged={'MC__gene_22054':'EXTENDED_ALIGNMENT_DISRUPTIONS','SP_Tet__g10591':'N_TERMINAL_COMPLETENESS_UNRESOLVED','SP_Hex__gene_Spor0461320.1':'EXON_STRUCTURE_AND_DISRUPTIONS_UNRESOLVED'}
old_screen={r['protein_id']:r for r in read(src/'11_priority_PPDK_RP_review.tsv')}
completeness=read(domain/'07_alignment_completeness.tsv');ledger=[];mapping_hashes={}
# Build exact mapping aliases only; never strip transcript suffixes to guess a gene.
needed=defaultdict(set)
for r in prior:
 if r['protein_id'] not in old_screen:raise SystemExit('Prior candidate absent from corrected screen')
 sp=old_screen[r['protein_id']]['species'];needed[sp].add(r['protein_id'])
resolved={}
for sp,ids in needed.items():
 path=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv')
 if not path.is_file():
  for pid in ids:resolved[pid]=set()
  continue
 mapping_hashes[str(path.relative_to(root))]=sha(path)
 for pid in ids:resolved[pid]=set()
 with path.open() as f:
  reader=csv.DictReader(f,delimiter='\t')
  if not {'gene_id','transcript_id','protein_id','cds_id'}.issubset(reader.fieldnames or []):raise SystemExit('Unexpected isoform mapping schema '+sp)
  for rr in reader:
   aliases={sp+'__'+rr[k] for k in ['gene_id','transcript_id','protein_id','cds_id'] if rr.get(k) and rr[k] not in ['NA','.']}
   for pid in ids & aliases:resolved[pid].add(rr['gene_id'])
for r in prior:
 pid=r['protein_id'];sp=old_screen[pid]['species'];genes=resolved[pid]
 cp=[x for x in completeness if x['query_id']==r['query_id']]
 status='UNRESOLVED_MODEL' if pid in flagged else 'SUPPORTED_MODEL'
 ledger.append(dict(protein_id=pid,species=sp,family=r['family'],wang_family='PPDK-RE' if r['family']=='PPDK_RP' else 'PPDK',length_aa=r['length_aa'],decision=status,reason=flagged.get(pid,'DOMAIN_AND_REFERENCE_ALIGNMENT_SUPPORT'),gene_id=next(iter(genes)) if len(genes)==1 else ';'.join(sorted(genes)),locus_mapping='UNAMBIGUOUS' if len(genes)==1 else 'UNRESOLVED_OR_AMBIGUOUS',min_reference_coverage=min(float(x['reference_position_coverage']) for x in cp),max_missing_N_reference_positions=max(int(x['missing_N_reference_positions']) for x in cp),read_validation='UNAVAILABLE_FOR_THREE_FLAGGED_LOCI' if pid in flagged else 'NOT_PERFORMED',final_count_status='PROVISIONAL_NOT_FINAL'))
table('09_PPDK_PPDK_RP_decision_ledger.tsv',ledger,list(ledger[0]))
counts=[]
for sp in sorted(needed):
 for fam in ['PPDK','PPDK_RP']:
  rr=[r for r in ledger if r['species']==sp and r['family']==fam]
  supported={r['gene_id'] for r in rr if r['decision']=='SUPPORTED_MODEL' and r['locus_mapping']=='UNAMBIGUOUS'}
  unresolved={r['gene_id'] for r in rr if r['decision']=='UNRESOLVED_MODEL' and r['locus_mapping']=='UNAMBIGUOUS'}
  counts.append(dict(species=sp,family=fam,candidate_records=len(rr),supported_distinct_mapped_loci=len(supported),unresolved_additional_mapped_loci=len(unresolved-supported),conflicting_status_loci=len(unresolved&supported),unresolved_mapping_records=sum(r['locus_mapping']!='UNAMBIGUOUS' for r in rr),count_status='PROVISIONAL_NOT_FINAL'))
table('10_PPDK_PPDK_RP_counts_PROVISIONAL.tsv',counts,list(counts[0]))
(out/'11_limitations.txt').write_text('Three PPDK regulatory-protein candidate loci showed incomplete coding models or alignment-implied coding disruptions. Suitable independent read evidence was reported unavailable for resolving these cases. Existing annotations were retained and these loci were reported separately as unresolved candidates, not confirmed pseudogenes or complete coding models. Supported-model counts do not imply experimentally validated function.\n')
(out/'LEDGER_PROVENANCE.json').write_text(json.dumps({'domain_release':os.path.relpath(domain,root),'domain_manifest_sha256':sha(domain/'SHA256SUMS.txt'),'audit_release':os.path.relpath(audit,root),'audit_manifest_sha256':sha(audit/'SHA256SUMS.txt'),'screen_manifest_sha256':sha(src/'SHA256SUMS.txt'),'mapping_sha256':mapping_hashes},indent=2)+'\n')

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
rows=[r for r in read(src/'09_specificity_evidence_REVIEW_ONLY.tsv') if r['category']=='V_PPASE']
if len(rows)!=48 or len({r['protein_id'] for r in rows})!=48:raise SystemExit('Expected 48 unique V_PPASE candidates; inspect screen release')
seqs=fasta(src/'05_candidate_sequences.faa');refs=fasta(src/'references/targets.faa')
refgroups={'V_PPASE':['P31414','Q56ZN6','Q9FWR2']}
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
(out/'README.txt').write_text("""PPDK/PPDK-RP evidence ledger and V-PPase domain/completeness review.
The ledger records the agreed evidence-based disposition, not experimentally proven function.
Three PPDK-RP loci remain unresolved; no pseudogene designation or sequence repair is made.
Suitable independent reads were reported unavailable by the investigator; not independently inventoried.
Only unambiguous mapping to a distinct gene contributes to provisional mapped-locus summaries.
V-PPase candidates remain PENDING. Pfam hits and alignment occupancy are descriptive.
Raw overlapping domain hits are retained. No automatic coverage cutoff or final gene acceptance.
Reference controls provide empirical domain architecture; missing domains require review.
C4/CAM functional recruitment and complete Wang-family coverage are not established by this step.
All existing releases and annotation files remain read-only.
""")
(out/'QC.json').write_text(json.dumps(dict(status='DOMAIN_REVIEW_COMPLETE_NOT_FINAL',candidates=len(rows),references=3,
 pfam_models=models,domain_hits=len(domains),reference_controls_without_hits=[q for q in query if q.startswith('REF_') and not byq[q]],final_counts=False),indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='PPDK_ledger_VPPASE_validation/'+str(p.relative_to(out)))
print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
