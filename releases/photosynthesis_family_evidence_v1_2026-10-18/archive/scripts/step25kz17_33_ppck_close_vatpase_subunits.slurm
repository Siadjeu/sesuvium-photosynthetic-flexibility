#!/usr/bin/env bash
#SBATCH --job-name=vatp_subunits
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=vatp_subunits_%j.out
#SBATCH --error=vatp_subunits_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
for tool in python3 diamond hmmsearch; do command -v "$tool"; done
python3 - <<'PY'
import gzip,csv,hashlib,json,os,re,subprocess,tarfile,shutil
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def release(step,env):
 paths=[root/os.environ[env]] if os.environ.get(env) else [p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
 if len(paths)!=1:raise SystemExit('Expected one release; set '+env+' to the intended release directory')
 p=paths[0].resolve(strict=True)
 for line in (p/'SHA256SUMS.txt').read_text().splitlines():
  h,n=line.split(maxsplit=1);f=(p/n).resolve()
  if not f.is_relative_to(p) or sha(f)!=h:raise SystemExit('Checksum mismatch: '+n)
 return p
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta(p):
 d={};key=None
 for line in p.read_text().splitlines():
  if line.startswith('>'):
   key=line[1:].split()[0]
   if key in d:raise SystemExit('Duplicate FASTA ID: '+key)
   d[key]=''
  elif line.strip():
   if key is None:raise SystemExit('Invalid FASTA')
   d[key]+=line.strip()
 return d
def writefa(p,d):p.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
closed=release('step25kz17.32','PPCK_LEDGER_RELEASE');screen=release('step25kz17.22','SCREEN_RELEASE')
if json.loads((closed/'PROVENANCE.json').read_text())['screen_manifest_sha256']!=sha(screen/'SHA256SUMS.txt'):raise SystemExit('PPCK/screen provenance mismatch')
ppck=read(closed/'01_PPCK_evidence_ledger.tsv');evidence=read(closed/'09_all_candidates_completeness_evidence.tsv')
assert len(ppck)==27 and len({(r['species'],r['gene_id']) for r in ppck})==27
assert len(evidence)==5121 and sum(r['status']=='REVIEW_TRIGGERED' for r in evidence)==762
assert not any('PROFILE_WITHIN' in r['triggers'] for r in evidence if r['status']!='ALREADY_SUPPORTED')
meta=read(screen/'02_reference_metadata.tsv');refmeta={r['accession']:r for r in meta}
vrefs={r['accession']:r['subunit'] for r in meta if r['assigned_category']=='V_ATPASE'}
subunits=['A','B','C','D','E','F','G','H','a','c',"c''",'d','e']
assert len(vrefs)==28 and set(vrefs.values())==set(subunits)
rows=[r for r in read(screen/'09_specificity_evidence_REVIEW_ONLY.tsv') if r['category']=='V_ATPASE']
ids={r['protein_id'] for r in rows};assert len(ids)==603 and len(rows)==970
seqs=fasta(screen/'05_candidate_sequences.faa');refs=fasta(screen/'references/target_and_background.faa')
assert set(vrefs)<=refs.keys() and set(refs)<=refmeta.keys()
# Locate the installed complete Pfam-A database; never choose between versions silently.
if os.environ.get('PFAM_HMM'):db=root/os.environ['PFAM_HMM']
else:
 paths=set();homes=[root/'databases',root/'00_databases',root.parent/'my_interproscan']
 if os.environ.get('INTERPROSCAN_HOME'):homes.append(Path(os.environ['INTERPROSCAN_HOME']))
 exe=shutil.which('interproscan.sh')
 if exe:homes.append(Path(exe).resolve().parent)
 for home in homes:
  if not home.is_dir():continue
  for folder,dirs,files in os.walk(home):
   dirs[:]=[d for d in dirs if not d.startswith('.') and d not in ['temp','tmp','work']]
   for n in files:
    if n.lower() in ['pfam-a.hmm','pfam_a.hmm','pfam-a.hmm.gz','pfam_a.hmm.gz']:paths.add((Path(folder)/n).resolve())
 if len(paths)!=1:raise SystemExit('Set PFAM_HMM to the intended installed Pfam-A ASCII .hmm or .hmm.gz file')
 db=next(iter(paths))
assert db.is_file()
out=base/'step25kz17.33'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
cpu=os.environ.get('SLURM_CPUS_PER_TASK','4');commands=[]
def run(args,log):
 commands.append(args)
 with (out/log).open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
# Closure records an agreed bounded endpoint, not exhaustive completeness or function.
for r in ppck:r['review_stage']='BOUNDED_REVIEW_CLOSED_WITH_LIMITATIONS'
table('01_PPCK_closed_ledger.tsv',ppck,list(ppck[0]))
shutil.copy2(closed/'02_PPCK_supported_counts.tsv',out/'02_PPCK_supported_counts.tsv')
(out/'03_PPCK_closure.txt').write_text('27 phylogenetically supported distinct loci retained. 5094 additional screened candidates checked. No additional candidate met the profile discovery threshold. 762 similarity-only flags retained as unresolved screening evidence; none promoted or categorically excluded. Profile training included the supported set; this was not independent validation. Whole-proteome/genome completeness, subtype assignments and biochemical function are not established. Bounded review closed by project decision.\n')
shutil.copy2(closed/'10_remaining_review_queue.tsv',out/'03_PPCK_unresolved_screening_flags.tsv')
table('04_VATPASE_original_screen_rows.tsv',rows,list(rows[0]))
table('05_reference_metadata.tsv',meta,list(meta[0]))
query={p:seqs[p] for p in sorted(ids)}
query.update({'REF_'+a:s for a,s in refs.items()})
writefa(out/'06_candidates_and_controls.faa',query)
writefa(out/'06_candidates.faa',{p:seqs[p] for p in sorted(ids)})
writefa(out/'06_reference_panel.faa',refs)
dbuse=db
if db.suffix=='.gz':
 dbuse=out/'work/Pfam-A.hmm'
 with gzip.open(db,'rb') as a,dbuse.open('wb') as b:shutil.copyfileobj(a,b)
models=ga=acc=0
with dbuse.open() as f:
 for l in f:
  if l.startswith('HMMER3/'):models+=1
  elif l.startswith('GA '):ga+=1
  elif l.startswith('ACC ') and 'PF' in l:acc+=1
assert models>1000 and models==ga==acc,'Expected complete Pfam-A with GA thresholds'
run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'07_pfam.domtbl'),str(dbuse),str(out/'06_candidates_and_controls.faa')],'07_hmmsearch.log')
domains=[];byq=defaultdict(list)
for l in (out/'07_pfam.domtbl').read_text().splitlines():
 if not l or l.startswith('#'):continue
 x=l.split(maxsplit=22)
 r=dict(protein_id=x[0],model=x[3],accession=x[4].split('.')[0],model_length=int(x[5]),independent_evalue=x[12],bitscore=float(x[13]),hmm_start=int(x[15]),hmm_end=int(x[16]),protein_start=int(x[17]),protein_end=int(x[18]),model_span=round((int(x[16])-int(x[15])+1)/int(x[5]),4))
 domains.append(r);byq[x[0]].append(r)
dfields=['protein_id','model','accession','model_length','independent_evalue','bitscore','hmm_start','hmm_end','protein_start','protein_end','model_span']
table('08_domain_evidence.tsv',domains,dfields)
table('09_reference_domain_baselines.tsv',[r for r in domains if r['protein_id'].startswith('REF_')],dfields)
run(['diamond','makedb','--threads',cpu,'--in',str(out/'06_reference_panel.faa'),'--db',str(out/'work/panel')],'10_makedb.log')
run(['diamond','blastp','--threads',cpu,'--very-sensitive','--query',str(out/'06_candidates.faa'),'--db',str(out/'work/panel'),'--evalue','1e-3','--max-target-seqs','0','--max-hsps','1','--outfmt','6','qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen','--out',str(out/'10_reference_hits.tsv')],'10_diamond.log')
hits=defaultdict(list)
for l in (out/'10_reference_hits.tsv').read_text().splitlines():
 x=l.split('\t');hits[x[0]].append(dict(ref=x[1],score=float(x[9]),evalue=float(x[8]),identity=float(x[2]),qcov=(abs(int(x[5])-int(x[4]))+1)/int(x[10]),rcov=(abs(int(x[7])-int(x[6]))+1)/int(x[11])))
# Exact mapping of annotated aliases; no suffix guessing and no sequence deduplication.
mapping={};maphash={}
for sp in sorted({p.split('__',1)[0] for p in ids}):
 path=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv')
 rr=read(path) if path.is_file() else [];maphash[sp]=sha(path) if path.is_file() else 'MISSING'
 aliases=defaultdict(set)
 for r in rr:
  for k in ['gene_id','transcript_id','protein_id','cds_id']:
   if r.get(k) and r.get('gene_id'):aliases[r[k]].add(r['gene_id'])
 for p in ids:
  if p.split('__',1)[0]==sp:mapping[p]=aliases[p.split('__',1)[1]]
# Reference domain signatures are descriptive. Shared ATPase domains cannot establish specificity.
models_by_ref={a:{d['accession'] for d in byq['REF_'+a]} for a in refs}
comparisons=[];ledger=[]
for p in sorted(ids):
 hs=sorted(hits[p],key=lambda h:(-h['score'],h['ref']));vh=[h for h in hs if h['ref'] in vrefs]
 best=vh[0] if vh else None;sub=vrefs[best['ref']] if best else ''
 candidate_models={d['accession'] for d in byq[p]}
 group_best={s:next((h for h in vh if vrefs[h['ref']]==s),None) for s in subunits}
 for s in subunits:
  h=group_best[s]
  if not h:continue
  baseline=models_by_ref[h['ref']];alternative=[z for z in hs if z['ref'] not in vrefs or vrefs[z['ref']]!=s]
  competitor=alternative[0] if alternative else None
  comparisons.append(dict(protein_id=p,subunit=s,reference=h['ref'],bitscore=h['score'],evalue=h['evalue'],identity=h['identity'],reference_span=round(h['rcov'],4),protein_span=round(h['qcov'],4),best_other_reference=competitor['ref'] if competitor else '',other_category=refmeta[competitor['ref']]['assigned_category'] if competitor else '',other_subunit=vrefs.get(competitor['ref'],'') if competitor else '',score_margin=round(h['score']-competitor['score'],2) if competitor else '',reference_models=';'.join(sorted(baseline)),candidate_models=';'.join(sorted(candidate_models)),missing_reference_models=';'.join(sorted(baseline-candidate_models)),additional_models=';'.join(sorted(candidate_models-baseline))))
 flags=[];alt=next((h for h in hs if h['ref'] not in vrefs or vrefs[h['ref']]!=sub),None) if best else None
 if not best:flags.append('NO_VATPASE_REFERENCE_HIT')
 else:
  baseline=models_by_ref[best['ref']]
  if not baseline:flags.append('REFERENCE_WITHOUT_GA_DOMAIN')
  if baseline-candidate_models:flags.append('MISSING_REFERENCE_DOMAIN_MODELS')
  if candidate_models-baseline:flags.append('ADDITIONAL_DOMAIN_MODELS')
  if best['rcov']<0.8 or best['qcov']<0.8:flags.append('PARTIAL_OR_LENGTH_MISMATCH')
  if alt and best['score']-alt['score']<20:flags.append('SUBUNIT_OR_COMPETITOR_AMBIGUITY')
  if not alt:flags.append('NO_COMPETITOR_REPORTED_IN_LIMITED_PANEL')
 if len(mapping[p])!=1:flags.append('LOCUS_MAPPING_UNRESOLVED')
 ledger.append(dict(protein_id=p,species=p.split('__',1)[0],gene_id=';'.join(sorted(mapping[p])),mapping_status='UNAMBIGUOUS' if len(mapping[p])==1 else 'UNRESOLVED',length_aa=len(seqs[p]),provisional_subunit=sub,best_VATPASE_reference=best['ref'] if best else '',best_other_reference=alt['ref'] if alt else '',score_margin=round(best['score']-alt['score'],2) if best and alt else '',flags=';'.join(flags),review_status='PRIORITY_REVIEW' if flags else 'CONCORDANT_SCREEN_DOMAIN_EVIDENCE_REVIEW_REQUIRED',final_decision='PENDING'))
# Explicitly mark multiple protein records mapping to one locus before any counting.
locus_records=Counter((r['species'],r['gene_id']) for r in ledger if r['mapping_status']=='UNAMBIGUOUS')
for r in ledger:
 if r['mapping_status']=='UNAMBIGUOUS' and locus_records[r['species'],r['gene_id']]>1:
  r['flags']+=';MULTIPLE_PROTEIN_RECORDS_AT_LOCUS';r['flags']=r['flags'].strip(';');r['review_status']='PRIORITY_REVIEW'
table('11_subunit_specific_evidence.tsv',comparisons,list(comparisons[0]))
table('12_VATPASE_candidate_ledger_REVIEW_ONLY.tsv',ledger,list(ledger[0]))
table('13_VATPASE_priority_review.tsv',[r for r in ledger if r['flags']],list(ledger[0]))
# Subunit-safe filenames retain case distinction even on case-insensitive systems.
(out/'subunits').mkdir()
for s in subunits:
 label=('V1_'+s) if s in list('ABCDEFGH') else ('V0_c_doubleprime' if s=="c''" else 'V0_'+s)
 group={r['protein_id']:seqs[r['protein_id']] for r in ledger if r['provisional_subunit']==s}
 group.update({'REF_'+a:refs[a] for a in vrefs if vrefs[a]==s})
 writefa(out/'subunits'/(label+'_PROVISIONAL.faa'),group)
summ=[]
for sp in sorted({r['species'] for r in ledger}):
 for s in subunits:
  rr=[r for r in ledger if r['species']==sp and r['provisional_subunit']==s]
  summ.append(dict(species=sp,provisional_subunit=s,candidate_records=len(rr),unambiguously_mapped_loci=len({r['gene_id'] for r in rr if r['mapping_status']=='UNAMBIGUOUS'}),priority_review_records=sum(bool(r['flags']) for r in rr),interpretation='SCREENING_COUNTS_NOT_ACCEPTED_GENE_COUNTS'))
table('14_subunit_screening_summary_NOT_FINAL.tsv',summ,list(summ[0]))
(out/'DATABASE_PROVENANCE.json').write_text(json.dumps(dict(source_sha256=sha(db),uncompressed_sha256=sha(dbuse),pfam_models=models,method='hmmsearch --cut_ga'),indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(ppck_manifest_sha256=sha(closed/'SHA256SUMS.txt'),screen_manifest_sha256=sha(screen/'SHA256SUMS.txt'),mapping_hashes=maphash),indent=2)+'\n')
(out/'QC.json').write_text(json.dumps(dict(status='PPCK_BOUNDED_REVIEW_CLOSED_VATPASE_REVIEW_REQUIRED',ppck_supported_loci=27,vatpase_candidates=len(ids),original_candidate_subunit_rows=len(rows),vatpase_references=len(vrefs),subunits=subunits,all_panel_references=len(refs),domain_hits=len(domains),references_without_GA_domains=[a for a in vrefs if not models_by_ref[a]],unresolved_locus_mappings=sum(r['mapping_status']!='UNAMBIGUOUS' for r in ledger),review_status_counts=dict(Counter(r['review_status'] for r in ledger)),final_VATPASE_counts=False),indent=2)+'\n')
(out/'README.txt').write_text('PPCK bounded review closed; limitations and 762 unresolved screening flags retained.\nV-ATPase validation covers 603 corrected-screen candidates and 28 V-type references across 13 case-sensitive subunits; all existing panel references are included as competitors and domain controls.\nCandidate/subunit evidence rows are not genes. Each candidate has one provisional best subunit; this is not an accepted assignment. Distinct mapped loci remain distinct even with identical sequences.\nPfam domain evidence uses installed GA thresholds. Overlapping and clan-related models are retained; do not interpret overlaps as extra physical domains. Reference domain sets are empirical baselines, not definitive required architectures. Missing domains and extra models are review flags, not exclusions.\nReference/protein spans are single-best-HSP spans, not summed coverage. The 80 percent span and 20-bit margin rules only prioritize review. No-competitor cases remain flagged.\nV-type/F-type related domains require specificity review. Subunits are evaluated separately, not aligned together; provisional FASTAs are provided for focused follow-up if needed. No new phylogeny in this step.\nScreening counts are not accepted V-ATPase totals. Whole-genome completeness, subcellular localization, physiological activity and C4/CAM specificity remain untested.\n')
for tool,args in [('diamond',['version']),('hmmsearch',['-h'])]:run([tool]+args,tool+'_version.txt')
(out/'COMMANDS.json').write_text(json.dumps([[v.replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for v in c] for c in commands],indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='PPCK_closed_VATPASE_subunits/'+str(p.relative_to(out)))
print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
