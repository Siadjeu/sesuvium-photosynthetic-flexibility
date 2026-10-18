#!/usr/bin/env bash
#SBATCH --job-name=ppt_combined
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=ppt_combined_%j.out
#SBATCH --error=ppt_combined_%j.err
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






src=release('step25kz17.65','PPT_PROVENANCE_RELEASE');cor=release('step25kz17.21','CORRECTED_INPUT_RELEASE')
prov=json.loads((src/'PROVENANCE.json').read_text());proteome=cor/'full_outputs/corrected_protein_candidates.faa'
need(sha(cor/'SHA256SUMS.txt')==prov['corrected_manifest_sha256'],'Corrected release changed')
need(sha(proteome)==prov['corrected_proteome_sha256'],'Corrected proteome changed')
refs,rheaders=fasta(src/'07_reference_panel.faa');panel=read(src/'05_reference_panel.tsv');roles={r['reference_id']:r for r in panel}
need(len(refs)==6 and {r['gene'] for r in panel}=={'PPT1','PPT2','TPT','GPT1','GPT2','XPT'},'Incomplete six-reference panel')
for r in panel:need(hashlib.sha256(refs[r['reference_id']].encode()).hexdigest()==r['sequence_sha256'],'Reference sequence changed')
inv=read(src/'01_project_PPT_inventory.tsv')
def evidence(suffix):
 rr=[r for r in inv if r['source_file'].endswith(suffix)];need(len(rr)==1,'Ambiguous or missing historical evidence '+suffix);r=rr[0]
 p=src/r['bundle_file'] if r['bundle_file'] else root/r['source_file']
 need(p.is_file() and sha(p)==r['sha256'],'Historical evidence changed '+suffix);return p
whole='/PPT_AMYLASE_whole_proteome_validation_primary_v1/'
refined='/PPT_AMYLASE_refined_domain_supported_primary_v1/'
hist=read(evidence(whole+'PPT.whole_proteome_classification.tsv'));refinedrows=read(evidence(refined+'PPT.refined_candidates.tsv'))
oldfa,_=fasta(evidence(whole+'PPT.retained_plus_references.fa'));refinedfa,_=fasta(evidence(refined+'PPT.refined_plus_references.fa'))
ogfa,_=fasta(evidence('/Orthogroup_Sequences/OG0000369.fa'))
need(len(hist)==283 and len(ogfa)==37,'Historical scope changed')
need(len({r['candidate_id'] for r in hist})==len(hist),'Duplicate historical candidate IDs')
current,_=fasta(proteome);aliases=defaultdict(set);mapping_hashes={}
for sp,h in prov['mapping_hashes'].items():
 mp=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv');need(sha(mp)==h,'Mapping changed '+sp);mapping_hashes[sp]=h
 for m in read(mp):
  for field in ['gene_id','protein_id','transcript_id']:
   v=m.get(field,'')
   if v and v not in {'.','NA'}:
    aliases[v].add((sp,m['gene_id']));aliases[sp+'__'+v].add((sp,m['gene_id']))
bygene=defaultdict(set)
for k in current:
 if len(aliases.get(k,set()))==1:bygene[next(iter(aliases[k]))].add(k)
# Recover only documented aliases, without sequence-only mapping or isoform guessing.
def resolve(k):
 mm=aliases.get(k,set())
 if len(mm)!=1:return '', 'AMBIGUOUS_OR_MISSING_LOCUS_MAPPING'
 poss=bygene[next(iter(mm))]
 if k in current:return k,'EXACT_IDENTIFIER'
 if len(poss)==1:return next(iter(poss)),'UNIQUE_DOCUMENTED_LOCUS_ALIAS'
 return '', 'AMBIGUOUS_OR_MISSING_CORRECTED_PROTEIN'
out=base/'step25kz17.66'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(n,rr,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):
 (out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in sorted(d.items())))
commands=[];cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
def run(cmd,log):
 commands.append(cmd)
 with (out/log).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
refaudit=[]
for acc in refs:
 r=roles[acc];oldids=[k for k in oldfa if k.startswith('REF__') and '__'+r['gene']+'__' in k]
 refaudit.append(dict(gene=r['gene'],current_accession=acc,old_reference_ids=';'.join(oldids),comparison='IDENTICAL_SEQUENCE' if len(oldids)==1 and norm(oldfa[oldids[0]])==norm(refs[acc]) else ('CHANGED_SEQUENCE' if oldids else 'MISSING_FROM_OLD_PANEL')))
table('01_reference_panel_comparison.tsv',refaudit,['gene','current_accession','old_reference_ids','comparison'])
refinedby={r['candidate_id']:r for r in refinedrows};history={};audit=[]
for k in sorted({r['candidate_id'] for r in hist}|set(ogfa)):
 corrected,status=resolve(k);known=[norm(d[k]) for d in [oldfa,refinedfa,ogfa] if k in d];seqstatus='EXACT_HISTORICAL_SEQUENCE' if corrected and known and all(x==norm(current[corrected]) for x in known) else ('HISTORICAL_SEQUENCE_CHANGED' if corrected and known else 'NO_UNALIGNED_HISTORICAL_SEQUENCE')
 if corrected:history.setdefault(corrected,[]).append(k)
 audit.append(dict(historical_id=k,corrected_protein_id=corrected,mapping_status=status,sequence_status=seqstatus,in_OG0000369=k in ogfa,prior_assignment=refinedby.get(k,{}).get('provisional_assignment',''),prior_evidence_tier=refinedby.get(k,{}).get('evidence_tier','')))
table('02_PPT_historical_reconciliation.tsv',audit,['historical_id','corrected_protein_id','mapping_status','sequence_status','in_OG0000369','prior_assignment','prior_evidence_tier'])
# Expanded panel and changed PPT2 require a uniform search of corrected proteins.
writefa('03_six_reference_panel.faa',refs)
print('Searching corrected proteome against six references',flush=True)
run(['diamond','makedb','--threads',cpu,'--in',str(out/'03_six_reference_panel.faa'),'--db',str(out/'work/panel')],'04_makedb.log')
run(['diamond','blastp','--threads',cpu,'--query',str(proteome),'--db',str(out/'work/panel'),'--out',str(out/'04_six_reference_hits.tsv'),'--very-sensitive','--evalue','1e-3','--max-target-seqs','0','--max-hsps','1','--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore'],'04_diamond.log')
hits=defaultdict(list)
for line in (out/'04_six_reference_hits.tsv').read_text().splitlines():
 x=line.split('\t');need(len(x)==12 and x[0] in current and x[1] in refs,'Unexpected DIAMOND record')
 hits[x[0]].append(dict(reference_id=x[1],gene=roles[x[1]]['gene'],role=roles[x[1]]['role'],pident=float(x[2]),qcov=(abs(int(x[7])-int(x[6]))+1)/int(x[4]),rcov=(abs(int(x[9])-int(x[8]))+1)/int(x[5]),evalue=float(x[10]),bitscore=float(x[11])))
ids=set(hits)|set(history);need(ids,'No candidates');need(len(ids)<=3000,'More than3000 candidates: stop for scope review before domain search')
proteins={k:norm(current[k]) for k in ids};need(all(valid(s) for s in proteins.values()),'Invalid candidate sequence')
writefa('05_current_candidate_proteins.faa',proteins);writefa('05_candidates_and_references.faa',{**proteins,**refs})
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

need(sha(dbuse)=='4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e','Pfam database differs from previously used release; set PFAM_HMM to the same database')
print('Checking candidate and reference domain architectures',flush=True)
run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'06_Pfam.domtbl'),str(dbuse),str(out/'05_candidates_and_references.faa')],'06_hmmsearch.log')
domains=defaultdict(list);dh=[]
for line in (out/'06_Pfam.domtbl').read_text().splitlines():
 if not line or line.startswith('#'):continue
 x=line.split(maxsplit=22);r=dict(sequence_id=x[0],model=x[3],accession=x[4].split('.')[0],model_length=int(x[5]),i_evalue=float(x[12]),score=float(x[13]),hmm_from=int(x[15]),hmm_to=int(x[16]),ali_from=int(x[17]),ali_to=int(x[18]));dh.append(r);domains[x[0]].append(r)
table('07_Pfam_domain_hits.tsv',dh,['sequence_id','model','accession','model_length','i_evalue','score','hmm_from','hmm_to','ali_from','ali_to'])
def models(k):return {h['accession'] for h in domains[k]}
def cov(k,m):
 hh=[h for h in domains[k] if h['accession']==m];covered=set()
 for h in hh:covered.update(range(h['hmm_from'],h['hmm_to']+1))
 return len(covered)/hh[0]['model_length'] if hh else 0.0
positive=[k for k in refs if roles[k]['role']=='PPT_POSITIVE'];shared=set.intersection(*(models(k) for k in positive));need(shared,'PPT positives share no GA-supported Pfam model: inspect outputs')
reference_models=set.union(*(models(k) for k in refs));cal=[]
for k in refs:
 for m in sorted(models(k)):cal.append(dict(reference_id=k,gene=roles[k]['gene'],model=m,model_coverage=round(cov(k,m),6),shared_by_PPT_positives=m in shared))
table('08_reference_domain_calibration.tsv',cal,['reference_id','gene','model','model_coverage','shared_by_PPT_positives'])
review=[]
for k in sorted(ids):
 hh=hits[k];best=max(hh,key=lambda h:h['bitscore']) if hh else None
 pp=[h for h in hh if h['role']=='PPT_POSITIVE'];cc=[h for h in hh if h['role']=='RELATED_COMPARATOR'];p=max(pp,key=lambda h:h['bitscore']) if pp else None;c=max(cc,key=lambda h:h['bitscore']) if cc else None
 category=('NO_CURRENT_PANEL_HIT' if not hh else 'COMPARATOR_ONLY' if not p else 'PPT_ONLY' if not c else 'PPT_FAVORED' if p['bitscore']>c['bitscore'] else 'COMPARATOR_FAVORED' if p['bitscore']<c['bitscore'] else 'PPT_COMPARATOR_TIE')
 margin=abs(p['bitscore']-c['bitscore'])/max(p['bitscore'],c['bitscore']) if p and c else None
 mm=aliases.get(k,set());sp,g=next(iter(mm)) if len(mm)==1 else ('','');flags=[]
 if len(mm)!=1:flags.append('UNRESOLVED_LOCUS_MAPPING')
 missing=shared-models(k)
 if missing:flags.append('MISSING_PPT_REFERENCE_SHARED_MODEL')
 if any(cov(k,m)<.8*min(cov(x,m) for x in positive) for m in shared):flags.append('LOW_SHARED_MODEL_COVERAGE_HEURISTIC')
 if best and best['rcov']<.8:flags.append('LOW_BEST_REFERENCE_COVERAGE_HEURISTIC')
 if models(k)-reference_models:flags.append('ADDITIONAL_DOMAIN_MODEL')
 if margin is not None and margin<.1:flags.append('SMALL_SCORE_MARGIN_HEURISTIC')
 if 'X' in proteins[k]:flags.append('AMBIGUOUS_RESIDUES')
 review.append(dict(protein_id=k,species=sp,gene_id=g,length_aa=len(proteins[k]),historical_ids=';'.join(history.get(k,[])),historical_OG0000369=any(a['in_OG0000369'] and a['corrected_protein_id']==k for a in audit),new_relative_to_historical_PPT_tables=k not in history,best_reference=best['reference_id'] if best else '',best_reference_gene=best['gene'] if best else '',best_evalue=best['evalue'] if best else '',best_pident=best['pident'] if best else '',best_query_coverage=round(best['qcov'],6) if best else '',best_reference_coverage=round(best['rcov'],6) if best else '',best_PPT_bitscore=p['bitscore'] if p else '',best_comparator_bitscore=c['bitscore'] if c else '',relative_score_margin=round(margin,6) if margin is not None else '',similarity_category=category,Pfam_models=';'.join(sorted(models(k))),shared_model_coverage=';'.join(m+':'+str(round(cov(k,m),6)) for m in sorted(shared)),review_flags=';'.join(flags) or 'NO_HEURISTIC_FLAG',evidence_status='REVIEW_ONLY',final_count='NA'))
table('09_combined_PPT_specificity_review.tsv',review,list(review[0]));priority=[r for r in review if r['similarity_category'] in {'PPT_ONLY','PPT_FAVORED','PPT_COMPARATOR_TIE'} or r['historical_OG0000369'] or ('SMALL_SCORE_MARGIN_HEURISTIC' in r['review_flags'])]
table('10_priority_review_queue.tsv',priority,list(review[0]));writefa('10_priority_candidates_with_references.faa',{**{r['protein_id']:proteins[r['protein_id']] for r in priority},**refs})
# Carry immutable counts, historical classifications, and references for review.
for n in ['source_01_integrated_locus_evidence.tsv','source_07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/n)
for p,n in [(evidence(whole+'PPT.whole_proteome_classification.tsv'),'historical_PPT_whole_proteome_classification.tsv'),(evidence(refined+'PPT.refined_candidates.tsv'),'historical_PPT_refined_candidates.tsv')]:shutil.copyfile(p,out/n)
shutil.copyfile(src/'05_reference_panel.tsv',out/'reference_panel_metadata.tsv')
coverage=read(src/'08_refreshed_family_coverage.tsv')
for r in coverage:
 if r['wang_family']=='PPT':r.update(current_status='CORRECTED_PROTEOME_SIX_REFERENCE_SCREEN_DOMAIN_REVIEW_ONLY',next_action='BOUNDED_ALIGNMENT_PHYLOGENETIC_SPECIFICITY_AND_COMPLETENESS_REVIEW',search_completeness='BOUNDED_PANEL_SCREEN_COMPLETED_EXHAUSTIVENESS_NOT_ESTABLISHED')
table('11_refreshed_family_coverage.tsv',coverage,list(coverage[0]))
run(['diamond','version'],'diamond_version.txt');run(['hmmsearch','-h'],'hmmer_version.txt')
def clean(s):return s.replace(str(dbuse),'<PFAM_HMM>').replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>')
for p in out.iterdir():
 if p.suffix in {'.log','.domtbl'}:p.write_text(clean(p.read_text()))
qc=dict(status='PPT_SCREEN_DOMAIN_REVIEW_READY',corrected_proteome_sequences=len(current),historical_table_candidates=len(hist),historical_OG0000369_candidates=len(ogfa),historical_reconciliation=dict(Counter(r['sequence_status'] for r in audit)),historical_unresolved=sum(not r['corrected_protein_id'] for r in audit),current_candidates=len(review),new_candidates=sum(r['new_relative_to_historical_PPT_tables'] for r in review),similarity_categories=dict(Counter(r['similarity_category'] for r in review)),priority_candidates=len(priority),shared_PPT_reference_models=sorted(shared),flagged_candidates=sum(r['review_flags']!='NO_HEURISTIC_FLAG' for r in review),unresolved_locus_mappings=sum(not r['gene_id'] for r in review),final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step65_manifest_sha256=sha(src/'SHA256SUMS.txt'),corrected_manifest_sha256=sha(cor/'SHA256SUMS.txt'),corrected_proteome_sha256=sha(proteome),mapping_hashes=mapping_hashes,pfam_sha256=sha(dbuse),commands=[[clean(x) for x in cmd] for cmd in commands]),indent=2)+'\n')
(out/'README.txt').write_text('PPT-only provenance reconciliation and expanded specificity/domain screen. Historical283 records retained,37 OG candidates tracked; no AMYLASE, aligned FASTA or count-table values admitted as historical candidates. Historical sequence agreement assessed only from unaligned PPT FASTAs. Missing historical sequence does not imply mismatch. Old4-reference screen lacks TPT/GPT2 and uses a different PPT2; uniform6-reference DIAMOND search of hash-verified corrected proteome performed. All panel hits at E<=1e-3 plus mapped historical candidates retained, including comparator-favored and no-current-hit entries. Pfam-A GA domain architecture screened for candidates and all6 references. Shared domain supports transporter family compatibility, not PPT specificity. Similarity winner and score margin are descriptive, not functional assignments. 80% coverage and10% margin are review heuristics, not acceptance thresholds. No tree or final counts generated. No protein promoted to counting ledger; all candidates REVIEW_ONLY. Previous integrated ledger unchanged. No independent read support or exhaustive genomic absence claims. Work databases and full proteome omitted from bundle; outputs checksummed.\n')
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='PPT_combined_specificity_domain_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle))

PY
