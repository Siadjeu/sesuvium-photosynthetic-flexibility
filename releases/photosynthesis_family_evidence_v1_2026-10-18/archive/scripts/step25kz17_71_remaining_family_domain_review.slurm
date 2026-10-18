#!/usr/bin/env bash
#SBATCH --job-name=remaining_batch_domains
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=remaining_batch_domains_%j.out
#SBATCH --error=remaining_batch_domains_%j.err
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










import urllib.parse
src=release('step25kz17.70','REMAINING_SCREEN_RELEASE');cor=release('step25kz17.21','CORRECTED_INPUT_RELEASE')
prov=json.loads((src/'PROVENANCE.json').read_text());pf=cor/'full_outputs/corrected_protein_candidates.faa'
need(sha(pf)==prov['corrected_proteome_sha256'],'Corrected protein provenance mismatch')
current,_=fasta(pf);refs,_=fasta(src/'05_combined_reference_panel.faa');proteins,_=fasta(src/'10_unique_candidate_proteins.faa')
for k,s in proteins.items():need(k in current and norm(s)==norm(current[k]),'Candidate provenance mismatch '+k)
roles=read(src/'07_family_reference_roles.tsv');metadata=read(src/'04_reference_metadata.tsv')
for r in metadata:need(hashlib.sha256(refs[r['reference_id']].encode()).hexdigest()==r['sequence_sha256'],'Reference changed')
alias=defaultdict(set)
for sp,h in prov['mapping_hashes'].items():
 p=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv');need(sha(p)==h,'Mapping changed '+sp)
 for m in read(p):
  for field in ['gene_id','protein_id','transcript_id']:
   v=m.get(field,'')
   if v and v not in {'.','NA'}:alias[v].add((sp,m['gene_id']));alias[sp+'__'+v].add((sp,m['gene_id']))
out=base/'step25kz17.71'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
for d in ['work','families','references']:(out/d).mkdir()
def table(n,rr,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in sorted(d.items())))
commands=[];cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
def run(cmd,log):
 commands.append(cmd)
 with (out/log).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
def ec_values(obj):
 found=set()
 if isinstance(obj,dict):
  for k,v in obj.items():
   if k=='ecNumbers':found.update(x.get('value','') for x in v)
   else:found.update(ec_values(v))
 elif isinstance(obj,list):
  for x in obj:found.update(ec_values(x))
 return found
# Broaden taxonomic scope, never silently relax the reviewed requirement.
query='(taxonomy_id:33090) AND (reviewed:true) AND (ec:1.1.1.27)'
url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode(dict(query=query,format='json',size=500))
ldh={};ldhmeta=[];retrieval_error=''
print('Recovering reviewed plant LDH comparator panel',flush=True)
try:
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Remaining-family-domain-review/1.0'}),timeout=60) as response:
  payload=response.read();link=response.headers.get('Link','')
 (out/'references/plant_LDH.json').write_bytes(payload)
 need('rel="next"' not in link,'Paginated comparator response; panel blocked')
 entries=json.loads(payload)['results'];need(len(entries)<=30,'More than30 LDH records; bounded panel blocked')
 for e in entries:
  desc=e.get('proteinDescription',{}).get('recommendedName',{}).get('fullName',{}).get('value','')
  if e.get('entryType')!='UniProtKB reviewed (Swiss-Prot)' or 'Viridiplantae' not in e.get('organism',{}).get('lineage',[]):continue
  if '1.1.1.27' not in ec_values(e.get('proteinDescription',{})) or not re.search(r'lactate dehydrogenase',desc,re.I):continue
  acc=e['primaryAccession'];seq=e['sequence']['value'];need(valid(seq),'Invalid comparator sequence');need(acc not in refs,'Unexpected existing comparator accession')
  ldh[acc]=seq;ldhmeta.append(dict(reference_id=acc,description=desc,taxon_id=e['organism']['taxonId'],sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),source_url='https://www.uniprot.org/uniprotkb/'+acc,role='MDH_COMPARATOR_REVIEWED_PLANT_LDH'))
except (OSError,ValueError,KeyError,SystemExit) as e:
 ldh={};ldhmeta=[];retrieval_error=str(e)
refs.update(ldh)
for k in ldh:roles.append(dict(wang_family='MDH',reference_id=k,role='COMPARATOR',subgroups='LDH_COMPARATOR',source_url='https://www.uniprot.org/uniprotkb/'+k))
table('01_LDH_reference_metadata.tsv',ldhmeta,['reference_id','description','taxon_id','sequence_sha256','source_url','role'])
(out/'01_LDH_retrieval.json').write_text(json.dumps(dict(url=url,status='READY' if ldh else 'BLOCKED_NO_VALIDATED_COMPARATOR',error=retrieval_error,accepted=len(ldh)),indent=2)+'\n')
mdhpos={r['reference_id'] for r in roles if r['wang_family']=='MDH' and r['role']=='POSITIVE'};need(mdhpos,'Missing MDH positives')
mdhhits=defaultdict(list)
def loadhits(path,allowed):
 for line in path.read_text().splitlines():
  x=line.split('\t');need(len(x)==12 and x[0] in current,'Invalid hit row')
  if x[1] not in allowed:continue
  mdhhits[x[0]].append(dict(reference_id=x[1],bitscore=float(x[11]),evalue=float(x[10]),qcov=(int(x[7])-int(x[6])+1)/int(x[4]),rcov=(int(x[9])-int(x[8])+1)/int(x[5])))
loadhits(src/'08_all_reference_hits.tsv',mdhpos)
if ldh:
 writefa('02_LDH_panel.faa',ldh)
 run(['diamond','makedb','--threads',cpu,'--in',str(out/'02_LDH_panel.faa'),'--db',str(out/'work/ldh')],'02_makedb.log')
 run(['diamond','blastp','--threads',cpu,'--query',str(pf),'--db',str(out/'work/ldh'),'--out',str(out/'02_LDH_hits.tsv'),'--very-sensitive','--evalue','1e-3','--max-target-seqs','0','--max-hsps','1','--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore'],'02_diamond.log')
 loadhits(out/'02_LDH_hits.tsv',set(ldh))
shared=set()
for r in read(src/'01_shared_table_records.tsv'):
 if 'MDH' in r['family_mentions'].split(';'):shared.update(filter(None,r['mapped_corrected_ids'].split(';')))
need(shared<=set(current),'Shared MDH ID not in corrected input')
mdhids=set(mdhhits)|shared;familyrows={p.name[:-11]:read(p) for p in (src/'families').glob('*_screen.tsv')}
mdhrows=[]
for k in sorted(mdhids):
 pp=[h for h in mdhhits[k] if h['reference_id'] in mdhpos];cc=[h for h in mdhhits[k] if h['reference_id'] in ldh]
 p=max(pp,key=lambda x:x['bitscore']) if pp else None;c=max(cc,key=lambda x:x['bitscore']) if cc else None
 cat='NO_PANEL_HIT' if not p and not c else 'COMPARATOR_ONLY' if not p else 'POSITIVE_ONLY' if not c else 'POSITIVE_FAVORED' if p['bitscore']>c['bitscore'] else 'COMPARATOR_FAVORED' if p['bitscore']<c['bitscore'] else 'TIED'
 mm=alias.get(k,set());sp,g=next(iter(mm)) if len(mm)==1 else ('','');flags=[]
 if not ldh:flags.append('MDH_COMPARATOR_PANEL_BLOCKED')
 if p and p['rcov']<.8:flags.append('LOW_REFERENCE_COVERAGE_HEURISTIC')
 if p and c and abs(p['bitscore']-c['bitscore'])/max(p['bitscore'],c['bitscore'])<.1:flags.append('SMALL_SCORE_MARGIN_HEURISTIC')
 mdhrows.append(dict(protein_id=k,species=sp,gene_id=g,best_positive=p['reference_id'] if p else '',best_comparator=c['reference_id'] if c else '',similarity_category=cat,review_flags=';'.join(flags) or 'NONE_AT_SIMILARITY_STAGE',shared_table_recovered=k in shared))
familyrows['MDH']=mdhrows
for k in mdhids:proteins[k]=current[k]
proteins={k:norm(s) for k,s in proteins.items()};need(len(proteins)<=5000,'Candidate union exceeds5000: stop before domain analysis');need(all(valid(s) for s in proteins.values()),'Invalid protein sequence')
# Deduplicate exact sequences, preserving all candidate/reference identities.
seqid={};seqs={};idmap={};idrows=[]
for kind,d in [('CANDIDATE',proteins),('REFERENCE',refs)]:
 for k,s in sorted(d.items()):
  s=norm(s)
  if s not in seqid:seqid[s]='D'+str(len(seqid)+1).zfill(6);seqs[seqid[s]]=s
  idmap[(kind,k)]=seqid[s];idrows.append(dict(kind=kind,original_id=k,domain_sequence_id=seqid[s],length_aa=len(s),sequence_sha256=hashlib.sha256(s.encode()).hexdigest()))
table('03_domain_sequence_map.tsv',idrows,list(idrows[0]));writefa('03_deduplicated_domain_input.faa',seqs);writefa('03_candidate_proteins.faa',proteins);writefa('03_reference_panel.faa',refs)
table('03_reference_roles.tsv',roles,list(roles[0]))
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
run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'06_Pfam.domtbl'),str(dbuse),str(out/'03_deduplicated_domain_input.faa')],'06_hmmsearch.log')
domains=defaultdict(list);dh=[]
for line in (out/'06_Pfam.domtbl').read_text().splitlines():
 if not line or line.startswith('#'):continue
 x=line.split(maxsplit=22);r=dict(sequence_id=x[0],model=x[3],accession=x[4].split('.')[0],model_length=int(x[5]),i_evalue=float(x[12]),score=float(x[13]),hmm_from=int(x[15]),hmm_to=int(x[16]),ali_from=int(x[17]),ali_to=int(x[18]));dh.append(r);domains[x[0]].append(r)
table('07_Pfam_domain_hits.tsv',dh,['sequence_id','model','accession','model_length','i_evalue','score','hmm_from','hmm_to','ali_from','ali_to'])

def hits(kind,k):return domains[idmap[(kind,k)]]
def models(kind,k):return {h['accession'] for h in hits(kind,k)}
def cov(kind,k,m):
 hh=[h for h in hits(kind,k) if h['accession']==m];covered=set()
 for h in hh:covered.update(range(h['hmm_from'],h['hmm_to']+1))
 return len(covered)/hh[0]['model_length'] if hh else 0.0
cal=[]
for r in roles:
 k=r['reference_id']
 for m in sorted(models('REFERENCE',k)):
  cal.append({**r,'model':m,'model_coverage':round(cov('REFERENCE',k,m),6)})
table('08_reference_domain_calibration.tsv',cal,list(roles[0])+['model','model_coverage'])
allrows=[];summary=[]
for fam,rows in sorted(familyrows.items()):
 panel={r['reference_id'] for r in roles if r['wang_family']==fam};pos={r['reference_id'] for r in roles if r['wang_family']==fam and r['role']=='POSITIVE'}
 panelmodels=set().union(*(models('REFERENCE',x) for x in panel));result=[]
 for r in rows:
  k=r['protein_id'];need(k in proteins,'Candidate missing');p=r.get('best_positive','');flags=[f for f in r.get('review_flags','').split(';') if f and f not in {'NONE_AT_SIMILARITY_STAGE','NO_HEURISTIC_FLAG'}]
  expected=models('REFERENCE',p) if p else set();missing=expected-models('CANDIDATE',k)
  if not hits('CANDIDATE',k):flags.append('NO_PFAM_GA_HIT')
  if not p:flags.append('NO_POSITIVE_SIMILARITY_HIT')
  elif not expected:flags.append('POSITIVE_REFERENCE_NO_PFAM_CALIBRATION')
  if missing:flags.append('MISSING_BEST_POSITIVE_MODEL_HEURISTIC')
  low=[m for m in expected if cov('CANDIDATE',k,m)<.8*cov('REFERENCE',p,m)]
  if low:flags.append('LOW_BEST_POSITIVE_MODEL_COVERAGE_HEURISTIC')
  # Additional overlapping models are recorded but not interpreted as extra domains.
  extra=[h for h in hits('CANDIDATE',k) if h['accession'] not in panelmodels]
  nonoverlap=[h for h in extra if not any(max(h['ali_from'],a['ali_from'])<=min(h['ali_to'],a['ali_to']) for a in hits('CANDIDATE',k) if a['accession'] in panelmodels)]
  if nonoverlap:flags.append('NONOVERLAPPING_ADDITIONAL_MODEL_HEURISTIC')
  if p and (len(proteins[k])<.6*len(refs[p]) or len(proteins[k])>1.5*len(refs[p])):flags.append('LENGTH_OUTLIER_VS_BEST_POSITIVE_HEURISTIC')
  if 'X' in proteins[k]:flags.append('AMBIGUOUS_RESIDUES')
  if not r.get('gene_id'):flags.append('UNRESOLVED_LOCUS_MAPPING')
  result.append(dict(wang_family=fam,protein_id=k,species=r.get('species',''),gene_id=r.get('gene_id',''),length_aa=len(proteins[k]),best_positive=p,best_comparator=r.get('best_comparator',''),similarity_category=r['similarity_category'],Pfam_models=';'.join(sorted(models('CANDIDATE',k))),best_positive_model_coverage=';'.join(m+':'+str(round(cov('CANDIDATE',k,m),6)) for m in sorted(expected)),missing_models=';'.join(sorted(missing)),additional_models=';'.join(sorted({h['accession'] for h in extra})),review_flags=';'.join(sorted(set(flags))) or 'NO_HEURISTIC_FLAG',count_status='REVIEW_ONLY'))
 table('families/'+fam+'_domain_review.tsv',result,list(result[0]));writefa('families/'+fam+'_candidates.faa',{r['protein_id']:proteins[r['protein_id']] for r in result});allrows.extend(result)
 summary.append(dict(wang_family=fam,candidate_records=len(result),flagged_records=sum(r['review_flags']!='NO_HEURISTIC_FLAG' for r in result),no_positive_hit=sum(not r['best_positive'] for r in result),status='DOMAIN_REVIEW_ONLY' if fam!='MDH' or ldh else 'DOMAIN_REVIEW_COMPARATOR_BLOCKED'))
table('09_combined_domain_review.tsv',allrows,list(allrows[0]));table('10_family_summary.tsv',summary,list(summary[0]))
queue=[r for r in allrows if r['review_flags']!='NO_HEURISTIC_FLAG'];table('10_exception_queue.tsv',queue,list(allrows[0]))
for n in ['source_01_integrated_locus_evidence.tsv','source_07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/n)
for n in ['04_reference_metadata.tsv','06_family_screening_scope.tsv','09_batch_screen_summary.tsv','01_shared_table_records.tsv']:shutil.copyfile(src/n,out/('source_'+n))
coverage=read(src/'11_refreshed_family_coverage.tsv')
for r in coverage:
 if r['wang_family'] in familyrows:r.update(current_status='BATCH_DOMAIN_REVIEW_ONLY',next_action='BOUNDED_FAMILY_SPECIFIC_DISCRIMINATION_AND_EXCEPTION_REVIEW',search_completeness='NOT_ESTABLISHED')
 if r['wang_family']=='MDH' and not ldh:r.update(current_status='DOMAIN_REVIEW_COMPARATOR_BLOCKED',next_action='RESOLVE_LDH_REFERENCE_GAP')
table('11_refreshed_family_coverage.tsv',coverage,list(coverage[0]))
run(['hmmsearch','-h'],'hmmer_version.txt');run(['diamond','version'],'diamond_version.txt')
def clean(s):return s.replace(str(dbuse),'<PFAM_HMM>').replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>')
for p in out.iterdir():
 if p.suffix in {'.log','.domtbl'}:p.write_text(clean(p.read_text()))
(out/'PROVENANCE.json').write_text(json.dumps(dict(step70_manifest_sha256=sha(src/'SHA256SUMS.txt'),corrected_manifest_sha256=sha(cor/'SHA256SUMS.txt'),corrected_proteome_sha256=sha(pf),mapping_hashes=prov['mapping_hashes'],pfam_sha256=sha(dbuse),commands=[[clean(x) for x in c] for c in commands]),indent=2)+'\n')
qc=dict(status='BATCH_DOMAIN_REVIEW_COMPLETE',unique_candidate_ids=len(proteins),reference_ids=len(refs),unique_domain_sequences=len(seqs),MDH_comparator_status='READY' if ldh else 'BLOCKED',MDH_candidates=len(mdhids),families=summary,integrated_records=len(read(src/'source_01_integrated_locus_evidence.tsv')),ledger_unchanged=sha(src/'source_01_integrated_locus_evidence.tsv')==sha(out/'source_01_integrated_locus_evidence.tsv'),final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'README.txt').write_text('One bounded full-Pfam GA pass over exact-sequence-deduplicated candidates and references. Existing seven-family similarity results reused. MDH positive hits reused from previous full-proteome search; mapped shared-table candidates retained. New reviewed Viridiplantae EC1.1.1.27 LDH comparators require matching name, taxonomic lineage and EC; complete JSON archived. Only new comparators searched against corrected proteome. Retrieval failure leaves MDH comparator blocked while domain checks proceed. MDH NAD/NADP groups retained in reference roles. Best-positive domain architecture is a descriptive comparator, not a mandatory universal family architecture. Domain model coverage is union of model positions, may be contributed by several fragments and does not prove intactness. Missing/low coverage, length0.6-1.5x and additional nonoverlapping models are review heuristics. Overlapping Pfam models are not independent domain evidence. No cofactor/substrate assignment from generic domains or score winners. MPC lacks independent non-MPC comparator. No trees, gene repair, functional C4/CAM assignment or count promotion. MEP/PPA stay definition-blocked. Existing ledger copied unchanged. Exhaustiveness and independent read validation not established. One pass only; review summaries before selecting any further family-specific analyses.\n')
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts);(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Remaining_family_domain_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle))
PY
