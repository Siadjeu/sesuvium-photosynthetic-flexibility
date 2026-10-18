#!/usr/bin/env bash
#SBATCH --job-name=remaining_batch_screen
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=remaining_batch_screen_%j.out
#SBATCH --error=remaining_batch_screen_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v diamond
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









import urllib.parse,copy
src=release('step25kz17.69','REMAINING_AUDIT_RELEASE');cor=release('step25kz17.21','CORRECTED_INPUT_RELEASE');prov=json.loads((src/'PROVENANCE.json').read_text())
pf=cor/'full_outputs/corrected_protein_candidates.faa';need(sha(pf)==prov['corrected_proteome_sha256'] and sha(cor/'SHA256SUMS.txt')==prov['corrected_manifest_sha256'],'Corrected input changed')
current,_=fasta(pf);families=['NAD-ME','NADP-ME','MDH','ASPAT','AlaAT','PEPCK','BASS2','MPC','MEP','PPA']
out=base/'step25kz17.70'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
for d in ['references','families','work']:(out/d).mkdir()
def table(n,rr,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in sorted(d.items())))
alias=defaultdict(set)
for sp,h in prov['mapping_hashes'].items():
 p=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv');need(sha(p)==h,'Mapping changed '+sp)
 for m in read(p):
  for fld in ['gene_id','protein_id','transcript_id']:
   v=m.get(fld,'')
   if v and v not in {'.','NA'}:alias[v].add((sp,m['gene_id']));alias[sp+'__'+v].add((sp,m['gene_id']))
bygene=defaultdict(set)
for k in current:
 if len(alias.get(k,set()))==1:bygene[next(iter(alias[k]))].add(k)
# Recover complete rows from exact shared-table locations found in step69.
mentions=read(src/'02_definition_and_provenance_mentions.tsv');requested=defaultdict(dict)
for m in mentions:
 if Path(m['source_file']).suffix.lower()=='.tsv':requested[m['source_file']][int(m['line_number'])]=m['families']
shared=[];sharedstatus=[];recovered=defaultdict(set);readbudget=512*1024**2;totalbytes=0
for rel,lines in sorted(requested.items(),key=lambda x:(0 if 'candidate' in x[0].lower() else 1,x[0])):
 p=root/rel
 if not p.is_file():sharedstatus.append(dict(source_file=rel,status='MISSING',sha256='',rows=0));continue
 size=p.stat().st_size
 if size>100*1024**2 or totalbytes+size>readbudget:sharedstatus.append(dict(source_file=rel,status='SIZE_BUDGET_SKIP',sha256='',rows=0));continue
 totalbytes+=size;n=0
 with p.open(errors='replace') as f:
  reader=csv.reader(f,delimiter='\t');fields=next(reader,[])
  for lineno,values in enumerate(reader,2):
   if lineno not in lines:continue
   if len(values)!=len(fields):continue
   row=dict(zip(fields,values));n+=1;ids=[]
   for field,value in row.items():
    if field.lower() in {'protein_id','candidate_id','gene_id','sequence_id','query_id','query','protein','gene','transcript_id'} and value and value not in {'NA','.','NONE'}:ids.append(value)
   sp=row.get('species',row.get('species_code',''));mapped=set();mapping_states=set()
   for ident in ids:
    mm=alias.get(ident,set())
    if sp in prov['mapping_hashes']:mm=alias.get(sp+'__'+ident,set()) or {x for x in mm if x[0]==sp}
    if len(mm)!=1:mapping_states.add('AMBIGUOUS_OR_MISSING_LOCUS');continue
    poss=bygene[next(iter(mm))]
    if ident in poss:mapped.add(ident)
    elif len(poss)==1:mapped.update(poss)
    else:mapping_states.add('AMBIGUOUS_OR_MISSING_CORRECTED_ISOFORM')
   # Row identity is not sequence equivalence or biochemical verification.
   fs=set(lines[lineno].split(';'))&set(families)
   if len(mapped)>1:mapping_states.add('MULTIPLE_PROTEINS_IN_ROW')
   for fam in fs:recovered[fam].update(mapped)
   shared.append(dict(source_file=rel,line_number=lineno,family_mentions=';'.join(sorted(fs)),mapped_corrected_ids=';'.join(sorted(mapped)),mapping_notes=';'.join(sorted(mapping_states)) or ('DOCUMENTED_ID_ONLY_NO_HISTORICAL_SEQUENCE_CHECK' if mapped else 'NO_IDENTIFIER_FIELD'),complete_row_json=json.dumps(row,ensure_ascii=False)))
 sharedstatus.append(dict(source_file=rel,status='READ_SELECTED_COMPLETE_ROWS',sha256=sha(p),rows=n))
table('01_shared_table_records.tsv',shared,['source_file','line_number','family_mentions','mapped_corrected_ids','mapping_notes','complete_row_json']);table('02_shared_table_provenance.tsv',sharedstatus,['source_file','status','sha256','rows'])
# Source-backed working screening scopes. None constitutes final Wang equivalence.
ecs={'NAD-ME':['1.1.1.39','1.1.1.38'],'NADP-ME':['1.1.1.40'],'MDH_NAD':['1.1.1.37'],'MDH_NADP':['1.1.1.82'],'ASPAT':['2.6.1.1'],'AlaAT':['2.6.1.2'],'AMINO_COMPARATOR':['2.6.1.4','2.6.1.42','2.6.1.5','2.6.1.44'],'LDH_COMPARATOR':['1.1.1.27'],'PEPCK':['4.1.1.49'],'PEPC_CONTEXT':['4.1.1.31']}
# A small reviewed human GTP-PEPCK panel supplies cofactor-class comparison, not plant positives.
queries={tag:('(organism_id:3702) AND (reviewed:true) AND ('+' OR '.join('ec:'+e for e in ee)+')',3702) for tag,ee in ecs.items()}
queries['PEPCK_GTP_COMPARATOR']=('(organism_id:9606) AND (reviewed:true) AND (ec:4.1.1.32)',9606)
queries['BASS']=('(organism_id:3702) AND (reviewed:true) AND ('+' OR '.join('gene_exact:BASS'+str(i) for i in range(1,7))+')',3702)
queries['MPC']=('(organism_id:3702) AND (reviewed:true) AND ('+' OR '.join('gene_exact:MPC'+str(i) for i in range(1,5))+')',3702)
refseq={};refmeta={};groups=defaultdict(set);downloads=[]
def ec_values(obj):
 found=set()
 if isinstance(obj,dict):
  for k,v in obj.items():
   if k=='ecNumbers':found.update(x.get('value','') for x in v)
   else:found.update(ec_values(v))
 elif isinstance(obj,list):
  for x in obj:found.update(ec_values(x))
 return found
print('Retrieving reviewed, metadata-checked reference records',flush=True)
for tag,(query,tax) in queries.items():
 saved=(dict(refseq),copy.deepcopy(refmeta),copy.deepcopy(groups))
 url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode({'query':query,'format':'json','size':500});accepted=0
 try:
  request=urllib.request.Request(url,headers={'User-Agent':'Remaining-family-reference-screen/1.0'})
  with urllib.request.urlopen(request,timeout=45) as response:
   payload=response.read();pagination=response.headers.get('Link','')
  (out/'references'/(tag+'.json')).write_bytes(payload);need('rel="next"' not in pagination,'Reference query requires pagination; blocked rather than truncated')
  entries=json.loads(payload).get('results',[]);need(len(entries)<=100,'Reference panel too large')
  for entry in entries:
   if entry.get('organism',{}).get('taxonId')!=tax or entry.get('entryType')!='UniProtKB reviewed (Swiss-Prot)':continue
   desc=entry.get('proteinDescription',{}).get('recommendedName',{}).get('fullName',{}).get('value','');genes={g.get('geneName',{}).get('value','').upper() for g in entry.get('genes',[])};ee=ec_values(entry.get('proteinDescription',{}));group=tag
   if tag in ecs and not(ee&set(ecs[tag])):continue
   if tag=='PEPCK_GTP_COMPARATOR' and '4.1.1.32' not in ee:continue
   if tag in {'NAD-ME','NADP-ME'} and not re.search(r'malic enzyme|malate dehydrogenase',desc,re.I):continue
   if tag.startswith('MDH') and not re.search(r'malate dehydrogenase',desc,re.I):continue
   if tag=='ASPAT' and not re.search(r'aspartate aminotransferase',desc,re.I):continue
   if tag=='AlaAT' and not re.search(r'^alanine aminotransferase',desc,re.I):continue
   if tag in {'PEPCK','PEPCK_GTP_COMPARATOR'} and not re.search(r'phosphoenolpyruvate carboxykinase',desc,re.I):continue
   if tag=='BASS':
    gg=genes&{'BASS'+str(i) for i in range(1,7)}
    if len(gg)!=1:continue
    group='BASS2' if 'BASS2' in gg else 'BASS_COMPARATOR'
    if group=='BASS2' and 'pyruvate' not in desc.lower():continue
   if tag=='MPC':
    if not(genes&{'MPC'+str(i) for i in range(1,5)}) or 'mitochondrial pyruvate carrier' not in desc.lower():continue
   acc=entry['primaryAccession'];sequence=entry['sequence']['value'];need(valid(sequence),'Invalid reference sequence')
   if acc in refseq:need(refseq[acc]==sequence,'Inconsistent reference sequence')
   refseq[acc]=sequence;groups[group].add(acc);accepted+=1
   if acc not in refmeta:refmeta[acc]=dict(reference_id=acc,genes=';'.join(sorted(genes)),description=desc,taxon_id=tax,ec_numbers=';'.join(sorted(ee)),length_aa=len(sequence),sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),entry_url='https://www.uniprot.org/uniprotkb/'+acc,query_groups=set())
   refmeta[acc]['query_groups'].add(group)
  downloads.append(dict(query_group=tag,status='READY' if accepted else 'NO_METADATA_VALIDATED_RECORD',accepted=accepted,url=url,error=''))
 except (OSError,ValueError,KeyError,SystemExit) as e:
  refseq,refmeta,groups=saved
  downloads.append(dict(query_group=tag,status='RETRIEVAL_OR_VALIDATION_BLOCKED',accepted=0,url=url,error=str(e)))
table('03_reference_retrieval_status.tsv',downloads,['query_group','status','accepted','url','error'])
rr=[]
for x in refmeta.values():rr.append({**x,'query_groups':';'.join(sorted(x['query_groups']))})
table('04_reference_metadata.tsv',rr,['reference_id','genes','description','taxon_id','ec_numbers','length_aa','sequence_sha256','entry_url','query_groups']);writefa('05_combined_reference_panel.faa',refseq)
# Positives/comparators are family-relative; same accession can have several documented activities.
def union(names):return set().union(*(groups[x] for x in names))
specs={
'NAD-ME':(['NAD-ME'],['NADP-ME']), 'NADP-ME':(['NADP-ME'],['NAD-ME']),
'MDH':(['MDH_NAD','MDH_NADP'],['LDH_COMPARATOR']),
'ASPAT':(['ASPAT'],['AlaAT','AMINO_COMPARATOR']), 'AlaAT':(['AlaAT'],['ASPAT','AMINO_COMPARATOR']),
'PEPCK':(['PEPCK'],['PEPCK_GTP_COMPARATOR']), 'BASS2':(['BASS2'],['BASS_COMPARATOR']), 'MPC':(['MPC'],[])}
scopes=[];panels={};roles=[]
for fam in families:
 if fam not in specs:scopes.append(dict(wang_family=fam,status='BLOCKED_UNRESOLVED_WANG_DEFINITION',positive_references=0,comparator_references=0,scope='NO_SCREEN',limitation='MEP_OR_PPA_ABBREVIATION_NOT_RESOLVED'));continue
 pg,cg=specs[fam];pos=union(pg);comp=union(cg)-pos
 # For MDH both cofactor groups are required; preserve subgroup metadata.
 missing=[g for g in pg if not groups[g]];enabled=bool(pos) and not missing and (bool(comp) or fam=='MPC')
 limitation='BOUNDED_ARABIDOPSIS_PANEL;WANG_COUNT_EQUIVALENCE_PENDING'
 if fam=='MPC':limitation+=';COMPONENT_PANEL_ONLY_NO_INDEPENDENT_NON_MPC_COMPARATOR'
 if fam=='MDH':limitation+=';NAD_AND_NADP_SUBGROUPS_RETAINED_SEPARATELY'
 if fam=='PEPCK':limitation+=';HUMAN_GTP_CLASS_COMPARATORS_NOT_PLANT_POSITIVES'
 scopes.append(dict(wang_family=fam,status='READY_FOR_REVIEW_ONLY_SCREEN' if enabled else 'BLOCKED_REFERENCE_PANEL_INCOMPLETE',positive_references=len(pos),comparator_references=len(comp),scope='WORKING_BIOCHEMICAL_FAMILY_NOT_C4_CAM_FUNCTION',limitation=limitation))
 if enabled:panels[fam]=(pos,comp)
 for acc in sorted(pos|comp):roles.append(dict(wang_family=fam,reference_id=acc,role='POSITIVE' if acc in pos else 'COMPARATOR',subgroups=';'.join(sorted(refmeta[acc]['query_groups'])),source_url=refmeta[acc]['entry_url']))
table('06_family_screening_scope.tsv',scopes,['wang_family','status','positive_references','comparator_references','scope','limitation']);table('07_family_reference_roles.tsv',roles,['wang_family','reference_id','role','subgroups','source_url'])
commands=[];cpu=os.environ.get('SLURM_CPUS_PER_TASK','4');allhits=defaultdict(list)
def run(cmd,log):
 commands.append(cmd)
 with (out/log).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
if panels:
 print('Running one combined DIAMOND screen',flush=True)
 run(['diamond','makedb','--threads',cpu,'--in',str(out/'05_combined_reference_panel.faa'),'--db',str(out/'work/panel')],'08_makedb.log')
 run(['diamond','blastp','--threads',cpu,'--query',str(pf),'--db',str(out/'work/panel'),'--out',str(out/'08_all_reference_hits.tsv'),'--very-sensitive','--evalue','1e-3','--max-target-seqs','0','--max-hsps','1','--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore'],'08_diamond.log')
 for line in (out/'08_all_reference_hits.tsv').read_text().splitlines():
  x=line.split('\t');need(len(x)==12 and x[0] in current and x[1] in refseq,'Invalid DIAMOND output')
  allhits[x[0]].append(dict(reference_id=x[1],pident=float(x[2]),qcov=(int(x[7])-int(x[6])+1)/int(x[4]),rcov=(int(x[9])-int(x[8])+1)/int(x[5]),evalue=float(x[10]),bitscore=float(x[11])))
summaries=[];unioncandidates={}
fields=['protein_id','species','gene_id','length_aa','shared_table_recovered','best_positive','best_positive_groups','best_positive_bitscore','best_comparator','best_comparator_bitscore','positive_query_coverage','positive_reference_coverage','positive_evalue','similarity_category','review_flags','count_status']
for fam in families:
 if fam not in panels:
  summaries.append(dict(wang_family=fam,screen_status='BLOCKED',screened_candidate_records='NA',positive_hit_records='NA',shared_table_recovered=len(recovered[fam]),similarity_categories=''));continue
 pos,comp=panels[fam];ids={k for k,hh in allhits.items() if any(h['reference_id'] in pos|comp for h in hh)}|recovered[fam];result=[]
 for k in sorted(ids):
  hh=allhits[k];pp=[h for h in hh if h['reference_id'] in pos];cc=[h for h in hh if h['reference_id'] in comp];p=max(pp,key=lambda x:x['bitscore']) if pp else None;c=max(cc,key=lambda x:x['bitscore']) if cc else None
  category='NO_PANEL_HIT' if not p and not c else 'COMPARATOR_ONLY' if not p else 'POSITIVE_ONLY' if not c else 'POSITIVE_FAVORED' if p['bitscore']>c['bitscore'] else 'COMPARATOR_FAVORED' if p['bitscore']<c['bitscore'] else 'TIED'
  mm=alias.get(k,set());sp,g=next(iter(mm)) if len(mm)==1 else ('','');flags=[]
  if len(mm)!=1:flags.append('UNRESOLVED_LOCUS_MAPPING')
  if p and p['rcov']<.8:flags.append('LOW_REFERENCE_COVERAGE_HEURISTIC')
  if p and c and abs(p['bitscore']-c['bitscore'])/max(p['bitscore'],c['bitscore'])<.1:flags.append('SMALL_SCORE_MARGIN_HEURISTIC')
  if fam=='MPC':flags.append('NO_NON_MPC_COMPARATOR')
  result.append(dict(protein_id=k,species=sp,gene_id=g,length_aa=len(current[k]),shared_table_recovered=k in recovered[fam],best_positive=p['reference_id'] if p else '',best_positive_groups=';'.join(sorted(refmeta[p['reference_id']]['query_groups'])) if p else '',best_positive_bitscore=p['bitscore'] if p else '',best_comparator=c['reference_id'] if c else '',best_comparator_bitscore=c['bitscore'] if c else '',positive_query_coverage=round(p['qcov'],6) if p else '',positive_reference_coverage=round(p['rcov'],6) if p else '',positive_evalue=p['evalue'] if p else '',similarity_category=category,review_flags=';'.join(flags) or 'NONE_AT_SIMILARITY_STAGE',count_status='REVIEW_ONLY'))
  unioncandidates[k]=current[k]
 table('families/'+fam+'_screen.tsv',result,fields);writefa('families/'+fam+'_candidates.faa',{k:current[k] for k in ids});writefa('families/'+fam+'_references.faa',{k:refseq[k] for k in pos|comp})
 summaries.append(dict(wang_family=fam,screen_status='SCREENED_REVIEW_ONLY',screened_candidate_records=len(result),positive_hit_records=sum(bool(x['best_positive']) for x in result),shared_table_recovered=len(recovered[fam]),similarity_categories=json.dumps(dict(Counter(x['similarity_category'] for x in result)),sort_keys=True)))
table('09_batch_screen_summary.tsv',summaries,list(summaries[0]));writefa('10_unique_candidate_proteins.faa',unioncandidates)
for n in ['source_01_integrated_locus_evidence.tsv','source_07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/n)
coverage=read(src/'08_refreshed_family_coverage.tsv')
for r in coverage:
 if r['wang_family'] in families:r.update(current_status='REFERENCE_SCREEN_REVIEW_ONLY' if r['wang_family'] in panels else 'BLOCKED_DEFINITION_OR_REFERENCE_PANEL',next_action='BATCH_DOMAIN_AND_COMPLETENESS_VALIDATION' if r['wang_family'] in panels else 'RESOLVE_DEFINITION_OR_REFERENCE_GAP',search_completeness='NOT_ESTABLISHED')
table('11_refreshed_family_coverage.tsv',coverage,list(coverage[0]));run(['diamond','version'],'diamond_version.txt')
(out/'SOURCE_DEFINITIONS.txt').write_text('Working family screening scopes, not final Wang-compatible counts. Original study: https://pmc.ncbi.nlm.nih.gov/articles/PMC10663116/ . Enzyme definitions: https://enzyme.expasy.org/EC/ followed by each EC in reference metadata. NAD-ME versus NADP-ME: https://pmc.ncbi.nlm.nih.gov/articles/PMC2287332/ . BASS specificity comparison: https://pmc.ncbi.nlm.nih.gov/articles/PMC5435425/ . MPC components: https://pmc.ncbi.nlm.nih.gov/articles/PMC6501077/ . Complete reviewed UniProt JSON query results and individual entry URLs are archived. Arabidopsis reference species is not project code AT. EC/name checks do not prove experimental support for every record; per-entry evidence remains in JSON. MEP and PPA remain blocked.\n')
(out/'README.txt').write_text('One shared-table recovery/reference preparation/combined DIAMOND job. Selected complete rows reread at step69 source line numbers; recovered IDs mapped through exact annotation aliases. Table mentions are provenance leads, not biological decisions; no historical sequence equality inferred. Reference records must be reviewed, expected taxon, correct EC or gene and relevant name. Missing references block the affected panel; other complete panels run. MPC is an explicitly limited component-only positive screen with no independent non-MPC comparator; no specificity conclusion. MDH NAD/NADP references remain distinct subgroups. PEPCK working ATP scope uses human GTP class comparators. Related biochemical activities can overlap: score winners are not substrate/cofactor proof. All hits at1e-3 retained with max-target-seqs0 and one HSP; shared-table candidates retained even without hits. Family tables may share proteins; unique union FASTA supports next batched domain validation. No domain checks/trees or count promotions in this step. Final Wang scope and C4/CAM function remain unestablished. Missing/blocked results NA, not zero.\n')
def clean(s):return s.replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>')
for p in out.glob('*.log'):p.write_text(clean(p.read_text()))
(out/'PROVENANCE.json').write_text(json.dumps(dict(step69_manifest_sha256=sha(src/'SHA256SUMS.txt'),corrected_proteome_sha256=sha(pf),mapping_hashes=prov['mapping_hashes'],reference_retrieval_time_UTC=datetime.now(timezone.utc).isoformat(),commands=[[clean(x) for x in cmd] for cmd in commands]),indent=2)+'\n')
qc=dict(status='BATCH_SCREEN_COMPLETE_WITH_EXPLICIT_BLOCKS',families_screened=list(panels),families_blocked=[x for x in families if x not in panels],reference_sequences=len(refseq),shared_rows=len(shared),unique_screen_candidates=len(unioncandidates),family_summary=summaries,integrated_records=len(read(src/'source_01_integrated_locus_evidence.tsv')),final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts);(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Remaining_family_reference_screen/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle))

PY
