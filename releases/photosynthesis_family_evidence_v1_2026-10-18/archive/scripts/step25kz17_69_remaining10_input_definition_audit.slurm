#!/usr/bin/env bash
#SBATCH --job-name=remaining10_audit
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=remaining10_audit_%j.out
#SBATCH --error=remaining10_audit_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
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








src=release('step25kz17.68','PPT_CLOSE_RELEASE');cor=release('step25kz17.21','CORRECTED_INPUT_RELEASE')
remaining=read(src/'08_families_without_integrated_ledger.tsv');families=['NAD-ME','NADP-ME','MDH','ASPAT','AlaAT','PEPCK','BASS2','MEP','MPC','PPA']
need({r['wang_family'] for r in remaining}==set(families),'Remaining family scope changed')
proteinfile=cor/'full_outputs/corrected_protein_candidates.faa';expected='c5643d9ce8ba9140cd994eb1518fa7f9cdc242443b3e66532b9ba3f8f039b568'
need(proteinfile.is_file() and sha(proteinfile)==expected,'Corrected protein provenance mismatch')
current,_=fasta(proteinfile);need(len(current)==287197,'Corrected protein record count changed')
out=base/'step25kz17.69'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
for d in ['collected','reconciled_candidates']:(out/d).mkdir()
def table(n,rr,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
# Narrow tokens for discovery. These are retrieval terms, never biological assignments.
terms={'NAD-ME':r'NAD[-_ ]?ME(?:[12])?', 'NADP-ME':r'NADP[-_ ]?ME(?:[1-4])?', 'MDH':r'(?:NADP?[-_ ])?MDH[0-9]?', 'ASPAT':r'ASPAT[0-9]?|ASP[-_ ]AT[0-9]?|ASPARTATE[-_ ]AMINOTRANSFERASE', 'AlaAT':r'ALAAT[0-9]?|ALA[-_ ]AT[0-9]?|ALANINE[-_ ]AMINOTRANSFERASE', 'PEPCK':r'PEPCK[0-9]?|PCK[12]|PHOSPHOENOLPYRUVATE[-_ ]CARBOXYKINASE', 'BASS2':r'BASS[-_ ]?2', 'MEP':r'MEP[0-9]?', 'MPC':r'MPC[1-4]?|MITOCHONDRIAL[-_ ]PYRUVATE[-_ ]CARRIER', 'PPA':r'PPA[1-9]?'}
patterns={f:re.compile(r'(?<![A-Z0-9])(?:'+p+r')(?![A-Z0-9])',re.I) for f,p in terms.items()}
def tags(s):return {f for f,p in patterns.items() if p.search(s)}
def exact_directory(s):return {f for f,p in terms.items() if re.fullmatch(p,s,re.I)}
questions={
'NAD-ME':'Confirm NAD-linked malic-enzyme scope and alpha/beta subunit treatment; distinguish NADP-ME and MDH.',
'NADP-ME':'Confirm NADP-linked malic-enzyme scope; distinguish NAD-ME and MDH; retain localization as separate evidence.',
'MDH':'Resolve whether Wang MDH scope includes NAD- and NADP-dependent forms; do not pool by acronym alone.',
'ASPAT':'Confirm aspartate aminotransferase definition and reference specificity against related aminotransferases.',
'AlaAT':'Confirm alanine aminotransferase definition and reference specificity against related aminotransferases.',
'PEPCK':'Confirm phosphoenolpyruvate carboxykinase definition; keep PEPC and PPCK distinct.',
'BASS2':'Confirm BASS2-specific scope and comparator coverage of other BASS family members.',
'MEP':'AMBIGUOUS_ABBREVIATION: recover explicit Wang/project definition and source-linked accessions before selecting positives.',
'MPC':'Confirm mitochondrial pyruvate carrier scope and component counting convention before pooling subunits.',
'PPA':'AMBIGUOUS_ABBREVIATION: recover explicit Wang/project definition; establish distinction from already reviewed V-PPase.'}
aliases=defaultdict(set);inputaudit=[];mappinghash={}
for mp in sorted((root/'03_harmonized/primary_v1').glob('*/*.all_isoform_mapping.tsv')):
 sp=mp.parent.name;mappinghash[sp]=sha(mp);mr=read(mp)
 for m in mr:
  for fld in ['gene_id','protein_id','transcript_id']:
   value=m.get(fld,'')
   if value and value not in {'.','NA'}:
    aliases[value].add((sp,m['gene_id']));aliases[sp+'__'+value].add((sp,m['gene_id']))
 inputaudit.append(dict(input_kind='ANNOTATION_MAPPING',species=sp,path=str(mp.relative_to(root)),status='PRESENT_HASHED',records=len(mr),sha256=mappinghash[sp]))
need(len(mappinghash)==8,'Expected eight annotation mapping files')
bygene=defaultdict(set)
for k in current:
 if len(aliases.get(k,set()))==1:bygene[next(iter(aliases[k]))].add(k)
manifest=root/'00_admin/manifests/input_sets.tsv';need(manifest.is_file(),'Missing input manifest')
for m in read(manifest):
 for fld in ['protein_all','cds_all','gff3','genome']:
  val=m.get(fld,'');p=root/val if val else None
  inputaudit.append(dict(input_kind=fld,species=m.get('input_set',''),path=(str(p.relative_to(root)) if p and p.is_relative_to(root) else '<EXTERNAL>/'+p.name if p else ''),status='PRESENT_NOT_HASHED' if p and p.is_file() else 'MISSING',records='',sha256='NOT_HASHED'))
inputaudit.append(dict(input_kind='CORRECTED_PROTEOME',species='ALL',path=str(proteinfile.relative_to(root)),status='PINNED_HASH_VERIFIED',records=len(current),sha256=expected))
table('00_input_availability.tsv',inputaudit,['input_kind','species','path','status','records','sha256'])
# Single metadata pass across project. Bound all content reads, copies and text excerpts.
inv=[];evidence=[];limits=[];fasta_jobs=[];scanned=0;content_bytes=0;copied=0;cap=500000;budget=1024**3;copybudget=60*1024**2
texttypes={'.tsv','.csv','.txt','.json','.yaml','.yml','.sh','.slurm','.py','.md'};fastatypes={'.faa','.fa','.fasta','.pep','.aa'}
def walkerror(e):limits.append(dict(path=str(e.filename).replace(str(root),'<PROJECT_ROOT>'),status='WALK_ERROR',detail=str(e).replace(str(root),'<PROJECT_ROOT>')))
print('Auditing ten families in one project scan',flush=True)
for folder,dirs,files in os.walk(root,followlinks=False,onerror=walkerror):
 dirs[:]=sorted(d for d in dirs if not d.startswith('.') and d not in {'node_modules','envs','miniforge3','miniconda3','__pycache__'} and (Path(folder)/d).resolve()!=out.parent)
 if scanned>=cap:limits.append(dict(path='PROJECT',status='FILE_ENTRY_CAP',detail=str(cap)));break
 for name in sorted(files):
  scanned+=1
  if scanned>cap:break
  p=Path(folder)/name;rel=str(p.relative_to(root));pathfamilies=tags(rel);basenamefamilies=tags(name)
  if p.is_symlink():
   if pathfamilies:limits.append(dict(path=rel,status='MATCHED_SYMLINK_NOT_FOLLOWED',detail=''))
   continue
  try:
   if not p.is_file():continue
   size=p.stat().st_size;suffix=p.suffix.lower();hint=pathfamilies
   if suffix in texttypes and size<=2*1024**2:
    if content_bytes+size<=budget:
     content_bytes+=size
     with p.open(errors='replace') as handle:
      for ln,line in enumerate(handle,1):
       tt=tags(line)
       if tt and len(evidence)<5000:
        evidence.append(dict(source_file=rel,line_number=ln,families=';'.join(sorted(tt)),text=line.rstrip()[:4000],interpretation='TEXT_MENTION_NOT_VALIDATED_DEFINITION'))
    else:limits.append(dict(path=rel,status='TEXT_BUDGET_SKIP',detail=str(size)))
   if not hint:continue
   copied_to='';digest=sha(p) if size<=10*1024**2 else 'NOT_HASHED_SIZE_LIMIT'
   if suffix in texttypes|fastatypes|{'.treefile','.nwk','.domtbl','.tbl','.log'} and size<=5*1024**2 and copied+size<=copybudget:
    copied_to='collected/F%05d_'%(len(inv)+1)+name;shutil.copyfile(p,out/copied_to);copied+=size
   # A shared directory token is insufficient to assign a FASTA to a family.
   parentfamilies=exact_directory(p.parent.name);assigned=basenamefamilies or parentfamilies
   isref=bool(re.search(r'reference|(?:^|[/_.-])refs?(?:[/_.-]|$)',rel,re.I))
   derived=bool(re.search(r'align|\.aln|codon|cds|nucleotide|genome|tree_input|plus_ref|with_ref',name,re.I))
   role='REFERENCE_FILE_UNVERIFIED' if isref else 'DERIVED_OR_MIXED_NOT_RECONCILED' if derived else 'CANDIDATE_FILE_UNVERIFIED'
   if suffix in fastatypes and size:
    if len(assigned)==1 and not derived and size<=10*1024**2:fasta_jobs.append((p,next(iter(assigned)),role))
    else:role='AMBIGUOUS_OR_DERIVED_FASTA_NOT_RECONCILED'
   inv.append(dict(source_file=rel,path_family_hints=';'.join(sorted(hint)),basename_family_hints=';'.join(sorted(basenamefamilies)),bytes=size,sha256=digest,file_role=role,bundle_file=copied_to))
  except (OSError,UnicodeError) as e:limits.append(dict(path=rel,status='READ_ERROR',detail=str(e).replace(str(root),'<PROJECT_ROOT>')))
if len(evidence)>=5000:limits.append(dict(path='PROJECT',status='TEXT_MATCH_CAP',detail='5000 excerpts; complete collected files remain available'))
table('01_existing_file_inventory.tsv',inv,['source_file','path_family_hints','basename_family_hints','bytes','sha256','file_role','bundle_file'])
table('02_definition_and_provenance_mentions.tsv',evidence,['source_file','line_number','families','text','interpretation'])
refrows=[];audit=[];ready=defaultdict(dict);parsed=0;parsebudget=100*1024**2
for p,fam,role in fasta_jobs:
 rel=str(p.relative_to(root));size=p.stat().st_size
 if parsed+size>parsebudget:limits.append(dict(path=rel,status='FASTA_PARSE_BUDGET_SKIP',detail=str(size)));continue
 parsed+=size
 try:d,headers=fasta(p)
 except (SystemExit,OSError) as e:limits.append(dict(path=rel,status='FASTA_PARSE_ERROR',detail=str(e)));continue
 for k,sequence in d.items():
  sequence=norm(sequence)
  if role=='REFERENCE_FILE_UNVERIFIED' or k.startswith(('REF__','sp|','tr|')):
   refrows.append(dict(wang_family=fam,source_file=rel,sequence_id=k,header=headers[k],length_aa=len(sequence),sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),sequence_format='PROTEIN_LIKE' if valid(sequence) else 'INVALID_OR_NUCLEOTIDE_LIKE',reference_role='UNVERIFIED_NO_AUTOMATIC_POSITIVE_ASSIGNMENT'));continue
  mm=aliases.get(k,set());corrected='';state='UNRESOLVED_IDENTIFIER'
  if not valid(sequence):state='INVALID_OR_NUCLEOTIDE_LIKE_SEQUENCE'
  elif len(mm)>1:state='AMBIGUOUS_LOCUS_MAPPING'
  elif len(mm)==1:
   poss=bygene[next(iter(mm))];exact={x for x in poss if norm(current[x])==sequence}
   if k in exact:corrected=k
   elif len(exact)==1:corrected=next(iter(exact))
   if corrected:state='EXACT_SEQUENCE_AND_DOCUMENTED_LOCUS';ready[fam][corrected]=current[corrected]
   elif len(exact)>1:state='AMBIGUOUS_CORRECTED_ISOFORM'
   else:state='SEQUENCE_MISMATCH' if poss else 'NO_CORRECTED_PROTEIN_FOR_LOCUS'
  sp,g=next(iter(mm)) if len(mm)==1 else ('','')
  audit.append(dict(wang_family=fam,source_file=rel,source_id=k,species=sp,gene_id=g,corrected_protein_id=corrected,status=state,interpretation='PROVENANCE_ONLY_FAMILY_MEMBERSHIP_UNVERIFIED'))
table('03_existing_reference_sequence_audit.tsv',refrows,['wang_family','source_file','sequence_id','header','length_aa','sequence_sha256','sequence_format','reference_role'])
table('04_candidate_input_reconciliation.tsv',audit,['wang_family','source_file','source_id','species','gene_id','corrected_protein_id','status','interpretation'])
summary=[];definitions=[]
for fam in families:
 ff=[x for x in inv if fam in x['path_family_hints'].split(';')];rr=[x for x in refrows if x['wang_family']==fam];aa=[x for x in audit if x['wang_family']==fam]
 state='NO_MATCHING_PATHS_IN_BOUNDED_SCAN' if not ff else 'EMPTY_MATCHES_ONLY' if not any(x['bytes'] for x in ff) else 'EXISTING_EVIDENCE_REQUIRES_SCOPE_REVIEW'
 gate='EXPLICIT_ABBREVIATION_DEFINITION_REQUIRED' if fam in {'MEP','PPA'} else 'REFERENCE_SCOPE_AND_COMPARATORS_REQUIRE_REVIEW'
 summary.append(dict(wang_family=fam,matching_paths=len(ff),nonempty_paths=sum(x['bytes']>0 for x in ff),reference_sequence_records=len(rr),unique_reference_sequences=len({x['sequence_sha256'] for x in rr}),candidate_source_records=len(aa),unique_reconciled_proteins=len(ready[fam]),definition_mentions=sum(fam in x['families'].split(';') for x in evidence),audit_status=state,blocking_review=gate,screen_ready=False,final_count='NA'))
 definitions.append(dict(wang_family=fam,retrieval_pattern=terms[fam],scope_question=questions[fam],Wang_definition_status='NOT_VERIFIED_FROM_ORIGINAL_METHODS',existing_reference_validation='PENDING_METADATA_AND_COMPARATOR_REVIEW',candidate_selection='EXACT_PROVENANCE_ONLY_NOT_FAMILY_ACCEPTANCE',batch_group='MALIC_ENZYMES' if fam in {'NAD-ME','NADP-ME'} else 'AMINOTRANSFERASES' if fam in {'ASPAT','AlaAT'} else 'FAMILY_SPECIFIC'))
 (out/'reconciled_candidates'/(fam+'.faa')).write_text(''.join('>'+k+'\n'+s+'\n' for k,s in sorted(ready[fam].items())))
table('05_ten_family_readiness.tsv',summary,list(summary[0]));table('06_reference_definition_review.tsv',definitions,list(definitions[0]));table('07_scan_limitations.tsv',limits,['path','status','detail'])
for name in ['01_integrated_locus_evidence.tsv','07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/name,out/('source_'+name))
coverage=read(src/'06_refreshed_family_coverage.tsv');bystatus={r['wang_family']:r for r in summary}
for r in coverage:
 if r['wang_family'] in bystatus:
  x=bystatus[r['wang_family']];r.update(current_status=x['audit_status'],next_action=x['blocking_review'],historical_matching_paths=str(x['matching_paths']),search_completeness='NOT_ESTABLISHED')
table('08_refreshed_family_coverage.tsv',coverage,list(coverage[0]))
(out/'README.txt').write_text('Ten-family input/reference-definition audit only; no search, domain scan, alignment, tree, gene-model repair or count promotion. Corrected proteome hash and count verified. Full project metadata pass<=500000 files; ordinary text<=2MiB with1GiB aggregate read budget and5000 bounded excerpts. Complete files<=5MiB collected under60MiB; other files inventoried. Raw FASTAs<=10MiB parsed within100MiB budget. Symlinks not followed. Exact family basename token or exact family directory required for FASTA attribution; shared-directory hints alone never assign candidates. Alignment/CDS/mixed FASTAs excluded from candidate reconciliation. Exact sequence plus documented locus mapping required; ambiguity retained. Reference headers and sequences inventoried, not promoted to trusted positive controls. Metadata in a filename is not proof of biological family. Generic filenames in complex family directories and text files>2MiB can be missed; review inventory/limitations. MEP/PPA definitions explicitly unresolved. Other definition questions are review prompts, not asserted Wang definitions. Raw reference sequence records may duplicate the same accession. Zero recovered inputs never means zero genes. Prior ledger/count state unchanged; final exhaustive counts pending. Batch compatible validation after one review of all ten readiness rows.\n')
qc=dict(status='TEN_FAMILY_AUDIT_READY_FOR_JOINT_REVIEW',families=len(families),files_examined=scanned,named_files=len(inv),content_bytes_read=content_bytes,parsed_fasta_bytes=parsed,scan_limitations=len(limits),family_summary=summary,integrated_records=len(read(src/'01_integrated_locus_evidence.tsv')),new_analysis_runs=0,final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');(out/'PROVENANCE.json').write_text(json.dumps(dict(step68_manifest_sha256=sha(src/'SHA256SUMS.txt'),corrected_manifest_sha256=sha(cor/'SHA256SUMS.txt'),corrected_proteome_sha256=expected,mapping_hashes=mappinghash,input_manifest_sha256=sha(manifest)),indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Remaining10_input_definition_audit/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle))

PY
