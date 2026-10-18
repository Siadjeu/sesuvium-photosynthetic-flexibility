#!/usr/bin/env bash
#SBATCH --job-name=ppt_provenance
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=ppt_provenance_%j.out
#SBATCH --error=ppt_provenance_%j.err
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





import urllib.parse
src=release('step25kz17.64','ALMT_CLOSE_RELEASE')
cor=release('step25kz17.21','CORRECTED_INPUT_RELEASE')
proteome=cor/'full_outputs/corrected_protein_candidates.faa'
expected='c5643d9ce8ba9140cd994eb1518fa7f9cdc242443b3e66532b9ba3f8f039b568'
need(proteome.is_file() and sha(proteome)==expected,'Corrected protein provenance mismatch')
ledger=read(src/'01_integrated_locus_evidence.tsv')
need(len(ledger)==1109 and not any(r['wang_family']=='PPT' for r in ledger),'Unexpected prior ledger')
out=base/'step25kz17.65'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False)
for d in ['collected','references']: (out/d).mkdir()
def table(n,rr,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
# Search metadata across the project, with bounded content reads. No directory symlinks.
pat=re.compile(r'(?<![A-Z0-9])PPT[12]?(?![A-Z0-9])|OG0000369',re.I)
inv=[];matches=[];issues=[];candidate_fastas=[];members=[];scanned=0;read_bytes=0;copied=0
max_entries=500000;content_budget=1024**3;copy_budget=40*1024**2
skip={'.git','__pycache__','.snakemake','.nextflow','node_modules','miniforge3','miniconda3','envs'}
def error(e):issues.append(dict(path=str(e.filename),status='WALK_ERROR',detail=str(e)))
print('Tracing PPT/OG0000369 paths and orthogroup membership',flush=True)
for folder,dirs,files in os.walk(root,followlinks=False,onerror=error):
 dirs[:]=sorted(d for d in dirs if not d.startswith('.') and d not in skip and (Path(folder)/d).resolve()!=out.parent)
 if scanned>=max_entries:
  issues.append(dict(path=str(Path(folder).relative_to(root)),status='ENTRY_CAP',detail=str(max_entries)));break
 for name in sorted(files):
  scanned+=1
  if scanned>max_entries:break
  p=Path(folder)/name;rel=str(p.relative_to(root))
  if p.is_symlink():continue
  try:
   if not p.is_file():continue
   size=p.stat().st_size;named=bool(pat.search(rel));suf=p.suffix.lower()
   texttype=suf in {'.tsv','.csv','.txt','.sh','.slurm','.py','.json','.yaml','.yml'}
   ogtable='orthogroups' in name.lower() and suf in {'.tsv','.txt','.csv'}
   content_ok=texttype and (size<=2*1024**2 or (ogtable and size<=512*1024**2))
   # Named text is preferred; other small text covers scripts/manifests and definitions.
   if content_ok and read_bytes+size<=content_budget:
    read_bytes+=size
    with p.open(errors='replace') as f:
     first=None
     for ln,line in enumerate(f,1):
      if first is None:first=line.rstrip('\n').split('\t')
      if re.search(r'OG0000369|PPT_OG|phosphoenolpyruvate.{0,3}phosphate',line,re.I):
       if len(matches)<10000:matches.append(dict(source_file=rel,line_number=ln,text=line.rstrip()[:6000]))
       if ogtable:
        cells=line.rstrip('\n').split('\t')
        if cells[0].rstrip(':')=='OG0000369' and len(cells)==len(first) and len(cells)>1:
         for sp,cell in zip(first[1:],cells[1:]):
          for k in re.split(r',\s*',cell):
           if k.strip():members.append(dict(source_file=rel,source_species=sp,source_id=k.strip(),sequence=''))
   elif content_ok:issues.append(dict(path=rel,status='CONTENT_BUDGET_SKIP',detail=str(size)))
   if not named:continue
   dest='';digest='NOT_HASHED_SIZE_LIMIT'
   if size<=10*1024**2:digest=sha(p)
   eligible=suf in {'.faa','.fa','.fasta','.tsv','.csv','.txt','.json','.sh','.slurm','.py','.log'}
   if eligible and size<=5*1024**2 and copied+size<=copy_budget:
    dest='collected/E%05d_'%(len(inv)+1)+name;shutil.copyfile(p,out/dest);copied+=size
   inv.append(dict(source_file=rel,bytes=size,sha256=digest,bundle_file=dest,status='EMPTY_FILE' if size==0 else ('COLLECTED' if dest else 'METADATA_ONLY')))
   # References are never treated as project candidates merely because named PPT.
   if size and size<=10*1024**2 and suf in {'.faa','.fa','.fasta'} and not re.search(r'reference|/refs?/',rel,re.I):
    candidate_fastas.append(p)
  except (OSError,UnicodeError) as e:issues.append(dict(path=rel,status='READ_ERROR',detail=str(e)))
for p in candidate_fastas:
 try:
  seqs,headers=fasta(p)
  for k,s in seqs.items():members.append(dict(source_file=str(p.relative_to(root)),source_species='',source_id=k,sequence=norm(s)))
 except (SystemExit,OSError) as e:issues.append(dict(path=str(p.relative_to(root)),status='FASTA_PARSE_ERROR',detail=str(e)))
table('01_project_PPT_inventory.tsv',inv,['source_file','bytes','sha256','bundle_file','status'])
table('02_provenance_text_matches.tsv',matches,['source_file','line_number','text'])
# Exact annotation aliases only; preserve collisions and all isoforms.
alias=defaultdict(set);mapping_hashes={}
for mp in sorted((root/'03_harmonized/primary_v1').glob('*/*.all_isoform_mapping.tsv')):
 sp=mp.parent.name;mapping_hashes[sp]=sha(mp)
 for m in read(mp):
  gene=m['gene_id']
  for field in ['gene_id','protein_id','transcript_id']:
   v=m.get(field,'')
   if v and v not in {'.','NA'}:
    alias[v].add((sp,gene));alias[sp+'__'+v].add((sp,gene))
need(len(mapping_hashes)==8,'Expected mappings for eight project species')
current,_=fasta(proteome);bygene=defaultdict(set)
for k in current:
 mm=alias.get(k,set())
 if len(mm)==1:bygene[next(iter(mm))].add(k)
reconciled=[];ready={}
for m in members:
 k=m['source_id'];hint=m['source_species'];mm=set(alias.get(k,set()))
 if hint in mapping_hashes:mm={x for x in mm if x[0]==hint} or set(alias.get(hint+'__'+k,set()))
 possible=set().union(*(bygene[x] for x in mm)) if mm else set()
 if k in current:possible={k} if k in possible else possible
 sequence=m['sequence'];exact={x for x in possible if norm(current[x])==sequence} if sequence else set()
 selected=next(iter(exact)) if len(exact)==1 else (next(iter(possible)) if len(possible)==1 and not sequence else '')
 status=('EXACT_SEQUENCE_AND_MAPPING' if sequence else 'DOCUMENTED_ID_MAPPING_ONLY') if selected and len(mm)==1 else ('AMBIGUOUS_MAPPING_OR_ISOFORM' if len(mm)>1 or len(possible)>1 else ('SEQUENCE_MISMATCH' if possible and sequence else 'UNRESOLVED_IDENTIFIER'))
 if selected and len(mm)==1:ready[selected]=current[selected]
 reconciled.append(dict(source_file=m['source_file'],source_species=hint,source_id=k,species=';'.join(sorted({x[0] for x in mm})),gene_id=';'.join(sorted({x[1] for x in mm})),corrected_protein_id=selected if len(mm)==1 else '',possible_corrected_ids=';'.join(sorted(possible)),status=status))
table('03_historical_candidate_reconciliation.tsv',reconciled,['source_file','source_species','source_id','species','gene_id','corrected_protein_id','possible_corrected_ids','status'])
(out/'04_reconciled_historical_candidates.faa').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in sorted(ready.items())))
# Retrieve independently annotated, reviewed reference records. No silent substitutions.
panel=[];ref_status=[]
for gene,role in [('PPT1','PPT_POSITIVE'),('PPT2','PPT_POSITIVE'),('TPT','RELATED_COMPARATOR'),('GPT1','RELATED_COMPARATOR'),('GPT2','RELATED_COMPARATOR'),('XPT','RELATED_COMPARATOR')]:
 query='(organism_id:3702) AND (reviewed:true) AND (gene_exact:'+gene+')'
 url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode({'query':query,'format':'json','size':100})
 try:
  request=urllib.request.Request(url,headers={'User-Agent':'PPT-provenance-review/1.0'})
  with urllib.request.urlopen(request,timeout=45) as response:
   payload=response.read();next_page=response.headers.get('Link','')
  data=json.loads(payload);(out/'references'/('query_'+gene+'.json')).write_bytes(payload)
  found=data.get('results',[]);need('rel="next"' not in next_page,'Unexpected reference pagination')
  selected=[]
  for entry in found:
   genes={x.get('geneName',{}).get('value','') for x in entry.get('genes',[])}
   desc=entry.get('proteinDescription',{}).get('recommendedName',{}).get('fullName',{}).get('value','')
   if gene not in genes or entry.get('organism',{}).get('taxonId')!=3702 or 'reviewed' not in entry.get('entryType','').lower() or 'unreviewed' in entry.get('entryType','').lower():continue
   if 'translocator' not in desc.lower():continue
   if role=='PPT_POSITIVE' and 'phosphoenolpyruvate' not in desc.lower():continue
   selected.append((entry,desc))
  if len(selected)!=1:raise ValueError('Expected one exact reviewed translocator record; found '+str(len(selected)))
  entry,desc=selected[0];acc=entry['primaryAccession'];seq=entry['sequence']['value'];need(valid(seq),'Invalid reference sequence')
  panel.append(dict(reference_id=acc,gene=gene,role=role,organism='Arabidopsis thaliana',taxon_id=3702,length_aa=len(seq),description=desc,source_url=url,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()))
  (out/'references'/(gene+'.faa')).write_text('>'+acc+' '+gene+' '+role+'\n'+seq+'\n')
  ref_status.append(dict(gene=gene,role=role,status='READY',detail=acc,url=url))
 except (OSError,ValueError,KeyError,SystemExit) as e:
  ref_status.append(dict(gene=gene,role=role,status='MISSING_OR_AMBIGUOUS',detail=str(e),url=url))
table('05_reference_panel.tsv',panel,['reference_id','gene','role','organism','taxon_id','length_aa','description','source_url','sequence_sha256'])
table('06_reference_retrieval_status.tsv',ref_status,['gene','role','status','detail','url'])
(out/'07_reference_panel.faa').write_text(''.join((out/'references'/(r['gene']+'.faa')).read_text() for r in panel))
# Carry ledgers unchanged; missing historical inputs never become zero counts.
shutil.copyfile(src/'01_integrated_locus_evidence.tsv',out/'source_01_integrated_locus_evidence.tsv')
shutil.copyfile(src/'07_final_Wang_counts_PENDING.tsv',out/'source_07_final_Wang_counts_PENDING.tsv')
coverage=read(src/'06_refreshed_family_coverage.tsv')
for r in coverage:
 if r['wang_family']=='PPT':r.update(current_status='HISTORICAL_CANDIDATES_RECONCILED_SCREEN_PENDING' if ready else 'NO_RECONCILED_HISTORICAL_CANDIDATES_SCREEN_PENDING',next_action='FULL_CORRECTED_PROTEOME_PPT_VERSUS_RELATED_TRANSLOCATORS_SCREEN_AFTER_PANEL_REVIEW',search_completeness='NOT_ESTABLISHED')
table('08_refreshed_family_coverage.tsv',coverage,list(coverage[0]))
table('09_scan_limitations.tsv',issues,['path','status','detail'])
qc=dict(status='PROVENANCE_REVIEW_READY',historical_source_records=len(members),unique_reconciled_proteins=len(ready),reconciliation_statuses=dict(Counter(r['status'] for r in reconciled)),named_files=len(inv),nonempty_named_files=sum(r['bytes']>0 for r in inv),project_files_examined=scanned,content_bytes_read=read_bytes,scan_limitations=len(issues),references_ready=len(panel),reference_panel_complete=len(panel)==6,screen_ready=len(panel)==6,integrated_records=len(ledger),final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step64_manifest_sha256=sha(src/'SHA256SUMS.txt'),corrected_manifest_sha256=sha(cor/'SHA256SUMS.txt'),corrected_proteome_sha256=expected,mapping_hashes=mapping_hashes,reference_retrieval_time_UTC=datetime.now(timezone.utc).isoformat(),definition_sources=['https://pubmed.ncbi.nlm.nih.gov/14617097/','https://pmc.ncbi.nlm.nih.gov/articles/PMC166879/']),indent=2)+'\n')
(out/'README.txt').write_text('PPT provenance and reference preparation only. No existing analysis rerun and no candidate promoted into counts. Project metadata scan capped at500000 files; content scan1GiB, ordinary text<=2MiB and orthogroup tables<=512MiB. Complete small named files copied within40MiB; path-only matches do not prove PPT identity. Symlink files/directories are not followed. Orthogroup recovery requires exact OG0000369 row in tabular Orthogroups files. Historical FASTA requires exact mapped identity and sequence agreement; ambiguous isoforms remain unresolved. Corrected-proteome full screen remains necessary even if no historical candidates recovered. Panel: reviewed Arabidopsis PPT1/PPT2 positives, TPT/GPT1/GPT2/XPT related comparators; retrieve JSON and verify taxon/exact primary gene/name. Missing references explicitly reported, never replaced automatically. Arabidopsis references are not project species AT. Six-reference panel is bounded, not exhaustive. Whole-project absence cannot be inferred from bounded scan. Source ledger unchanged; final counts pending.\n')
files=sorted(p for p in out.rglob('*') if p.is_file())
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='PPT_provenance_reference_preparation/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('RELEASE='+str(out));print('BUNDLE='+str(bundle))

PY
