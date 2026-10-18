#!/usr/bin/env bash
#SBATCH --job-name=family_coverage
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:30:00
#SBATCH --output=family_coverage_%j.out
#SBATCH --error=family_coverage_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from the project directory}"
python3 - <<'PY'
import csv,gzip,hashlib,json,os,re,tarfile
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families'
if not base.is_dir():raise SystemExit('Missing 08_focused_families')
out=base/'validation/step25kz17.17'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False)
families=['BETA_CA','PEPC','PEPC_KINASE','NAD_ME','NADP_ME','MDH','ASPAT','ALAAT','PEPCK','PPDK','PPDK_RP','BASS2','MEP','DIT2','DIC','DTC','MPC','ALMT','TDT','PPA','PPT','V_ATPASE','V_PPASE']
priority=['PPDK','PPDK_RP','PEPC_KINASE','V_ATPASE','V_PPASE','TDT_DTC_DIC','TDT','DTC','DIC']
species=['AT','BV','MC','PO','SLP_Hap1','SP_Hex','SP_Tet','SS']
manifest=[];issues=[];headers=[];statusrows=[];fastarows=[];rules=[];inventory=[];taken=set()
limit=200_000_000;total=0
def table(name,rows,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
def digest(b):return hashlib.sha256(b).hexdigest()
def relative(p):return str(p.relative_to(root))
def norm(s):return re.sub('[^A-Z0-9]+','_',s.upper())
def matches(path,fam):
 return bool(re.search(r'(^|_)'+re.escape(fam)+r'(_|$)',norm(path)))
def family_hint(path):
 return ';'.join(f for f in families+['TDT_DTC_DIC'] if matches(path,f))
def source_id(raw):
 # Recognize exact project identifiers and verified family/species wrappers.
 x=raw.split()[0] if raw else ''
 if '|' in x:
  wrapper,tail=x.split('|',1);sp=tail.split('__',1)[0]
  if wrapper==sp or wrapper.endswith('__'+sp):x=tail
 sp=x.split('__',1)[0]
 return (x,sp) if sp in species and '__' in x else (x,'UNRESOLVED_OR_REFERENCE')
def collect(p,reason):
 global total
 rel=relative(p)
 if rel in taken:return
 taken.add(rel)
 try:
  size=p.stat().st_size
  # Do not pull whole proteomes, databases or enormous historical hit tables.
  if size>8_000_000 or total+size>limit:
   manifest.append(dict(path=rel,bytes=size,sha256='',collection_status='SKIPPED_SIZE_LIMIT',reason=reason));return
  data=p.read_bytes();total+=len(data)
  # Text snapshots redact absolute cluster paths; originals stay unchanged.
  text=data.decode('utf-8')
  text=re.sub(r'/(?:dss|fs|home|scratch|work|lustre|gpfs)/[^\s\t"\x27]+','<ABS_PATH>',text)
  q=out/'sources'/rel;q.parent.mkdir(parents=True,exist_ok=True);q.write_text(text)
  manifest.append(dict(path=rel,bytes=len(data),sha256=digest(data),collection_status='COLLECTED_TEXT_REDACTED',reason=reason))
  lines=text.splitlines()
  if p.suffix in ['.tsv','.csv']:
   delimiter='\t' if p.suffix=='.tsv' else ','
   reader=csv.DictReader(lines,delimiter=delimiter);fields=reader.fieldnames or []
   headers.append(dict(path=rel,columns=' | '.join(fields)))
   keys=[k for k in fields if any(w in k.lower() for w in ['status','class','reason','provenance','flag'])]
   counts=Counter()
   for row in reader:
    for k in keys:counts[(k,str(row.get(k,'')))]+=1
   statusrows.extend(dict(path=rel,column=k,value=v,rows=n) for (k,v),n in sorted(counts.items()))
  if p.suffix in ['.fa','.faa','.fasta']:
   records=[];name=None;seq=[]
   for line in lines:
    if line.startswith('>'):
     if name is not None:records.append((name,''.join(seq)))
     name=line[1:];seq=[]
    elif name is not None:seq.append(line.strip())
   if name is not None:records.append((name,''.join(seq)))
   seen=Counter(h.split()[0] for h,s in records)
   for h,s in records:
    pid,sp=source_id(h)
    fastarows.append(dict(path=rel,source_header=h,canonical_id=pid,species=sp,
      length_aa=len(s),sequence_sha256=digest(s.encode()),id_occurrences=seen[h.split()[0]],
      role='REFERENCE_TO_CURATE' if 'reference' in rel.lower() else 'CANDIDATE_NOT_FINAL'))
  if p.suffix in ['.py','.sh','.slurm','.md','.txt','.yaml','.yml','.json']:
   pat=re.compile(r'length|min_len|max_len|coverage|qcov|identity|pident|evalue|e-value|bitscore|pfam|interpro|domain|threshold|PPDK|ATPASE|TDT|DTC|DIC',re.I)
   hit=[i for i,l in enumerate(lines) if pat.search(l)]
   context=set(j for i in hit for j in range(max(0,i-2),min(len(lines),i+3)))
   rules.extend(dict(path=rel,line_number=i+1,text=lines[i]) for i in sorted(context))
 except (OSError,UnicodeError,csv.Error) as e:
  issues.append(dict(path=rel,issue=type(e).__name__,detail=str(e).replace(str(root),'<PROJECT_ROOT>')))
# Complete checkpoints necessary to reconstruct inclusion/exclusion and reference scope.
for fam in priority:
 for suffix in ['candidate_classification.tsv','primary_candidates.faa','sensitivity_candidates.faa']:
  p=base/'validation/step25kz17.2/release'/fam/suffix
  if p.is_file():collect(p,'PRIMARY_SENSITIVITY_FILTER_AUDIT')
 for suffix in ['candidate_validation.tsv',fam+'.accepted.proteins.fa']:
  p=base/'validation/step25kn/results'/fam/suffix
  if p.is_file():collect(p,'LEGACY_VALIDATION_AUDIT')
 for ext in ['.faa','.fa','.fasta']:
  p=base/'validation/step25kn/references'/(fam+ext)
  if p.is_file():collect(p,'REFERENCE_SCOPE_AUDIT')
# Discover definitions and additional family evidence under focused-family module.
# Skip previous audit packages to avoid recursively collecting copied sources.
for folder,dirs,files in os.walk(base,followlinks=False):
 dirs[:]=sorted(d for d in dirs if d not in ['step25kz17.16','step25kz17.17','.git','sources','source_checkpoints'] and not d.startswith('.'))
 for name in sorted(files):
  p=Path(folder)/name
  if p.is_symlink():continue
  rel=relative(p);low=rel.lower();hint=family_hint(rel)
  if p.suffix.lower() not in ['.tsv','.csv','.fa','.faa','.fasta','.py','.sh','.slurm','.md','.txt','.yaml','.yml','.json']:continue
  iscode=p.suffix in ['.py','.sh','.slurm']
  relevant=bool(hint) or any(k in low for k in ['master_candidate','family_def','family_panel','reference_inventory','reference_manifest','family_config','family_registry','all_species_input_inventory'])
  if relevant:
   inventory.append(dict(path=rel,bytes=p.stat().st_size,family_name_hint=hint))
  nameddef=any(k in low for k in ['family_def','family_panel','family_config','reference_inventory','reference_manifest','family_registry'])
  focal=any('/'+step+'/' in rel for step in ['step25kn','step25kz17.1','step25kz17.2'])
  ref=('reference' in low and p.suffix in ['.faa','.fa','.fasta','.tsv'])
  # Select small complete definitions and focal scripts, rather than all historical searches.
  if nameddef or (focal and iscode) or (relevant and ref) or (
      relevant and any(k in low for k in ['master_candidate_inventory','all_species_input_inventory'])):
   collect(p,'DEFINITION_REFERENCE_OR_METHOD')
# Coverage here measures located evidence only, never biological absence.
coverage=[]
for fam in families:
 paths=[r for r in inventory if matches(r['path'],fam)]
 refs=[r for r in paths if 'reference' in r['path'].lower()]
 coverage.append(dict(family=fam,matching_paths=len(paths),reference_named_paths=len(refs),
   definition_status='PENDING',search_completeness='PENDING',
   next_action='INSPECT_LOCATED_EVIDENCE' if paths else 'LOCATE_EVIDENCE_OR_DESIGN_TARGETED_SEARCH'))
table('01_family_evidence_coverage.tsv',coverage,['family','matching_paths','reference_named_paths','definition_status','search_completeness','next_action'])
table('02_source_manifest.tsv',manifest,['path','bytes','sha256','collection_status','reason'])
table('03_table_schemas.tsv',headers,['path','columns'])
table('04_filter_status_summary.tsv',statusrows,['path','column','value','rows'])
table('05_sequence_and_reference_inventory.tsv',fastarows,['path','source_header','canonical_id','species','length_aa','sequence_sha256','id_occurrences','role'])
table('06_filter_rule_context.tsv',rules,['path','line_number','text'])
table('07_available_family_paths.tsv',inventory,['path','bytes','family_name_hint'])
table('08_collection_issues.tsv',issues,['path','issue','detail'])
review=[]
for r in fastarows:
 if r['role']=='REFERENCE_TO_CURATE':
  review.append(dict(source=r['path'],reference_header=r['source_header'],length_aa=r['length_aa'],
   proposed_Wang_category='PENDING',subunit_or_subfamily='PENDING',curation_evidence='',
   decision='PENDING'))
table('09_reference_definition_review.tsv',review,['source','reference_header','length_aa','proposed_Wang_category','subunit_or_subfamily','curation_evidence','decision'])
(out/'README.txt').write_text("""Family coverage and definition evidence audit.
No searches, alignments or trees are run. Existing files are read only.
01 is file-discovery coverage, NOT evidence that all genes were searched.
02 records original byte hashes. Source snapshots redact absolute cluster paths;
snapshot checksums are in SHA256SUMS.txt and therefore differ when redaction occurs.
04 summarizes existing filter/status labels without treating PASS as final membership.
05 retains distinct identifiers even when protein sequences are identical.
06 gives filter context; complete selected scripts are also included.
07 lists additional paths for follow-up; file-name hints are not family assignments.
09 is a reference curation template; no family identity is inferred from names alone.

Review order:
1. PPDK primary vs sensitivity vs legacy sets: identify exactly why candidates were lost.
2. V_ATPASE reference subunits and length filters: evaluate subunit-specific requirements.
3. PPDK_RP references and rejected candidates: do not interpret empty output as zero.
4. TDT_DTC_DIC references and candidates: split only with family-specific evidence.
5. Remaining Wang categories: distinguish missing search from uncollected evidence.
Do not force counts to match Wang. Final locus mapping and functional C4/CAM annotation
remain separate tasks. This audit uses the focused-family directory; absence of files
here does not exclude useful annotations elsewhere in the project.
Collection limits: 8 MB per source, 200 MB total. Skipped sources are explicit.
""")
qc={'status':'AUDIT_COMPLETE_WITH_COLLECTION_ISSUES' if issues or any(r['collection_status'].startswith('SKIPPED') for r in manifest) else 'AUDIT_COMPLETE',
 'source_files':len(manifest),'snapshot_bytes':total,'sequence_records':len(fastarows),
 'issues':len(issues),'analyses_rerun':False,'final_counts_produced':False}
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file())
(out/'SHA256SUMS.txt').write_text(''.join(digest(p.read_bytes())+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='Family_coverage_review')
print(json.dumps(qc,indent=2));print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
