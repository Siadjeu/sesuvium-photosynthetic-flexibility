#!/usr/bin/env bash
#SBATCH --job-name=proteome_audit
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=proteome_audit_%j.out
#SBATCH --error=proteome_audit_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from the project directory}"
python3 - <<'PY'
import csv,gzip,hashlib,itertools,json,os,re,tarfile
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
source=root/os.environ.get('PROTEOME','08_focused_families/validation/step25kz16.7/release/work/all_species_primary_proteins.faa')
if not source.is_file():raise SystemExit('Primary proteome not found')
out=base/'step25kz17.20'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False)
species={'AT','BV','MC','PO','SLP_Hap1','SP_Hex','SP_Tet','SS'}
aa=set('ACDEFGHIKLMNPQRSTVWYBXZJUO*')
code=dict(zip((''.join(x) for x in itertools.product('TCAG',repeat=3)),
 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def digest(s):return hashlib.sha256(s.encode()).hexdigest()
def fasta(p):
 h=None;seq=[]
 with p.open() as f:
  for line in f:
   if line.startswith('>'):
    if h is not None:yield h,''.join(seq)
    h=line[1:].strip();seq=[]
   elif h is not None:seq.append(''.join(line.split()))
   elif line.strip():raise ValueError('Sequence before first FASTA header')
 if h is not None:yield h,''.join(seq)
def trimstop(s):return s[:-1] if s.endswith('*') else s
def table(name,rows,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
  w.writeheader();w.writerows(rows)
rows=[];suspects={};index=defaultdict(list);idrows=defaultdict(list);counts=Counter()
for ordinal,(header,original) in enumerate(fasta(source),1):
 if not header:raise SystemExit('Empty FASTA header')
 pid=header.split()[0];sp=pid.split('__',1)[0];s=original.upper()
 chars=set(s);nt=sum(s.count(x) for x in 'ACGTN')/len(s) if s else 0
 flags=[]
 if not s:flags.append('EMPTY_SEQUENCE')
 if sp not in species or '__' not in pid:flags.append('UNEXPECTED_ID_OR_SPECIES')
 if chars-set(aa):flags.append('NON_PROTEIN_CHARACTERS')
 if '*' in trimstop(s):flags.append('INTERNAL_STOP')
 if s.endswith('*'):flags.append('TERMINAL_STOP')
 if s and chars<=set('ACGTN'):flags.append('NUCLEOTIDE_ALPHABET_ONLY')
 elif len(s)>=90 and nt>=.95:flags.append('HIGH_NUCLEOTIDE_CHARACTER_FRACTION')
 # Composition is a flag, not proof: some genuine short proteins can use this alphabet.
 suspected=any(x in flags for x in ['NUCLEOTIDE_ALPHABET_ONLY','HIGH_NUCLEOTIDE_CHARACTER_FRACTION'])
 row=dict(record_number=ordinal,record_id=pid,species=sp,length=len(s),
  alphabet=''.join(sorted(chars)),ACGTN_fraction=round(nt,6),
  internal_stops=trimstop(s).count('*'),sequence_sha256=digest(s),flags=';'.join(flags) or 'NO_FLAGS')
 rows.append(row);idrows[pid].append(ordinal)
 if suspected:suspects[ordinal]=(pid,sp,s)
 elif s and not chars-set(aa) and '*' not in trimstop(s) and sp in species:
  index[(sp,digest(trimstop(s)))].append((ordinal,pid,trimstop(s)))
 counts[(sp,'TOTAL')]+=1
 for flag in flags:counts[(sp,flag)]+=1
table('01_all_record_audit.tsv',rows,list(rows[0]) if rows else ['record_number','record_id','flags'])
table('02_species_flag_summary.tsv',[dict(species=s,metric=m,records=n) for (s,m),n in sorted(counts.items())],['species','metric','records'])
table('03_duplicate_identifiers.tsv',[dict(record_id=p,record_numbers=';'.join(map(str,v)),occurrences=len(v)) for p,v in sorted(idrows.items()) if len(v)>1],['record_id','record_numbers','occurrences'])
# Frame +1 standard-code translation only. No translation is inserted into the input.
# Same-species matching avoids equating identical paralogs with the same locus.
matches=[];susp_ids=set();matching_ids=set()
for ordinal,(pid,sp,s) in suspects.items():
 susp_ids.add(pid);status='NOT_TRANSLATABLE_AS_FRAME1_CDS';translated=''
 if len(s)%3==0 and set(s)<=set('ACGTN'):
  translated=''.join(code.get(s[i:i+3],'X') for i in range(0,len(s),3))
  core=trimstop(translated)
  found=[x for x in index.get((sp,digest(core)),[]) if x[2]==core] if '*' not in core else []
  status='EXACT_TRANSLATION_MATCH' if found else ('INTERNAL_STOP_IN_FRAME1' if '*' in core else 'NO_EXACT_PROTEIN_MATCH')
 else:found=[]
 for recnum,other,seq in found:matching_ids.add(other)
 matches.append(dict(record_number=ordinal,suspect_id=pid,species=sp,nucleotide_like_length=len(s),
  translation_length=len(trimstop(translated)),translation_internal_stops=trimstop(translated).count('*'),
  translation_status=status,matching_protein_ids=';'.join(x[1] for x in found),
  matching_record_numbers=';'.join(str(x[0]) for x in found),
  same_stem_matches=';'.join(x[1] for x in found if re.sub(r'\.t\d+$','',x[1])==re.sub(r'\.t\d+$','',pid)),
  locus_status='REQUIRES_SOURCE_OR_ANNOTATION_MAPPING'))
table('04_translation_matches.tsv',matches,['record_number','suspect_id','species','nucleotide_like_length','translation_length','translation_internal_stops','translation_status','matching_protein_ids','matching_record_numbers','same_stem_matches','locus_status'])
# Export only a small illustrative sequence set; complete record and match tables are retained.
sample_ord={x['record_number'] for x in matches[:20]}
for x in matches:
 if x['suspect_id']=='SS__g2993':sample_ord.add(x['record_number'])
sample_ids={pid for row in matches if row['record_number'] in sample_ord for pid in row['matching_protein_ids'].split(';') if pid}
with (out/'05_example_records.fasta').open('w') as f:
 for i,(h,s) in enumerate(fasta(source),1):
  if i in sample_ord or h.split()[0] in sample_ids:f.write('>'+str(i)+'|'+h.split()[0]+'\n'+s+'\n')
# Metadata source tracing: exact cell matches only, retaining original field names.
metadata=base/'step25kz16.7/release/01_all_species_protein_metadata.tsv'
meta_hits=[];meta_fields=[];meta_coverage=set();metadata_status='MISSING'
if metadata.is_file():
 with metadata.open() as f:
  reader=csv.DictReader(f,delimiter='\t');meta_fields=reader.fieldnames or [];metadata_status='READ'
  targets=susp_ids|matching_ids
  for row in reader:
   hit=targets.intersection(v for v in row.values() if isinstance(v,str))
   if hit:
    meta_coverage.update(hit)
    meta_hits.append(dict(matched_ids=';'.join(sorted(hit)),**row))
 table('06_matching_source_metadata.tsv',meta_hits,['matched_ids']+meta_fields)
table('07_metadata_match_coverage.tsv',[dict(record_id=p,metadata_exact_cell_match='YES' if p in meta_coverage else 'NO') for p in sorted(susp_ids|matching_ids)],['record_id','metadata_exact_cell_match'])
# Collect original input inventory and producer scripts; no recursive project scan.
provenance=[dict(source='PRIMARY_PROTEOME',relative_path=str(source.relative_to(root)) if source.is_relative_to(root) else '<EXTERNAL_INPUT>',sha256=sha(source))]
collect=[metadata,base/'step25kz16.7/release/00_all_species_input_inventory.tsv']
collect+=sorted((base/'step25kz16.7/scripts').glob('*'))
def redact(t):return re.sub(r'/(?:dss|fs|home|scratch|work|lustre|gpfs)/[^\s\t"\x27]+','<ABS_PATH>',t)
for p in collect:
 if not p.is_file():continue
 provenance.append(dict(source=p.name,relative_path=str(p.relative_to(root)),sha256=sha(p)))
 if p==metadata:continue # filtered rows above provide complete relevant records
 if p.stat().st_size<=2_000_000:
  q=out/'source_context'/p.relative_to(root);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(redact(p.read_text()))
# Redact machine-specific paths in extracted metadata while retaining other source attributes.
p=out/'06_matching_source_metadata.tsv'
if p.exists():p.write_text(redact(p.read_text()))
table('08_input_provenance.tsv',provenance,['source','relative_path','sha256'])
qc=dict(status='AUDIT_COMPLETE_REVIEW_REQUIRED',total_records=len(rows),nucleotide_like_records=len(suspects),
 exact_translation_matches=sum(x['translation_status']=='EXACT_TRANSLATION_MATCH' for x in matches),
 duplicate_ids=sum(len(v)>1 for v in idrows.values()),metadata_status=metadata_status,
 suspect_ids_with_metadata=sum(p in meta_coverage for p in susp_ids),
 original_input_modified=False,corrected_proteome_created=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'README.txt').write_text("""Primary-proteome integrity audit. Existing files are never modified.
01 reports every record, including ordinal indices to preserve duplicate IDs.
Nucleotide alphabet/composition flags are screening criteria, not universal proof
of DNA: amino-acid sequences can use the same letters. No records are deleted.
04 tests frame +1, standard genetic code, allowing N -> X and one terminal stop.
Exact translation matches support nucleotide/protein duplication; they do not
prove a unique locus mapping when multiple proteins share a sequence.
Same-stem matches are hints only. No automatic suffix stripping for final counts.
No match does not establish that the sequence is a protein or unusable CDS.
06 traces exact metadata-cell ID matches. Missing matches are explicitly in 07.
Producer scripts and input inventory are included where present. Source sequence
files are not automatically opened or copied based on metadata paths.
Absolute cluster paths in context/metadata are redacted; originals stay unchanged.
No corrected proteome is released until source provenance and locus mappings
are reviewed. Do not resume family-count finalization from composition alone.
All audit tables are included; only illustrative FASTA records are bundled.
""")
files=sorted(p for p in out.rglob('*') if p.is_file())
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='Primary_proteome_audit')
print(json.dumps(qc,indent=2));print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
