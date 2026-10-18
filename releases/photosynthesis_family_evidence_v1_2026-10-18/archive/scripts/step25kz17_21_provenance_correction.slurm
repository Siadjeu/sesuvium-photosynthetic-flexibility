#!/usr/bin/env bash
#SBATCH --job-name=protein_correct
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=protein_correct_%j.out
#SBATCH --error=protein_correct_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from project directory}"
python3 - <<'PY'
import csv,hashlib,json,os,re,tarfile
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
choices=[root/os.environ['AUDIT_RELEASE']] if os.environ.get('AUDIT_RELEASE') else sorted(p for p in (base/'step25kz17.20').glob('release_*') if p.is_dir() and (p/'QC.json').is_file())
if len(choices)!=1:raise SystemExit('Set AUDIT_RELEASE to intended step17.20 release')
audit=choices[0]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
for l in (audit/'SHA256SUMS.txt').read_text().splitlines():
 h,n=l.split(maxsplit=1)
 if sha(audit/n)!=h:raise SystemExit('Audit checksum mismatch: '+n)
old=base/'step25kz16.7/release/work/all_species_primary_proteins.faa'
meta=base/'step25kz16.7/release/01_all_species_protein_metadata.tsv'
for label,p in [('PRIMARY_PROTEOME',old),('01_all_species_protein_metadata.tsv',meta)]:
 expected=[r['sha256'] for r in read(audit/'08_input_provenance.tsv') if r['source']==label]
 if len(expected)!=1 or sha(p)!=expected[0]:raise SystemExit('Input changed since audit: '+label)
out=base/'step25kz17.21'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'full_outputs').mkdir()
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rr)
def fasta(p):
 h=None;s=[]
 with p.open() as f:
  for l in f:
   if l.startswith('>'):
    if h is not None:yield h,''.join(s)
    h=l[1:].strip();s=[]
   elif h is not None:s.append(''.join(l.split()))
 if h is not None:yield h,''.join(s)
oldseq={}
for h,s in fasta(old):
 pid=h.split()[0]
 if pid in oldseq:raise SystemExit('Duplicate combined ID')
 oldseq[pid]=s
metadata=read(meta)
if len({r['tree_id'] for r in metadata})!=len(metadata) or {r['tree_id'] for r in metadata}!=set(oldseq):raise SystemExit('Metadata/FASTA IDs do not correspond one-to-one')
badsource='06_positive_selection/sequences/positive_selection_sequence_recovery_primary_v1/species_sequences/SS.candidate.cds.fa'
removed=[r for r in metadata if r['source_file']==badsource]
if len(removed)!=1073 or any(r['species']!='SS' for r in removed):raise SystemExit('Confirmed CDS source inventory changed; expected 1073 SS records')
keep=[r for r in metadata if r['source_file']!=badsource]
groups=defaultdict(list)
for r in keep:groups[r['source_file']].append(r)
sourcehash=[];ledger=[];problems=[];stats=Counter()
aa=set('ACDEFGHIKLMNPQRSTVWYBXZJUO')
rawout=out/'full_outputs/source_verified_records.fasta'
cleanout=out/'full_outputs/corrected_protein_candidates.faa'
qout=out/'full_outputs/quarantined_protein_source_records.fasta'
with rawout.open('w') as raw,cleanout.open('w') as clean,qout.open('w') as quarantine:
 for source_name,rr in sorted(groups.items()):
  p=Path(source_name)
  if p.is_absolute() or '..' in p.parts:raise SystemExit('Expected project-relative source path')
  p=root/p
  if not p.is_file():raise SystemExit('Missing source protein file: '+source_name)
  sourcehash.append(dict(source_file=source_name,sha256=sha(p)))
  wanted={r['original_id'] for r in rr};found={}
  for h,s in fasta(p):
   pid=h.split()[0]
   if pid in wanted:
    if pid in found:raise SystemExit('Duplicate source ID: '+pid)
    found[pid]=(h,s)
  if wanted-set(found):raise SystemExit('Metadata IDs absent from source: '+source_name)
  for r in rr:
   pid=r['tree_id'];h,s=found[r['original_id']];prior=oldseq[pid]
   # Reproduce historical normalization ONLY to verify provenance.
   reproduced=re.sub('[^A-Za-z]','',s)
   if reproduced!=prior:raise SystemExit('Source does not reproduce old record: '+pid)
   raw.write('>'+pid+'\n'+s+'\n')
   normalized=s.upper()
   terminal=normalized.endswith('*')
   if terminal:normalized=normalized[:-1] # one terminal stop only; internal stops are retained in raw/quarantine
   flags=[]
   if not normalized:flags.append('EMPTY')
   if '*' in normalized:flags.append('INTERNAL_STOP')
   if set(normalized)-aa:flags.append('INVALID_PROTEIN_SYMBOLS')
   frac=sum(normalized.count(x) for x in 'ACGTN')/len(normalized) if normalized else 0
   composition=bool(normalized and (set(normalized)<=set('ACGTN') or (len(normalized)>=90 and frac>=.95)))
   # AT composition flags retained explicitly for review, not classified as CDS.
   action='QUARANTINE' if flags else 'RETAIN'
   if action=='QUARANTINE':
    quarantine.write('>'+pid+'\n'+s+'\n');problems.append(dict(record_id=pid,reason=';'.join(flags)))
   else:clean.write('>'+pid+'\n'+normalized+'\n')
   stats[(r['species'],action)]+=1
   if composition:stats[(r['species'],'COMPOSITION_REVIEW')]+=1
   ledger.append(dict(record_id=pid,species=r['species'],original_id=r['original_id'],source_file=source_name,
    old_length=len(prior),source_length=len(s),output_length=len(normalized),
    action=action,terminal_stop_removed='YES' if terminal else 'NO',
    uppercased='YES' if s!=s.upper() else 'NO',composition_review='YES' if composition else 'NO',
    flags=';'.join(flags),same_as_old_protein='YES' if normalized==prior else 'NO',
    source_sequence_sha256=hashlib.sha256(s.encode()).hexdigest(),
    output_sequence_sha256=hashlib.sha256(normalized.encode()).hexdigest()))
with (out/'full_outputs/removed_confirmed_CDS.fna').open('w') as f:
 for r in removed:f.write('>'+r['tree_id']+'\n'+oldseq[r['tree_id']]+'\n')
table('01_removed_CDS_provenance.tsv',removed,list(metadata[0]))
table('02_retained_record_ledger.tsv',ledger,list(ledger[0]))
table('03_species_summary.tsv',[dict(species=s,action=a,records=n) for (s,a),n in sorted(stats.items())],['species','action','records'])
table('04_source_protein_checksums.tsv',sourcehash,['source_file','sha256'])
table('05_quarantine_reasons.tsv',problems,['record_id','reason'])
table('06_composition_review.tsv',[r for r in ledger if r['composition_review']=='YES'],list(ledger[0]))
# Keep unmatched CDS and its original audit evidence available without inserting translations.
trans=read(audit/'04_translation_matches.tsv');rmids={r['tree_id'] for r in removed}
unmatched=[r for r in trans if r['suspect_id'] in rmids and r['translation_status']!='EXACT_TRANSLATION_MATCH']
table('07_unmatched_CDS_review.tsv',unmatched,list(trans[0]))
with (out/'07_unmatched_CDS.fna').open('w') as f:
 for r in unmatched:f.write('>'+r['suspect_id']+'\n'+oldseq[r['suspect_id']]+'\n')
# Provenance-verified source count is distinct from number passing protein-symbol QC.
assert len(keep)==287197
retained=sum(r['action']=='RETAIN' for r in ledger)
assert retained+len(problems)==len(keep)
outputids=[h.split()[0] for h,s in fasta(cleanout)]
assert len(outputids)==retained and len(set(outputids))==retained and not set(outputids)&rmids
table('08_full_output_checksums.tsv',[dict(file=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted((out/'full_outputs').iterdir())],['file','bytes','sha256'])
qc=dict(status='PROVENANCE_CORRECTION_COMPLETE_REVIEW_REQUIRED',input_records=len(metadata),
 removed_confirmed_CDS=len(removed),source_verified_protein_records=len(keep),
 retained_protein_candidates=retained,quarantined=len(problems),
 composition_flagged=sum(r['composition_review']=='YES' for r in ledger),
 retained_sequences_different_from_old=sum(r['action']=='RETAIN' and r['same_as_old_protein']=='NO' for r in ledger),
 final_gene_count=False,approved_for_downstream_searches=False,old_inputs_modified=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'README.txt').write_text("""Provenance-based correction of step16.7 combined input.
Only records attributed to the confirmed SS candidate.cds.fa source are removed.
Every retained original identifier is looked up in its metadata-declared protein file.
The original sequence must reproduce the old record under its historical letter-only
normalization, or the script stops. Source checksums record the inputs used.
source_verified_records.fasta preserves source sequence characters.
corrected_protein_candidates.faa uppercases and removes at most one terminal stop.
Internal stops, invalid symbols and empty proteins are quarantined, never silently
stripped. The complete ledger records each action and sequence checksum.
Composition flags, including AT N-rich proteins, are retained for explicit review.
No identical protein sequences or transcript suffixes are collapsed. No CDS translation
is added. 287197 is a source record count, not an approved nonredundant locus count.
Unmatched CDS remains separate. Existing pipelines and inputs are not overwritten.
Do not automatically launch searches: review QC, quarantine, AT flags and source
selection before approval. Source provenance is not validation of biological annotation.
Full FASTAs stay in full_outputs on cluster. Their hashes are in 08; upload bundle
contains review tables and unmatched CDS but omits the large full_outputs directory.
""")
files=sorted(p for p in out.rglob('*') if p.is_file() and 'full_outputs' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Provenance_correction/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
