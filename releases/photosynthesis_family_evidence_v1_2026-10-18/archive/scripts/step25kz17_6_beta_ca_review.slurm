#!/usr/bin/env bash
#SBATCH --job-name=beta_ca_review
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=beta_ca_review_%j.out
#SBATCH --error=beta_ca_review_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
: "${REVIEW_TSV:?Export REVIEW_TSV pointing to the corrected manual-review TSV}"
command -v python3 >/dev/null
command -v mafft >/dev/null
export BETA_SOURCE="$PROJECT_ROOT/08_focused_families/validation/step25kz17.5/release/BETA_CA"
export BETA_REF="${BETA_REF:-$PROJECT_ROOT/08_focused_families/validation/step25kn/references/BETA_CA.faa}"
export BETA_OUT="$PROJECT_ROOT/08_focused_families/validation/step25kz17.6/release_${SLURM_JOB_ID:-manual}_$(date -u +%Y%m%dT%H%M%SZ)"
python3 - <<'PY'
import os,csv,json,hashlib
from pathlib import Path
from collections import Counter
source=Path(os.environ['BETA_SOURCE'])
review=Path(os.environ['REVIEW_TSV'])
ref=Path(os.environ['BETA_REF'])
out=Path(os.environ['BETA_OUT'])
paths={'classification':source/'05_BETA_CA_candidate_classification.tsv',
       'hmmer':source/'03_BETA_CA_HMMER_hits.tbl',
       'candidates':source/'04_BETA_CA_rescue_candidates.faa',
       'review':review,'references':ref}
for k,p in paths.items():
 if not p.is_file() or p.stat().st_size==0: raise SystemExit('Missing or empty input: '+k)
 if p.suffix.lower()=='.xlsx': raise SystemExit('Use the corrected TSV, not an XLSX workbook.')
def load_tsv(p):
 with p.open(encoding='utf-8-sig',newline='') as f:
  r=csv.DictReader(f,delimiter='\t')
  if not r.fieldnames or 'protein_id' not in r.fieldnames: raise SystemExit('Invalid TSV header: protein_id required')
  rows=list(r)
  if any(None in x for x in rows):raise SystemExit('Malformed TSV row')
  ids=[r['protein_id'] for r in rows]
  if len(ids)!=len(set(ids)) or any(not x for x in ids):raise SystemExit('Duplicate or empty protein identifiers')
  return r.fieldnames,rows
def fasta(p):
 records={};key=None;chunks=[]
 def save():
  if key is None:return
  seq=''.join(chunks).upper()
  if not seq:raise SystemExit('Empty FASTA sequence: '+key)
  records[key]=seq
 with p.open() as f:
  for line in f:
   if line.startswith('>'):
    save();key=line[1:].split()[0]
    if key in records:raise SystemExit('Duplicate FASTA ID: '+key)
    chunks=[]
   elif line.strip():
    if key is None:raise SystemExit('Sequence before FASTA header')
    chunks.append(''.join(line.split()))
 save();return records
def write_table(p,fields,rows):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def write_fasta(p,records):
 with p.open('w') as f:
  for key,seq in records:
   f.write('>'+key+'\n')
   for start in range(0,len(seq),80):f.write(seq[start:start+80]+'\n')
fields,base=load_tsv(paths['classification'])
_,review_rows=load_tsv(review)
reviews={r['protein_id']:r for r in review_rows}
seqs=fasta(paths['candidates']);refs=fasta(ref)
ids={r['protein_id'] for r in base}
if ids!=set(reviews):raise SystemExit('Review/classification identifier sets differ; stopping without output.')
if ids!=set(seqs):raise SystemExit('FASTA/classification identifier sets differ; stopping without output.')
hits={}
for line in paths['hmmer'].read_text().splitlines():
 if not line.strip() or line.startswith('#'):continue
 t=line.split()
 if len(t)<18:raise SystemExit('Unexpected HMMER tblout format')
 # HMMER --tblout full-sequence E-value, bit score, bias.
 try:h={'evalue':float(t[4]),'score':float(t[5]),'bias':float(t[6]),'query':t[2]}
 except ValueError:raise SystemExit('Non-numeric HMMER full-sequence fields')
 if t[0] in hits:raise SystemExit('Multiple HMMER hits for one target; explicit multi-profile policy required')
 hits[t[0]]=h
if not ids.issubset(hits):raise SystemExit('One or more candidates lack a raw HMMER record')
allowed={'PENDING','ACCEPT_FORMAL_BETA_CA','REVIEW_MANUAL','EXCLUDE_NON_BETA_CA','EXCLUDE_FRAGMENT'}
review_fields=['reviewer_decision','catalytic_region_check','truncation_check','reviewer_reason']
extra=['hmmer_full_sequence_bias','hmmer_review_flag','sequence_review_flag','formal_export_status','review_action']
fields=fields+ [k for k in review_fields+extra if k not in fields]
result=[];focus=[];selected=[];eligible=[]
for row in base:
 pid=row['protein_id'];rv=reviews[pid];seq=seqs[pid];h=hits[pid]
 if len(seq)!=int(row['length_aa']):raise SystemExit('Sequence length mismatch: '+pid)
 for k in review_fields:row[k]=(rv.get(k) or '').strip()
 decision=row['reviewer_decision'] or 'PENDING'
 # Also support the original worksheet export with decision labels placed in reason.
 if decision=='PENDING' and row['reviewer_reason'] in allowed-{'PENDING'}:
  decision=row['reviewer_reason'];row['reviewer_reason']='Reviewer label moved from reviewer_reason to reviewer_decision.'
 if decision not in allowed:raise SystemExit('Unknown decision for '+pid+': '+decision)
 row['reviewer_decision']=decision
 for k in ['catalytic_region_check','truncation_check']:
  row[k]=row[k] or 'PENDING'
 row['hmmer_hit']='YES';row['hmmer_score']=format(h['score'],'.8g');row['hmmer_evalue']=format(h['evalue'],'.8g')
 row['hmmer_full_sequence_bias']=format(h['bias'],'.8g')
 # A review flag only; no automatic family acceptance or exclusion.
 weak=h['evalue']>1e-5
 row['hmmer_review_flag']='REVIEW_EVALUE_GT_1E-5' if weak else 'EVALUE_LE_1E-5_NOT_FORMAL_VALIDATION'
 flags=[]
 if '*' in seq:flags.append('STOP_SYMBOL_PRESENT')
 if '-' in seq or '.' in seq:flags.append('GAP_SYMBOL_IN_SOURCE')
 if 'X' in seq:flags.append('AMBIGUOUS_X_PRESENT')
 if any(c not in 'ACDEFGHIKLMNPQRSTVWYBXZJUO*-.' for c in seq):flags.append('UNEXPECTED_SEQUENCE_SYMBOL')
 row['sequence_review_flag']=';'.join(flags) or 'NO_BASIC_SEQUENCE_FLAG'
 actions=[]
 if decision=='PENDING':actions.append('RECORD_DECISION')
 if decision=='REVIEW_MANUAL':actions.append('RESOLVE_MANUAL_REVIEW')
 if decision=='ACCEPT_FORMAL_BETA_CA':
  selected.append(pid)
  if row['catalytic_region_check']!='PASS':actions.append('CATALYTIC_REGION_CHECK_REQUIRED')
  if row['truncation_check']!='PASS':actions.append('TRUNCATION_CHECK_REQUIRED')
  if weak:actions.append('RESOLVE_WEAK_HMMER_SUPPORT')
  if flags:actions.append('RESOLVE_SEQUENCE_FLAGS')
 excluded=decision.startswith('EXCLUDE_')
 ready=decision=='ACCEPT_FORMAL_BETA_CA' and not actions
 row['formal_export_status']='ELIGIBLE_FOR_SEPARATE_EXPORT' if ready else ('EXCLUDED_BY_REVIEWER' if excluded else 'WITHHELD_PENDING_REVIEW')
 if ready:eligible.append(pid)
 row['review_action']=';'.join(actions) or ('NO_EXPORT_REVIEWER_EXCLUSION' if excluded else 'READY_FOR_SEPARATE_EXPORT')
 if (actions or weak or flags) and not excluded:focus.append(pid)
 result.append(row)
out.mkdir(parents=True,exist_ok=False)
write_table(out/'07_BETA_CA_manual_review_corrected.tsv',fields,result)
write_table(out/'08_BETA_CA_review_queue.tsv',fields,[r for r in result if r['protein_id'] in set(focus)])
write_fasta(out/'09_BETA_CA_all_candidates_with_references.faa',
 [(pid,seqs[pid]) for pid in sorted(seqs)]+ [('REF__'+pid,seq) for pid,seq in sorted(refs.items())])
# Keep the five historically unresolved/weak candidates easy to inspect.
priority={'BV__Bv4_100840_nphg','SLP_Hap1__SP10G01360','SP_Hex__gene_Spor0152910.1','SS__g13232.t1','PO__evm.TU.LG23.56'}
write_fasta(out/'10_BETA_CA_priority_candidates_with_references.faa',
 [(pid,seqs[pid]) for pid in sorted(priority & ids)]+ [('REF__'+pid,seq) for pid,seq in sorted(refs.items())])
write_table(out/'11_reference_id_map.tsv',['alignment_id','original_reference_id'],
 [{'alignment_id':'REF__'+pid,'original_reference_id':pid} for pid in sorted(refs)])
qc={'status':'REVIEW_PREPARATION_COMPLETE_NO_FORMAL_EXPORT','candidates':len(ids),'references':len(refs),
 'reviewer_selected':len(selected),'pending_decisions':sum(r['reviewer_decision']=='PENDING' for r in result),
 'eligible_for_separate_export':len(eligible),'hmmer_evalue_gt_1e_5':sum(hits[i]['evalue']>1e-5 for i in ids),
 'rescue_status_policy':'Original labels retained as historical screen labels, not recomputed formal evidence.',
 'acceptance_policy':'No formal FASTA or tree is exported by this preparation job.',
 'check_instructions':'Use PASS only after actually inspecting the catalytic region and truncation. Otherwise retain PENDING or enter FAIL.',
 'hmmer_flag_policy':'1e-5 is a diagnostic review threshold, not a validated universal family cutoff.',
 'sequence_policy':'Candidate and reference sequences preserved; no automatic trimming, stop removal or deduplication.',
 'input_sha256':{k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()}}
(out/'REVIEW_QC.json').write_text(json.dumps(qc,indent=2)+'\n')
print('Candidate records:',len(ids));print('Reviewer selections:',len(selected));print('Pending decisions:',qc['pending_decisions'])
print('Formal export performed: NO')
PY
mafft --thread "${SLURM_CPUS_PER_TASK:-4}" --auto \
  "$BETA_OUT/09_BETA_CA_all_candidates_with_references.faa" \
  > "$BETA_OUT/12_BETA_CA_review_alignment.faa" \
  2> "$BETA_OUT/12_mafft.log"
python3 - <<'PY'
import os,hashlib,json
from pathlib import Path
out=Path(os.environ['BETA_OUT'])
p=out/'12_BETA_CA_review_alignment.faa'
seqs={};key=None
for l in p.read_text().splitlines():
 if l.startswith('>'):
  key=l[1:].split()[0]
  if key in seqs:raise SystemExit('Duplicate alignment ID')
  seqs[key]=''
 elif l.strip():
  if key is None:raise SystemExit('Invalid alignment')
  seqs[key]+=l.strip()
expected={l[1:].split()[0] for l in (out/'09_BETA_CA_all_candidates_with_references.faa').read_text().splitlines() if l.startswith('>')}
if set(seqs)!=expected or len({len(s) for s in seqs.values()})!=1:raise SystemExit('Alignment validation failed')
q=json.loads((out/'REVIEW_QC.json').read_text());q['alignment_records']=len(seqs);q['alignment_columns']=len(next(iter(seqs.values())))
(out/'REVIEW_QC.json').write_text(json.dumps(q,indent=2)+'\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.iterdir()):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n')
(out/'COMPLETE_REVIEW_PREPARATION.txt').write_text('Review preparation completed. No formal acceptance/export performed.\n')
print('Review preparation and alignment checks completed.')
PY
printf 'Output: %s\n' "$BETA_OUT"
