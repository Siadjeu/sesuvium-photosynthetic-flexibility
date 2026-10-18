#!/usr/bin/env bash
#SBATCH --job-name=vatp_close
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=vatp_close_%j.out
#SBATCH --error=vatp_close_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
for tool in python3; do command -v "$tool"; done
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
src=release('step25kz17.35','VATPASE_CDS_RELEASE')
bounded=release('step25kz17.34','VATPASE_REVIEW_RELEASE')
assert json.loads((src/'PROVENANCE.json').read_text())['source_manifest_sha256']==sha(bounded/'SHA256SUMS.txt'),'Source chain mismatch'
rows=read(src/'prior_01_reconciled_evidence_ledger.tsv');checks=read(src/'01_protein_CDS_consistency.tsv')
assert len(rows)==603 and len({r['protein_id'] for r in rows})==603
assert all(r['mapping_status']=='UNAMBIGUOUS' for r in rows)
assert len({(r['species'],r['gene_id']) for r in rows})==603
assert len(checks)==23 and len({r['protein_id'] for r in checks})==23
assert all(r['status']=='EXACT_MATCH_TERMINAL_STOP_PRESENT' and r['internal_stops']=='0' and r['ambiguous_codons']=='0' and r['length_mod3']=='0' for r in checks),'Unexpected consistency result; review before closure'
expected={'CONCORDANT_EVIDENCE_REVIEW_REQUIRED':85,'CONCORDANT_WITH_LIMITED_PANEL_CAVEAT':208,'STRUCTURE_COMPLETENESS_EXCEPTION':23,'F_TYPE_FAVOURED_NOT_COUNTED_AS_VATPASE':287}
assert dict(Counter(r['evidence_category'] for r in rows))==expected
assert {r['protein_id'] for r in rows if r['evidence_category']=='STRUCTURE_COMPLETENESS_EXCEPTION'}=={r['protein_id'] for r in checks}
# Explicit reviewed caveats; no automatic promotion from translation consistency.
notes={
'AT__gene-LOC130807904':'Subunit E reference-relative C-terminal loss; AT CDS recovered through protein_id XP_057529280.1; biological completeness unresolved.',
'BV__Bv5_140060_npkj':'Subunit c-double-prime internal reference gaps; biological completeness unresolved. Source terminal-stop formatting normalized.',
'MC__gene_11085':'Near-complete subunit F reference occupancy with an insertion; length flag does not establish truncation.',
'MC__gene_21307':'Subunit G core with long C-terminal extension and DUF761 match; unusual architecture or annotation fusion unresolved.',
'MC__gene_21309':'Subunit G core with C-terminal extension; model architecture unresolved.',
'MC__gene_24691':'Subunit c-double-prime internal gaps and insertion; domain evidence retained, completeness unresolved.',
'MC__gene_5426':'Complete subunit D reference occupancy with insertion; structural interpretation unresolved.',
'MC__gene_8424':'Partial subunit B architecture with long N-terminal extension and Thioredoxin_8 matches; possible annotation fusion unresolved.',
'PO__evm.TU.LG02.1420':'Subunit c C-terminal reference loss; incomplete model evidence retained.',
'PO__evm.TU.LG04.1643':'Complete subunit c reference occupancy with N-terminal extension; model architecture unresolved.',
'PO__evm.TU.LG06.705':'Partial subunit c region with long extension and Steroid_dh match; possible fusion unresolved.',
'PO__evm.TU.LG08.1443':'Subunit c internal reference gap and C-terminal extension; overlapping domain matches are not independent gene copies.',
'PO__evm.TU.LG10.1244':'Subunit c terminal/internal reference losses; completeness unresolved.',
'PO__evm.TU.LG10.623':'Complete subunit B core with upstream MFS_5 domain; possible fusion unresolved; distinct locus retained.',
'PO__evm.TU.LG16.1221':'Strongly partial subunit a evidence; functional completeness not established.',
'PO__evm.TU.LG16.551':'Strongly partial subunit a evidence with large reference gap; functional completeness not established.',
'PO__evm.TU.LG17.635':'Complete subunit B core with upstream MFS_5 domain; possible fusion unresolved; distinct locus retained.',
'PO__evm.TU.LG18.169':'Strongly partial subunit a evidence; functional completeness not established.',
'SP_Hex__gene_Spor0365890.1':'Subunit E internal gap, insertion and extension; split profile hits retained without interpreting them as copies.',
'SP_Hex__gene_Spor0371460.1':'Incomplete subunit B model with substantial N-terminal reference loss and missing reference domain evidence.',
'SP_Hex__gene_Spor0389790.1':'Subunit a evidence with GATase match and internal architecture differences; possible fusion unresolved.',
'SP_Hex__gene_Spor0475540.1':'Subunit G alignment occupancy conflicts with partial domain-model coverage; extensions and DUF761 match; architecture unresolved.',
'SP_Hex__gene_Spor0482030.1':'Short subunit F fragment with low domain coverage and reference loss; functional completeness not established.'}
assert set(notes)=={r['protein_id'] for r in checks}
out=base/'step25kz17.36'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False)
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
labels={'CONCORDANT_EVIDENCE_REVIEW_REQUIRED':'CONCORDANT_NO_PRIORITY_MODEL_FLAG','CONCORDANT_WITH_LIMITED_PANEL_CAVEAT':'CONCORDANT_LIMITED_PANEL_CAVEAT','STRUCTURE_COMPLETENESS_EXCEPTION':'MODEL_EXCEPTION_RETAINED_WITH_CAVEAT','F_TYPE_FAVOURED_NOT_COUNTED_AS_VATPASE':'F_TYPE_FAVOURED_OUTSIDE_VATPASE_SET'}
ledger=[]
for r in rows:
 p=r['protein_id'];cat=labels[r['evidence_category']]
 note=notes.get(p,'')
 if cat=='CONCORDANT_LIMITED_PANEL_CAVEAT':note='No alternative reference hit in limited panel; not a model-defect flag. Specificity remains panel-dependent.'
 elif cat=='CONCORDANT_NO_PRIORITY_MODEL_FLAG':note='Concordant reference similarity, domain and coverage evidence; complete biochemical function and exhaustive membership not established.'
 elif cat=='F_TYPE_FAVOURED_OUTSIDE_VATPASE_SET':
  assert r['best_other_reference'] in ['P92549','P83483'] and float(r['score_margin'])<=-166
  note='Much stronger F-type ATP synthase reference match; not counted as V-ATPase. No locus deletion or definitive F-type functional annotation.'
 ledger.append(dict(protein_id=p,species=r['species'],gene_id=r['gene_id'],subunit=r['provisional_subunit'],subunit_status='PROVISIONAL_EVIDENCE_ASSIGNMENT',evidence_category=cat,original_evidence_category=r['evidence_category'],original_flags=r['flags'],best_VATPASE_reference=r['best_VATPASE_reference'],best_other_reference=r['best_other_reference'],score_margin=r['score_margin'],retained_VATPASE_candidate='NO' if cat=='F_TYPE_FAVOURED_OUTSIDE_VATPASE_SET' else 'YES',protein_CDS_status='EXACT_MATCH_TERMINAL_STOP_PRESENT' if p in notes else 'NOT_TESTED_IN_EXCEPTION_CHECK',bounded_review_status='CLOSED_WITH_LIMITATIONS',model_caveat=note,confirmed_complete_functional_gene='NOT_ESTABLISHED',whole_genome_completeness='NOT_ESTABLISHED'))
table('01_VATPASE_closed_evidence_ledger.tsv',ledger,list(ledger[0]))
table('02_VATPASE_exception_decisions.tsv',[r for r in ledger if r['protein_id'] in notes],list(ledger[0]))
table('03_F_type_favoured_outside_VATPASE.tsv',[r for r in ledger if r['retained_VATPASE_candidate']=='NO'],list(ledger[0]))
cats=list(labels.values());species=sorted({r['species'] for r in ledger});subs=['A','B','C','D','E','F','G','H','a','c',"c''",'d','e']
summary=[];subsummary=[]
for sp in species:
 rr=[r for r in ledger if r['species']==sp];c=Counter(r['evidence_category'] for r in rr)
 summary.append(dict(species=sp,**{k:c[k] for k in cats},retained_candidates=sum(r['retained_VATPASE_candidate']=='YES' for r in rr),confirmed_functional_count='NOT_ESTABLISHED'))
 for sub in subs:
  cc=Counter(r['evidence_category'] for r in rr if r['subunit']==sub)
  subsummary.append(dict(species=sp,subunit=sub,**{k:cc[k] for k in cats},count_type='EVIDENCE_CATEGORIES_NOT_FINAL_GENE_COUNTS'))
table('04_species_evidence_counts.tsv',summary,list(summary[0]));table('05_subunit_evidence_counts.tsv',subsummary,list(subsummary[0]))
for name in ['01_protein_CDS_consistency.tsv','02_AT_CDS_recovery_evidence.tsv','prior_05_exception_alignment_completeness.tsv','prior_06_exception_and_reference_domains.tsv','prior_07_exception_locus_summary.tsv']:
 shutil.copy2(src/name,out/('support_'+name))
# Carry forward the closed PPCK results unchanged for continuity.
for p in bounded.glob('PPCK_*'):
 if p.is_file():shutil.copy2(p,out/p.name)
(out/'LIMITATIONS.txt').write_text('Bounded V-ATPase annotation review closed by project decision. All 603 records map to distinct annotated loci. Counts preserve 85 concordant, 208 limited-panel-caveat, 23 model-exception and 287 F-type-favoured categories.\nThe retained 316 candidates comprise 293 without detected priority model problems plus 23 unresolved model exceptions. They are not 316 confirmed complete functional genes. No blanket acceptance of the 23 exceptions, pseudogene designation, model repair, locus merging or sequence deduplication.\nAll 23 exception CDSs translate exactly to the corresponding proteins allowing one terminal stop; no internal stops, ambiguous codons or incomplete codons detected. This establishes internal consistency, not independent validation of gene structures, fusions or completeness.\nThe AT extraction issue was resolved by a documented protein_id; the BV source terminal-stop formatting difference was recorded.\nSubunit assignments remain evidence-based and provisional. Shared domains and overlapping Pfam matches do not establish specificity or additional copies. Candidate ascertainment was limited to prior screened proteins; missing/divergent genes, localization, enzyme activity and C4/CAM-specific roles were not established.\nNo new searches, trees, alignments or genome/read validation were performed during closure. Earlier releases are unchanged.\n')
(out/'QC.json').write_text(json.dumps(dict(status='VATPASE_BOUNDED_REVIEW_CLOSED_WITH_LIMITATIONS',distinct_loci=603,evidence_category_counts=dict(Counter(r['evidence_category'] for r in ledger)),retained_candidates=316,exceptions_consistent_with_CDS=23,final_functional_counts=False,models_modified=False),indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(cds_release_manifest_sha256=sha(src/'SHA256SUMS.txt'),bounded_release_manifest_sha256=sha(bounded/'SHA256SUMS.txt'),closure_basis='User-authorized bounded review closure; explicit per-locus caveats retained'),indent=2)+'\n')
files=sorted(p for p in out.iterdir() if p.is_file())
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+p.name+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='VATPASE_closed_ledger')
print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
