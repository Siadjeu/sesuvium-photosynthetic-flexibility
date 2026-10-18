#!/usr/bin/env bash
#SBATCH --job-name=beta_residue_review
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=beta_residue_review_%j.out
#SBATCH --error=beta_residue_review_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
: "${CA_RELEASE:?Export CA_RELEASE as the completed Step25kz17.7 release directory}"
python3 - <<'PY'
import os,csv,json,hashlib,tarfile
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
source=Path(os.environ['CA_RELEASE']).resolve(strict=True)
inputs={n:source/n for n in ['02_CA_family_specific_review.tsv','03_domain_evidence.tsv','08_beta_only_review_alignment.faa']}
for p in inputs.values():
 if not p.is_file():raise SystemExit('Missing required input: '+p.name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read_tsv(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,fields,rows):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
rows=read_tsv(inputs['02_CA_family_specific_review.tsv']);domains=read_tsv(inputs['03_domain_evidence.tsv'])
if len({r['protein_id'] for r in rows})!=len(rows):raise SystemExit('Duplicate review identifier')
a={};key=None
for line in inputs['08_beta_only_review_alignment.faa'].read_text().splitlines():
 if line.startswith('>'):
  key=line[1:].split()[0]
  if key in a:raise SystemExit('Duplicate alignment identifier')
  a[key]=''
 elif line.strip():
  if key is None:raise SystemExit('Invalid FASTA alignment')
  a[key]+=line.strip().upper()
if not a or len(set(map(len,a.values())))!=1:raise SystemExit('Alignment lengths differ')
refs=[k for k in a if k.startswith('REF__')]
expected_refs={'P27140','P42737','Q94CE3','Q94CE4','Q9C6F5','Q9ZUC2'}
if len(refs)!=6 or {k.split('|')[1] for k in refs}!=expected_refs:raise SystemExit('Unexpected beta reference set')
anchor=next(k for k in refs if '|P27140|' in k)
positions=[167,169,171,227,230];expected='CDRHC';mapping={};n=0
for i,c in enumerate(a[anchor]):
 if c!='-':n+=1;mapping[n]=i
cols=[mapping[p] for p in positions]
for k in refs:
 if ''.join(a[k][i] for i in cols)!=expected:raise SystemExit('Reference residue mapping disagrees: '+k)
# Dense reference-alignment columns are only a diagnostic core approximation.
# Five of six references must contain an amino acid at a core column.
core=[i for i in range(len(a[anchor])) if sum(a[k][i]!='-' for k in refs)>=5]
if not core:raise SystemExit('No diagnostic reference core found')
primary=[r for r in rows if r['screen_class']=='BETA_SUPPORTED_REVIEW']
partial=[r for r in rows if r['screen_class']=='BETA_PARTIAL_REVIEW']
if set(a)-set(refs)!={r['protein_id'] for r in primary}:raise SystemExit('Alignment/review candidate sets differ')
basefields=list(rows[0])
fields=basefields+['anchor_residues','reference_residue_check','diagnostic_core_columns','covered_core_columns',
 'core_coverage_fraction','missing_core_at_N_end','missing_core_at_C_end','longest_internal_core_gap',
 'best_domain_hmm_from','best_domain_hmm_to','best_domain_model_length',
 'identical_protein_group','identical_protein_group_size','review_recommendation']
for pos in positions:fields+=['residue_at_BCA1_'+str(pos),'window_at_BCA1_'+str(pos)]
groups={}
for r in primary:groups.setdefault(a[r['protein_id']].replace('-',''),[]).append(r['protein_id'])
group_ids={seq:'SEQGROUP_'+str(i).zfill(3) for i,seq in enumerate(sorted(groups),1)}
review=[]
for old in primary:
 r=dict(old);pid=r['protein_id'];s=a[pid]
 if len(s.replace('-',''))!=int(r['length_aa']):raise SystemExit('Alignment length disagrees with review: '+pid)
 residues=''.join(s[i] for i in cols);r['anchor_residues']=residues
 r['reference_residue_check']='MATCH_ALL_FIVE' if residues==expected else 'REVIEW_GAP_AT_SITE' if '-' in residues else 'REVIEW_SUBSTITUTION_AT_SITE'
 covered=[j for j,i in enumerate(core) if s[i]!='-']
 r['diagnostic_core_columns']=len(core);r['covered_core_columns']=len(covered)
 r['core_coverage_fraction']=format(len(covered)/len(core),'.6f')
 r['missing_core_at_N_end']=covered[0] if covered else len(core)
 r['missing_core_at_C_end']=len(core)-1-covered[-1] if covered else len(core)
 longest=run=0
 if covered:
  for j in range(covered[0],covered[-1]+1):
   run=run+1 if s[core[j]]=='-' else 0;longest=max(longest,run)
 r['longest_internal_core_gap']=longest
 ds=[d for d in domains if d['protein_id']==pid and d['family']=='BETA_CA' and float(d['independent_evalue'])<=1e-5]
 if not ds:raise SystemExit('Missing significant beta domain: '+pid)
 best=max(ds,key=lambda d:(float(d['model_coverage']),float(d['domain_score'])))
 for dest,src in [('best_domain_hmm_from','hmm_from'),('best_domain_hmm_to','hmm_to'),('best_domain_model_length','model_length')]:r[dest]=best[src]
 seq=s.replace('-','');r['identical_protein_group']=group_ids[seq];r['identical_protein_group_size']=len(groups[seq])
 flags=[]
 if residues!=expected:flags.append('INSPECT_REFERENCE_RESIDUE_SITES')
 if r['missing_core_at_N_end']>0:flags.append('INSPECT_N_TERMINAL_CORE_COVERAGE')
 if r['missing_core_at_C_end']>0:flags.append('INSPECT_C_TERMINAL_CORE_COVERAGE')
 if longest>=10:flags.append('INSPECT_INTERNAL_CORE_GAP_GE_10')
 if r['sequence_flags']!='NONE':flags.append('INSPECT_SEQUENCE_FLAGS')
 r['review_recommendation']=';'.join(flags) if flags else 'NO_DIAGNOSTIC_FLAG_REVIEW_COMPLETENESS'
 for pos,i in zip(positions,cols):
  r['residue_at_BCA1_'+str(pos)]=s[i]
  # Five alignment columns centered on the mapped site, including gaps.
  r['window_at_BCA1_'+str(pos)]=s[max(0,i-2):i+3]
 # Preserve existing reviewer inputs exactly. Diagnostics never become manual PASS.
 review.append(r)
stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
out=root/'08_focused_families/validation/step25kz17.8'/('release_'+os.environ.get('SLURM_JOB_ID','manual')+'_'+stamp)
out.mkdir(parents=True,exist_ok=False)
table(out/'01_BETA_CA_residue_completeness_review.tsv',fields,review)
partial_rows=[]
for old in partial:
 r=dict(old);r['review_recommendation']='PARTIAL_BETA_SUPPORT_HOLD_OUTSIDE_PRIMARY_SET';partial_rows.append(r)
table(out/'02_BETA_CA_partial_hits_review.tsv',basefields+['review_recommendation'],partial_rows)
table(out/'03_reference_coordinate_map.tsv',['reference','reference_position','alignment_column','expected_residue'],
 [{'reference':anchor,'reference_position':pos,'alignment_column':i+1,'expected_residue':res} for pos,i,res in zip(positions,cols,expected)])
counts=Counter(r['reference_residue_check'] for r in review)
species=Counter((r['species'],r['reference_residue_check']) for r in review)
table(out/'04_species_residue_check_counts.tsv',['species','reference_residue_check','candidate_records'],
 [{'species':s,'reference_residue_check':c,'candidate_records':v} for (s,c),v in sorted(species.items())])
# Benchmark concordance does not influence assignments or acceptance.
bench=[]
for sp,wang in [('PO',11),('BV',3)]:
 n=sum(r['species']==sp for r in primary)
 bench.append({'species':sp,'wang_S6_user_supplied_beta_CA':wang,'provisional_supported_candidates':n,'difference':n-wang,'interpretation':'COUNT_BENCHMARK_ONLY_NOT_SEQUENCE_VALIDATION'})
table(out/'05_Wang_S6_count_benchmark.tsv',list(bench[0]),bench)
notes={
 'status':'COMPLETE_RESIDUE_COMPLETENESS_DIAGNOSTICS_REVIEW_ONLY',
 'primary_beta_candidates':len(review),'partial_beta_candidates':len(partial),
 'reference_residue_checks':dict(counts),'diagnostic_core_columns':len(core),
 'candidates_with_diagnostic_flags':sum(r['review_recommendation']!='NO_DIAGNOSTIC_FLAG_REVIEW_COMPLETENESS' for r in review),
 'formal_export_performed':False,'tree_built':False,
 'input_sha256':{k:sha(p) for k,p in inputs.items()},
 'interpretation':[
 'MATCH_ALL_FIVE means alignment residues match all six references; it does not prove catalytic activity.',
 'Core columns have residues in at least five of six Arabidopsis references. This is a diagnostic alignment core, not an experimentally defined domain.',
 'Terminal gaps, long internal gaps and reference residues guide inspection; they are not automatic gene completeness calls.',
 'All existing manual decisions/checks are preserved; no new PASS or formal acceptance is assigned.',
 'Identical proteins retain separate identifiers. Sequence equality does not establish distinct loci.',
 'Five-character windows are alignment columns centered on a mapped BCA1 position; they can contain gaps.',
 'Wang S6 counts were supplied by the user and are external benchmarks, not fitting targets.',
 'Review fields: new_reviewer_decision, new_catalytic_region_check, new_truncation_check, new_reviewer_reason.',
 'Use PASS for a manual check only after inspection. Keep uncertain candidates REVIEW_MANUAL or PENDING.',
 'Recommended decision vocabulary: ACCEPT_FORMAL_BETA_CA, REVIEW_MANUAL, EXCLUDE_NON_BETA_CA, EXCLUDE_FRAGMENT, PENDING.'
 ]}
(out/'RESIDUE_REVIEW_QC.json').write_text(json.dumps(notes,indent=2)+'\n')
(out/'COMPLETE.txt').write_text('Diagnostic review prepared. No formal export or tree built.\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.iterdir()):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+p.name+'\n')
archive=out/'BETA_CA_residue_review_bundle.tar.gz'
with tarfile.open(archive,'w:gz') as t:
 for p in sorted(out.iterdir()):
  if p.is_file() and p!=archive:t.add(p,arcname='BETA_CA_residue_review/'+p.name)
print('Reference-site checks:',dict(counts))
print('Output: <PROJECT_ROOT>/'+out.relative_to(root).as_posix())
print('Upload BETA_CA_residue_review_bundle.tar.gz')
PY
