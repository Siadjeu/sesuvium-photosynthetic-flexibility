#!/usr/bin/env bash
#SBATCH --job-name=wang_ledger
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=wang_ledger_%j.out
#SBATCH --error=wang_ledger_%j.err
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
vatp=release('step25kz17.36','VATPASE_CLOSED_RELEASE')
vpp=release('step25kz17.30','VPPASE_LEDGER_RELEASE')
ppdk=release('step25kz17.27','PPDK_LEDGER_RELEASE')
coverage=release('step25kz17.17','COVERAGE_RELEASE')
prior=release('step25kz17.16','WANG_RELEASE')
species=['AT','BV','MC','PO','SLP_Hap1','SP_Hex','SP_Tet','SS']
families=['beta-CA','PEPC','PPCK','NAD-ME','NADP-ME','MDH','ASPAT','AlaAT','PEPCK','PPDK','PPDK-RE','BASS2','MEP','DIT2','DIC','DTC','MPC','ALMT','TDT','PPA','PPT','V_ATPase','V_Ppase']
lookup=dict(zip(['BETA_CA','PEPC','PEPC_KINASE','NAD_ME','NADP_ME','MDH','ASPAT','ALAAT','PEPCK','PPDK','PPDK_RP','BASS2','MEP','DIT2','DIC','DTC','MPC','ALMT','TDT','PPA','PPT','V_ATPASE','V_PPASE'],families))
out=base/'step25kz17.37'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'sources').mkdir()
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
ledger=[];sourceinfo=[]
def archive(p,label):
 dest=out/'sources'/(label+'__'+p.name);shutil.copy2(p,dest)
 sourceinfo.append(dict(source_label=label,file=p.name,sha256=sha(p),release_manifest_sha256=sha(p.parent/'SHA256SUMS.txt'),archived_file=str(dest.relative_to(out))))
 return label+'__'+p.name
def add(fam,r,cat,tier,source,sub='',caveat=''):
 assert r['species'] in species and r['gene_id'] and ';' not in r['gene_id']
 ledger.append(dict(wang_family=fam,species=r['species'],gene_id=r['gene_id'],protein_id=r['protein_id'],subunit=sub,evidence_category=cat,count_tier=tier,source_file=source,model_or_scope_limitation=caveat,exhaustive_count='NOT_ESTABLISHED',C4_CAM_specific_function='NOT_ESTABLISHED'))
p=ppdk/'09_PPDK_PPDK_RP_decision_ledger.tsv';source=archive(p,'PPDK')
rr=read(p);assert len(rr)==32
for r in rr:
 assert r['locus_mapping']=='UNAMBIGUOUS' and r['decision'] in ['SUPPORTED_MODEL','UNRESOLVED_MODEL']
 add('PPDK' if r['family']=='PPDK' else 'PPDK-RE',r,r['decision'],'SUPPORTED' if r['decision']=='SUPPORTED_MODEL' else 'UNRESOLVED',source,caveat=r['reason']+'; '+r['read_validation'])
p=vpp/'11_VPPASE_decision_ledger.tsv';source=archive(p,'VPPASE');rr=read(p);assert len(rr)==48
for r in rr:
 assert r['mapping_status']=='UNAMBIGUOUS'
 tiers={'SUPPORTED_NO_PRIORITY_MODEL_FLAG':'SUPPORTED','FAMILY_SUPPORTED_MODEL_CAVEAT':'MODEL_CAVEAT','UNRESOLVED_FRAGMENT':'UNRESOLVED'}
 add('V_Ppase',r,r['decision'],tiers[r['decision']],source,caveat=r['review_flag']+'; suspected_region='+r.get('suspected_region','')+'; '+r['read_validation'])
p=vatp/'PPCK_01_PPCK_closed_ledger.tsv';source=archive(p,'PPCK');rr=read(p);assert len(rr)==27
for r in rr:
 assert r['decision']=='PHYLOGENETICALLY_SUPPORTED_PPCK_LOCUS'
 add('PPCK',r,r['decision'],'SUPPORTED',source,caveat='Bounded screened-candidate check; profile calibration was in-sample; PPCK1/2 subtype and function unresolved.')
p=vatp/'01_VATPASE_closed_evidence_ledger.tsv';source=archive(p,'VATPASE');rr=read(p);assert len(rr)==603
for r in rr:
 tiers={'CONCORDANT_NO_PRIORITY_MODEL_FLAG':'SUPPORTED','CONCORDANT_LIMITED_PANEL_CAVEAT':'LIMITED_PANEL_CAVEAT','MODEL_EXCEPTION_RETAINED_WITH_CAVEAT':'MODEL_CAVEAT','F_TYPE_FAVOURED_OUTSIDE_VATPASE_SET':'OUTSIDE_FAMILY_SET'}
 add('V_ATPase',r,r['evidence_category'],tiers[r['evidence_category']],source,sub=r['subunit'],caveat=r['model_caveat'])
assert len(ledger)==710
keys=[(r['wang_family'],r['species'],r['gene_id']) for r in ledger];assert len(keys)==len(set(keys)),'Duplicate family/locus entries'
for p,label in [(ppdk/'11_limitations.txt','PPDK'),(vpp/'13_VPPASE_limitations.txt','VPPASE'),(vatp/'LIMITATIONS.txt','VATPASE'),(vatp/'PPCK_03_PPCK_closure.txt','PPCK'),(coverage/'01_family_evidence_coverage.tsv','COVERAGE'),(coverage/'09_reference_definition_review.tsv','COVERAGE')]:archive(p,label)
# Preserve the prior framework unchanged, including BETA-CA/PEPC evidence and benchmark inputs.
for name in ['04_gene_evidence.tsv','05_gene_count_review.tsv','07_count_completeness_review.tsv','08_final_gene_counts_PENDING.tsv','09_PO_BV_benchmark_comparison.tsv','11_issues.tsv']:
 p=prior/name
 if p.is_file():archive(p,'PRIOR_WANG')
table('01_integrated_locus_evidence.tsv',ledger,list(ledger[0]))
tiers=['SUPPORTED','LIMITED_PANEL_CAVEAT','MODEL_CAVEAT','UNRESOLVED','OUTSIDE_FAMILY_SET']
active={r['wang_family'] for r in ledger};counts=[]
for sp in species:
 for fam in families:
  rr=[r for r in ledger if r['species']==sp and r['wang_family']==fam];c=Counter(r['count_tier'] for r in rr)
  counts.append(dict(species=sp,wang_family=fam,**{k:c[k] if fam in active else 'NA' for k in tiers},records_integrated=len(rr) if fam in active else 'NA',coverage_status='RECENT_LEDGER_INTEGRATED' if fam in active else 'PRIOR_EVIDENCE_RECONCILIATION_REQUIRED',final_gene_count='NA'))
table('02_species_family_evidence_counts.tsv',counts,list(counts[0]))
# Each matrix is a separate evidence layer; never sum unlike layers into final counts.
for tier in tiers:
 matrix=[dict(species=sp,**{fam:next(r[tier] for r in counts if r['species']==sp and r['wang_family']==fam) for fam in families}) for sp in species]
 table('03_'+tier+'_matrix.tsv',matrix,['species']+families)
table('04_final_Wang_counts_PENDING.tsv',[dict(species=sp,**{f:'NA' for f in families}) for sp in species],['species']+families)
old={lookup[r['family']]:r for r in read(coverage/'01_family_evidence_coverage.tsv') if r['family'] in lookup}
prior_ev=read(prior/'04_gene_evidence.tsv') if (prior/'04_gene_evidence.tsv').is_file() else []
refresh=[]
for fam in families:
 o=old.get(fam,{});n=sum(r.get('wang_category')==fam for r in prior_ev)
 if fam in active:state='RECENT_LEDGER_INTEGRATED_WITH_LIMITATIONS';action='RETAIN_CATEGORY_COUNTS_NO_ANALYSIS_RERUN';priority=4
 elif fam in ['beta-CA','PEPC']:state='PRIOR_FAMILY_ANALYSES_REQUIRE_LOCUS_LEDGER_RECONCILIATION';action='RECONCILE_EXISTING_FAMILY_SPECIFIC_RELEASES_AND_LOCUS_DECISIONS';priority=1
 elif fam in ['DIC','DTC','TDT']:state='POOLED_FAMILY_DEFINITION_UNRESOLVED';action='SEPARATE_POOLED_TDT_DTC_DIC_EVIDENCE_WITH_CURATED_FAMILY_DEFINITIONS';priority=2
 elif int(o.get('matching_paths','0'))>0:state='HISTORICAL_EVIDENCE_NOT_YET_INTEGRATED';action='REVIEW_EXISTING_DEFINITIONS_AND_EVIDENCE_BEFORE_NEW_SEARCH';priority=2
 else:state='NO_LOCATED_EVIDENCE_IN_PRIOR_COVERAGE_AUDIT';action='DEFINE_SCOPE_AND_LOCATE_EXISTING_INPUTS_BEFORE_TARGETED_SEARCH';priority=3
 refresh.append(dict(wang_family=fam,current_status=state,integrated_records=sum(r['wang_family']==fam for r in ledger) if fam in active else 'NA',prior_evidence_rows=n,historical_matching_paths=o.get('matching_paths','NA'),historical_reference_paths=o.get('reference_named_paths','NA'),search_completeness='NOT_ESTABLISHED',priority=priority,next_action=action))
table('05_refreshed_family_coverage.tsv',refresh,list(refresh[0]))
table('06_unresolved_family_priorities.tsv',sorted([r for r in refresh if r['wang_family'] not in active],key=lambda r:(r['priority'],families.index(r['wang_family']))),list(refresh[0]))
# Detect cross-family locus reuse without silently merging or discarding it.
bylocus=defaultdict(list)
for r in ledger:
 if r['count_tier']!='OUTSIDE_FAMILY_SET':bylocus[(r['species'],r['gene_id'])].append(r)
conflicts=[dict(species=k[0],gene_id=k[1],families=';'.join(sorted({r['wang_family'] for r in rr})),status='REVIEW_CROSS_FAMILY_MEMBERSHIP') for k,rr in bylocus.items() if len({r['wang_family'] for r in rr})>1]
table('07_cross_family_locus_flags.tsv',conflicts,['species','gene_id','families','status'])
table('08_source_manifest.tsv',sourceinfo,list(sourceinfo[0]))
(out/'README.txt').write_text('Integration of five recent family ledgers into the 23-column Wang framework. 710 family/locus records are retained including outside-family records. Source evidence categories and limitations are preserved. SUPPORTED means supported at the scope of the source analysis, not proven biochemical function, exhaustive family membership, or uniform evidence strength across families.\nV-ATPase 85 concordant, 208 limited-panel caveat, 23 model caveats and 287 outside-set records remain separate. V-PPase fragments and suspected split regions remain unresolved; no merging. PPDK-RP corresponds to Wang PPDK-RE. PPCK is distinct from PEPC. Subunit case is preserved.\nPrior BETA-CA, PEPC and other evidence is archived unchanged; no counts are inferred from old FASTA sizes or tree tips. Refreshed coverage distinguishes completed integration from earlier analyses requiring locus-ledger reconciliation. The historical coverage inventory is not a new whole-project search; absent evidence is not biological absence.\nAll final-count cells remain NA. Zeros in evidence layers apply only to integrated ledger categories, not genome-wide absence. No single grand total across categories/families is asserted. No new biological analysis, source editing, sequence deduplication, or replacement of prior decisions.\nNext: reconcile existing BETA-CA and PEPC locus decisions, then review unresolved transporter definitions and remaining family coverage. Preserve distinct loci even with identical proteins; preserve unresolved functional and completeness questions.\n')
(out/'QC.json').write_text(json.dumps(dict(status='RECENT_LEDGERS_INTEGRATED_FINAL_COUNTS_PENDING',species=8,wang_families=23,integrated_families=sorted(active),integrated_records=len(ledger),cross_family_locus_flags=len(conflicts),analysis_rerun=False,final_counts=False),indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps({label:sha(p/'SHA256SUMS.txt') for label,p in [('VATPASE_17_36',vatp),('VPPASE_17_30',vpp),('PPDK_17_27',ppdk),('COVERAGE_17_17',coverage),('WANG_17_16',prior)]},indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file())
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='Integrated_family_ledgers')
print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
