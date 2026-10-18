#!/usr/bin/env bash
#SBATCH --job-name=pepc_close_wang
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=pepc_close_wang_%j.out
#SBATCH --error=pepc_close_wang_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
python3 - <<'PY'
import csv,gzip,hashlib,json,os,re,shutil,subprocess,tarfile
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
src=release('step25kz17.43','PEPC_TREE_RELEASE');prior=release('step25kz17.39','INTEGRATED_BETA_RELEASE')
need(sha(src/'05_prior_integrated_ledger_UNCHANGED.tsv')==sha(prior/'01_integrated_locus_evidence.tsv'),'Integrated ledger provenance mismatch')
ledger=read(prior/'01_integrated_locus_evidence.tsv');pepc=read(src/'03_PEPC_evidence_ledger.tsv');coverage=read(prior/'09_refreshed_family_coverage.tsv')
need(len(ledger)==820 and len(pepc)==55 and not any(r['wang_family']=='PEPC' for r in ledger),'Unexpected prior ledger state')
expected={'PEPC_DOMAIN_SUPPORTED_BROAD_COVERAGE':51,'PEPC_FAMILY_SUPPORTED_MODEL_CAVEAT':2,'EXCLUDED_HISTORICAL_CDS_RECORD':2}
need(dict(Counter(r['evidence_category'] for r in pepc))==expected,'Unexpected PEPC evidence categories')
# Parse the generated unquoted safe-ID Newick trees and record unrooted splits.
def parse(s):
 tok=re.findall(r'[^\s(),:;]+|[(),:;]',s);i=0
 def node():
  nonlocal i
  ch=[];label=''
  if tok[i]=='(':
   i+=1;ch.append(node())
   while tok[i]==',':i+=1;ch.append(node())
   need(tok[i]==')','Invalid Newick');i+=1
   if tok[i] not in ':,);':label=tok[i];i+=1
  else:label=tok[i];i+=1
  if tok[i]==':':i+=2
  tips=set().union(*(x[2] for x in ch)) if ch else {label}
  return ch,label,tips
 result=node();need(tok[i]==';' and i==len(tok)-1,'Unexpected trailing Newick data');return result
def nodes(n):
 yield n
 for c in n[0]:yield from nodes(c)
def splits(t,common):
 out={}
 for n in nodes(t):
  a=n[2]&common;b=common-a
  if min(len(a),len(b))<2:continue
  key=tuple(sorted(a)) if (len(a),sorted(a))<(len(b),sorted(b)) else tuple(sorted(b))
  out.setdefault(key,[]).append(n[1])
 return out
tiprows=read(src/'01_tree_tip_map.tsv');tipmap={r['tip_id']:r for r in tiprows};trees={m:parse((src/m/'PEPC_corrected.treefile').read_text()) for m in ['main','sensitivity']}
for mode,t in trees.items():
 expectedtips={r['tip_id'] for r in tiprows if r[mode+'_tree']=='True'}
 need(t[2]==expectedtips,'Tree tip mismatch')
common=trees['main'][2];sa=splits(trees['main'],common);sb=splits(trees['sensitivity'],common)
ppc4=[r['tip_id'] for r in tiprows if 'Q8GVE8' in r['original_id']];need(len(ppc4)==1,'PPC4 control missing')
# Test the exact reviewed PPC4-associated split, regardless of Newick root placement.
broad_split={r['protein_id'] for r in pepc if r['prior_domain_category']=='ORDERED_SPLIT_BROAD_COVERAGE'}
caveats={r['protein_id'] for r in pepc if r['count_tier']=='MODEL_CAVEAT'}
support=[]
for mode,t in trees.items():
 wanted={r['tip_id'] for r in tiprows if r['original_id'] in broad_split|(caveats if mode=='sensitivity' else set())}|set(ppc4)
 found=[n[1] for n in nodes(t) if n[2]==wanted or t[2]-n[2]==wanted]
 need(len(found)==1 and found[0]=='100/100','Reviewed PPC4-associated split changed; do not auto-close')
 support.append(dict(analysis=mode,candidate_loci=len(wanted)-1,SH_aLRT_UFBoot=found[0],tip_ids=';'.join(sorted(wanted))))
out=base/'step25kz17.44'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'pooled_sources').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
closed=[]
for r in pepc:
 r=dict(r);r['phylogenetic_review']='BOUNDED_REVIEW_CLOSED_WITH_LIMITATIONS' if r['count_tier']!='OUTSIDE_FAMILY_SET' else 'EXCLUDED_INVALID_HISTORICAL_RECORD'
 r['reference_group']='PPC4_ASSOCIATED_SUPPORTED_SPLIT' if r['protein_id'] in broad_split|caveats else 'OTHER_PEPC_REFERENCE_SIDE' if r['count_tier']=='SUPPORTED' else 'NOT_APPLICABLE'
 closed.append(r)
 ledger.append(dict(wang_family='PEPC',species=r['species'],gene_id=r['gene_id'],protein_id=r['protein_id'],subunit='',evidence_category=r['evidence_category'],count_tier=r['count_tier'],source_file='step25kz17.44/02_PEPC_closed_evidence_ledger.tsv',model_or_scope_limitation=r['limitation']+'; '+r['phylogenetic_review']+'; C4/CAM function and whole-genome completeness not established',exhaustive_count='NOT_ESTABLISHED',C4_CAM_specific_function='NOT_ESTABLISHED'))
need(len(ledger)==875 and len({(r['wang_family'],r['species'],r['gene_id']) for r in ledger})==875,'Duplicate family/species/locus keys')
table('01_integrated_locus_evidence.tsv',ledger,list(ledger[0]));table('02_PEPC_closed_evidence_ledger.tsv',closed,list(closed[0]));table('03_PPC4_associated_split_support.tsv',support,list(support[0]))
diffs=[dict(split_side=';'.join(k),main_support=';'.join(sa.get(k,[])),sensitivity_support=';'.join(sb.get(k,[])),presence='SHARED' if k in sa and k in sb else 'MAIN_ONLY' if k in sa else 'SENSITIVITY_ONLY') for k in sorted(sa.keys()|sb.keys())]
table('04_common_tip_split_comparison.tsv',diffs,['split_side','main_support','sensitivity_support','presence'])
families=[r['wang_family'] for r in coverage];species=sorted({r['species'] for r in ledger});tiers=['SUPPORTED','MODEL_CAVEAT','LIMITED_PANEL_CAVEAT','UNRESOLVED','OUTSIDE_FAMILY_SET']
need(len(set(families))==23 and len(species)==8,'Expected 23 families and eight species')
counts=Counter((r['wang_family'],r['species'],r['count_tier']) for r in ledger);integrated={r['wang_family'] for r in ledger}
summary=[]
for family in families:
 for sp in species:
  summary.append(dict(wang_family=family,species=sp,**{t:counts[(family,sp,t)] if family in integrated else 'NA' for t in tiers},records_integrated=sum(counts[(family,sp,t)] for t in tiers) if family in integrated else 'NA',final_gene_count='NA',whole_genome_completeness='NOT_ESTABLISHED'))
table('05_species_family_evidence_counts.tsv',summary,list(summary[0]))
for r in coverage:
 if r['wang_family']=='PEPC':r.update(current_status='BOUNDED_REVIEW_CLOSED_LEDGER_INTEGRATED',integrated_records=55,priority=4,next_action='RETAIN_51_SUPPORTED_2_MODEL_CAVEAT_2_EXCLUDED; COMPLETENESS_NOT_ESTABLISHED')
 elif r['wang_family']=='beta-CA':r.update(priority=4,next_action='RETAIN_EXISTING_SUPPORTED_CAVEAT_PARTIAL_AND_OUTSIDE_CATEGORIES')
 elif r['wang_family'] in {'TDT','DTC','DIC'}:r.update(priority=1,next_action='REVIEW_POOLED_SOURCE_IDENTIFIERS_AND_CURATE_SEPARATE_FAMILY_DEFINITIONS; DO_NOT_TRIPLE_COUNT')
table('06_refreshed_family_coverage.tsv',coverage,list(coverage[0]));shutil.copyfile(src/'06_final_Wang_counts_PENDING.tsv',out/'07_final_Wang_counts_PENDING.tsv')
# Detect retained loci assigned across families; retain exclusions only as historical evidence.
bylocus=defaultdict(set)
for r in ledger:
 if r['count_tier']!='OUTSIDE_FAMILY_SET':bylocus[(r['species'],r['gene_id'])].add(r['wang_family'])
conflicts=[dict(species=k[0],gene_id=k[1],families=';'.join(sorted(v)),decision='REVIEW_NO_AUTOMATIC_REASSIGNMENT') for k,v in bylocus.items() if len(v)>1]
table('08_cross_family_flags.tsv',conflicts,['species','gene_id','families','decision'])
# Collect small, explicitly scoped pooled sources. Never turn one pooled set into
# three accepted families; no biological screening is run at this step.
scopes=[base/'step25kn/results/TDT_DTC_DIC',base/'step25kz17.2/release/TDT_DTC_DIC']
paths=set()
for p in [base/'step25kn/references/TDT_DTC_DIC.faa',base/'step25kn/phylogeny_inputs/TDT_DTC_DIC.faa']:
 if p.is_file():paths.add(p)
for scope in scopes:
 if scope.is_dir():
  for p in scope.rglob('*'):
   if p.is_file() and len(p.relative_to(scope).parts)<=3 and p.suffix.lower() in {'.tsv','.csv','.faa','.fa','.fasta','.json','.txt'}:paths.add(p)
manifest=[];issues=[];headers=[];total=0
for i,p in enumerate(sorted(paths),1):
 size=p.stat().st_size
 if size>5_000_000 or total+size>20_000_000 or len(manifest)>=50:
  issues.append(dict(source=str(p.relative_to(root)),issue='BOUNDED_COLLECTION_SIZE_OR_COUNT_LIMIT'));continue
 total+=size;name=f'S{i:03d}_'+p.name;target=out/'pooled_sources'/name;original_hash=sha(p)
 text=p.read_text().replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>');target.write_text(text)
 manifest.append(dict(source=str(p.relative_to(root)),original_sha256=original_hash,snapshot='pooled_sources/'+name,snapshot_sha256=sha(target),bytes=size))
 if p.suffix.lower() in {'.faa','.fa','.fasta'}:
  for line in text.splitlines():
   if line.startswith('>'):headers.append(dict(source='pooled_sources/'+name,sequence_id=line[1:].split()[0],header=line[1:],family_assignment='UNRESOLVED_POOLED_EVIDENCE'))
 need(sha(p)==original_hash,'Pooled source changed during collection')
if not manifest:issues.append(dict(source='TDT_DTC_DIC bounded scopes',issue='NO_SOURCES_LOCATED_NOT_EVIDENCE_OF_ABSENCE'))
table('09_pooled_source_manifest.tsv',manifest,['source','original_sha256','snapshot','snapshot_sha256','bytes']);table('10_pooled_FASTA_headers.tsv',headers,['source','sequence_id','header','family_assignment']);table('11_pooled_collection_issues.tsv',issues,['source','issue'])
# Pin the evidence used to close the tree review in this small release.
for mode in ['main','sensitivity']:
 shutil.copyfile(src/mode/'PEPC_corrected.treefile',out/(mode+'_PEPC_corrected.treefile'))
shutil.copyfile(src/'01_tree_tip_map.tsv',out/'12_PEPC_tree_tip_map.tsv')
(out/'PROVENANCE.json').write_text(json.dumps(dict(PEPC_tree_manifest_sha256=sha(src/'SHA256SUMS.txt'),prior_integrated_manifest_sha256=sha(prior/'SHA256SUMS.txt')),indent=2)+'\n')
qc=dict(status='PEPC_BOUNDED_REVIEW_CLOSED_INTEGRATED_WITH_LIMITATIONS',integrated_records=len(ledger),integrated_families=len(integrated),species_family_cells=len(summary),PEPC_categories=expected,common_tips=len(common),shared_splits=len(sa.keys()&sb.keys()),main_only_splits=len(sa.keys()-sb.keys()),sensitivity_only_splits=len(sb.keys()-sa.keys()),cross_family_flags=len(conflicts),pooled_sources=len(manifest),pooled_collection_issues=len(issues),next_priority='TDT_DTC_DIC_DEFINITION_AND_SOURCE_REVIEW',analysis_rerun=False,final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'README.txt').write_text('Bounded PEPC review closed using corrected domain, annotation and phylogenetic evidence. 51 broadly supported loci; two model-caveat loci retained separately; two historical CDS-derived records excluded. PPC4-associated split verified at 100/100 in both trees. Reference grouping is not a biochemical, localization, or C4/CAM-specific function assignment. Unrooted common-tip split comparison is descriptive; separate alignments/masks and arbitrary branching among identical sequences can affect topology. Prior 820 rows preserved; 55 PEPC rows appended. Counts are evidence-category summaries, not exhaustive final Wang counts. Unintegrated families retain NA rather than zero. Cross-family flags require review. Coverage carries forward prior bounded searches; no fresh project-wide inventory is claimed. Next priority: separate pooled TDT/DTC/DIC definitions, references and candidate evidence before validation. Collected files are historical inputs only, not validated membership and never three copies of one family count. Missing/empty sources do not indicate biological absence. No analyses rerun.\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='PEPC_closed_Wang_coverage')
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
