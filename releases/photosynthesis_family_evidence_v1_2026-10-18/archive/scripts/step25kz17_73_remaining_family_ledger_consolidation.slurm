#!/usr/bin/env bash
#SBATCH --job-name=remaining_consolidate
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=remaining_consolidate_%j.out
#SBATCH --error=remaining_consolidate_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v iqtree2 || command -v iqtree
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












import signal
src=release('step25kz17.72','REMAINING_TREE_RELEASE');source_manifest=sha(src/'SHA256SUMS.txt')
out=base/'step25kz17.73'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
def table(name,rows,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
roles=read(src/'03_reference_roles.tsv');membership=read(src/'01_all_membership.tsv');reviews=read(src/'09_combined_domain_review.tsv');old=read(src/'source_01_integrated_locus_evidence.tsv');oldfields=list(old[0]);groupqc={x['group']:dict(x) for x in json.loads((src/'QC.json').read_text())['groups']}
for group in groupqc:shutil.copytree(src/group,out/group)
iq=shutil.which('iqtree2') or shutil.which('iqtree');need(iq,'IQ-TREE missing');cpu=os.environ.get('SLURM_CPUS_PER_TASK','4');commands=[]
version=subprocess.run([iq,'--version'],capture_output=True,text=True,timeout=15,check=True).stdout
(out/'iqtree_version.txt').write_text(version)
need(version.splitlines()[0]==(src/'iqtree_version.txt').read_text().splitlines()[0],'IQ-TREE version differs from checkpoint producer; activate original environment')
d=out/'Aminotransferase';resume=dict(status='NOT_NEEDED',error='',source_checkpoint_sha256='',resume_message_observed=False)
if groupqc['Aminotransferase']['status']!='TREE_COMPLETE_REVIEW_ONLY':
 need((d/'main.ckp.gz').is_file(),'Aminotransferase checkpoint missing; refusing fresh restart')
 resume['source_checkpoint_sha256']=sha(d/'main.ckp.gz');need(sha(d/'03_main_core.faa')==sha(src/'Aminotransferase/03_main_core.faa'),'Alignment changed')
 cmd=[iq,'-s',str(d/'03_main_core.faa'),'-st','AA','-m','LG+G4','-nt',cpu,'-seed','17272','-bb','1000','-alrt','1000','-bnni','-keep-ident','-pre',str(d/'main')];commands.append(cmd)
 print('Resuming copied aminotransferase checkpoint; 60-minute cap; no redo',flush=True)
 with (d/'resume_stdout.log').open('w') as a,(d/'resume_stderr.log').open('w') as b:
  p=subprocess.Popen(cmd,stdout=a,stderr=b,start_new_session=True)
  try:
   rc=p.wait(timeout=3600);resume['status']='COMPLETED' if rc==0 else 'BLOCKED_EXIT_'+str(rc)
  except subprocess.TimeoutExpired:
   os.killpg(p.pid,signal.SIGTERM)
   try:p.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
   resume['status']='TIME_LIMIT_CHECKPOINT_RETAINED'
 log=(d/'resume_stdout.log').read_text(errors='replace')+(d/'resume_stderr.log').read_text(errors='replace')
 resume['resume_message_observed']=bool(re.search(r'checkpoint|resum',log,re.I))
 if resume['status']=='COMPLETED' and not resume['resume_message_observed']:resume['status']='BLOCKED_RESUME_NOT_CONFIRMED'
 if resume['status']=='COMPLETED':
  need((d/'main.treefile').is_file() and (d/'main.iqtree').is_file(),'Completed tree report missing');groupqc['Aminotransferase']['status']='TREE_COMPLETE_REVIEW_ONLY';groupqc['Aminotransferase']['error']=''
 else:groupqc['Aminotransferase']['status']=resume['status']
(d/'QC.json').write_text(json.dumps(groupqc['Aminotransferase'],indent=2)+'\n');(out/'00_resume_status.json').write_text(json.dumps(resume,indent=2)+'\n')
# Parse unrooted tree edge bipartitions; no assumed root or nearest-tip assignment.
def parse(s):
 tokens=re.findall(r'[(),:;]|[^(),:;\s]+',s);i=0;nodes=[];leaves=[]
 def node():
  nonlocal i
  kids=[];label=''
  if tokens[i]=='(':
   i+=1;kids.append(node())
   while tokens[i]==',':i+=1;kids.append(node())
   need(tokens[i]==')','Malformed Newick');i+=1
  if tokens[i] not in [':',',',')',';']:label=tokens[i];i+=1
  if tokens[i]==':':i+=2
  tips=set().union(*(n[0] for n in kids)) if kids else {label}
  if not kids:leaves.append(label)
  n=(tips,label);nodes.append(n);return n
 alltips=node()[0];need(len(leaves)==len(set(leaves)),'Duplicate tree tips');return alltips,nodes
families=['NAD-ME','NADP-ME','MDH','ASPAT','AlaAT','BASS2','PEPCK','MPC'];groups={'Malic_enzyme':['NAD-ME','NADP-ME'],'Aminotransferase':['ASPAT','AlaAT'],'MDH_LDH':['MDH'],'BASS':['BASS2'],'PEPCK':['PEPCK'],'MPC':['MPC']}
associations=defaultdict(list);edge_rows=[]
for group,fams in groups.items():
 if groupqc[group]['status']!='TREE_COMPLETE_REVIEW_ONLY':continue
 rr=[x for x in membership if x['group']==group];by={x['tip_id']:x for x in rr};tips,nodes=parse((out/group/'main.treefile').read_text());need(tips=={x['tip_id'] for x in rr if x['alignment_category']=='MAIN_CORE'},'Tree membership mismatch '+group)
 refids={k for k in tips if by[k]['kind']=='REFERENCE'}
 for fam in fams:
  positive={k for k in refids if any(r['wang_family']==fam and r['role']=='POSITIVE' and r['reference_id']==by[k]['original_id'] for r in roles)}
  need(positive,'No positive anchors '+fam);eligible={}
  for tt,label in nodes:
   sup=label.split('/')
   if len(sup)!=2:continue
   try:sh,uf=map(float,sup)
   except ValueError:continue
   for side in [tt,tips-tt]:
    anchors=side&refids
    if not anchors or not anchors<=positive or not side-refids:continue
    key=frozenset(side)
    eligible[key]=(sh,uf)
  # For each anchor, the largest reference-exclusive side is retained before support assessment.
  selected=set()
  for anchor in positive:
   choices=[side for side in eligible if anchor in side]
   if choices:selected.add(max(choices,key=lambda x:(len(x),tuple(sorted(x)))))
  for j,side in enumerate(sorted(selected,key=lambda x:tuple(sorted(x))),1):
   sh,uf=eligible[side];supported=sh>=80 and uf>=95;edgeid=group+'_'+fam+'_'+str(j)
   edge_rows.append(dict(group=group,wang_family=fam,edge_id=edgeid,reference_accessions=';'.join(sorted(by[k]['original_id'] for k in side&refids)),contains_all_family_anchors=positive<=side,candidate_proteins=len(side-refids),SH_aLRT=sh,UFBoot=uf,support_category='SUPPORTED_BOUNDED_ASSOCIATION' if supported else 'SUPPORT_LIMITED_ASSOCIATION'))
   for k in side-refids:associations[(fam,by[k]['original_id'])].append((edgeid,sh,uf,supported))
table('02_reference_edge_associations.tsv',edge_rows,['group','wang_family','edge_id','reference_accessions','contains_all_family_anchors','candidate_proteins','SH_aLRT','UFBoot','support_category'])
mem={(r['group'],r['original_id']):r for r in membership if r['kind']=='CANDIDATE'};new=[];proteinrows=[]
for r in reviews:
 fam=r['wang_family'];k=r['protein_id'];group=next(g for g,ff in groups.items() if fam in ff);m=mem[(group,k)];aa=associations.get((fam,k),[]);supported=any(a[3] for a in aa);flags=m['architecture_flags'];cat=m['alignment_category']
 limitation='FIXED_MODEL_BOUNDED_REFERENCE_PANEL;NO_INDEPENDENT_READ_VALIDATION;EXHAUSTIVENESS_NOT_ESTABLISHED'
 if fam in ['MPC','PEPCK']:
  limitation+=';NO_DOMAIN_COMPATIBLE_INDEPENDENT_COMPARATOR_TREE'
  if cat=='MAIN_CORE' and r['best_positive']:tier='LIMITED_PANEL_CAVEAT';evidence='DOMAIN_AND_REFERENCE_ALIGNMENT_COMPATIBLE'
  elif cat=='PARTIAL_OR_ARCHITECTURE' and r['best_positive']:tier='PARTIAL';evidence='PARTIAL_OR_ARCHITECTURE_REVIEW'
  else:tier='UNRESOLVED';evidence='OUTSIDE_ALIGNMENT_NOT_EXCLUDED'
 elif cat=='MAIN_CORE' and supported:
  tier='MODEL_CAVEAT' if flags!='NONE' else 'SUPPORTED';evidence='BOUNDED_REFERENCE_PHYLOGENETIC_ASSOCIATION'
 elif cat=='MAIN_CORE' and aa:tier='UNRESOLVED';evidence='PHYLOGENETIC_SUPPORT_LIMITED'
 elif cat=='MAIN_CORE' and groupqc[group]['status']=='TREE_COMPLETE_REVIEW_ONLY':tier='OUTSIDE_FAMILY_SET';evidence='OUTSIDE_REFERENCE_EXCLUSIVE_GROUP_NOT_PROVEN_NONMEMBER'
 elif cat=='PARTIAL_OR_ARCHITECTURE':tier='PARTIAL';evidence='CORE_PARTIAL_OR_ARCHITECTURE_NOT_PHYLOGENETICALLY_ASSIGNED'
 else:tier='UNRESOLVED';evidence='NO_COMPLETED_FAMILY_ASSIGNMENT'
 if flags!='NONE':limitation+=';'+flags
 if not r['species'] or not r['gene_id']:tier='UNRESOLVED';evidence='LOCUS_MAPPING_UNRESOLVED'
 proteinrows.append(dict(wang_family=fam,protein_id=k,species=r['species'],gene_id=r['gene_id'],alignment_category=cat,edge_ids=';'.join(a[0] for a in aa),count_tier=tier,evidence_category=evidence,limitation=limitation))
# Explicit mutually exclusive family conflict gate. Shared screening rows are not added twice as supported counts.
conflicts=[]
for fa,fb in [('NAD-ME','NADP-ME'),('ASPAT','AlaAT')]:
 pa={(r['species'],r['gene_id']) for r in proteinrows if r['wang_family']==fa and r['count_tier'] in {'SUPPORTED','MODEL_CAVEAT','LIMITED_PANEL_CAVEAT'}}
 pb={(r['species'],r['gene_id']) for r in proteinrows if r['wang_family']==fb and r['count_tier'] in {'SUPPORTED','MODEL_CAVEAT','LIMITED_PANEL_CAVEAT'}}
 for sp,g in sorted(pa&pb):
  conflicts.append(dict(species=sp,gene_id=g,families=fa+';'+fb))
  for r in proteinrows:
   if r['species']==sp and r['gene_id']==g and r['wang_family'] in {fa,fb}:r.update(count_tier='UNRESOLVED',evidence_category='CROSS_FAMILY_LOCUS_CONFLICT')
table('03_protein_decisions.tsv',proteinrows,list(proteinrows[0]));table('04_cross_family_conflicts.tsv',conflicts,['species','gene_id','families'])
loci=defaultdict(list)
for r in proteinrows:loci[(r['wang_family'],r['species'],r['gene_id'] or 'UNMAPPED:'+r['protein_id'])].append(r)
for (fam,sp,g),rr in sorted(loci.items()):
 tiers={r['count_tier'] for r in rr};tier=next(iter(tiers)) if len(tiers)==1 else 'UNRESOLVED'
 new.append(dict(wang_family=fam,species=sp,gene_id=g,protein_id=';'.join(sorted({r['protein_id'] for r in rr})),subunit='',evidence_category=';'.join(sorted({r['evidence_category'] for r in rr})) if len(tiers)==1 else 'WITHIN_LOCUS_EVIDENCE_CONFLICT',count_tier=tier,source_file='03_protein_decisions.tsv',model_or_scope_limitation=';'.join(sorted({r['limitation'] for r in rr})),exhaustive_count='NOT_ESTABLISHED',C4_CAM_specific_function='NOT_ESTABLISHED'))
need(not(set(families)&{r['wang_family'] for r in old}),'Target families already integrated; stop to avoid duplicate ledger insertion')
ledger=old+new;table('01_integrated_locus_evidence.tsv',ledger,oldfields);table('05_new_family_locus_ledger.tsv',new,oldfields)
counts=Counter((r['wang_family'],r['species'],r['count_tier']) for r in new);table('06_new_family_species_evidence_counts.tsv',[dict(wang_family=f,species=sp,count_tier=t,evidence_loci=n,formal_Wang_count='NA') for (f,sp,t),n in sorted(counts.items())],['wang_family','species','count_tier','evidence_loci','formal_Wang_count'])
coverage=read(src/'11_refreshed_family_coverage.tsv')
for r in coverage:
 if r['wang_family'] in families:
  fam=r['wang_family'];group=next(g for g,ff in groups.items() if fam in ff);r.update(current_status='BOUNDED_EVIDENCE_LEDGER_INTEGRATED',integrated_records=sum(x['wang_family']==fam for x in new),next_action='RETAIN_LIMITATIONS_AND_RESOLVE_REMAINING_DEFINITIONS',search_completeness='NOT_ESTABLISHED')
  if groupqc[group]['status'] not in {'TREE_COMPLETE_REVIEW_ONLY','ALIGNED_LIMITED_PANEL_NO_SPECIFICITY_TREE'}:r.update(current_status='LEDGER_INTEGRATED_WITH_UNFINISHED_PHYLOGENY',next_action='REVIEW_RETAINED_CHECKPOINT_STATUS')
table('08_refreshed_family_coverage.tsv',coverage,list(coverage[0]))
allspecies=sorted({r['species'] for r in ledger if r['species']});pending=[dict(wang_family=f,species=sp,evidence_records=sum(r['wang_family']==f and r['species']==sp for r in ledger),formal_Wang_count='NA',status='FINAL_SCOPE_AND_COMPLETENESS_NOT_ESTABLISHED') for f in [r['wang_family'] for r in coverage] for sp in allspecies];table('07_final_Wang_counts_PENDING.tsv',pending,list(pending[0]))
for n in ['source_01_integrated_locus_evidence.tsv','source_07_final_Wang_counts_PENDING.tsv','03_reference_roles.tsv','09_combined_domain_review.tsv','01_all_membership.tsv']:shutil.copyfile(src/n,out/('prior_'+n))
def clean(s):return s.replace(str(out),'<RELEASE>').replace(str(src),'<SOURCE_RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>')
for p in out.rglob('*'):
 if p.is_file() and p.suffix in {'.log','.iqtree','.txt'}:p.write_text(clean(p.read_text(errors='replace')))
(out/'PROVENANCE.json').write_text(json.dumps(dict(step72_manifest_sha256=source_manifest,upstream=json.loads((src/'PROVENANCE.json').read_text()),resume=resume,commands=[[clean(x) for x in c] for c in commands]),indent=2)+'\n')
(out/'README.txt').write_text('Combined checkpoint resume and evidence-ledger consolidation. Previous release left unchanged; aminotransferase directory copied including unmodified checkpoint and identical alignment. Same IQ-TREE version and analysis options, no redo; paths relocated to new release.60min resume cap. Completed other trees reused. Unrooted edge sides containing target-positive anchors and no other reference tips identified; largest eligible side per anchor selected before support gate80SH-aLRT/95UFBoot. This is bounded panel association, not validated substrate/cofactor biochemistry or C4/CAM function. Not all target-family references need be monophyletic; subset-anchor groups are explicit. No bootstrap-threshold chasing or sensitivity reruns. MAIN_CORE means domain-alignment eligibility, not independently verified full gene model. Flags preserved; partial/core-ineligible candidates remain PARTIAL or UNRESOLVED, never discarded from lack of tree placement. Outside reference group is OUTSIDE_FAMILY_SET, not a definitive nonmember claim. PEPCK/MPC use LIMITED_PANEL_CAVEAT, never supported counts without comparators. One record per family/species/locus, transcript IDs retained; conflicting locus assignments unresolved. Paired-family supported conflicts blocked. Old ledger rows retained verbatim in value/order. New evidence tiers are provisional; formal Wang counts all NA. MEP/PPA remain definition-blocked; unsupported function/exhaustiveness not inferred. Source checkpoint and copied outputs may retain software-embedded paths; human-readable logs sanitized. IQ-TREE checkpoint documentation: https://iqtree.github.io/doc/Command-Reference .\n')
verify(src);need(sha(src/'SHA256SUMS.txt')==source_manifest,'Source release changed during run')
qc=dict(status='COMBINED_LEDGER_CONSOLIDATION_COMPLETE',resume=resume,group_statuses=groupqc,prior_ledger_records=len(old),new_locus_records=len(new),integrated_records=len(ledger),integrated_families=len({r['wang_family'] for r in ledger}),protein_evidence_records=len(proteinrows),unique_candidate_proteins=len({r['protein_id'] for r in proteinrows}),new_tiers=dict(Counter(r['count_tier'] for r in new)),cross_family_conflicts=len(conflicts),final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Remaining_family_ledger_consolidation/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle))
PY
