#!/usr/bin/env bash
#SBATCH --job-name=beta_ca_trees
#SBATCH --partition=all_cpu.p
#SBATCH --array=0-3%2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=beta_ca_trees_%A_%a.out
#SBATCH --error=beta_ca_trees_%A_%a.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from the project directory}"
python3 - <<'PY'
import csv,hashlib,json,os,re,shutil,subprocess,tarfile
from pathlib import Path
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
requirements=['QC.json','SHA256SUMS.txt','02_main_tip_map.tsv','03_sensitivity_tip_map.tsv','main/alignment.faa','main/diagnostic_reference_core.faa','domain_sensitivity/alignment.faa','domain_sensitivity/diagnostic_reference_core.faa']
if os.environ.get('ALIGNMENT_RELEASE'):choices=[Path(os.environ['ALIGNMENT_RELEASE']).resolve(strict=True)]
else:choices=[p for p in (base/'step25kz17.12').glob('release*') if p.is_dir() and all((p/f).is_file() for f in requirements)]
if len(choices)!=1:raise SystemExit('Set ALIGNMENT_RELEASE to one completed Step17.12 release; found '+str(len(choices)))
source=choices[0]
for name in requirements:
 if not (source/name).is_file() or not (source/name).stat().st_size:raise SystemExit('Missing source: '+name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
checks={}
for line in (source/'SHA256SUMS.txt').read_text().splitlines():
 h,name=line.split(maxsplit=1);checks[name.lstrip('*')]=h
for name in requirements:
 if name=='SHA256SUMS.txt':continue
 if name not in checks or sha(source/name)!=checks[name]:raise SystemExit('Source checksum mismatch: '+name)
qc=json.loads((source/'QC.json').read_text())
if qc['proposed_loci']!=45 or any(not x['reference_sites_consistent'] or x['inspection_flags'] for x in qc['alignments'].values()):raise SystemExit('Alignment review differs from the reviewed checkpoint')
scenarios=[('main_full','main','alignment.faa',51),('main_core','main','diagnostic_reference_core.faa',51),('domains_full','domain_sensitivity','alignment.faa',52),('domains_core','domain_sensitivity','diagnostic_reference_core.faa',52)]
i=int(os.environ.get('SLURM_ARRAY_TASK_ID','0'))
if not 0<=i<4:raise SystemExit('Array index must be 0..3')
name,group,filename,n=scenarios[i]
run_id=os.environ.get('TREE_RUN_ID') or os.environ.get('SLURM_ARRAY_JOB_ID') or os.environ.get('SLURM_JOB_ID')
if not run_id or not re.fullmatch(r'[A-Za-z0-9_-]+',run_id):raise SystemExit('Use SLURM, or set a safe TREE_RUN_ID')
out=base/'step25kz17.13'/('release_'+run_id)/name
out.mkdir(parents=True,exist_ok=True)
exe=os.environ.get('IQTREE_BIN') or next((x for x in ['iqtree2','iqtree3','iqtree'] if shutil.which(x)),None)
if not exe or not shutil.which(exe):raise SystemExit('IQ-TREE not found. Activate its environment or set IQTREE_BIN.')
version=subprocess.run([exe,'--version'],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True).stdout
if not re.search(r'version\s+[23]\.',version,re.I):raise SystemExit('This script requires IQ-TREE 2 or 3; inspect --version')
seqs={};key=None
for line in (source/group/filename).read_text().splitlines():
 if line.startswith('>'):
  key=line[1:].split()[0]
  if key in seqs:raise SystemExit('Duplicate input ID')
  seqs[key]=''
 elif line.strip():seqs[key]+=line.strip().upper()
if len(seqs)!=n or len({len(s) for s in seqs.values()})!=1:raise SystemExit('Unexpected alignment dimensions')
if any(not s.replace('-','') or any(c not in 'ACDEFGHIKLMNPQRSTVWYBXZJUO-' for c in s) for s in seqs.values()):raise SystemExit('Invalid protein alignment')
# Stable short labels avoid any Newick punctuation changes. A complete mapping is bundled.
ordered=sorted(seqs);safe={f'T{j:03d}':k for j,k in enumerate(ordered,1)}
input_text=''.join('>'+k+'\n'+seqs[v]+'\n' for k,v in safe.items())
threads=str(int(os.environ.get('SLURM_CPUS_PER_TASK','4')));seed=251713+i
args=[exe,'-s','input.faa','-st','AA','-m','MFP','-bb','1000','-bnni','-alrt','1000','-nt',threads,'-mem','6G','-seed',str(seed),'-keep-ident','-pre','tree']
state={'scenario':name,'source_alignment_sha256':sha(source/group/filename),'input_sha256':hashlib.sha256(input_text.encode()).hexdigest(),'seed':seed,'command_options':args[1:],'iqtree_version':version,'tips':n,'columns':len(next(iter(seqs.values()))),'rooting':'UNROOTED_NO_BIOLOGICAL_OUTGROUP_SPECIFIED','status':'RESEARCH_TREE_NOT_FINAL_FAMILY_COUNTS'}
state_path=out/'RUN_CONFIG.json'
if state_path.exists():
 if json.loads(state_path.read_text())!=state:raise SystemExit('Existing run configuration differs; use a new TREE_RUN_ID')
else:
 state_path.write_text(json.dumps(state,indent=2)+'\n');(out/'input.faa').write_text(input_text)
 if not (out/'input.faa').is_file():raise SystemExit('Input creation failed')
if sha(out/'input.faa')!=state['input_sha256']:raise SystemExit('Existing input altered')
mapfile='02_main_tip_map.tsv' if group=='main' else '03_sensitivity_tip_map.tsv'
with (source/mapfile).open() as f:maps={r['tip_id']:r for r in csv.DictReader(f,delimiter='\t')}
if set(maps)!=set(seqs):raise SystemExit('Tip mapping mismatch')
with (out/'TIP_LABEL_MAP.tsv').open('w',newline='') as f:
 fields=['safe_tip_id']+list(next(iter(maps.values())))
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
 for k,v in safe.items():w.writerow({'safe_tip_id':k,**maps[v]})
shutil.copyfile(source/group/'column_occupancy.tsv',out/'source_alignment_column_occupancy.tsv')
# Completed matching runs are verified and repackaged; interrupted runs resume via IQ-TREE's checkpoint.
if not (out/'COMPLETE.txt').is_file():
 with (out/'console.txt').open('a') as f:subprocess.run(args,cwd=out,stdout=f,stderr=subprocess.STDOUT,check=True)
for filename2 in ['tree.treefile','tree.iqtree','tree.log','tree.contree']:
 if not (out/filename2).is_file() or not (out/filename2).stat().st_size:raise SystemExit('Missing IQ-TREE output: '+filename2)
def label_tree(filename2):
 text=(out/filename2).read_text()
 tips=re.findall(r'(?<=[(,])\s*(T\d+)(?=[:),;\s])',text)
 if len(tips)!=n or set(tips)!=set(safe):raise SystemExit('Tree tip verification failed: '+filename2)
 # Quote all full labels in Newick; doubled apostrophes are the Newick escape convention.
 text=re.sub(r'(?<=[(,])\s*(T\d+)(?=[:),;\s])',lambda m:"'"+safe[m.group(1)].replace("'","''")+"'",text)
 (out/(filename2+'.labeled.nwk')).write_text(text)
label_tree('tree.treefile');label_tree('tree.contree')
report=(out/'tree.iqtree').read_text()
models=re.findall(r'^.*Best-fit model.*$',report,re.M)
(out/'TREE_QC.json').write_text(json.dumps({'status':'COMPLETE_TIPS_VERIFIED','scenario':name,'tips':n,'alignment_columns':state['columns'],'model_report_lines':models,'support_order':'SH-aLRT / ultrafast bootstrap; verify tree.iqtree','rooting':'UNROOTED','notes':['Identical sequence labels are retained; their relative branching is not resolved by sequence evidence.','AT domain tips represent one locus, not two genes.','Core filtering is a diagnostic sensitivity analysis, not a validated universal trimming rule.','Do not compare raw likelihoods between different alignments.']},indent=2)+'\n')
(out/'README.txt').write_text('Main and domain-sensitivity trees each have full and diagnostic-core variants. No extra column filtering or realignment is applied.\nCore columns were selected in Step17.12 by occupancy in at least 5 of 6 references. source_alignment_column_occupancy.tsv maps that selection in the full alignment.\nThe full tree is a control for the core result. Interpret robust shared clades; inspect disagreements before biological conclusions.\nNo alpha-CA outgroup is introduced. Newick serialization does not constitute biological rooting.\nIdentical tips are preserved; AT domain fragments are not additional loci.\nThe checkpoint remains on the cluster for restart. TREE_RUN_ID can identify the same release on resubmission; configuration must match.\n')
(out/'COMPLETE.txt').write_text('TREE_RUN_COMPLETE_TIPS_VERIFIED\n')
# Small review bundle excludes the resumable binary checkpoint and large bootstrap samples.
files=[p for p in sorted(out.iterdir()) if p.is_file() and not p.name.endswith(('.tar.gz','.ckp.gz','.ufboot','.boottrees')) and p.name!='SHA256SUMS.txt']
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in files:f.write(sha(p)+'  '+p.name+'\n')
bundle=out/('BETA_CA_tree_'+name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='BETA_CA_tree_'+name+'/'+p.name)
print('Bundle relative to PROJECT_ROOT:',bundle.relative_to(root))
PY
