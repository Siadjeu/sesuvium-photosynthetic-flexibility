#!/usr/bin/env bash
#SBATCH --job-name=beta_ca_compare
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=beta_ca_compare_%j.out
#SBATCH --error=beta_ca_compare_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
python3 - <<'PY'
import csv,hashlib,html,itertools,json,os,re,tarfile
from pathlib import Path
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
names=['main_full','main_core','domains_full','domains_core']
if os.environ.get('TREE_RELEASE'):choices=[Path(os.environ['TREE_RELEASE']).resolve(strict=True)]
else:choices=[p for p in (base/'step25kz17.13').glob('release*') if all((p/n/'COMPLETE.txt').is_file() for n in names)]
if len(choices)!=1:raise SystemExit('Set TREE_RELEASE to the directory containing all four scenarios; found '+str(len(choices)))
source=choices[0]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,rows,fields):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def parse(s):
 tokens=re.findall(r'[^(),:;\s]+|[(),:;]',s);i=0;nodes=[]
 def sub():
  nonlocal i
  children=[]
  if tokens[i]=='(':
   i+=1;children.append(sub())
   while tokens[i]==',':i+=1;children.append(sub())
   if tokens[i]!=')':raise ValueError('Invalid Newick')
   i+=1
  label=''
  if tokens[i] not in ':,);':label=tokens[i];i+=1
  length=0.0
  if tokens[i]==':':length=float(tokens[i+1]);i+=2
  leaves=set().union(*(c['leaves'] for c in children)) if children else {label}
  n=dict(label=label,length=length,children=children,leaves=leaves);nodes.append(n);return n
 tree=sub()
 if tokens[i:]!=[';']:raise ValueError('Unexpected Newick ending')
 return tree,nodes
def split_key(a,universe):
 b=universe-a
 if min(len(a),len(b))<2:return None
 return tuple(sorted(a)) if (len(a),tuple(sorted(a)))<(len(b),tuple(sorted(b))) else tuple(sorted(b))
def support(label):
 if not re.fullmatch(r'[0-9.]+/[0-9.]+',label):return ('NA','NA')
 return tuple(float(x) for x in label.split('/'))
trees={};inputs={};report=[]
for name in names:
 d=source/name
 for line in (d/'SHA256SUMS.txt').read_text().splitlines():
  h,f=line.split(maxsplit=1);p=d/f
  if not p.is_file() or sha(p)!=h:raise SystemExit('Checksum failed: '+name+'/'+f)
 for f in ['tree.treefile','TIP_LABEL_MAP.tsv','tree.iqtree']:inputs[name+'/'+f]=sha(d/f)
 mapping={r['safe_tip_id']:r for r in read(d/'TIP_LABEL_MAP.tsv')}
 tree,nodes=parse((d/'tree.treefile').read_text())
 leafnodes=[n for n in nodes if not n['children']]
 if len(leafnodes)!=len(mapping) or tree['leaves']!=set(mapping):raise SystemExit('Tree tips do not match map')
 for node in nodes:node['leaves']={mapping[k]['tip_id'] for k in node['leaves']}
 splits={}
 for node in nodes:
  key=split_key(node['leaves'],tree['leaves']);node['split']=key
  if key:
   if key in splits:raise SystemExit('Duplicate split; expected unrooted degree-three serialization')
   splits[key]=node
 text=(d/'tree.iqtree').read_text();qc=json.loads((d/'TREE_QC.json').read_text())
 models=re.findall(r'Best-fit model according to BIC:\s*(\S+)',text)
 info=re.findall(r'Number of parsimony informative sites:\s*(\d+)',text)
 warnings=[l.strip() for l in text.splitlines() if 'WARNING:' in l]
 report.append(dict(scenario=name,tips=len(mapping),sites=qc['alignment_columns'],informative_sites=info[0] if info else 'NA',model=models[0] if models else 'SEE_REPORT',warnings=' | '.join(warnings)))
 trees[name]=dict(tree=tree,nodes=nodes,mapping=mapping,splits=splits)
out=base/'step25kz17.14'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
out.mkdir(parents=True,exist_ok=False)
table(out/'01_tree_summary.tsv',report,list(report[0]))
comparison=[];branches=[]
for prefix in ['main','domains']:
 a,b=trees[prefix+'_full'],trees[prefix+'_core']
 if a['tree']['leaves']!=b['tree']['leaves']:raise SystemExit('Paired trees have different tip sets')
 sa,sb=a['splits'],b['splits'];shared=set(sa)&set(sb)
 comparison.append(dict(comparison=prefix+'_full_vs_core',tips=len(a['tree']['leaves']),full_splits=len(sa),core_splits=len(sb),shared_splits=len(shared),full_only=len(set(sa)-set(sb)),core_only=len(set(sb)-set(sa)),unweighted_RF=len(set(sa)^set(sb)),note='Includes weak and near-zero branches; topology differences are not all strong conflicts'))
 for num,key in enumerate(sorted(set(sa)|set(sb)),1):
  fa,fb=support(sa[key]['label']) if key in sa else ('NA','NA');ca,cb=support(sb[key]['label']) if key in sb else ('NA','NA')
  branches.append(dict(comparison=prefix,split_id=prefix+'_'+str(num),side_A=';'.join(key),side_B=';'.join(sorted(a['tree']['leaves']-set(key))),presence='SHARED' if key in shared else 'FULL_ONLY' if key in sa else 'CORE_ONLY',full_SH_aLRT=fa,full_UFBoot=fb,core_SH_aLRT=ca,core_UFBoot=cb,full_length=sa[key]['length'] if key in sa else 'NA',core_length=sb[key]['length'] if key in sb else 'NA'))
 a['shared']=b['shared']=shared
table(out/'02_topology_comparison.tsv',comparison,list(comparison[0]))
table(out/'03_branch_support_comparison.tsv',branches,list(branches[0]))
# Exact focal splits; absence means the exact grouping was not recovered.
focal=[]
for name,t in trees.items():
 targets={'SLP_PO_group':{'SLP_Hap1__SP3G00376','PO__evm.TU.LG08.1094','PO__evm.TU.LG22.838'}}
 if name.startswith('domains'):targets['AT_domain_pair']={'AT__gene-LOC130820426__domain_1','AT__gene-LOC130820426__domain_2'}
 for label,tips in targets.items():
  node=t['splits'].get(split_key(tips,t['tree']['leaves']));a,b=support(node['label']) if node else ('NA','NA')
  focal.append(dict(scenario=name,group=label,members=';'.join(sorted(tips)),present=node is not None,SH_aLRT=a,UFBoot=b))
table(out/'04_focal_groups.tsv',focal,list(focal[0]))
colors={'AT':'#7b3294','BV':'#52616b','MC':'#008080','PO':'#d2691e','SLP_Hap1':'#0072b2','SP_Hex':'#2e7d32','SP_Tet':'#8a7800','SS':'#c62828'}
def esc(s):return html.escape(str(s),quote=True)
def draw_pair(prefix):
 width=2500;height=1700;svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">','<rect width="100%" height="100%" fill="white"/>']
 def text(x,y,s,size=13,color='#222',weight='normal'):
  svg.append(f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{size}" fill="{color}" font-weight="{weight}">{esc(s)}</text>')
 def line(x1,y1,x2,y2,color='#555',stroke=1.4):svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{stroke}"/>')
 text(30,35,'BETA-CA | '+('Main isoform trees' if prefix=='main' else 'AT domain sensitivity trees'),24,weight='bold')
 text(30,63,'Unrooted trees shown in their stored orientation; the left-hand junction does not identify an ancestor.',15)
 text(30,87,'Orange: split absent from the paired tree. Grey: shared split. Branch labels: SH-aLRT / UFBoot.',15)
 text(30,111,'Different leaf orders reflect tree topology. Identical sequences remain separate loci; AT domains represent one locus.',15)
 # Use the same branch-length scale in both panels.
 maxdepth=0
 for suffix in ['full','core']:
  t=trees[prefix+'_'+suffix];counter=[0]
  def layout(n,depth):
   n['x']=depth
   if not n['children']:n['y']=180+counter[0]*25;counter[0]+=1
   else:
    for c in n['children']:layout(c,depth+c['length'])
    n['y']=sum(c['y'] for c in n['children'])/len(n['children'])
  layout(t['tree'],0);maxdepth=max(maxdepth,max(n['x'] for n in t['nodes']))
 scale=610/maxdepth
 for panel,suffix in enumerate(['full','core']):
  t=trees[prefix+'_'+suffix];offset=35+panel*1240
  text(offset,150,suffix.upper()+' ALIGNMENT',18,weight='bold')
  for n in t['nodes']:
   x=offset+n['x']*scale;y=n['y']
   if n['children']:
    ys=[c['y'] for c in n['children']];line(x,min(ys),x,max(ys))
    for c in n['children']:
     xx=offset+c['x']*scale;color='#d55e00' if c['split'] and c['split'] not in t['shared'] else '#555'
     line(x,c['y'],xx,c['y'],color)
     if c['children']:text(xx+2,c['y']-4,c['label'],9,color)
   else:
    label=t['mapping'][n['label']]['tip_id'];species=label.split('__')[0]
    text(x+7,y+4,label,12,colors.get(species,'#111'))
  y=1530;line(offset,y,offset+.2*scale,y,'#222',2);text(offset,y+23,'0.2 substitutions/site',12)
 text(30,1600,'Tip colours: '+', '.join(colors)+'; reference sequences in black.',14)
 text(30,1625,'These are review figures. Shared topology does not imply strong support, and coloured changes need not be well-supported conflicts.',14)
 text(30,1650,'The core dataset is a sensitivity analysis using reference occupancy, not an independently validated trimming procedure.',14)
 svg.append('</svg>');(out/(prefix+'_full_vs_core.svg')).write_text('\n'.join(svg))
for prefix in ['main','domains']:draw_pair(prefix)
(out/'README.txt').write_text('Tree comparison and review figures only; no inference rerun or final gene-count decision.\nSVG files are vector figures readable in a browser or Inkscape. Supports are written for every internal branch.\nSplit comparison is unrooted and uses full sequence identifiers, never only the safe T labels.\nRF comparisons are restricted to matching tip sets: main full/core and domain full/core. No direct RF between different tip sets is reported.\nExact split absence does not establish a strongly supported conflict. Identical-sequence branching and near-zero edges remain unresolved.\nRectangular tree orientation is for display only; no outgroup or biological root is inferred.\nAT domain fragments remain one locus. The 45-locus tree set is not a final family census.\n')
(out/'QC.json').write_text(json.dumps({'status':'COMPLETE_COMPARISON','input_hashes':inputs,'comparisons':comparison,'trees':report},indent=2)+'\n')
(out/'COMPLETE.txt').write_text('COMPARISON_COMPLETE\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.iterdir()):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+p.name+'\n')
bundle=out/'BETA_CA_tree_comparison_bundle.tar.gz'
with tarfile.open(bundle,'w:gz') as t:
 for p in sorted(out.iterdir()):
  if p.is_file() and p!=bundle:t.add(p,arcname='BETA_CA_tree_comparison/'+p.name)
print('Bundle relative to PROJECT_ROOT:',bundle.relative_to(root))
PY
