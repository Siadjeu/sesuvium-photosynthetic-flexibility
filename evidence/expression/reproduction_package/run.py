#!/usr/bin/env python3
import csv,json,hashlib,shutil,tarfile,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
P=Path(__file__).resolve().parent
plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.size':9})
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def verify(p):
 for l in (p/'SHA256SUMS.txt').read_text().splitlines():
  h,n=l.split(None,1);f=p/n.lstrip('*');assert f.resolve().is_relative_to(p.resolve());assert sha(f)==h,n
root=Path(os.environ['PROJECT_ROOT']).resolve()
if os.environ.get('PRESENTATION_RELEASE'):src=Path(os.environ['PRESENTATION_RELEASE']).resolve()
else:
 candidates=[p for p in (root/'08_focused_families/downstream_supported/presentation').glob('release_*') if p.is_dir() and (p/'ALL_TREE_ID_CORRESPONDENCE.tsv').is_file() and (p/'log2FC').is_dir()]
 assert len(candidates)==1,'Set PRESENTATION_RELEASE explicitly: expected exactly one previous presentation release'
 src=candidates[0]
verify(src)
for l in (P/'PACKAGE_SHA256SUMS.txt').read_text().splitlines():
 h,n=l.split(None,1);assert sha(P/n)==h,n
out=root/'08_focused_families/downstream_supported/presentation'/('mean_zscore_v4_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
out.mkdir(parents=True,exist_ok=False)
# Selective export deliberately excludes the old log2FC folder and its plotting scripts.
shutil.copytree(src/'trees',out/'trees')
shutil.copyfile(src/'ALL_TREE_ID_CORRESPONDENCE.tsv',out/'ALL_TREE_ID_CORRESPONDENCE.tsv')
shutil.copytree(P,out/'reproduction_package')
ledger=read(P/'data/01_supported_locus_evidence.tsv');qc=[];manifest=[];members=[]
colors={'leaf':'#4daf4a','stem':'#377eb8','cotyledon':'#984ea3','drought':'#e41a1c','well_watered':'#377eb8','08H':'#fee8c8','16H':'#fdbb84','24H':'#e34a33','DAY':'#fdd835','NIGHT':'#3949ab'}
cell=.30 # inches: identical square cells in PO-PEPC and every other panel
for sp in ['PO','SS']:
 annotations=read(P/f'data/{sp}_sample_annotation.tsv');counts=read(P/f'data/{sp}_supported_normalized_counts.tsv');byid={x['gene_id']:x for x in counts};groups={}
 fields=['tissue','condition','time','day_night'] if sp=='PO' else ['tissue','day_night']
 for x in annotations:
  key=tuple(x[k] for k in fields);groups.setdefault(key,[]).append(x['sample_id'])
 order={'leaf':0,'stem':1,'cotyledon':1,'well_watered':0,'drought':1,'08H':0,'16H':1,'24H':2,'DAY':0,'NIGHT':1}
 keys=sorted(groups,key=lambda k:tuple(order[v] for v in k));assert all(len(groups[k])==3 for k in keys)
 samples=[s for k in keys for s in groups[k]];assert len(samples)==len(set(samples));assert set(samples)==set(counts[0])-{'gene_id'}
 names=['_'.join(k) for k in keys]
 for k,n in zip(keys,names):
  for sample in groups[k]:members.append(dict(species=sp,group=n,sample_id=sample,replicate_count=3))
 for fam in sorted({x['wang_family'] for x in ledger}):
  loci=sorted([x for x in ledger if x['species']==sp and x['wang_family']==fam],key=lambda x:x['gene_id'])
  if not loci:continue
  means=[];zs=[];mr=[];zr=[]
  for x in loci:
   row=byid[x['harmonized_gene_id']];raw=np.array([float(row[s]) for s in samples]);assert np.isfinite(raw).all() and (raw>=0).all()
   m=np.array([np.mean([float(row[s]) for s in groups[k]]) for k in keys]);sd=m.std(ddof=0);z=(m-m.mean())/sd if sd>0 else np.full(len(m),np.nan)
   means.append(m);zs.append(z);mr.append(dict(gene_id=x['gene_id'],**dict(zip(names,m))));zr.append(dict(gene_id=x['gene_id'],**{n:float(v) if np.isfinite(v) else 'NA' for n,v in zip(names,z)}))
   qc.append(dict(family=fam,species=sp,gene_id=x['gene_id'],groups=len(keys),replicates_per_group=3,status='PLOTTED' if sd>0 else 'ZERO_VARIANCE_MASKED'))
  dest=out/'expression_mean_zscore'/fam;dest.mkdir(parents=True,exist_ok=True);write(dest/f'{sp}_group_means.tsv',mr);write(dest/f'{sp}_group_mean_zscores.tsv',zr)
  # Explicit physical axes geometry, no tight-layout or tight bounding box resizing.
  nr,nc=len(loci),len(keys);left=2.5;bottom=2.2;top=2.25;right=1.45
  w=left+nc*cell+right;h=bottom+nr*cell+top;fig=plt.figure(figsize=(w,h))
  ax=fig.add_axes([left/w,bottom/h,nc*cell/w,nr*cell/h]);cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#d9d9d9')
  im=ax.imshow(np.ma.masked_invalid(zs),aspect='equal',interpolation='nearest',vmin=-2.5,vmax=2.5,cmap=cmap)
  display=lambda value: {'DAY':'day','NIGHT':'night','well_watered':'well-watered'}.get(value,value)
  sample_labels=['_'.join(display(v) for v in k) for k in keys]
  ax.set_yticks(range(nr),[sp+' | '+x['gene_id'] for x in loci]);ax.set_xticks(range(nc),sample_labels,rotation=60,ha='right',fontsize=8)
  ax.set_xticks(np.arange(nc+1)-.5,minor=True);ax.set_yticks(np.arange(nr+1)-.5,minor=True)
  ax.grid(which='minor',color='black',linestyle='-',linewidth=.45)
  ax.tick_params(which='minor',bottom=False,left=False)
  bar=fig.add_axes([left/w,(bottom+nr*cell+.08)/h,nc*cell/w,len(fields)*.14/h])
  for i,field in enumerate(fields):
   for j,k in enumerate(keys):bar.add_patch(plt.Rectangle((j,i),1,1,facecolor=colors[k[i]],edgecolor='black',linewidth=.45))
  bar.set_xlim(0,nc);bar.set_ylim(len(fields),0);bar.set_xticks([]);bar.set_yticks(np.arange(len(fields))+.5,[{'tissue':'Tissue','condition':'Condition','time':'Time','day_night':'Period'}[f] for f in fields],fontsize=8)
  cax=fig.add_axes([(left+nc*cell+.22)/w,bottom/h,.16/w,min(nr*cell,1.8)/h]);fig.colorbar(im,cax=cax,label='Row z-score')
# Figure-coordinate legend uses physical inches for consistent row/column alignment.
  legend_rows=[('Tissue',['leaf','stem'] if sp=='PO' else ['leaf','cotyledon'])]
  if sp=='PO':legend_rows += [('Condition',['drought','well_watered']),('Time',['08H','16H','24H'])]
  legend_rows += [('Period',['DAY','NIGHT'])]
  for ri,(caption,values) in enumerate(legend_rows):
   yy=h-.65-ri*.28
   fig.text(.55/w,yy/h,caption+':',ha='left',va='center',fontsize=9)
   for ci,value in enumerate(values):
    xx=1.65+ci*1.10
    fig.add_artist(plt.Rectangle((xx/w,(yy-.065)/h),.18/w,.13/h,transform=fig.transFigure,facecolor=colors[value],edgecolor='black',linewidth=.6))
    fig.text((xx+.24)/w,yy/h,display(value),ha='left',va='center',fontsize=9)
  fig.suptitle(f'{fam} | {sp} | mean of 3 replicates',y=1-.12/h,fontsize=11)
  fig.text(.12/w,.15/h,'Z-score across group means; grey = zero variance.\nColours limited to ±2.5; square cells = 0.30 inches.',fontsize=8)
  prefix=dest/f'{sp}_{fam}_mean_zscore'
  fig.canvas.draw();bbox=ax.get_window_extent();assert abs(bbox.width/nc-bbox.height/nr)<1e-6
  for ext in ['svg','pdf','png']:fig.savefig(str(prefix)+'.'+ext,dpi=200)
  plt.close(fig);manifest.append(dict(family=fam,species=sp,genes=nr,groups=nc,cell_width_in=cell,cell_height_in=cell,prefix=str(prefix.relative_to(out))))
write(out/'GROUP_MEMBERSHIP.tsv',members);write(out/'EXPRESSION_QC.tsv',qc);write(out/'FIGURE_INDEX.tsv',manifest)
assert len(manifest)==36 and len(qc)==226
shutil.copyfile(P/'HANDOVER.md',out/'WORKFLOW_HANDOVER.md')
(out/'QC.json').write_text(json.dumps(dict(status='COMPLETE',panels=len(manifest),expression_loci=len(qc),zero_variance_rows=sum(x['status']=='ZERO_VARIANCE_MASKED' for x in qc),source_release=src.name,source_manifest_sha256=sha(src/'SHA256SUMS.txt'),log2FC_folder_present=False,previous_release_modified=False,cell_inches=cell),indent=2)+'\n')
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in sorted(out.rglob('*')) if p.is_file() and p!=out/'SHA256SUMS.txt'))
verify(out);verify(src)
with tarfile.open(str(out)+'_bundle.tar.gz','w:gz') as t:t.add(out,arcname=out.name)
print((out/'QC.json').read_text());print('Upload: '+str(out)+'_bundle.tar.gz')
