#!/usr/bin/env python3
"""Independent-cohort SS transcript/metabolite presentation. No cross-omics correlations."""
import argparse,csv,hashlib,json,itertools,shutil,sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser(); p.add_argument('--out',required=True); a=p.parse_args()
base=Path(__file__).resolve().parent; out=Path(a.out); out.mkdir(parents=True,exist_ok=False)
(out/'figures').mkdir(); (out/'tables').mkdir(); shutil.copytree(base/'data',out/'inputs')
def save(df,name,index=False): df.to_csv(out/'tables'/name,sep='\t',index=index,na_rep='NA')
with (base/'data/metabolites.csv').open(encoding='utf-8-sig') as f: raw=list(csv.reader(f))
names=[s.strip() for s in raw[0][1:]]; ids=[r[0].strip() for r in raw[1:]]
x=np.array([[float(v) for v in r[1:]] for r in raw[1:]])
assert np.isfinite(x).all() and (x>=0).all(), 'Missing/negative abundance requires explicit handling'
sel=[i for i,s in enumerate(ids) if s in ['C1','C2','C3'] or s in [f'SS{j}' for j in range(1,7)] or s in [f'SSS{j}' for j in range(1,7)]]
samples=[ids[i] for i in sel]; tissue=np.array(['cotyledon' if s.startswith('C') else 'leaf' for s in samples]); xss=x[sel]
assert len(set(samples))==9 and sum(tissue=='cotyledon')==3 and sum(tissue=='leaf')==6
# Explicit conservative exclusion of amino acids and related amino compounds/derivatives.
terms=['alanine','asparagin','aspartic','cysteine','glutamic','glutamine','glycine','homoserine','isoleucine','leucine','lysine','methionine','ornithine','phenylalanine','proline','serine','threonine','tryptophan','tyrosine','valine','dihydroxyphenylalanine','pyroglutamic','piperidinecarboxylic','adipic_acid,2-amino','putrescine','phenethylamine','ethanolamine']
audit=[]; keep=[]
for j,n in enumerate(names):
 reasons=[]
 if any(t in n.lower() for t in terms): reasons.append('amino_acid_or_related_compound')
 if names.count(n)>1: reasons.append('duplicate_name_ambiguous')
 if any(k!=j and np.array_equal(x[:,j],x[:,k]) for k in range(len(names))): reasons.append('identical_column_ambiguous_identity')
 audit.append({'column':j+2,'name':n,'decision':'EXCLUDE' if reasons else 'RETAIN','reason':';'.join(reasons)})
 if not reasons: keep.append(j)
save(pd.DataFrame(audit),'column_audit.tsv')
save(pd.DataFrame({'sample_id':samples,'tissue':tissue,'period':'day','independent_biological_replicate':True}),'metabolite_samples.tsv')
save(pd.DataFrame(xss[:,keep],index=samples,columns=[names[j] for j in keep]),'retained_metabolites.tsv',True)
qc=[]
for i,s in enumerate(samples):
 vals,n=np.unique(xss[i],return_counts=True); k=np.argmax(n)
 qc.append({'sample':s,'zero_count':int(sum(xss[i]==0)),'most_repeated_value':vals[k],'repetitions':int(n[k]),'interpretation':'possible_preprocessing_floor_unconfirmed'})
save(pd.DataFrame(qc),'repeated_values_QC.tsv')
# Predefined small non-amino-acid panel; contextual organic acids are not all CAM-specific.
target=['Malic acid','Pyruvic acid','Citric acid','Fumaric acid','Succinic acid','Glutaric acid, 2-oxo-']
assert all(n in [names[j] for j in keep] for n in target)
y=xss[:,[names.index(n) for n in target]]; labels=['Malate','Pyruvate','Citrate','Fumarate','Succinate','2-Oxoglutarate']
save(pd.DataFrame(y,index=samples,columns=target),'target_metabolites.tsv',True)
# Exhaustive 3-versus-6 label permutations of absolute arithmetic-mean difference (84).
assign=list(itertools.combinations(range(9),3)); contrasts=[]
for k,n in enumerate(target):
 v=y[:,k]; c=v[tissue=='cotyledon']; l=v[tissue=='leaf']; obs=abs(l.mean()-c.mean())
 null=[]
 for idx in assign:
  mask=np.zeros(9,bool); mask[list(idx)]=True; null.append(abs(v[~mask].mean()-v[mask].mean()))
 contrasts.append(dict(metabolite=n,cotyledon_mean=c.mean(),leaf_mean=l.mean(),leaf_minus_cotyledon=l.mean()-c.mean(),ratio_leaf_cotyledon=l.mean()/c.mean() if c.mean()>0 else np.nan,p_exact=np.mean(np.array(null)>=obs-1e-10)))
def bh(ps):
 ps=np.asarray(ps); order=np.argsort(ps); q=np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1]; ans=np.empty(len(ps)); ans[order]=np.minimum(q,1); return ans
ct=pd.DataFrame(contrasts); ct['q_BH']=bh(ct.p_exact); save(ct,'tissue_comparison_exploratory.tsv')
# Rank residualization controls tissue intercept. Exact within-tissue permutations: 3!*6!=4320.
r=np.apply_along_axis(rankdata,0,y); residual=r.copy()
for t in ['cotyledon','leaf']: residual[tissue==t]-=residual[tissue==t].mean(axis=0)
def corr(z):
 norms=np.linalg.norm(z,axis=0); return np.divide(z.T@z,np.outer(norms,norms),out=np.full((z.shape[1],z.shape[1]),np.nan),where=np.outer(norms,norms)>1e-12)
pooled=corr(r-r.mean(axis=0)); adjusted=corr(residual)
ci=np.where(tissue=='cotyledon')[0]; li=np.where(tissue=='leaf')[0]
perms=[]
for c in itertools.permutations(ci):
 for l in itertools.permutations(li):
  order=np.arange(9); order[ci]=c; order[li]=l; perms.append(order)
perms=np.array(perms); pairs=[]
for i in range(6):
 for j in range(i+1,6):
  den=np.linalg.norm(residual[:,i])*np.linalg.norm(residual[:,j]); pv=np.nan
  if den>1e-12:
   null=(residual[perms,j]@residual[:,i])/den; pv=float(np.mean(abs(null)>=abs(adjusted[i,j])-1e-12))
  pairs.append(dict(metabolite_1=target[i],metabolite_2=target[j],pooled_spearman=pooled[i,j],tissue_adjusted_rank_correlation=adjusted[i,j],p_exact_within_tissue=pv))
pt=pd.DataFrame(pairs); pt['q_BH']=np.nan; ok=pt.p_exact_within_tissue.notna(); pt.loc[ok,'q_BH']=bh(pt.loc[ok,'p_exact_within_tissue']); save(pt,'metabolite_correlations_exploratory.tsv')
counts=pd.read_csv(base/'data/SS_supported_normalized_counts.tsv',sep='\t',index_col=0)
meta=pd.read_csv(base/'data/SS_sample_annotation.tsv',sep='\t'); ledger=pd.read_csv(base/'data/01_supported_locus_evidence.tsv',sep='\t'); ledger=ledger[ledger.species=='SS']
assert counts.index.is_unique and np.isfinite(counts.to_numpy()).all() and (counts.to_numpy()>=0).all()
groups=[('leaf','DAY'),('cotyledon','DAY'),('cotyledon','NIGHT')]; means=[]
for t,d in groups:
 cols=meta.loc[(meta.tissue==t)&(meta.day_night==d),'sample_id'].tolist(); assert len(cols)==3
 means.append(counts[cols].mean(axis=1))
m=pd.concat(means,axis=1); m.columns=['leaf_day','cotyledon_day','cotyledon_night']; sd=m.std(axis=1,ddof=0); z=m.sub(m.mean(axis=1),axis=0).div(sd.replace(0,np.nan),axis=0)
save(m,'expression_group_means.tsv',True); save(z,'expression_row_zscore.tsv',True)
plt.rcParams.update({'font.size':9,'svg.fonttype':'none','pdf.fonttype':42})
colors=['#B67842','#268577']; cmap=plt.get_cmap('RdBu_r').copy(); cmap.set_bad('#cccccc')
for fam, rows in ledger.groupby('wang_family',sort=True):
 genes=rows.harmonized_gene_id.drop_duplicates().tolist(); missing=sorted(set(genes)-set(z.index)); assert not missing,(fam,missing)
 genes=sorted(genes); n=len(genes); height=max(9,n*.30+3)
 fig=plt.figure(figsize=(16,height)); gs=fig.add_gridspec(2,3,width_ratios=[1.35,1.35,1.4],height_ratios=[1,1],wspace=.6,hspace=.65)
 ax=fig.add_subplot(gs[:,0]); im=ax.imshow(z.loc[genes],vmin=-2.5,vmax=2.5,cmap=cmap,aspect='equal')
 ax.set_xticks(range(3),['Leaf\nday','Cotyledon\nday','Cotyledon\nnight'],rotation=45,ha='right'); ax.set_yticks(range(n),[g.replace('SS__','') for g in genes],fontsize=8)
 ax.set_xticks(np.arange(-.5,3,1),minor=True); ax.set_yticks(np.arange(-.5,n,1),minor=True); ax.grid(which='minor',color='black',lw=.45); ax.tick_params(which='minor',length=0)
 ax.set_title('A  '+fam+' expression\nMean of 3 RNA-seq replicates',loc='left'); fig.colorbar(im,ax=ax,shrink=.35,label='Row z-score (all 3 groups)')
 sub=gs[:,1].subgridspec(3,2,hspace=.65,wspace=.65)
 for k in range(6):
  ax=fig.add_subplot(sub[k//2,k%2])
  for g,t in enumerate(['cotyledon','leaf']):
   v=y[tissue==t,k]; ax.scatter(g+np.linspace(-.12,.12,len(v)),v,color=colors[g],s=25,zorder=3); ax.errorbar(g,v.mean(),yerr=v.std(ddof=1),fmt='_',color='black',capsize=4,ms=15)
  ax.set_xticks([0,1],['Cotyledon','Leaf'],rotation=35,ha='right'); ax.set_title(('B  ' if k==0 else '')+labels[k]); ax.set_ylabel('Relative abundance'); ax.set_xlim(-.5,1.5)
  ax.spines[['top','right']].set_visible(False)
  q = float(ct.loc[k, 'q_BH'])
  mark = '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'ns'
  ax.plot(
      [0, 0, 1, 1], [0.86, 0.90, 0.90, 0.86],
      transform=ax.get_xaxis_transform(),
      color='black', linewidth=0.8
  )
  ax.text(
      0.5, 0.92, f'{mark}  (BH-adjusted p={q:.3f})',
      transform=ax.get_xaxis_transform(),
      ha='center', va='bottom', fontsize=8
  )
  lower, upper = ax.get_ylim()
  ax.set_ylim(lower, lower + (upper - lower) * 1.35)
 for slot,mat,title in [(gs[0,2],pooled,'C  Pooled Spearman'),(gs[1,2],adjusted,'D  Tissue-adjusted rank correlation')]:
  ax=fig.add_subplot(slot); im=ax.imshow(mat,vmin=-1,vmax=1,cmap=cmap); ax.set_xticks(range(6),labels,rotation=60,ha='right',fontsize=8); ax.set_yticks(range(6),labels,fontsize=8); ax.set_title(title,fontsize=10,loc='left')
  for i in range(6):
   for j in range(6): ax.text(j,i,'NA' if not np.isfinite(mat[i,j]) else f'{mat[i,j]:.2f}',ha='center',va='center',fontsize=7,color='white' if abs(mat[i,j])>.65 else 'black')
  fig.colorbar(im,ax=ax,shrink=.6,label='Correlation')
 fig.suptitle('SS | '+fam+' expression and daytime metabolite context',fontsize=15,y=.99)
 fig.text(.5,.015,'Independent cohorts; no gene–metabolite sample pairing. Metabolomics: cotyledon n=3, leaf n=6.\nDots: biological replicates; bars: mean ± SD. Correlations exploratory; see permutation tests and FDR tables.\nShared metabolite panel provides pathway context, not proof of family-specific function or CAM flux.',ha='center',fontsize=9)
 fig.subplots_adjust(top=.91,bottom=.14,left=.1,right=.97)
 name=fam.replace('/','_')
 for ext in ['pdf','svg','png']: fig.savefig(out/'figures'/f'SS_{name}_combined.{ext}',dpi=180)
 plt.close(fig)
methods='''Expression: arithmetic means of three DESeq2-normalized biological replicates in each of three groups; row z-score across all three means, population SD (ddof=0), no log transformation. Zero-variance genes are grey; display clipped at ±2.5, source tables unclipped. Metabolites: relative abundances, untransformed replicate dot plots, mean ± sample SD; no pseudocount or inferred detection-limit replacement. Independent cohorts preclude sample-level gene–metabolite correlations. Amino acids/related compounds and ambiguous duplicated columns excluded with audit. Six predefined organic acids used without selection based on observed significance. Tissue contrasts: absolute difference of arithmetic means, all 84 group-label allocations, two-sided exact randomization test under exchangeability, BH adjustment across six tests. Pooled metabolite correlations: Spearman, descriptive only. Tissue-adjusted correlations: Pearson correlations of global ranks after subtracting within-tissue rank means; exact permutation of one variable within each tissue (4320 permutations), two-sided absolute statistic; BH across estimable pairs. These exploratory results are limited by n=9, unconfirmed preprocessing/detection floors, and tissue/cohort effects. No night metabolomics, flux inference, or gene-function validation. Same metabolite panel repeated for each family is contextual, not a family-specific target assignment. No RNA-seq differential testing is performed. Metabolite cohort label permutation does not remove unrecorded batch effects. Correlation adjustment removes tissue intercepts only, not all confounding.\n'''
(out/'METHODS.txt').write_text(methods)
qc={'status':'COMPLETE_EXPLORATORY','RNA_genes':len(z),'families':ledger.wang_family.nunique(),'metabolite_samples':len(samples),'target_metabolites':target,'excluded_columns':sum(a['decision']=='EXCLUDE' for a in audit),'gene_metabolite_correlation_performed':False,'preprocessing_floor_unresolved':True,'python':sys.version.split()[0],'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__}
(out/'QC.json').write_text(json.dumps(qc,indent=2)); shutil.copy2(__file__,out/'run.py')
with (out/'SHA256SUMS.txt').open('w') as f:
 for path in sorted(out.rglob('*')):
  if path.is_file() and path.name!='SHA256SUMS.txt': f.write(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+str(path.relative_to(out))+'\n')
print(json.dumps(qc,indent=2)); print('RESULTS:',out)
