#!/usr/bin/env python3
import argparse,hashlib,json,shutil,sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--out',required=True);args=p.parse_args()
base=Path(__file__).resolve().parent; data=base/'data'; out=Path(args.out)
for n,h in json.loads((data/'input_hashes.json').read_text()).items():
 assert hashlib.sha256((data/n).read_bytes()).hexdigest()==h, 'Input changed: '+n
out.mkdir(parents=True,exist_ok=False); (out/'figures').mkdir();(out/'tables').mkdir();shutil.copytree(data,out/'inputs')
def write(d,n): d.to_csv(out/'tables'/n,sep='\t',index=False,na_rep='NA')
ledger=pd.read_csv(data/'01_supported_locus_evidence.tsv',sep='\t');ledger=ledger[ledger.species=='SS']
d=pd.read_csv(data/'04_DE_locus_contrast_evidence.tsv',sep='\t');d=d[(d.species=='SS')&(d.contrast_id=='SS_LEAF_VS_COTYLEDON_DAY')].copy()
assert set(d.numerator_group)=={'leaves_day'} and set(d.denominator_group)=={'cotyledons_day'}
assert not d.harmonized_gene_id.duplicated().any() and not ledger.harmonized_gene_id.duplicated().any()
assert set(d.harmonized_gene_id)==set(ledger.harmonized_gene_id)
assert d.set_index('harmonized_gene_id').wang_family.sort_index().equals(ledger.set_index('harmonized_gene_id').wang_family.sort_index())
for c in ['log2FoldChange_apeglm','padj_within_contrast','assessable','primary_DE']:d[c]=pd.to_numeric(d[c],errors='coerce')
# Preserve the upstream DE flag; p-adjustment alone is NOT substituted for the upstream decision.
d['display_status']=np.where(d.assessable!=1,'NOT_ASSESSABLE',np.where(d.primary_DE==1,'UPSTREAM_PRIMARY_DE','NOT_PRIMARY_DE'))
met=pd.read_csv(data/'tissue_comparison_exploratory.tsv',sep='\t')
assert met.metabolite.is_unique and (met[['leaf_mean','cotyledon_mean']]>0).all().all()
met['log2_leaf_cotyledon_mean_ratio']=np.log2(met.leaf_mean/met.cotyledon_mean)
met['metabolite_FDR_lt_005']=met.q_BH<.05
met['floor_caution']=np.where(met.metabolite=='Malic acid','cotyledon_values_near_possible_preprocessing_floor','preprocessing_not_fully_documented')
# Family-level enzyme reactions; no claim that each locus is active in this compartment/pathway.
relations=[
 ('MDH','Malic acid','substrate_in_written_reaction; reversible','malate + NAD+ = oxaloacetate + NADH + H+','1.1.1.37'),
 ('NAD-ME','Malic acid','substrate_in_written_reaction','malate + NAD+ = pyruvate + CO2 + NADH','1.1.1.39'),
 ('NAD-ME','Pyruvic acid','product_in_written_reaction','malate + NAD+ = pyruvate + CO2 + NADH','1.1.1.39'),
 ('NADP-ME','Malic acid','substrate_in_written_reaction','malate + NADP+ = pyruvate + CO2 + NADPH','1.1.1.40'),
 ('NADP-ME','Pyruvic acid','product_in_written_reaction','malate + NADP+ = pyruvate + CO2 + NADPH','1.1.1.40'),
 ('PPDK','Pyruvic acid','substrate_in_written_reaction; reversible','ATP + pyruvate + phosphate = AMP + phosphoenolpyruvate + diphosphate','2.7.9.1'),
 ('ASPAT','Glutaric acid, 2-oxo-','substrate_in_written_reaction; reversible','aspartate + 2-oxoglutarate = oxaloacetate + glutamate','2.6.1.1'),
 ('AlaAT','Glutaric acid, 2-oxo-','substrate_in_written_reaction; reversible','alanine + 2-oxoglutarate = pyruvate + glutamate','2.6.1.2'),
 ('AlaAT','Pyruvic acid','product_in_written_reaction; reversible','alanine + 2-oxoglutarate = pyruvate + glutamate','2.6.1.2')]
r=pd.DataFrame(relations,columns=['wang_family','metabolite','role','reaction','EC'])
r['source_url']=r.EC.map(lambda x:'https://iubmb.qmul.ac.uk/enzyme/EC'+x.replace('.','/')+'.html')
r['scope']='family-level reaction; no locus-specific activity/localization validation'
write(r,'reaction_mapping.tsv');write(d,'SS_daytime_DE_preserved.tsv');write(met,'metabolite_effects.tsv')
pairs=d.merge(r,on='wang_family',how='inner').merge(met,on='metabolite',how='left',validate='many_to_one')
pairs['direction_pattern']=np.where(pairs.assessable!=1,'NOT_ASSESSABLE',np.where(pairs.log2FoldChange_apeglm*pairs.log2_leaf_cotyledon_mean_ratio>0,'SAME_SIGN',np.where(pairs.log2FoldChange_apeglm*pairs.log2_leaf_cotyledon_mean_ratio<0,'OPPOSITE_SIGN','ZERO_OR_UNDEFINED')))
pairs['joint_threshold_support']=(pairs.assessable==1)&(pairs.primary_DE==1)&pairs.metabolite_FDR_lt_005
pairs['interpretation']='descriptive sign comparison; no correlation, causal inference, or expected substrate/product direction'
write(pairs,'gene_metabolite_concordance.tsv')
coverage=d.groupby('wang_family').agg(supported=('gene_id','size'),assessable=('assessable','sum'),primary_DE=('primary_DE','sum')).reset_index()
coverage['mapping_status']=np.where(coverage.wang_family.isin(r.wang_family),'DIRECT_REACTION_PANEL','CONTEXT_ONLY_NO_DIRECT_MAPPING_IN_THIS_PANEL');write(coverage,'family_coverage.tsv')
labels={'Malic acid':'Malate','Pyruvic acid':'Pyruvate','Citric acid':'Citrate','Fumaric acid':'Fumarate','Succinic acid':'Succinate','Glutaric acid, 2-oxo-':'2-Oxoglutarate'}
plt.rcParams.update({'font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
for fam,g in d.groupby('wang_family',sort=True):
 total_supported = len(g)
 g = g.loc[
     (g['assessable'] == 1)
     & np.isfinite(g['log2FoldChange_apeglm'])
 ].copy()
 if g.empty:
  continue
 g=g.sort_values('harmonized_gene_id'); relevant=r[r.wang_family==fam]; direct=len(relevant)>0
 mm=met[met.metabolite.isin(relevant.metabolite)] if direct else met
 fig,(ax,bx)=plt.subplots(1,2,figsize=(12,max(5.5,len(g)*.32+2.7)),gridspec_kw={'width_ratios':[1.35,1]})
 for j,(_,row) in enumerate(g.iterrows()):
  fc=row.log2FoldChange_apeglm
  valid=row.assessable==1 and np.isfinite(fc)
  if not valid: ax.text(0,j,'NA',ha='center',va='center',color='#777777');continue
  color='#b44335' if fc>0 else '#297da4'; primary=row.primary_DE==1
  ax.plot([0,fc],[j,j],color=color,alpha=.5)
  ax.scatter(fc,j,facecolor=color if primary else 'white',edgecolor=color,s=45,zorder=3)
  ax.annotate('q='+format(row.padj_within_contrast,'.2g'),(fc,j),xytext=(7 if fc >= 0 else -7,0),textcoords='offset points',ha='left' if fc >= 0 else 'right',va='center',fontsize=8,annotation_clip=False)
 ax.set_yticks(range(len(g)),g.gene_id); ax.invert_yaxis(); ax.axvline(0,color='gray',lw=.7)
 vals=g.loc[g.assessable==1,'log2FoldChange_apeglm'].dropna(); limit=max(1,float(vals.abs().max()) if len(vals) else 1); ax.set_xlim(-limit*1.85,limit*1.85)
 ax.set_xlabel('Gene expression: shrunken log2FC\nleaf / cotyledon, daytime');ax.set_title(f'A  {fam}: {len(g)} of {total_supported} supported genes assessable',loc='left',fontsize=10)
 for j,(_,row) in enumerate(mm.iterrows()):
  fc=row.log2_leaf_cotyledon_mean_ratio;color='#b44335' if fc>0 else '#297da4';signif=row.q_BH<.05
  bx.barh(j,fc,color=color,alpha=.8 if signif else .35)
  bx.annotate((' *' if signif else ' ns')+'; q='+format(row.q_BH,'.3f'),(fc,j),xytext=(5 if fc>=0 else -5,0),textcoords='offset points',ha='left' if fc>=0 else 'right',va='center',fontsize=9)
 bx.set_yticks(range(len(mm)),[labels[v] for v in mm.metabolite]); bx.invert_yaxis();bx.axvline(0,color='gray',lw=.7)
 lim=max(1,float(mm.log2_leaf_cotyledon_mean_ratio.abs().max()));bx.set_xlim(-lim*1.65,lim*1.65)
 bx.set_xlabel('Metabolites: log2 ratio of arithmetic means\nleaf / cotyledon, daytime');bx.set_title('B  '+('Reaction-associated metabolites' if direct else 'Metabolic context only'),loc='left')
 for axis in [ax,bx]:axis.spines[['top','right']].set_visible(False)
 fig.suptitle('SS daytime transcript–metabolite comparison | '+fam,fontsize=14)
 fig.text(.5,.025,'Unassessable genes omitted from plots; all genes retained in evidence tables. Filled dots: primary DE; open: not primary DE.\nMetabolite *: BH q<0.05; ns: q≥0.05. Gene and metabolite q-values come from separate test sets.\nFamily-level reaction mapping does not establish gene-specific activity, flux or localization. Malate ratio may reflect a preprocessing floor.',ha='center',fontsize=8)
 fig.subplots_adjust(left=.12,right=.96,bottom=.25,top=.86,wspace=.75)
 for ext in ['pdf','svg','png']:fig.savefig(out/'figures'/('SS_'+fam.replace('/','_')+'_concordance.'+ext),dpi=180)
 plt.close(fig)
methods='''Existing SS_LEAF_VS_COTYLEDON_DAY results were reused without refitting or re-adjusting p-values. Numerator leaves_day; denominator cotyledons_day. Expression effects are existing apeglm-shrunken log2 fold changes; original within-contrast adjusted p-values, assessability and primary_DE calls are retained (primary_DE may include an effect-size criterion). Missing/unassessable results are not zeros. Metabolite effects are log2 ratios of arithmetic relative-abundance means, not DESeq2 estimates. Existing exact-permutation/BH q-values across six metabolites were reused. No new significance tests or combined p-values were calculated. Nine family-level enzyme/metabolite links were defined from IUBMB reaction entries; direction is only as written, not inferred in vivo. MDH cofactor/subcellular specificity is not established by this broad family label. Amino acids remain excluded from measurements; names in reaction equations are context only. Other families are explicitly context-only: no automatic substrate assignment to transporters or regulators. Same/opposite sign is descriptive and is not evidence of causal regulation. Separate plants and only two overlapping conditions preclude sample-level gene/metabolite correlation. No sequences, RNA-seq alignment or differential-expression fitting rerun. Metabolite preprocessing/floors remain unresolved, and no nighttime metabolomics or flux data are available. Existing z-score figures remain unchanged.\n'''
(out/'METHODS.txt').write_text(methods)
qc=dict(status='COMPLETE_DESCRIPTIVE',genes=len(d),assessable=int(sum(d.assessable==1)),families=int(d.wang_family.nunique()),direct_family_mappings=int(r.wang_family.nunique()),mapped_gene_metabolite_rows=len(pairs),new_DE_tests=False,gene_metabolite_correlations=False,python=sys.version.split()[0],numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__)
(out/'QC.json').write_text(json.dumps(qc,indent=2));shutil.copy2(__file__,out/'run.py')
with (out/'SHA256SUMS.txt').open('w') as f:
 for path in sorted(out.rglob('*')):
  if path.is_file() and path.name!='SHA256SUMS.txt':f.write(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+str(path.relative_to(out))+'\n')
print(json.dumps(qc,indent=2));print('RESULTS:',out)
