#!/usr/bin/env python3
"""Integrate existing evidence; Python 3 standard library only. No analyses rerun."""
import argparse,csv,gzip,hashlib,json,math,shutil
from pathlib import Path
from collections import defaultdict,Counter
P=argparse.ArgumentParser();P.add_argument('--source-release',type=Path,required=True);P.add_argument('--output',type=Path,required=True);a=P.parse_args();root=a.source_release.resolve();out=a.output
assert not out.exists(),'Output exists; choose a new directory'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
for line in (root/'SHA256SUMS.txt').read_text().splitlines():
 digest,name=line.split(maxsplit=1);p=(root/name.lstrip('*')).resolve();assert p.is_relative_to(root);assert sha(p)==digest,'Checksum mismatch: '+name

def rows(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f: yield from csv.DictReader(f,delimiter='\t')
def find(name):
 hits=list((root/'sources').rglob(name));assert len(hits)==1,(name,len(hits));return hits[0]
def write(name,data,fields=None):
 data=list(data);fields=fields or list(data[0]);
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
def finite(x):
 try:return math.isfinite(float(x))
 except (ValueError,TypeError):return False
def unique(data,key):
 d={}
 for x in data:
  k=key(x);assert k not in d,('Duplicate key',k);d[k]=x
 return d
species='AT BV MC PO SLP_Hap1 SP_Hex SP_Tet SS'.split()
families='beta-CA PEPC PPCK NAD-ME NADP-ME MDH ASPAT AlaAT PEPCK PPDK PPDK-RE BASS2 MEP DIT2 DIC DTC MPC ALMT TDT PPA PPT V_ATPase V_Ppase'.split()
ledger=list(rows(root/'01_supported_locus_ledger.tsv'));assert len(ledger)==657
unique(ledger,lambda r:(r['wang_family'],r['species'],r['gene_id']))
assert all(r['count_tier']=='SUPPORTED' for r in ledger),Counter(r['count_tier'] for r in ledger)
maps={}
for sp in species:
 for r in rows(find(sp+'.representative_gene_map.tsv')):
  k=(sp,r['original_gene_id']);assert k not in maps;maps[k]=r
og={}
for r in rows(find('Orthogroups.tsv')):
 for sp in species:
  for gene in filter(None,r[sp].split(', ')):
   assert (sp,gene) not in og;og[sp,gene]=r['Orthogroup']
wgd=unique(rows(find('all_species_gene_duplication_classes.tsv')),lambda r:(r['species'],r['gene_id']))
cafe=unique(rows(find('all_species_terminal_model_comparison.tsv')),lambda r:(r['terminal_species'],r['family_id']))
ps=defaultdict(list)
for r in rows(find('absrel_functional_integration.tsv')):ps[r['foreground_species'],r['candidate_gene_id']].append(r)
rnamap=unique(rows(find('orthofinder_to_annotation_gene_map.tsv')),lambda r:(r['species'],r['orthofinder_gene_id']))
contrasts=[r for r in rows(find('CONTRAST_REGISTRY.tsv')) if r['status']=='AUTHORIZED_AFTER_VIGNETTE_QC'];assert len(contrasts)==14
wanted={}
for r in ledger:
 m=maps.get((r['species'],r['gene_id']));
 if m:wanted[r['species'],m['harmonized_gene_id']]=True
de={}
for c in contrasts:
 for r in rows(find(c['contrast_id']+'.deseq2.tsv.gz')):
  assert r['contrast_id']==c['contrast_id'] and r['species_code']==c['species_code']
  if (r['species_code'],r['gene_id']) in wanted:
   k=(r['species_code'],r['gene_id'],r['contrast_id']);assert k not in de;de[k]=r
result=[];events=[];pslong=[];cafelong=[];issues=[]
for r in ledger:
 sp=r['species'];gene=r['gene_id'];fam=r['wang_family'];base={'wang_family':fam,'species':sp,'gene_id':gene};m=maps.get((sp,gene),{});hid=m.get('harmonized_gene_id','NA');k=(sp,hid);group=og.get(k,'NA');w=wgd.get(k,{});c=cafe.get((sp,group),{});tests=ps.get(k,[])
 if w and w['orthogroup_assigned']=='YES':assert w['orthogroup']==group,(k,'OG disagreement')
 valid=[x for x in tests if x['status']=='VALID_ABSREL_RESULT' and finite(x['global_bh_q_value'])]
 for x in valid:assert (float(x['global_bh_q_value'])<.05)==(x['significant_global_fdr_0.05']=='YES')
 for x in tests:pslong.append(dict(base,**x))
 if c:cafelong.append(dict(base,**c))
 cp=c.get('primary_significance') in ['SIGNIFICANT','NOT_SIGNIFICANT'];cs=c.get('sensitivity_significance') in ['SIGNIFICANT','NOT_SIGNIFICANT']
 robust=c.get('robust_significant_same_direction')=='YES' and c.get('primary_direction')=='EXPANSION'
 if robust:assert cp and cs and c['sensitivity_direction']=='EXPANSION'
 wd=w.get('duplication_class','NA');wa=wd in ['SINGLETON','WGD_SEGMENTAL','DISPERSED','PROXIMAL','TANDEM']
 rr=dict(r,harmonized_gene_id=hid,selected_transcript_id=m.get('selected_transcript_id','NA'),selected_protein_id=m.get('selected_protein_id','NA'),id_mapping='EXACT_ORIGINAL_GENE' if m else 'UNMAPPED',orthogroup=group,duplication_class=wd,wgd_assessed=int(wa),wgd_segmental=int(wd=='WGD_SEGMENTAL') if wa else 'NA',cafe_model_present=int(bool(c)),cafe_primary_branch_tested=int(cp),cafe_primary_expansion=int(c.get('primary_significance')=='SIGNIFICANT' and c.get('primary_direction')=='EXPANSION') if cp else 'NA',cafe_both_models_branch_tested=int(cp and cs),robust_expansion=int(robust) if cp and cs else 'NA',ps_valid_tests=len(valid),ps_assessed=int(bool(valid)),positive_selection=int(any(x['significant_global_fdr_0.05']=='YES' for x in valid)) if valid else 'NA')
 geneevents=[]
 for ct in [x for x in contrasts if x['species_code']==sp]:
  d=de.get((sp,hid,ct['contrast_id']),{});rm=rnamap.get(k,{});assert not d or rm.get('mapping_status')=='EXACT_GENE',(k,'RNA mapping missing')
  assess=bool(d) and finite(d.get('padj_within_contrast')) and finite(d.get('log2FoldChange_apeglm'));globalass=bool(d) and finite(d.get('padj_global_14_contrasts')) and finite(d.get('log2FoldChange_apeglm'))
  positive=assess and float(d['padj_within_contrast'])<.05 and abs(float(d['log2FoldChange_apeglm']))>=1
  e=dict(base,contrast_id=ct['contrast_id'],priority=ct['priority'],interpretation=ct['interpretation'],numerator_group=ct['numerator_group'],denominator_group=ct['denominator_group'],harmonized_gene_id=hid,gff_gene_id=rm.get('gff_gene_id','NA'),result_present=int(bool(d)),assessable=int(assess),primary_DE=int(positive) if assess else 'NA',global_sensitivity_DE=int(float(d['padj_global_14_contrasts'])<.05 and abs(float(d['log2FoldChange_apeglm']))>=1) if globalass else 'NA',direction=('UP' if float(d['log2FoldChange_apeglm'])>0 else 'DOWN') if positive else 'NA')
  for col in ['baseMean','pvalue','padj_within_contrast','log2FoldChange_apeglm','padj_global_14_contrasts','independent_filter_flag','cooks_outlier_or_undefined_flag']:e[col]=d.get(col,'NA')
  events.append(e);geneevents.append(e)
 rr.update(DE_contrasts=len(geneevents),DE_assessable_contrasts=sum(x['assessable'] for x in geneevents),DE_positive_contrasts=sum(x['primary_DE']==1 for x in geneevents),DE_any=int(any(x['primary_DE']==1 for x in geneevents)) if any(x['assessable'] for x in geneevents) else 'NA',DE_all_contrasts_assessable=int(all(x['assessable'] for x in geneevents)) if geneevents else 'NA')
 result.append(rr)
 if not m or group=='NA' or not wa:issues.append(dict(base,mapping=rr['id_mapping'],orthogroup=group,duplication_class=wd))
# Categories are observations with explicit assessment denominators; missing is not negative.
def summary(rs,sp,fam):
 blocked=fam in ['MEP','PPA'];n=len(rs);z={'wang_family':fam,'species':sp,'scope':'DEFINITION_BLOCKED' if blocked else 'SUPPORTED_ONLY_NOT_EXHAUSTIVE','supported_loci':'NA' if blocked else n}
 def metric(name,assess,positive,scope=True):
  eligible=[r for r in rs if assess(r)];z[name+'_assessed']='NA' if blocked or not scope else len(eligible);z[name+'_positive']='NA' if blocked or not scope or (n>0 and not eligible) else sum(positive(r) for r in eligible)
 metric('WGD_segmental',lambda r:r['wgd_assessed'],lambda r:r['wgd_segmental']==1)
 metric('CAFE_primary_expansion',lambda r:r['cafe_primary_branch_tested'],lambda r:r['cafe_primary_expansion']==1)
 metric('CAFE_robust_expansion',lambda r:r['cafe_both_models_branch_tested'],lambda r:r['robust_expansion']==1)
 z['CAFE_model_present']='NA' if blocked else sum(r['cafe_model_present'] for r in rs)
 metric('positive_selection',lambda r:r['ps_assessed'],lambda r:r['positive_selection']==1,sp not in ['BV','SLP_Hap1'])
 metric('DE_any_contrast',lambda r:r['DE_assessable_contrasts']>0,lambda r:r['DE_any']==1,sp in ['PO','SS'])
 z['DE_all_contrasts_assessed']='NA' if blocked or sp not in ['PO','SS'] else sum(r['DE_all_contrasts_assessable']==1 for r in rs)
 return z
summaries=[summary([r for r in result if r['species']==sp and r['wang_family']==fam],sp,fam) for fam in families for sp in species]
# Species totals deduplicate loci across families explicitly.
spresults=[]
for sp in species:
 dedup={r['gene_id']:r for r in result if r['species']==sp};spresults.append(summary(list(dedup.values()),sp,'ALL_UNIQUE_LOCI'))
desummary=[]
for fam in families:
 for ct in contrasts:
  es=[e for e in events if e['wang_family']==fam and e['contrast_id']==ct['contrast_id']];ass=sum(e['assessable'] for e in es);blocked=fam in ['MEP','PPA']
  desummary.append(dict(wang_family=fam,species=ct['species_code'],contrast_id=ct['contrast_id'],priority=ct['priority'],supported_loci='NA' if blocked else len(es),assessed='NA' if blocked else ass,DE_positive='NA' if blocked or (es and not ass) else sum(e['primary_DE']==1 for e in es),UP='NA' if blocked or (es and not ass) else sum(e['direction']=='UP' for e in es),DOWN='NA' if blocked or (es and not ass) else sum(e['direction']=='DOWN' for e in es)))
out.mkdir(parents=True)
write('01_supported_locus_evidence.tsv',result);write('02_family_species_summary.tsv',summaries);write('03_species_summary.tsv',spresults);write('04_DE_locus_contrast_evidence.tsv',events);write('05_DE_family_contrast_summary.tsv',desummary);write('06_positive_selection_tests.tsv',pslong);write('07_CAFE_model_evidence.tsv',cafelong);write('08_mapping_exceptions.tsv',issues,['wang_family','species','gene_id','mapping','orthogroup','duplication_class'])
long=[]
for r in result:
 for category,col in [('SUPPORTED','count_tier'),('WGD_SEGMENTAL','wgd_segmental'),('PRIMARY_EXPANDED_OG_MEMBER','cafe_primary_expansion'),('ROBUST_EXPANDED_OG_MEMBER','robust_expansion'),('POSITIVE_SELECTION','positive_selection'),('DE_ANY_CONTRAST','DE_any')]:
  if r[col] in [1,'SUPPORTED']:long.append({k:r[k] for k in ['wang_family','species','gene_id','protein_id','harmonized_gene_id']}|{'evidence_category':category})
write('09_family_species_gene_lists.tsv',long)
for metric in ['supported_loci','WGD_segmental_positive','WGD_segmental_assessed','CAFE_primary_expansion_positive','CAFE_primary_expansion_assessed','CAFE_robust_expansion_positive','positive_selection_positive','positive_selection_assessed','DE_any_contrast_positive','DE_any_contrast_assessed']:
 matrix=[]
 for sp in species:
  row={'species':sp}
  for fam in families:row[fam]=next(r[metric] for r in summaries if r['species']==sp and r['wang_family']==fam)
  matrix.append(row)
 write('matrix_'+metric+'.tsv',matrix)
qc=dict(status='COMPLETE_WITH_EXPLICIT_ASSESSMENT_LIMITS',supported_records=len(result),unique_species_loci=len(set((r['species'],r['gene_id']) for r in result)),mapped_loci=sum(r['id_mapping']=='EXACT_ORIGINAL_GENE' for r in result),mapping_exceptions=len(issues),DE_gene_contrast_records=len(events),DE_result_rows=sum(e['result_present'] for e in events),DE_assessable_rows=sum(e['assessable'] for e in events),new_analyses=False,source_ledger_sha256=sha(root/'01_supported_locus_ledger.tsv'))
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');shutil.copyfile(root/'02_full_source_manifest.tsv',out/'SOURCE_MANIFEST.tsv');shutil.copyfile(root/'SHA256SUMS.txt',out/'SOURCE_BUNDLE_SHA256SUMS.txt')
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+p.name+'\n' for p in sorted(out.iterdir()) if p.is_file()))
print(json.dumps(qc,indent=2));print(json.dumps(spresults,indent=2))
