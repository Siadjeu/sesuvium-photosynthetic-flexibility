#!/usr/bin/env python3
"""PEPC representative phylogenies and iTOL export; no third-party Python modules."""
import argparse,collections,csv,hashlib,json,math,os,re,shutil,subprocess,tarfile
from pathlib import Path
from datetime import datetime,timezone
AA=set('ACDEFGHIKLMNPQRSTVWY'); POS=(519,780,890)
SPECIES={'AT':'Amaranthus tricolor','BV':'Beta vulgaris','MC':'Mesembryanthemum crystallinum','PO':'Portulaca oleracea','SLP_Hap1':'Selenicereus megalanthus','SP_Hex':'Sesuvium portulacastrum (hexaploid)','SP_Tet':'Sesuvium portulacastrum (tetraploid)','SS':'Sesuvium sesuvioides'}
PALETTE={'AT':'#D55E00','BV':'#0072B2','MC':'#009E73','PO':'#CC79A7','SLP_Hap1':'#E69F00','SP_Hex':'#56B4E9','SP_Tet':'#8C6BB1','SS':'#CC3311','REFERENCE':'#777777'}
SOURCE_SS='08_focused_families/downstream_supported/pepc_reference_review/ss_annotation_20835561_20261022T205826998329Z'
SOURCE_REC='08_focused_families/downstream_supported/pepc_reference_review/reconcile_20834956_20261022T203501981696Z'
META='reference_metadata/00_step25kz12_pepc_36_recovered_metadata.tsv'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def table(p,rows,fields=None):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def readtsv(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta(p):
 d={};k=None
 for line in p.read_text().splitlines():
  if line.startswith('>'):
   k=line[1:].split()[0];need(k not in d,'Duplicate FASTA ID: '+k);d[k]=''
  elif line.strip():
   need(k is not None,'FASTA sequence before header');d[k]+=line.strip().upper()
 return d
def writefa(p,d):p.write_text(''.join('>'+k+'\n'+s+'\n' for k,s in d.items()))
def ungap(s):return s.replace('-','').replace('.','')
def checked(folder,names):
 m=folder/'SHA256SUMS.txt';expected={}
 for line in m.read_text().splitlines():
  h,n=line.split(maxsplit=1);expected[n.lstrip('*').removeprefix('./')]=h
 result={str(m):sha(m)}
 for name in names:
  p=folder/name;need(p.is_file(),'Missing input: '+name)
  h=sha(p);need(expected.get(name)==h,'Checksum mismatch: '+name);result[str(p)]=h
 return result
# Parser for the safe, unquoted IDs emitted by this job and IQ-TREE support labels.
def parse_newick(text):
 tokens=re.findall(r'[^\s(),:;]+|[(),:;]',text.strip());i=0
 def node():
  nonlocal i
  children=[];name='';length=''
  need(i<len(tokens),'Truncated Newick')
  if tokens[i]=='(':
   i+=1;children.append(node())
   while tokens[i]==',':i+=1;children.append(node())
   need(tokens[i]==')','Invalid Newick closing bracket');i+=1
   if i<len(tokens) and tokens[i] not in [':',',',')',';']:name=tokens[i];i+=1
  else:
   name=tokens[i];need(name not in [':',',',')',';'],'Invalid Newick tip');i+=1
  if i<len(tokens) and tokens[i]==':':
   i+=1;length=tokens[i];float(length);i+=1
  return {'name':name,'length':length,'children':children}
 tree=node();need(i==len(tokens)-1 and tokens[i]==';','Unexpected Newick tail');return tree

def leaves(t):
 return [t['name']] if not t['children'] else [x for child in t['children'] for x in leaves(child)]
def splits(t):
 universe=set(leaves(t));result=set()
 def visit(n):
  s=set(leaves(n));other=universe-s
  if len(s)>1 and len(other)>1:
   a,b=tuple(sorted(s)),tuple(sorted(other));result.add(min(a,b))
  for child in n['children']:visit(child)
 visit(t);return result

def dataset_text(path,label,ordered,texts,colors=None,bold=False):
 lines=['DATASET_TEXT','SEPARATOR TAB','DATASET_LABEL\t'+label,'COLOR\t#222222','MARGIN\t12','ALIGN_TO_TREE\t1','SIZE_FACTOR\t1','SHOW_LABELS\t1','LABEL_ROTATION\t0','DATA']
 for k in ordered:
  if k in texts:lines.append('\t'.join([k,texts[k],'-1',(colors or {}).get(k,'#222222'),'bold' if bold else 'normal','1','0']))
 path.write_text('\n'.join(lines)+'\n')

def export(out,prep):
 metadata={r['tip_id']:r for r in readtsv(out/'tip_metadata.tsv')};wins={(r['tip_id'],int(r['maize_position'])):r for r in readtsv(out/'maize_windows.tsv')}
 trees={};orders={}
 for mode in ['main','sensitivity']:
  p=out/'trees'/(mode+'.treefile');t=parse_newick(p.read_text());order=leaves(t)
  need(len(order)==len(metadata) and set(order)==set(metadata),'Tree tip mismatch: '+mode);trees[mode]=t;orders[mode]=order
  folder=out/'itol'/mode;folder.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,folder/'01_PEPC.treefile')
  (folder/'02_labels.txt').write_text('LABELS\nSEPARATOR TAB\nDATA\n'+''.join(k+'\t'+metadata[k]['display_label']+'\n' for k in order))
  groups=list(PALETTE);lines=['DATASET_COLORSTRIP','SEPARATOR TAB','DATASET_LABEL\tSpecies / reference','COLOR\t#777777','COLOR_BRANCHES\t0','STRIP_WIDTH\t18','MARGIN\t5','BORDER_WIDTH\t0.5','BORDER_COLOR\t#FFFFFF','COMPLETE_BORDER\t1','LEGEND_TITLE\tSpecies codes','LEGEND_SHAPES\t'+'\t'.join('1' for _ in groups),'LEGEND_COLORS\t'+'\t'.join(PALETTE[x] for x in groups),'LEGEND_LABELS\t'+'\t'.join(groups),'DATA']
  lines += [k+'\t'+PALETTE[metadata[k]['species_code']]+'\t'+metadata[k]['species'] for k in order]
  (folder/'03_species_strip.txt').write_text('\n'.join(lines)+'\n')
  ordered_rows=[]
  for p in POS:
   texts={k:wins[(k,p)]['window'] for k in order};colors={k:'#222222' if wins[(k,p)]['full_window_covered']=='1' else '#999999' for k in order}
   dataset_text(folder/f'04_window_{p}.txt',f'{p} (maize {p-2}-{p+2})',order,texts,colors)
  dataset_text(folder/'05_historical_reference_subtypes.txt','Historical reference subtype (not inferred)',order,{k:metadata[k]['historical_subtype'] for k in order if metadata[k]['historical_subtype'] and metadata[k]['historical_subtype']!='UNASSIGNED'})
  for i,k in enumerate(order,1):
   ordered_rows.append(dict(tree_order=i,tip_id=k,label=metadata[k]['display_label'],represented_loci=metadata[k]['represented_loci'],window_519=wins[(k,519)]['window'],window_780=wins[(k,780)]['window'],window_890=wins[(k,890)]['window']))
  table(folder/'06_windows_in_tree_order.tsv',ordered_rows)
  (folder/'07_windows_in_tree_order.txt').write_text('tip_id\t519\t780\t890\n'+''.join(k+'\t'+'\t'.join(wins[(k,p)]['window'] for p in POS)+'\n' for k in order))
 a,b=splits(trees['main']),splits(trees['sensitivity']);common=a&b
 comparison=dict(main_nontrivial_splits=len(a),sensitivity_nontrivial_splits=len(b),shared_nontrivial_splits=len(common),RF_symmetric_difference=len(a^b),normalized_RF=len(a^b)/(len(a)+len(b)) if a or b else 0)
 (out/'TREE_SENSITIVITY.json').write_text(json.dumps(comparison,indent=2)+'\n')
 table(out/'TREE_SPLIT_COMPARISON.tsv',[dict(split_side=';'.join(s),main=int(s in a),sensitivity=int(s in b)) for s in sorted(a|b)],['split_side','main','sensitivity'])
 qc=dict(status='TREE_AND_ITOL_EXPORT_COMPLETE_REVIEW_TOPOLOGY',**prep,all_tree_tip_sets_verified=True,rooting='UNROOTED_NO_BIOLOGICAL_OUTGROUP_ASSIGNED',subtype_assignments_to_supported_genes='NOT_PERFORMED',reference_species_conflict_retained=True,existing_locus_counts_modified=False,tree_sensitivity=comparison)
 (out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');return qc

