#!/usr/bin/env bash
#SBATCH --job-name=beta_ca_itol
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=beta_ca_itol_%j.out
#SBATCH --error=beta_ca_itol_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from the project directory}"
python3 - <<'PY'
import csv,hashlib,json,os,re,shutil,zipfile
from pathlib import Path
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
names=['main_full','main_core','domains_full','domains_core']
if os.environ.get('TREE_RELEASE'):choices=[Path(os.environ['TREE_RELEASE']).resolve(strict=True)]
else:choices=[p for p in (base/'step25kz17.13').glob('release*') if all((p/n/'COMPLETE.txt').is_file() for n in names)]
if len(choices)!=1:raise SystemExit('Set TREE_RELEASE to the directory containing the four completed tree scenarios; found '+str(len(choices)))
source=choices[0]
out=base/'step25kz17.15'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
out.mkdir(parents=True,exist_ok=False)
colors={'AT':'#7b3294','BV':'#52616b','MC':'#008080','PO':'#d2691e','SLP_Hap1':'#0072b2','SP_Hex':'#2e7d32','SP_Tet':'#8a7800','SS':'#c62828','REFERENCE':'#111111'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,rr,fields):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def strip(p,title,palette,data):
 lines=['DATASET_COLORSTRIP','SEPARATOR TAB','DATASET_LABEL\t'+title,'COLOR\t#555555','COLOR_BRANCHES\t0','STRIP_WIDTH\t20','MARGIN\t3','SHOW_STRIP_LABELS\t0','LEGEND_TITLE\t'+title,'LEGEND_SHAPES\t'+'\t'.join('1' for _ in palette),'LEGEND_COLORS\t'+'\t'.join(palette.values()),'LEGEND_LABELS\t'+'\t'.join(palette),'DATA']
 p.write_text('\n'.join(lines+['\t'.join(r) for r in data])+'\n')
summary=[]
for name in names:
 src=source/name;d=out/name;d.mkdir()
 for l in (src/'SHA256SUMS.txt').read_text().splitlines():
  h,f=l.split(maxsplit=1)
  if sha(src/f)!=h:raise SystemExit('Source checksum mismatch: '+name+'/'+f)
 mapping=read(src/'TIP_LABEL_MAP.tsv');byid={r['safe_tip_id']:r for r in mapping}
 tree=(src/'tree.treefile').read_text().strip();leaves=re.findall(r'(?<=[(,])\s*(T\d+)(?=[:),;\s])',tree)
 if len(leaves)!=len(byid) or set(leaves)!=set(byid):raise SystemExit('Tree mapping mismatch')
 # Numeric Newick support = UFBoot; both measures also retained as named NHX metadata.
 # Stable internal IDs link the support table to iTOL nodes.
 support_rows=[]
 def convert(m):
  a,b,length=m.groups();node='N'+str(len(support_rows)+1).zfill(3)
  if not (0<=float(a)<=100 and 0<=float(b)<=100):raise ValueError('Invalid support value')
  support_rows.append(dict(node_id=node,SH_aLRT=a,UFBoot=b,branch_length=length))
  return ')'+b+':'+length+'[&&NHX:name='+node+':SH_aLRT='+a+':UFBoot='+b+']'
 converted=re.sub(r'\)([0-9.]+)/([0-9.]+):([0-9.eE+-]+)',convert,tree)
 if len(support_rows)!=len(byid)-3:raise SystemExit('Unexpected number of internal support labels')
 (d/'01_tree_for_iTOL.nwk').write_text(converted+'\n')
 (d/'original_IQTREE_dual_support.nwk').write_text(tree+'\n')
 labels=['LABELS','SEPARATOR TAB','DATA'];styles=['TREE_COLORS','SEPARATOR TAB','DATA'];species=[];roles=[]
 role_colors={'Standard candidate':'#dddddd','AT alternative isoform':'#cc79a7','AT domain fragment - one locus':'#cc79a7','SLP reviewed candidate':'#56b4e9','Reference':'#444444'}
 for r in mapping:
  sid=r['safe_tip_id'];tip=r['tip_id'];sp='REFERENCE' if r['role']=='REFERENCE' else tip.split('__')[0]
  if sp not in colors:raise SystemExit('Unknown species code: '+sp)
  labels.append(sid+'\t'+tip)
  special='LOC130820426' in tip or tip=='SLP_Hap1__SP3G00376'
  styles.append('\t'.join([sid,'label',colors[sp],'bold' if special else 'normal','1']))
  if special:styles.append('\t'.join([sid,'label_background','#fff2cc']))
  species.append((sid,colors[sp],sp))
  role='Reference' if sp=='REFERENCE' else 'AT domain fragment - one locus' if '__domain_' in tip else 'AT alternative isoform' if '__isoform_' in tip else 'SLP reviewed candidate' if tip=='SLP_Hap1__SP3G00376' else 'Standard candidate'
  roles.append((sid,role_colors[role],role))
 (d/'02_tip_labels.txt').write_text('\n'.join(labels)+'\n')
 (d/'03_label_colors.txt').write_text('\n'.join(styles)+'\n')
 strip(d/'04_species_strip.txt','Species code',colors,species)
 present={r[2] for r in roles};palette={k:v for k,v in role_colors.items() if k in present}
 strip(d/'05_sequence_role_strip.txt','Sequence representation',palette,roles)
 table(d/'06_branch_support.tsv',support_rows,['node_id','SH_aLRT','UFBoot','branch_length'])
 shutil.copyfile(src/'TIP_LABEL_MAP.tsv',d/'07_tip_identity_map.tsv')
 # Ensure every annotation ID points to an existing leaf.
 assert len(species)==len(byid) and {x[0] for x in species}==set(byid)
 # Only support metadata changed; preserve exact leaf IDs and every branch-length token.
 assert re.findall(r':([0-9.eE+-]+)(?=[,);\[])',converted)==re.findall(r':([0-9.eE+-]+)(?=[,);\[])',tree)
 summary.append(dict(scenario=name,tips=len(byid),internal_support_records=len(support_rows),source_tree_sha256=sha(src/'tree.treefile'),export_tree_sha256=sha(d/'01_tree_for_iTOL.nwk')))
(out/'UPLOAD_INSTRUCTIONS.txt').write_text('''BETA-CA iTOL upload package

1. Extract this ZIP on your computer. Start with the main_full folder.
2. At https://itol.embl.de upload 01_tree_for_iTOL.nwk as a new tree.
3. Once displayed, drag these files from the SAME scenario folder onto the tree:
   02_tip_labels.txt
   03_label_colors.txt
   04_species_strip.txt
   05_sequence_role_strip.txt
   Alternatively use Upload annotation files in the Datasets tab.
4. Choose rectangular or circular layout. In rectangular view, align tip labels
   and increase vertical spacing to prevent overlap. Keep branch lengths enabled
   when presenting a phylogram and show the scale bar.
5. Under Advanced, use Bootstraps / metadata. Numeric bootstrap labels are
   UFBoot, not SH-aLRT. Named NHX fields UFBoot and SH_aLRT are also included.
   Select the desired data source explicitly. Start with symbols rather than
   every support number to reduce clutter; state any display threshold in the legend.
   All values remain available in 06_branch_support.tsv, keyed to NHX node IDs.
6. Repeat with main_core as a separate sensitivity tree. The two domains_* folders
   provide supplementary checks. Never mix annotation files between folders:
   safe T identifiers are local to each scenario.
7. Export SVG or PDF from iTOL for review, subject to your account's export options.

Interpretation
The primary phylogenetic set has 45 candidate loci plus six references. The
278-aa AT isoform is used only for this analysis. The domain trees replace that
isoform with two fragments of one AT locus, not two genes. Identical sequences
retain separate locus IDs. Species codes follow the project; REFERENCE denotes
Arabidopsis reference sequences, whereas AT denotes Amaranthus tricolor.
Trees remain biologically unrooted. The stored display junction is not an
inferred ancestor. Do not choose alpha-CA as an outgroup for these exports.
Orange/blue role strips and highlighted labels describe sequence representation,
not confidence or biological acceptance. Final family counts remain unresolved.

Provenance
The export preserves topology and branch lengths exactly. UFBoot is copied into
the numeric internal-label field; both original measures are recorded in NHX
metadata and a TSV. original_IQTREE_dual_support.nwk preserves the source tree.
No inference, rerooting, sequence filtering or gene-count update is performed.
Official formats: https://itol.embl.de/help.cgi
https://itol.embl.de/help/labels_template.txt
https://itol.embl.de/help/colors_styles_template.txt
https://itol.embl.de/help/dataset_color_strip_template.txt
''')
table(out/'EXPORT_QC.tsv',summary,list(summary[0]))
(out/'COMPLETE.txt').write_text('ITOL_EXPORT_COMPLETE_LOCAL_VALIDATION_ONLY\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+str(p.relative_to(out))+'\n')
bundle=out/'BETA_CA_iTOL_upload_package.zip'
with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p!=bundle:z.write(p,'BETA_CA_iTOL/'+str(p.relative_to(out)))
print('Package relative to PROJECT_ROOT:',bundle.relative_to(root))
PY
