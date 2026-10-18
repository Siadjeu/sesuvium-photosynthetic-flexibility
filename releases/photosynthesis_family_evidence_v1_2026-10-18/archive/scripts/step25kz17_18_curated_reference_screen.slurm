#!/usr/bin/env bash
#SBATCH --job-name=curated_rescue
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=curated_rescue_%j.out
#SBATCH --error=curated_rescue_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from the project directory}"
command -v python3
command -v diamond
python3 - <<'PY'
import csv,hashlib,json,os,re,subprocess,tarfile
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
refs=base/'step25kn/references'
proteome=root/os.environ.get('PROTEOME','08_focused_families/validation/step25kz16.7/release/work/all_species_primary_proteins.faa')
if not proteome.is_file():raise SystemExit('Primary proteome not found; set PROTEOME explicitly')
out=base/'step25kz17.18'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir();(out/'references').mkdir()
cpus=os.environ.get('SLURM_CPUS_PER_TASK','4')
species={'AT','BV','MC','PO','SLP_Hap1','SP_Hex','SP_Tet','SS'}
def fasta(p):
 h=None;s=[]
 with p.open() as f:
  for l in f:
   if l.startswith('>'):
    if h is not None:yield h,''.join(s)
    h=l[1:].strip();s=[]
   elif h is not None:s.append(l.strip())
 if h is not None:yield h,''.join(s)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rr)
source=[];records={};audit=[]
target_accessions={'O23404':'PPDK','O49562':'PPDK_RP','Q9MAC9':'PPDK_RP',
 'P31414':'V_PPASE','Q56ZN6':'V_PPASE','Q9FWR2':'V_PPASE',
 'Q93VK0':'PPCK','Q9SPK4':'PPCK'}
for name in ['PPDK','PPDK_RP','V_ATPASE','V_PPASE','PEPC']:
 p=refs/(name+'.faa')
 if not p.is_file():raise SystemExit('Required reference file missing: '+name)
 source.append(dict(source=str(p.relative_to(root)),sha256=sha(p)))
 for header,seq in fasta(p):
  bits=header.split()[0].split('|')
  if len(bits)<3:raise SystemExit('Unexpected reference identifier: '+header)
  acc=bits[1]
  if acc in records and records[acc]['sequence']!=seq:raise SystemExit('Conflicting sequences for '+acc)
  cat=target_accessions.get(acc,'BACKGROUND_REFERENCE')
  sub=''
  if 'V-type proton ATPase' in header:
   cat='V_ATPASE'
   m=re.search(r'GN=VHA-([^ ]+)',header)
   if not m:raise SystemExit('V-ATPase subunit annotation absent: '+acc)
   sub=re.sub(r'\d+$','',m.group(1)) # A and a remain distinct; c and c-double-prime remain distinct.
  if acc in target_accessions:
   expected={'PPDK':'Pyruvate, phosphate dikinase 1','PPDK_RP':'dikinase regulatory protein',
    'V_PPASE':'Pyrophosphate-energized','PPCK':'Phosphoenolpyruvate carboxylase kinase'}[cat]
   if expected not in header:raise SystemExit('Accession/header mismatch: '+acc)
  records[acc]=dict(accession=acc,category=cat,subunit=sub,header=header,sequence=seq)
  audit.append(dict(source_file=name+'.faa',accession=acc,assigned_category=cat,subunit=sub,header=header))
if not set(target_accessions).issubset(records):raise SystemExit('Missing expected curated accessions')
counts=Counter(r['category'] for r in records.values())
if counts['V_ATPASE']!=28:raise SystemExit('V-ATPase reference inventory changed: expected 28; review before running')
table('01_reference_curation.tsv',audit,['source_file','accession','assigned_category','subunit','header'])
table('02_reference_counts.tsv',[dict(category=k,reference_count=v) for k,v in sorted(counts.items())],['category','reference_count'])
table('03_reference_source_checksums.tsv',source,['source','sha256'])
# Search targets only. Background references retained for future explicit specificity tests.
targets={a:r for a,r in records.items() if r['category']!='BACKGROUND_REFERENCE'}
for cat in sorted(counts):
 with (out/'references'/(cat+'.faa')).open('w') as f:
  for acc,r in sorted(records.items()):
   if r['category']==cat:f.write('>'+acc+' '+r['header']+'\n'+r['sequence']+'\n')
query=out/'references/curated_targets.faa'
with query.open('w') as f:
 for acc,r in sorted(targets.items()):f.write('>'+acc+'\n'+r['sequence']+'\n')
table('04_unresolved_categories.tsv',[dict(category=x,status='REFERENCE_REQUIRED_NOT_SCREENED') for x in ['TDT','DTC','DIC']],['category','status'])
# Validate unique primary IDs before searching, avoiding accidental isoform or reference inclusion.
lengths={};spcounts=Counter()
for header,seq in fasta(proteome):
 pid=header.split()[0];sp=pid.split('__',1)[0]
 if sp not in species or '__' not in pid:raise SystemExit('Unexpected primary proteome ID: '+pid)
 if pid in lengths:raise SystemExit('Duplicate primary protein ID: '+pid)
 lengths[pid]=len(seq);spcounts[sp]+=1
if set(spcounts)!=species:raise SystemExit('Primary proteome must contain all eight species')
table('05_primary_proteome_inventory.tsv',[dict(species=s,records=spcounts[s]) for s in sorted(species)],['species','records'])
(out/'PROTEOME_SHA256.txt').write_text(sha(proteome)+'  primary_proteome_input\n')
commands=[]
def run(args,log):
 commands.append(args)
 with (out/log).open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
db=out/'work/primary'
run(['diamond','version'],'diamond_version.txt')
run(['diamond','makedb','--in',str(proteome),'--db',str(db),'--threads',cpus],'06_makedb.log')
hit=out/'07_curated_reference_hits.tsv'
run(['diamond','blastp','--query',str(query),'--db',str(db),'--out',str(hit),
 '--threads',cpus,'--sensitive','--evalue','1e-5','--max-target-seqs','0',
 '--outfmt','6','qseqid','sseqid','pident','length','qstart','qend','sstart','send',
 'evalue','bitscore','qlen','slen'],'07_diamond.log')
fields=['reference','protein_id','pident','alignment_length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen']
best={};queryhits=Counter()
with hit.open() as f:
 for r in csv.DictReader(f,fieldnames=fields,delimiter='\t'):
  if r['reference'] not in targets or r['protein_id'] not in lengths:raise SystemExit('Unrecognized search ID')
  ref=targets[r['reference']];pid=r['protein_id'];key=(ref['category'],ref['subunit'],pid)
  queryhits[r['reference']]+=1
  entry=dict(category=ref['category'],subunit=ref['subunit'],species=pid.split('__',1)[0],
   protein_id=pid,protein_length=lengths[pid],best_reference=r['reference'],evalue=r['evalue'],
   bitscore=r['bitscore'],identity_percent=r['pident'],
   reference_span_coverage=round((abs(int(r['qend'])-int(r['qstart']))+1)/int(r['qlen']),4),
   candidate_span_coverage=round((abs(int(r['send'])-int(r['sstart']))+1)/int(r['slen']),4),
   status='SIMILARITY_CANDIDATE_REVIEW_ONLY')
  if key not in best or float(entry['bitscore'])>float(best[key]['bitscore']):best[key]=entry
rows=[best[k] for k in sorted(best)]
table('08_family_subunit_candidates.tsv',rows,['category','subunit','species','protein_id','protein_length','best_reference','evalue','bitscore','identity_percent','reference_span_coverage','candidate_span_coverage','status'])
byid=defaultdict(set)
for r in rows:byid[r['protein_id']].add(r['category']+(':'+r['subunit'] if r['subunit'] else ''))
table('09_ambiguous_assignments.tsv',[dict(protein_id=p,hit_categories=';'.join(sorted(c)),status='REVIEW_MULTIPLE_REFERENCE_GROUPS') for p,c in sorted(byid.items()) if len(c)>1],['protein_id','hit_categories','status'])
table('10_reference_search_recovery.tsv',[dict(reference=a,category=r['category'],subunit=r['subunit'],hit_rows=queryhits[a]) for a,r in sorted(targets.items())],['reference','category','subunit','hit_rows'])
with (out/'11_candidate_sequences.faa').open('w') as f:
 for h,s in fasta(proteome):
  if h.split()[0] in byid:f.write('>'+h.split()[0]+'\n'+s+'\n')
table('12_screen_candidate_counts_NOT_FINAL.tsv',[dict(category=c,species=s,candidate_ids=len({r['protein_id'] for r in rows if r['category']==c and r['species']==s})) for c in sorted({r['category'] for r in targets.values()}) for s in sorted(species)],['category','species','candidate_ids'])
# Store reproducible commands without machine-specific absolute paths.
(out/'COMMANDS.json').write_text(json.dumps([[a.replace(str(root),'<PROJECT_ROOT>') for a in cmd] for cmd in commands],indent=2)+'\n')
(out/'README.txt').write_text("""Curated reference screening, not final gene assignment.
PPDK O23404; PPDK regulators O49562/Q9MAC9; PPCK Q93VK0/Q9SPK4;
membrane PPases P31414/Q56ZN6/Q9FWR2.
28 V-type ATPase references selected from explicit source-header annotations and
separated by VHA subunit; A/a and C/c case distinctions are preserved.
Curation uses the supplied sequences and annotations, not a fresh external database.
Other references are retained as BACKGROUND_REFERENCE; not used to define targets.
DIAMOND searches all eight primary proteomes with E <= 1e-5 and no length filter.
This is a permissive recovery screen: no hits or zero candidates is not proven absence.
08 reports best single alignment per candidate/category/subunit. Coverage is alignment
coordinate span, not domain coverage or proof of full-length function.
09 flags multi-group hits; it does not resolve them by highest score.
Do not sum overlapping subunit/category rows as gene counts.
Next validate domain architecture and specificity, including relevant background
families, then reconcile candidates with old retained/rejected lists and gene loci.
TDT/DTC/DIC await independent curated references. No combined-family HMM is built.
Existing BETA-CA and PEPC releases are not modified.
The work database stays on the cluster and is omitted from this review bundle.
""")
qc=dict(status='SCREEN_COMPLETE_REVIEW_ONLY',reference_counts=dict(counts),
 target_references=len(targets),unique_candidate_ids=len(byid),candidate_group_rows=len(rows),
 multi_group_candidates=sum(len(c)>1 for c in byid.values()),final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Curated_reference_screen/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
