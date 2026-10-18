#!/usr/bin/env bash
#SBATCH --job-name=vatp_excepts
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=vatp_excepts_%j.out
#SBATCH --error=vatp_excepts_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
for tool in python3 mafft; do command -v "$tool"; done
python3 - <<'PY'
import gzip,csv,hashlib,json,os,re,subprocess,tarfile,shutil
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def release(step,env):
 paths=[root/os.environ[env]] if os.environ.get(env) else [p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
 if len(paths)!=1:raise SystemExit('Expected one release; set '+env+' to the intended release directory')
 p=paths[0].resolve(strict=True)
 for line in (p/'SHA256SUMS.txt').read_text().splitlines():
  h,n=line.split(maxsplit=1);f=(p/n).resolve()
  if not f.is_relative_to(p) or sha(f)!=h:raise SystemExit('Checksum mismatch: '+n)
 return p
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta(p):
 d={};key=None
 for line in p.read_text().splitlines():
  if line.startswith('>'):
   key=line[1:].split()[0]
   if key in d:raise SystemExit('Duplicate FASTA ID: '+key)
   d[key]=''
  elif line.strip():
   if key is None:raise SystemExit('Invalid FASTA')
   d[key]+=line.strip()
 return d
def writefa(p,d):p.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
src=release('step25kz17.33','VATPASE_RELEASE')
rows=read(src/'12_VATPASE_candidate_ledger_REVIEW_ONLY.tsv')
meta={r['accession']:r for r in read(src/'05_reference_metadata.tsv')}
assert len(rows)==603 and len({(r['species'],r['gene_id']) for r in rows})==603
for r in rows:
 f=r['flags'];other=r['best_other_reference']
 if 'SUBUNIT_OR_COMPETITOR_AMBIGUITY' in f:
  assert other in ['P92549','P83483'] and float(r['score_margin'])<=-166
  cat='F_TYPE_FAVOURED_NOT_COUNTED_AS_VATPASE'
 elif not f:cat='CONCORDANT_EVIDENCE_REVIEW_REQUIRED'
 elif f=='NO_COMPETITOR_REPORTED_IN_LIMITED_PANEL':cat='CONCORDANT_WITH_LIMITED_PANEL_CAVEAT'
 else:cat='STRUCTURE_COMPLETENESS_EXCEPTION'
 r['evidence_category']=cat;r['final_decision']='PENDING'
counts=Counter(r['evidence_category'] for r in rows)
assert sorted(counts.values())==[23,85,208,287]
exceptions=[r for r in rows if r['evidence_category']=='STRUCTURE_COMPLETENESS_EXCEPTION']
seqs=fasta(src/'06_candidates_and_controls.faa');refs=fasta(src/'06_reference_panel.faa')
out=base/'step25kz17.34'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False)
def table(name,rr,fields):
 p=out/name;p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
table('01_reconciled_evidence_ledger.tsv',rows,list(rows[0]))
table('02_exception_queue.tsv',exceptions,list(rows[0]))
for name in ['01_PPCK_closed_ledger.tsv','02_PPCK_supported_counts.tsv','03_PPCK_closure.txt']:
 shutil.copy2(src/name,out/('PPCK_'+name))
summary=[]
for sp in sorted({r['species'] for r in rows}):
 for sub in sorted({r['provisional_subunit'] for r in rows}):
  for cat in sorted(counts):
   n=sum(r['species']==sp and r['provisional_subunit']==sub and r['evidence_category']==cat for r in rows)
   summary.append(dict(species=sp,subunit=sub,evidence_category=cat,loci=n,final_count='NO'))
table('03_evidence_category_counts_NOT_FINAL.tsv',summary,list(summary[0]))
commands=[];cpu=os.environ.get('SLURM_CPUS_PER_TASK','4');metrics=[]
(out/'alignments').mkdir()
for sub in sorted({r['provisional_subunit'] for r in exceptions}):
 label='V1_'+sub if sub.isupper() else ('V0_c_doubleprime' if sub=="c''" else 'V0_'+sub)
 rr=[r for r in rows if r['provisional_subunit']==sub and r['evidence_category']!='F_TYPE_FAVOURED_NOT_COUNTED_AS_VATPASE']
 aa=[a for a,m in meta.items() if m['assigned_category']=='V_ATPASE' and m['subunit']==sub]
 group={r['protein_id']:seqs[r['protein_id']] for r in rr};group.update({'REF_'+a:refs[a] for a in aa})
 inp=out/'alignments'/(label+'_input.faa');alnp=out/'alignments'/(label+'_alignment.faa');writefa(inp,group)
 cmd=['mafft','--auto','--thread',cpu,str(inp)];commands.append(cmd)
 with alnp.open('w') as f,(out/'alignments'/(label+'_mafft.log')).open('w') as log:subprocess.run(cmd,stdout=f,stderr=log,check=True)
 aln=fasta(alnp)
 assert set(aln)==set(group) and len({len(s) for s in aln.values()})==1
 assert all(aln[k].replace('-','').upper()==s.upper() for k,s in group.items())
 for r in rr:
  q=r['protein_id'];cs=aln[q]
  for a in aa:
   rs=aln['REF_'+a];pos=[i for i,c in enumerate(rs) if c!='-'];occ=[j for j,i in enumerate(pos) if cs[i]!='-'];pair=[i for i in pos if cs[i]!='-']
   gap=longest=0
   for i in pos:
    gap=gap+1 if cs[i]=='-' else 0;longest=max(longest,gap)
   metrics.append(dict(protein_id=q,subunit=sub,evidence_category=r['evidence_category'],reference=a,reference_length=len(pos),protein_length=len(group[q]),reference_occupancy=round(len(occ)/len(pos),4),paired_identity=round(sum(cs[i].upper()==rs[i].upper() for i in pair)/len(pair),4) if pair else 0,missing_N_reference=occ[0] if occ else len(pos),missing_C_reference=len(pos)-1-occ[-1] if occ else len(pos),longest_reference_gap=longest,N_extension=sum(c!='-' for c in cs[:pos[0]]),C_extension=sum(c!='-' for c in cs[pos[-1]+1:]),insertion_residues=sum(cs[i]!='-' and rs[i]=='-' for i in range(pos[0],pos[-1]+1))))
table('04_alignment_completeness_all_context.tsv',metrics,list(metrics[0]))
exids={r['protein_id'] for r in exceptions}
table('05_exception_alignment_completeness.tsv',[r for r in metrics if r['protein_id'] in exids],list(metrics[0]))
dom=read(src/'08_domain_evidence.tsv')
table('06_exception_and_reference_domains.tsv',[r for r in dom if r['protein_id'] in exids or r['protein_id'].startswith('REF_')],list(dom[0]))
# Collect existing annotation evidence, without model repair or genomic alignment.
def op(p):return gzip.open(p,'rt') if p.suffix=='.gz' else p.open()
def records(p):
 with op(p) as f:
  h=None;s=[]
  for l in f:
   if l.startswith('>'):
    if h is not None:yield h,''.join(s)
    h=l[1:].split()[0];s=[]
   elif l.strip():s.append(l.strip())
  if h is not None:yield h,''.join(s)
manifest=root/'00_admin/manifests/input_sets.tsv';inputs=read(manifest)
provenance=[];status=[];loci=[]
oldhash=json.loads((src/'PROVENANCE.json').read_text())['mapping_hashes']
for sp in sorted({r['species'] for r in exceptions}):
 folder=Path('annotation')/sp;(out/folder).mkdir(parents=True)
 mp=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv')
 assert sha(mp)==oldhash[sp],'Mapping changed since step17.33'
 mappings=read(mp);genes={r['gene_id'] for r in exceptions if r['species']==sp}
 selected=[m for m in mappings if m['gene_id'] in genes]
 assert {m['gene_id'] for m in selected}==genes
 table(str(folder/'01_all_isoform_mapping.tsv'),selected,list(selected[0]))
 aliases=set();windows=[]
 for g in sorted(genes):
  mm=[m for m in selected if m['gene_id']==g];locations={(m['seqid'],m['strand']) for m in mm}
  assert len(locations)==1,'Inconsistent locus coordinates'
  start=min(int(m['start']) for m in mm);end=max(int(m['end']) for m in mm)
  loci.append(dict(species=sp,gene_id=g,seqid=mm[0]['seqid'],strand=mm[0]['strand'],start=start,end=end,isoform_mapping_rows=len(mm)))
  windows.append((mm[0]['seqid'],max(1,start-5000),end+5000))
  for m in mm:
   for k in ['gene_id','transcript_id','protein_id','cds_id']:
    if m.get(k) and m[k] not in ['.','NA']:aliases.add(m[k])
 aliases|={sp+'__'+x for x in list(aliases)}
 configs=[m for m in inputs if m['input_set']==sp];assert len(configs)==1
 conf=configs[0]
 for kind,key in [('gff','gff3'),('protein','protein_all'),('cds','cds_all')]:
  path=root/conf[key]
  provenance.append(dict(species=sp,kind=kind,source=conf[key],sha256=sha(path) if path.is_file() else 'MISSING'))
  if not path.is_file():status.append(dict(species=sp,item=kind,status='MISSING',records=0));continue
  n=0
  if kind=='gff':
   with op(path) as f,(out/folder/'02_annotation_windows.gff3').open('w') as w:
    w.write('##gff-version 3\n')
    for l in f:
     if l.startswith('##FASTA'):break
     if l.startswith('#'):continue
     c=l.rstrip('\n').split('\t')
     if len(c)==9 and any(c[0]==chrom and int(c[4])>=a and int(c[3])<=b for chrom,a,b in windows):w.write(l.rstrip('\n')+'\n');n+=1
  else:
   found={}
   for h,s in records(path):
    if h in aliases:
     assert h not in found,'Duplicate source FASTA ID'
     found[h]=s
   writefa(out/folder/('03_'+kind+'_isoforms.fasta'),found);n=len(found)
  status.append(dict(species=sp,item=kind,status='COLLECTED' if n else 'NO_EXACT_MATCH',records=n))
 provenance.append(dict(species=sp,kind='mapping',source=str(mp.relative_to(root)),sha256=sha(mp)))
table('07_exception_locus_summary.tsv',loci,list(loci[0]))
table('08_annotation_collection_status.tsv',status,['species','item','status','records'])
table('09_annotation_provenance.tsv',provenance,['species','kind','source','sha256'])
review=[dict(protein_id=r['protein_id'],species=r['species'],gene_id=r['gene_id'],subunit=r['provisional_subunit'],original_flags=r['flags'],alignment_review='PENDING',annotation_review='PENDING',family_decision='PENDING',model_limitation='',reason='') for r in exceptions]
table('10_manual_review_template.tsv',review,list(review[0]))
(out/'README.txt').write_text('Bounded V-ATPase review: 85 concordant records, 208 concordant with limited-panel caveat, 287 strongly F-type-favoured cross-hits not counted as V-ATPase, and 23 structural/completeness exceptions. Categories are not final accepted gene counts.\nAlignments contain only non-F-type-favoured candidates of the same provisional subunit plus its V-type references; unchanged distinct loci are retained. No joint alignment of different subunits, no tree, no new search or model repair.\nMAFFT occupancy and terminal gaps are descriptive, not proof of biological truncation. Inspect all isoforms and domain evidence, including overlapping Pfam models, before interpretation.\nGFF windows retain original coordinates and include neighboring genes within 5 kb; neighboring features are context, not assigned to the target. Protein/CDS extraction uses exact mapped aliases. Missing files or unmatched IDs are explicit collection limitations. No read support, genome translation validation, or experimental function is inferred.\nPPCK closure is copied unchanged. Prior releases are not altered.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),input_manifest_sha256=sha(manifest)),indent=2)+'\n')
(out/'QC.json').write_text(json.dumps(dict(status='BOUNDED_ALIGNMENT_ANNOTATION_REVIEW_COMPLETE',evidence_categories=dict(counts),exception_loci=len(exceptions),aligned_subunits=sorted({r['provisional_subunit'] for r in exceptions}),annotation_collection_issues=[r for r in status if r['status']!='COLLECTED'],final_counts=False),indent=2)+'\n')
with (out/'mafft_version.txt').open('w') as f:subprocess.run(['mafft','--version'],stdout=f,stderr=subprocess.STDOUT,check=True)
(out/'COMMANDS.json').write_text(json.dumps([[v.replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for v in c] for c in commands],indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file())
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='VATPASE_bounded_review')
print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
