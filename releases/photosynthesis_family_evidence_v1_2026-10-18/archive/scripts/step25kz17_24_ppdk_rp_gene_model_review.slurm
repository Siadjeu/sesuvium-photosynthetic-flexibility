#!/usr/bin/env bash
#SBATCH --job-name=ppdk_rp_models
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=ppdk_rp_models_%j.out
#SBATCH --error=ppdk_rp_models_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from your project directory}"
python3 - <<'PY'
import csv,gzip,hashlib,json,os,re,tarfile
from pathlib import Path
from datetime import datetime,timezone
from urllib.parse import unquote
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
out=root/'08_focused_families/validation/step25kz17.24'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+stamp)
out.mkdir(parents=True,exist_ok=False)
targets={'MC':'gene_22054','SP_Tet':'g10591','SP_Hex':'gene_Spor0461320.1'}
manifest=root/'00_admin/manifests/input_sets.tsv'
with manifest.open() as f:
 manifest_rows=list(csv.DictReader(f,delimiter='\t'))
config={}
for sp,gene in targets.items():
 rr=[r for r in manifest_rows if r['input_set']==sp]
 if len(rr)!=1:raise SystemExit('Ambiguous or missing manifest input_set: '+sp)
 r=rr[0]
 config[sp]=([gene],r['gff3'],r['genome'],r['protein_all'],r['cds_all'])
base=root/'08_focused_families/validation/step25kz17.23'
choices=[root/os.environ['DOMAIN_RELEASE']] if os.environ.get('DOMAIN_RELEASE') else [p for p in base.glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
if len(choices)!=1:raise SystemExit('Set DOMAIN_RELEASE to the intended step25kz17.23 release')
src=choices[0]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
for line in (src/'SHA256SUMS.txt').read_text().splitlines():
 digest,name=line.split(maxsplit=1)
 p=(src/name).resolve()
 if not p.is_relative_to(src.resolve()):raise SystemExit('Invalid checksum path')
 if sha(p)!=digest:raise SystemExit('Domain-release checksum mismatch: '+name)
with (src/'01_sequence_map.tsv').open() as f:
 baseline=[r for r in csv.DictReader(f,delimiter='\t') if r['protein_id'] in {sp+'__'+g for sp,g in targets.items()}]
if len(baseline)!=3:raise SystemExit('Expected three target records in domain release')
import shutil
checkpoint=out/'domain_checkpoint';checkpoint.mkdir()
for name in ['01_sequence_map.tsv','04_domain_evidence.tsv','07_alignment_completeness.tsv','02_candidates_and_references.faa','QC.json']:
 shutil.copy2(src/name,checkpoint/name)
def op(p):return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)
def table(p,fields,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def records(p):
 with op(p) as f:
  h=None;s=[]
  for l in f:
   if l.startswith('>'):
    if h is not None:yield h,''.join(s)
    h=l[1:].strip();s=[]
   elif l.strip():s.append(l.strip())
  if h is not None:yield h,''.join(s)
def write_seq(f,h,s):
 f.write('>'+h+'\n')
 for i in range(0,len(s),80):f.write(s[i:i+80]+'\n')
def attr(s):
 d={}
 for a in s.split(';'):
  if '=' in a:
   k,v=a.split('=',1);d[k]=[unquote(x) for x in v.split(',')]
 return d
status=[];sources=[];matches=[];intervals=[]
for sp,(genes,gff,genome,protein,cds) in config.items():
 d=out/sp;d.mkdir()
 mapping='03_harmonized/primary_v1/'+sp+'/'+sp+'.all_isoform_mapping.tsv'
 for kind,path in [('mapping',mapping),('gff',gff),('genome',genome),('protein',protein),('cds',cds)]:
  p=root/path;sources.append(dict(species=sp,kind=kind,relative_path=path,exists=p.is_file(),size_bytes=p.stat().st_size if p.is_file() else 0))
 if not (root/mapping).is_file():
  status.append(dict(species=sp,item='mapping',status='MISSING',count=0));continue
 with op(root/mapping) as f:
  reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames
  if not fields or not {'gene_id','seqid','start','end'}.issubset(fields):raise SystemExit('Unexpected mapping schema: '+mapping)
  rows=[r for r in reader if r.get('gene_id') in genes]
 table(d/'01_isoform_mapping.tsv',fields,rows)
 aliases=set(genes);windows=[]
 for gene in genes:
  rr=[r for r in rows if r['gene_id']==gene]
  status.append(dict(species=sp,item=gene,status='MAPPED' if rr else 'UNMAPPED',count=len(rr)))
  groups={}
  for r in rr:
   for k in ['gene_id','transcript_id','protein_id','cds_id']:
    if r.get(k) and r[k] not in ['NA','.']:aliases.add(r[k])
   groups.setdefault(r['seqid'],[]).append(r)
  for seqid,rs in groups.items():
   start=max(0,min(int(r['start']) for r in rs)-1-5000);end=max(int(r['end']) for r in rs)+5000
   windows.append(dict(species=sp,gene_id=gene,seqid=seqid,start_0based=start,end_exclusive_requested=end))
 # Collect all overlapping features, retaining original genomic coordinates.
 if (root/gff).is_file():
  selected=[]
  with op(root/gff) as f:
   for line in f:
    if line.startswith('##FASTA'):break
    c=line.rstrip('\n').split('\t')
    if len(c)==9 and not line.startswith('#'):
     if any(c[0]==v['seqid'] and int(c[4])>v['start_0based'] and int(c[3])<=v['end_exclusive_requested'] for v in windows):selected.append(c)
  changed=True
  while changed:
   old=set(aliases)
   for c in selected:
    aa=attr(c[8]);ids=set(sum([aa.get(k,[]) for k in ['ID','Parent','protein_id','transcript_id']],[]))
    if ids & aliases:aliases.update(ids)
   changed=aliases!=old
  n=0
  with op(root/gff) as f,(d/'02_locus_annotations.gff3').open('w') as w:
   w.write('##gff-version 3\n')
   for l in f:
    if l.startswith('##FASTA'):break
    if l.startswith('#'):continue
    c=l.rstrip('\n').split('\t')
    if len(c)!=9:continue
    try:a,b=int(c[3]),int(c[4])
    except ValueError:continue
    if any(c[0]==v['seqid'] and b>v['start_0based'] and a<=v['end_exclusive_requested'] for v in windows):
     w.write(l if l.endswith('\n') else l+'\n');n+=1
     aa=attr(c[8]);ids=set(sum([aa.get(k,[]) for k in ['ID','Parent','protein_id','transcript_id']],[]))
     # Only expand exact target-linked aliases, never neighboring-gene aliases.
     if ids & aliases:aliases.update(ids)
  status.append(dict(species=sp,item='GFF_features_in_windows',status='FOUND' if n else 'NO_COORDINATE_MATCH',count=n))
 else:status.append(dict(species=sp,item='GFF',status='MISSING',count=0))
 aliases |= {sp+'__'+x for x in list(aliases)}
 for kind,path in [('protein',protein),('cds',cds)]:
  n=0;seen=set()
  if (root/path).is_file():
   with (d/('03_'+kind+'_isoforms.fasta')).open('w') as w:
    for h,s in records(root/path):
     tokens={h.split()[0]}|set(re.findall(r'\[(?:protein_id|transcript_id|gene|locus_tag)=([^\]]+)\]',h))
     hit=tokens & aliases
     if hit:
      if h.split()[0] in seen:raise RuntimeError('Duplicate matched FASTA identifier: '+sp+' '+kind)
      seen.add(h.split()[0]);write_seq(w,h,s);n+=1
      matches.append(dict(species=sp,kind=kind,record_id=h.split()[0],matched_aliases=','.join(sorted(hit)),length=len(s)))
  status.append(dict(species=sp,item=kind,status='FOUND_REVIEW_COMPLETENESS' if n else 'MISSING_OR_NO_EXACT_ID_MATCH',count=n))
 # Stream the genome and retain only requested windows; do not load whole chromosomes.
 seqid=None;pos=0;chunks={i:[] for i in range(len(windows))};lengths={}
 if (root/genome).is_file():
  with op(root/genome) as f:
   for l in f:
    if l.startswith('>'):
     seqid=l[1:].split()[0];pos=0
     if seqid in {v['seqid'] for v in windows}:
      if seqid in lengths:raise RuntimeError('Duplicate target genomic sequence ID')
      lengths[seqid]=0
     continue
    s=''.join(l.split());n=len(s)
    for i,v in enumerate(windows):
     if seqid==v['seqid']:
      a=max(pos,v['start_0based']);b=min(pos+n,v['end_exclusive_requested'])
      if a<b:chunks[i].append(s[a-pos:b-pos])
    pos+=n
    if seqid in lengths:lengths[seqid]=pos
  with (d/'04_genomic_windows.fna').open('w') as w:
   for i,v in enumerate(windows):
    s=''.join(chunks[i]);end=min(v['end_exclusive_requested'],lengths.get(v['seqid'],0))
    record=dict(v,end_exclusive_actual=end,extracted_bp=len(s),status='EXTRACTED' if s else 'NO_SEQID_OR_INTERVAL_MATCH')
    intervals.append(record)
    if s:write_seq(w,v['gene_id']+'|'+v['seqid']+':'+str(v['start_0based'])+'-'+str(end)+'|zero_based_half_open|genomic_plus_strand',s)
 else:
  status.append(dict(species=sp,item='genome',status='MISSING',count=0))
 (d/'05_exact_aliases.txt').write_text('\n'.join(sorted(aliases))+'\n')
table(out/'SOURCE_FILES.tsv',['species','kind','relative_path','exists','size_bytes'],sources)
table(out/'COLLECTION_STATUS.tsv',['species','item','status','count'],status)
table(out/'SEQUENCE_MATCHES.tsv',['species','kind','record_id','matched_aliases','length'],matches)
table(out/'GENOMIC_INTERVALS.tsv',['species','gene_id','seqid','start_0based','end_exclusive_requested','end_exclusive_actual','extracted_bp','status'],intervals)
(out/'README.txt').write_text('Targeted evidence collection only. No gene repair, isoform replacement, acceptance or counting.\nGFF coordinates remain original 1-based inclusive. FASTA windows use genomic plus orientation; interval table is 0-based half-open.\nAnnotations include neighboring features in the windows. No sequence-ID renaming or fuzzy matching is performed.\nFOUND does not certify that every expected isoform was recovered: compare mapping and sequence-match tables.\nMissing sources or exact matches are reported; COMPLETE means collection finished, not biological validation.\nSource paths are relative to PROJECT_ROOT. Full input genomes are not copied.\n')
table(out/'TARGETS.tsv',['species','gene_id','review_status'],[dict(species=sp,gene_id=g,review_status='PENDING_GENE_MODEL_REVIEW') for sp,g in targets.items()])
(out/'QC.json').write_text(json.dumps({'status':'COLLECTION_FINISHED_NOT_BIOLOGICAL_VALIDATION','targets':3,'final_counts':False,'missing_or_unresolved':[r for r in status if r['count']==0]},indent=2)+'\n')
(out/'INPUT_PROVENANCE.json').write_text(json.dumps({'domain_release':os.path.relpath(src,root),'domain_manifest_sha256':sha(src/'SHA256SUMS.txt'),'input_manifest_sha256':sha(manifest),'flanking_bp':5000},indent=2)+'\n')
(out/'COMPLETE.txt').write_text('COLLECTION_FINISHED_REVIEW_STATUS_TABLES\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(out))+'\n')
bundle=out/'PPDK_RP_gene_model_bundle.tar.gz'
with tarfile.open(bundle,'w:gz') as t:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p!=bundle:t.add(p,arcname='PPDK_RP_gene_model_review/'+str(p.relative_to(out)))
print('Collection finished. Review COLLECTION_STATUS.tsv and SEQUENCE_MATCHES.tsv.')
print('Bundle relative to PROJECT_ROOT:',bundle.relative_to(root))
PY
