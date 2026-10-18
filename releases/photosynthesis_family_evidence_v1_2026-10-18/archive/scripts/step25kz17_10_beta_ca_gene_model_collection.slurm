#!/usr/bin/env bash
#SBATCH --job-name=beta_ca_models
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=beta_ca_models_%j.out
#SBATCH --error=beta_ca_models_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from your project directory}"
python3 - <<'PY'
import csv,gzip,hashlib,json,os,re,tarfile
from pathlib import Path
from datetime import datetime,timezone
from urllib.parse import unquote
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
out=root/'08_focused_families/validation/step25kz17.10'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+stamp)
out.mkdir(parents=True,exist_ok=False)
slp='01_inputs/raw/SLP/Dahong_haplotype_resolved_v1/figshare_27999314_v1/'
config={
 'AT':(['gene-LOC130820426'],'01_inputs/standardized/AT/annotation/AT.gff3.gz','01_inputs/standardized/AT/genome/AT.genome.fa.gz','01_inputs/standardized/AT/protein/AT.proteins.all.fa.gz','01_inputs/standardized/AT/cds/AT.cds.all.fa.gz'),
 'SP_Hex':(['gene_Spor0140320.1','gene_Spor0527280.1'],'01_inputs/standardized/SP_Hex/annotation/SP_Hex.gff3.gz','01_inputs/standardized/SP_Hex/genome/SP_Hex.genome.fa.gz','01_inputs/standardized/SP_Hex/protein/SP_Hex.proteins.all.fa.gz','01_inputs/standardized/SP_Hex/cds/SP_Hex.cds.all.fa.gz'),
 'SLP_Hap1':(['SP3G00376'],'01_inputs/standardized/SLP/Hap1/annotation/SLP_Hap1.cr_normalized.gff3.gz','01_inputs/standardized/SLP/Hap1/genome/SLP_Hap1.genome.oriseqid.fa.gz',slp+'hap1_prot.fa.gz',slp+'hap1.CDS.fa.gz')}
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
  reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=[r for r in reader if r.get('gene_id') in genes]
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
(out/'COMPLETE.txt').write_text('COLLECTION_FINISHED_REVIEW_STATUS_TABLES\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(out))+'\n')
bundle=out/'BETA_CA_gene_model_bundle.tar.gz'
with tarfile.open(bundle,'w:gz') as t:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p!=bundle:t.add(p,arcname='BETA_CA_gene_model_review/'+str(p.relative_to(out)))
print('Collection finished. Review COLLECTION_STATUS.tsv and SEQUENCE_MATCHES.tsv.')
print('Bundle relative to PROJECT_ROOT:',bundle.relative_to(root))
PY
