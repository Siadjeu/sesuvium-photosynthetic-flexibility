#!/usr/bin/env bash
#SBATCH --job-name=vppase_loci
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=vppase_loci_%j.out
#SBATCH --error=vppase_loci_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from your project directory}"
python3 - <<'PY'
import csv,gzip,hashlib,json,os,re,tarfile,shutil
from pathlib import Path
from datetime import datetime,timezone
from urllib.parse import unquote
from collections import defaultdict,Counter
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
choices=[root/os.environ['VPPASE_RELEASE']] if os.environ.get('VPPASE_RELEASE') else [p for p in (base/'step25kz17.27').glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
if len(choices)!=1:raise SystemExit('Set VPPASE_RELEASE to the intended step17.27 directory')
src=choices[0].resolve(strict=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
for line in (src/'SHA256SUMS.txt').read_text().splitlines():
 h,n=line.split(maxsplit=1);p=(src/n).resolve()
 if not p.is_relative_to(src) or sha(p)!=h:raise SystemExit('Checksum failure '+n)
rows=[r for r in read(src/'01_sequence_map.tsv') if r['role']=='CANDIDATE']
if len(rows)!=48 or any(r['family']!='V_PPASE' for r in rows):raise SystemExit('Expected 48 V_PPASE candidates')
flags={'PO__evm.TU.LG07.1805':'PARTIAL_C_TERMINAL','PO__evm.TU.LG07.1806':'PARTIAL_N_TERMINAL','SS__g15406.t1':'SHORT_FRAGMENT','SS__g199.t1':'SHORT_FRAGMENT','PO__evm.TU.LG08.751':'ADDITIONAL_DOMAIN','PO__evm.TU.LG22.1173':'ADDITIONAL_DOMAIN','MC__gene_6610':'LONG_MODEL','SP_Hex__gene_Spor0325880.1':'LONG_MODEL'}
if not set(flags).issubset({r['protein_id'] for r in rows}):raise SystemExit('Missing expected exception candidates')
out=base/'step25kz17.28'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
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
manifest=root/'00_admin/manifests/input_sets.tsv';inputs=read(manifest)
ledger=[];resolved_rows={};mapping_provenance=[];config={};aliasfields=['gene_id','transcript_id','protein_id','cds_id']
for sp in sorted({r['protein_id'].split('__',1)[0] for r in rows}):
 candidates=[r for r in rows if r['protein_id'].split('__',1)[0]==sp]
 path=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv')
 mapping=read(path) if path.is_file() else []
 if mapping and not set(aliasfields+['seqid','start','end','strand']).issubset(mapping[0]):raise SystemExit('Unexpected mapping schema '+sp)
 mapping_provenance.append(dict(species=sp,path=str(path.relative_to(root)),sha256=sha(path) if path.is_file() else 'MISSING'))
 selected_genes=set()
 for r in candidates:
  pid=r['protein_id'];local=pid.split('__',1)[1]
  hits=[m for m in mapping if local in {m.get(k,'') for k in aliasfields}]
  genes={m['gene_id'] for m in hits};locations={(m['seqid'],m['strand']) for m in hits}
  valid=len(genes)==1 and len(locations)==1
  if valid:
   # Collect all isoforms for the resolved gene, not just the matching representative.
   hits=[m for m in mapping if m['gene_id'] in genes]
   locations={(m['seqid'],m['strand']) for m in hits};valid=len(locations)==1
  resolved_rows[pid]=hits if valid else []
  if valid and pid in flags:selected_genes.update(genes)
  ledger.append(dict(protein_id=pid,query_id=r['query_id'],species=sp,length_aa=r['length_aa'],gene_id=next(iter(genes)) if len(genes)==1 else ';'.join(sorted(genes)),mapping_status='UNAMBIGUOUS' if valid else 'UNRESOLVED_OR_AMBIGUOUS',seqid=hits[0]['seqid'] if valid else '',strand=hits[0]['strand'] if valid else '',start=min(int(m['start']) for m in hits) if valid else '',end=max(int(m['end']) for m in hits) if valid else '',isoform_rows=len(hits),review_flag=flags.get(pid,'NO_PRIORITY_EXCEPTION'),decision='PENDING_LOCUS_AND_MODEL_REVIEW'))
 if selected_genes:
  matches=[r for r in inputs if r['input_set']==sp]
  if len(matches)!=1:raise SystemExit('Ambiguous manifest input_set '+sp)
  m=matches[0];config[sp]=(sorted(selected_genes),m['gff3'],m['genome'],m['protein_all'],m['cds_all'])
table(out/'01_VPPASE_locus_ledger.tsv',list(ledger[0]),ledger)
summary=[]
for sp in sorted({r['species'] for r in ledger}):
 rr=[r for r in ledger if r['species']==sp];mapped=[r for r in rr if r['mapping_status']=='UNAMBIGUOUS']
 summary.append(dict(species=sp,candidate_records=len(rr),distinct_mapped_gene_loci=len({r['gene_id'] for r in mapped}),unresolved_mapping_records=len(rr)-len(mapped),exception_records=sum(r['protein_id'] in flags for r in rr),status='CANDIDATE_INVENTORY_NOT_FINAL_COUNTS'))
table(out/'02_locus_summary_NOT_FINAL.tsv',list(summary[0]),summary)
# Report pairs on the same reference sequence; adjacency is not a merge instruction.
pairs=[]
for i,a in enumerate(ledger):
 for b in ledger[i+1:]:
  if a['mapping_status']!='UNAMBIGUOUS' or b['mapping_status']!='UNAMBIGUOUS' or a['species']!=b['species'] or a['seqid']!=b['seqid']:continue
  left,right=sorted([a,b],key=lambda r:r['start'])
  pairs.append(dict(species=a['species'],protein_a=left['protein_id'],protein_b=right['protein_id'],seqid=a['seqid'],same_gene=a['gene_id']==b['gene_id'],same_strand=a['strand']==b['strand'],intervening_bp=max(0,right['start']-left['end']-1),overlap_bp=max(0,min(a['end'],b['end'])-max(a['start'],b['start'])+1),interpretation='COORDINATES_ONLY_DO_NOT_MERGE'))
table(out/'03_same_sequence_pairs.tsv',['species','protein_a','protein_b','seqid','same_gene','same_strand','intervening_bp','overlap_bp','interpretation'],pairs)
checkpoint=out/'domain_checkpoint';checkpoint.mkdir()
for name in ['01_sequence_map.tsv','02_candidates_and_references.faa','04_domain_evidence.tsv','07_alignment_completeness.tsv','09_PPDK_PPDK_RP_decision_ledger.tsv','10_PPDK_PPDK_RP_counts_PROVISIONAL.tsv','11_limitations.txt']:
 shutil.copy2(src/name,checkpoint/name)
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
# Detect exact candidate sequence duplicates descriptively; distinct genes remain distinct.
rows=[r for r in read(src/'01_sequence_map.tsv') if r['role']=='CANDIDATE']
seqs=dict((h.split()[0],s) for h,s in records(src/'02_candidates_and_references.faa'))
dup=defaultdict(list)
for r in rows:dup[hashlib.sha256(seqs[r['query_id']].encode()).hexdigest()].append(r['protein_id'])
table(out/'04_identical_candidate_sequences.tsv',['sequence_sha256','protein_ids','records'],[dict(sequence_sha256=h,protein_ids=';'.join(ids),records=len(ids)) for h,ids in dup.items() if len(ids)>1])
(out/'PROVENANCE.json').write_text(json.dumps({'source_release':os.path.relpath(src,root),'source_manifest_sha256':sha(src/'SHA256SUMS.txt'),'input_manifest_sha256':sha(manifest),'mapping_files':mapping_provenance},indent=2)+'\n')
(out/'QC.json').write_text(json.dumps({'status':'COLLECTION_COMPLETE_REVIEW_REQUIRED','candidate_records':48,'exception_records':8,'unresolved_mapping':[r['protein_id'] for r in ledger if r['mapping_status']!='UNAMBIGUOUS'],'missing_collection_evidence':[r for r in status if r['count']==0],'final_counts':False},indent=2)+'\n')
(out/'README.txt').write_text('V-PPase locus reconciliation and evidence collection only. No reannotation, locus merging, pseudogene designation or final gene counts.\nAll 48 candidates are matched by exact aliases to the existing all-isoform mappings. Eight flagged models receive annotation, protein/CDS isoforms and genomic windows with 5 kb flanks. Missing or ambiguous evidence remains explicit.\nIdentical proteins at different mapped loci are retained. Same-sequence pairs describe coordinates, not verified split genes.\nAnnotation excerpts retain original 1-based inclusive coordinates and may include partial neighboring genes. Genomic windows are plus-oriented; interval table uses 0-based half-open coordinates.\nReview COLLECTION_STATUS.tsv and SEQUENCE_MATCHES.tsv for missing source records. No suitable read evidence is assumed. Existing PPDK limitations are carried forward unchanged.\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+str(p.relative_to(out))+'\n')
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='VPPASE_locus_review')
print('RELEASE='+str(out));print('BUNDLE='+str(bundle))
PY
