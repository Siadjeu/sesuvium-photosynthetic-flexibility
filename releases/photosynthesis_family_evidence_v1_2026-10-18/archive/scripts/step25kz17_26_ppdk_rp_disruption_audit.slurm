#!/usr/bin/env bash
#SBATCH --job-name=ppdk_rp_audit
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=ppdk_rp_audit_%j.out
#SBATCH --error=ppdk_rp_audit_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
python3 - <<'PY'
import csv,hashlib,json,os,re,tarfile,shutil
from pathlib import Path
from datetime import datetime,timezone
from urllib.parse import unquote
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True);base=root/'08_focused_families/validation'
choices=[root/os.environ['ALIGNMENT_RELEASE']] if os.environ.get('ALIGNMENT_RELEASE') else [p for p in (base/'step25kz17.25').glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
if len(choices)!=1:raise SystemExit('Set ALIGNMENT_RELEASE to the intended step17.25 release')
src=choices[0].resolve(strict=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def tab(name,fields,rows):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def attr(s):return dict((k,unquote(v)) for k,v in (x.split('=',1) for x in s.split(';') if '=' in x))
def fasta(p):return ''.join(l.strip() for l in p.read_text().splitlines() if not l.startswith('>')).upper()
def features(p):
 for l in p.read_text().splitlines():
  c=l.split('\t')
  if not l.startswith('#') and len(c)==9:yield c,attr(c[8])
for l in (src/'SHA256SUMS.txt').read_text().splitlines():
 h,n=l.split(maxsplit=1);p=(src/n).resolve()
 if not p.is_relative_to(src) or sha(p)!=h:raise SystemExit('Checksum failure: '+n)
out=base/'step25kz17.26'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
review=read(src/'05_alignment_review.tsv');intervals=read(src/'GENOMIC_INTERVALS.tsv');events=[];junctions=[];checks=[];exons=[];upstream=[]
for sp in ['MC','SP_Tet','SP_Hex']:
 d=src/sp;v=next(x for x in intervals if x['species']==sp);offset=int(v['start_0based']);seq=fasta(d/'window.fna');m=read(d/'source_01_isoform_mapping.tsv')[0]
 def get(a,b):
  if a<offset+1 or b>offset+len(seq) or b<a:return 'OUTSIDE_WINDOW'
  return seq[a-offset-1:b-offset]
 old=[(int(c[3]),int(c[4])) for c,a in features(d/'source_02_locus_annotations.gff3') if c[2]=='CDS' and m['transcript_id'] in a.get('Parent','').split(',')];old.sort()
 if m['strand']!='+':raise SystemExit('This focused audit expects the three known plus-strand loci')
 def junction_rows(parts,model,q):
  for (_,left),(right,_) in zip(parts,parts[1:]):
   a,b=left+1,right-1
   junctions.append(dict(species=sp,model=model,query=q,intron_start=a,intron_end=b,length=b-a+1,donor=get(a,a+1),acceptor=get(b-1,b),annotated_junction=(left,right) in [(x[1],y[0]) for x,y in zip(old,old[1:])],donor_context=get(a-10,a+11),acceptor_context=get(b-11,b+10)))
 junction_rows(old,'ANNOTATED',m['transcript_id'])
 selected={r['hit_id']:r for r in review if r['species']==sp and r['review_status']=='TARGET_OVERLAP_REVIEW'}
 for c,a in features(d/'04_miniprot.genomic.gff3'):
  if c[2]=='CDS' and a.get('Parent') in selected:
   hit=a['Parent'];start,end=int(c[3]),int(c[4]);exons.append(dict(species=sp,hit=hit,query=selected[hit]['query_id'],start=start,end=end,phase=c[7],exact_annotated_CDS=(start,end) in old,reported_frameshifts=a.get('Frameshift','NA'),reported_stops=a.get('StopCodon','NA')))
 for hit,r in selected.items():junction_rows(sorted((x['start'],x['end']) for x in exons if x['species']==sp and x['hit']==hit),hit,r['query_id'])
 # Associate each embedded PAF with its immediately following mRNA feature.
 pending=None
 for line in (d/'03_miniprot.local.gff3').read_text().splitlines():
  if line.startswith('##PAF\t'):pending=line.split('\t')[1:];continue
  c=line.split('\t')
  if len(c)!=9 or c[2]!='mRNA':continue
  a=attr(c[8]);hit=a['ID'];p=pending;pending=None
  if hit not in selected:continue
  if p is None or p[4]!='+':raise SystemExit('Missing PAF or unsupported strand')
  r=selected[hit]
  if p[0]!=r['query_id'] or int(p[7])+offset+1!=int(c[3])+offset:raise SystemExit('PAF/GFF association mismatch')
  tags={x.split(':',2)[0]:x.split(':',2)[2] for x in p[12:]};cg=tags['cg'];ops=re.findall(r'(\d+)([MIDFGNUV])',cg)
  if ''.join(n+op for n,op in ops)!=cg:raise SystemExit('Unsupported CIGAR '+cg)
  pos=int(p[7]);qpos=int(p[2]);nf=ns=0
  def event(kind,positions,q,op):
   coords=[offset+x+1 for x in positions];lo,hi=min(coords),max(coords)
   events.append(dict(species=sp,hit=hit,query=p[0],query_role=r['query_role'],event=kind,genomic_positions=','.join(map(str,coords)),start=lo,end=hi,query_position_1based=q+1,operation=op,bases=''.join(seq[x] for x in positions),context_start=max(offset+1,lo-20),context=get(max(offset+1,lo-20),min(offset+len(seq),hi+20))))
  for n,op in ops:
   n=int(n)
   if op in 'MD':
    for j in range(n):
     positions=list(range(pos+3*j,pos+3*j+3))
     if ''.join(seq[x] for x in positions) in ['TAA','TAG','TGA']:
      event('INTERNAL_STOP',positions,qpos+(j if op=='M' else 0),op);ns+=1
    pos+=3*n
    if op=='M':qpos+=n
   elif op=='I':qpos+=n
   elif op in 'FG':
    event('FRAMESHIFT',list(range(pos,pos+n)),qpos,op);nf+=1;pos+=n
    if op=='G':qpos+=1
   elif op=='N':pos+=n
   else:
    phase=1 if op=='U' else 2
    positions=list(range(pos,pos+phase))+list(range(pos+n-(3-phase),pos+n))
    if ''.join(seq[x] for x in positions) in ['TAA','TAG','TGA']:event('SPLIT_INTERNAL_STOP',positions,qpos,op);ns+=1
    pos+=n;qpos+=1
  checks.append(dict(species=sp,hit=hit,query=p[0],parsed_frameshifts=nf,reported_frameshifts=tags.get('fs','NA'),parsed_stops=ns,reported_stops=tags.get('st','NA'),target_consumption_ok=pos==int(p[8]),query_consumption_ok=qpos==int(p[3]),event_counts_match=str(nf)==tags.get('fs') and str(ns)==tags.get('st')))
 # Upstream context is supplied for inspection; no ORF is declared a rescued exon.
 a=offset+1;b=old[0][0]-1
 (out/(sp+'_upstream.fna')).write_text('>'+sp+'|'+v['seqid']+':'+str(a)+'-'+str(b)+'|plus_strand\n'+get(a,b)+'\n')
 upstream.append(dict(species=sp,seqid=v['seqid'],start=a,end=b,length=b-a+1,N_bases=get(a,b).count('N'),interpretation='UNTESTED_UPSTREAM_CONTEXT_NOT_ABSENCE_EVIDENCE'))
 shutil.copy2(d/'04_miniprot.genomic.gff3',out/(sp+'_alignment.genomic.gff3'))
 shutil.copy2(d/'source_02_locus_annotations.gff3',out/(sp+'_annotation.gff3'))
# Always emit headers, including when a category is empty.
tab('01_disruption_events.tsv',['species','hit','query','query_role','event','genomic_positions','start','end','query_position_1based','operation','bases','context_start','context'],events)
tab('02_splice_junctions.tsv',['species','model','query','intron_start','intron_end','length','donor','acceptor','annotated_junction','donor_context','acceptor_context'],junctions)
tab('03_exon_comparison.tsv',['species','hit','query','start','end','phase','exact_annotated_CDS','reported_frameshifts','reported_stops'],exons)
tab('04_parser_checks.tsv',['species','hit','query','parsed_frameshifts','reported_frameshifts','parsed_stops','reported_stops','target_consumption_ok','query_consumption_ok','event_counts_match'],checks)
tab('05_upstream_context.tsv',['species','seqid','start','end','length','N_bases','interpretation'],upstream)
# BED intervals for checking existing read alignments in a genome viewer.
with (out/'06_disruption_sites.bed').open('w') as f:
 for sp,a,b,kind in sorted({(r['species'],r['start'],r['end'],r['event']) for r in events}):
  chrom=next(v['seqid'] for v in intervals if v['species']==sp);f.write(f'{chrom}\t{a-1}\t{b}\t{sp}_{kind}\t0\t+\n')
bad=[r for r in checks if not all(r[k] for k in ['target_consumption_ok','query_consumption_ok','event_counts_match'])]
(out/'QC.json').write_text(json.dumps({'status':'AUDIT_COMPLETE' if not bad else 'PARSER_REVIEW_REQUIRED','alignments_audited':len(checks),'event_rows':len(events),'parser_discrepancies':bad,'final_counts':False},indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps({'source_release':os.path.relpath(src,root),'source_manifest_sha256':sha(src/'SHA256SUMS.txt'),'method':'miniprot protein CIGAR plus genomic window; no realignment'},indent=2)+'\n')
(out/'README.txt').write_text('Coordinates are original genomic 1-based inclusive, except BED uses 0-based half-open.\nOnly target-overlapping plus-strand hits are audited; flagged queries remain controls.\nCIGAR events describe alignment-implied disruptions, not validated variants or pseudogenes.\nTerminal stop codons outside the PAF aligned interval are excluded. Split codons have explicit genomic_positions; BED spans their full interval.\nCheck 04_parser_checks.tsv before interpreting exact event positions. Any discrepancy requires parser review.\nSplice dinucleotides and contexts are descriptive; canonical motifs do not validate an exon.\nUpstream FASTA is existing 5 kb context only; no new search or gene rescue is claimed.\nNo BAM/CRAM evidence is inspected and no gene count or annotation is changed.\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.iterdir()):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+p.name+'\n')
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='PPDK_RP_disruption_audit')
print('Upload:',bundle)
if bad:print('WARNING: parser discrepancies; see QC.json before biological interpretation')
PY
