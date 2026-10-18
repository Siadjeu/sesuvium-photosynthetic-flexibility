#!/usr/bin/env bash
#SBATCH --job-name=ppdk_rp_genome
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=ppdk_rp_genome_%j.out
#SBATCH --error=ppdk_rp_genome_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
command -v python3
command -v miniprot || { echo 'ERROR: activate an environment containing miniprot.' >&2; exit 1; }
python3 - <<'PY'
import csv,hashlib,json,os,shutil,subprocess,tarfile
from pathlib import Path
from datetime import datetime,timezone
from urllib.parse import unquote
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
choices=[root/os.environ['MODEL_RELEASE']] if os.environ.get('MODEL_RELEASE') else [p for p in (base/'step25kz17.24').glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
if len(choices)!=1:raise SystemExit('Set MODEL_RELEASE to the intended step25kz17.24 release directory')
src=choices[0].resolve(strict=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,fields,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def fasta(p):
 d={};h=None;s=[]
 for l in p.read_text().splitlines():
  if l.startswith('>'):
   if h is not None:d[h]=''.join(s)
   h=l[1:].split()[0]
   if h in d:raise ValueError('Duplicate FASTA identifier '+h)
   s=[]
  elif l.strip():s.append(l.strip())
 if h is not None:d[h]=''.join(s)
 return d
def attrs(x):return {k:unquote(v) for k,v in (a.split('=',1) for a in x.split(';') if '=' in a)}
def features(p):
 for l in p.read_text().splitlines():
  if l.startswith('#'):continue
  c=l.split('\t')
  if len(c)==9:yield c,attrs(c[8])
def overlap(a,b,c,d):return max(0,min(b,d)-max(a,c)+1)
def union(intervals):
 merged=[]
 for a,b in sorted(intervals):
  if merged and a<=merged[-1][1]+1:merged[-1][1]=max(b,merged[-1][1])
  else:merged.append([a,b])
 return merged
for l in (src/'SHA256SUMS.txt').read_text().splitlines():
 h,n=l.split(maxsplit=1);p=(src/n).resolve()
 if not p.is_relative_to(src) or sha(p)!=h:raise SystemExit('Input checksum failure: '+n)
targets={'MC':'gene_22054','SP_Tet':'g10591','SP_Hex':'gene_Spor0461320.1'}
checkpoint=src/'domain_checkpoint'
seqs=fasta(checkpoint/'02_candidates_and_references.faa')
queries=[r for r in read(checkpoint/'01_sequence_map.tsv') if r['family']=='PPDK_RP']
if len(queries)!=19:raise SystemExit('Expected 17 candidate and 2 reference PPDK-RP queries')
flagged={sp+'__'+g for sp,g in targets.items()}
for r in queries:
 r['alignment_role']='FLAGGED_MODEL_CONTROL' if r['protein_id'] in flagged else ('REFERENCE' if r['role']=='REFERENCE' else 'SUPPORTED_CANDIDATE')
 if r['query_id'] not in seqs or len(seqs[r['query_id']])!=int(r['length_aa']):raise SystemExit('Query FASTA/map mismatch')
version=subprocess.check_output(['miniprot','--version'],text=True).strip()
out=base/'step25kz17.25'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False)
qfile=out/'01_PPDK_RP_queries.faa'
qfile.write_text(''.join('>'+r['query_id']+'\n'+seqs[r['query_id']]+'\n' for r in queries))
table(out/'02_query_map.tsv',list(queries[0]),queries)
qmap={r['query_id']:r for r in queries};intervals=read(src/'GENOMIC_INTERVALS.tsv');summary=[];counts=[];commands=[]
for sp,gene in targets.items():
 d=out/sp;d.mkdir();source=src/sp
 for name in ['01_isoform_mapping.tsv','02_locus_annotations.gff3','03_protein_isoforms.fasta','03_cds_isoforms.fasta','04_genomic_windows.fna']:
  shutil.copy2(source/name,d/('source_'+name))
 vv=[v for v in intervals if v['species']==sp and v['gene_id']==gene]
 if len(vv)!=1 or vv[0]['status']!='EXTRACTED':raise SystemExit('Missing genomic window: '+sp)
 v=vv[0];offset=int(v['start_0based']);genomes=fasta(source/'04_genomic_windows.fna')
 if len(genomes)!=1:raise SystemExit('Expected one genomic window per species')
 genome=next(iter(genomes.values()));length=len(genome)
 if length!=int(v['extracted_bp']):raise SystemExit('Window length mismatch')
 window=d/'window.fna';window.write_text('>WINDOW\n'+genome+'\n')
 mapping=read(source/'01_isoform_mapping.tsv')
 if len(mapping)!=1 or mapping[0]['gene_id']!=gene:raise SystemExit('Unexpected locus mapping')
 m=mapping[0];old=list(features(source/'02_locus_annotations.gff3'))
 target_cds=union([(int(c[3]),int(c[4])) for c,a in old if c[2]=='CDS' and m['transcript_id'] in a.get('Parent','').split(',')])
 if not target_cds:raise SystemExit('No annotated target CDS')
 genes=[(a.get('ID',''),int(c[3]),int(c[4])) for c,a in old if c[2]=='gene']
 cmd=['miniprot','-t',os.environ.get('SLURM_CPUS_PER_TASK','4'),'-G','20000','--gff','--aln','--outs=0.5',str(window),str(qfile)]
 commands.append({'species':sp,'argv':['<WINDOW>' if x==str(window) else '<QUERY_FASTA>' if x==str(qfile) else x for x in cmd]})
 raw=d/'03_miniprot.local.gff3'
 with raw.open('w') as f,(d/'03_miniprot.log').open('w') as log:subprocess.run(cmd,stdout=f,stderr=log,check=True)
 fs=list(features(raw));mrnas=[(c,a) for c,a in fs if c[2]=='mRNA'];mapped=set()
 with (d/'04_miniprot.genomic.gff3').open('w') as f:
  f.write('##gff-version 3\n')
  for c,a in fs:
   c=c.copy()
   if c[0]!='WINDOW':raise SystemExit('Unexpected target sequence')
   c[0]=v['seqid'];c[3]=str(int(c[3])+offset);c[4]=str(int(c[4])+offset);f.write('\t'.join(c)+'\n')
 for c,a in mrnas:
  target=a.get('Target','').split()
  if len(target)<3 or target[0] not in qmap:raise SystemExit('Unexpected miniprot Target attribute')
  q=target[0];mapped.add(q);hit=a['ID'];start=int(c[3])+offset;end=int(c[4])+offset
  cds=union([(int(cc[3])+offset,int(cc[4])+offset) for cc,aa in fs if cc[2]=='CDS' and hit in aa.get('Parent','').split(',')])
  shared=sum(overlap(x,y,u,w) for x,y in cds for u,w in target_cds)
  others=sorted({g for g,u,w in genes if g!=gene and any(overlap(x,y,u,w) for x,y in cds)})
  row=dict(species=sp,target_gene=gene,hit_id=hit,query_id=q,protein_id=qmap[q]['protein_id'],query_role=qmap[q]['alignment_role'],query_length=len(seqs[q]),query_start=int(target[1]),query_end=int(target[2]),query_coverage=round((int(target[2])-int(target[1])+1)/len(seqs[q]),4),seqid=v['seqid'],start=start,end=end,strand=c[6],same_strand=c[6]==m['strand'],target_CDS_overlap_bp=shared,other_gene_span_overlap=';'.join(others),identity=a.get('Identity','NA'),frameshift=a.get('Frameshift','NA'),stop_codon=a.get('StopCodon','NA'),extends_left_bp=max(0,min(x for x,y in target_cds)-min((x for x,y in cds),default=start)),extends_right_bp=max(0,max((y for x,y in cds),default=end)-max(y for x,y in target_cds)),window_edge_distance_bp=min(int(c[3])-1,length-int(c[4])),review_status='TARGET_OVERLAP_REVIEW' if shared and c[6]==m['strand'] else 'OFF_TARGET_OR_OPPOSITE_STRAND')
  summary.append(row)
 for q in qmap:counts.append(dict(species=sp,query_id=q,status='ALIGNMENT_REPORTED' if q in mapped else 'NO_REPORTED_ALIGNMENT'))
fields=['species','target_gene','hit_id','query_id','protein_id','query_role','query_length','query_start','query_end','query_coverage','seqid','start','end','strand','same_strand','target_CDS_overlap_bp','other_gene_span_overlap','identity','frameshift','stop_codon','extends_left_bp','extends_right_bp','window_edge_distance_bp','review_status']
table(out/'05_alignment_review.tsv',fields,summary)
table(out/'06_query_status.tsv',['species','query_id','status'],counts)
shutil.copy2(src/'GENOMIC_INTERVALS.tsv',out/'GENOMIC_INTERVALS.tsv')
(out/'PROVENANCE.json').write_text(json.dumps({'source_release':os.path.relpath(src,root),'source_manifest_sha256':sha(src/'SHA256SUMS.txt'),'miniprot_version':version,'commands':commands},indent=2)+'\n')
(out/'QC.json').write_text(json.dumps({'status':'ALIGNMENTS_COMPLETE_REVIEW_REQUIRED','queries':len(queries),'windows':3,'reported_alignments':len(summary),'target_overlapping_alignments':sum(r['review_status']=='TARGET_OVERLAP_REVIEW' for r in summary),'final_counts':False},indent=2)+'\n')
(out/'README.txt').write_text('Exploratory protein-to-genome alignments; no annotation replacement or final counts.\n14 supported candidate queries plus 2 reference queries and 3 flagged-model controls.\nAll queries searched against each window. Flagged controls are not independent rescue evidence.\nLocal GFF coordinates refer to WINDOW, genomic GFF coordinates to original chromosomes (1-based inclusive).\nOriginal annotation excerpts can contain only part of neighboring genes. Neighbor overlap column means predicted CDS overlaps a neighboring gene span, not necessarily its CDS.\nTarget overlap uses CDS intervals, not the intron-containing alignment span.\nExtension sizes are descriptive and meaningful only for target-associated alignments; they do not establish gene repair.\nQuery coverage is aligned query span, not identical-residue coverage. Missing GFF attributes are NA, not zero.\nMaximum intron is 20 kb; windows have only 5 kb flanks. No hit is not absence; edge hits may require larger windows.\nRaw residue alignments and logs are retained. Proposed models may reflect paralogs or alignment artifacts.\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+str(p.relative_to(out))+'\n')
# Exclude the checksum file itself from its own manifest.
p=out/'SHA256SUMS.txt';p.write_text(''.join(l for l in p.read_text().splitlines(True) if not l.rstrip().endswith('  SHA256SUMS.txt')))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='PPDK_RP_genome_alignment')
print('Release:',out)
print('Upload bundle:',bundle)
PY
