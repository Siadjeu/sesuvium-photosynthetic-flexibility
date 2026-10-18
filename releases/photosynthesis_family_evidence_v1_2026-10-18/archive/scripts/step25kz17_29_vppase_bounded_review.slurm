#!/usr/bin/env bash
#SBATCH --job-name=vppase_bounded
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=vppase_bounded_%j.out
#SBATCH --error=vppase_bounded_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
command -v python3
command -v miniprot
python3 - <<'PY'
import csv,hashlib,json,os,subprocess,tarfile,shutil,itertools
from pathlib import Path
from datetime import datetime,timezone
from urllib.parse import unquote
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True);base=root/'08_focused_families/validation'
choices=[root/os.environ['LOCUS_RELEASE']] if os.environ.get('LOCUS_RELEASE') else [p for p in (base/'step25kz17.28').glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
if len(choices)!=1:raise SystemExit('Set LOCUS_RELEASE to intended step17.28 directory')
src=choices[0].resolve(strict=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta(p):
 d={};h=None;s=[]
 for l in p.read_text().splitlines():
  if l.startswith('>'):
   if h is not None:d[h]=''.join(s).upper()
   h=l[1:].split()[0];s=[]
   if h in d:raise ValueError('Duplicate FASTA ID '+h)
  elif l.strip():s.append(l.strip())
 if h is not None:d[h]=''.join(s).upper()
 return d
def attrs(s):return dict((k,unquote(v)) for k,v in (a.split('=',1) for a in s.split(';') if '=' in a))
def feats(p):
 for l in p.read_text().splitlines():
  c=l.split('\t')
  if not l.startswith('#') and len(c)==9:yield c,attrs(c[8])
def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
code=dict(zip(map(''.join,itertools.product('TCAG',repeat=3)),'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
def translate(s):return ''.join(code.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
def overlap(a,b,c,d):return max(0,min(b,d)-max(a,c)+1)
for l in (src/'SHA256SUMS.txt').read_text().splitlines():
 h,n=l.split(maxsplit=1);p=(src/n).resolve()
 if not p.is_relative_to(src) or sha(p)!=h:raise SystemExit('Checksum failure '+n)
ledger=read(src/'01_VPPASE_locus_ledger.tsv');targets=[r for r in ledger if r['review_flag']!='NO_PRIORITY_EXCEPTION']
if len(targets)!=8 or any(r['mapping_status']!='UNAMBIGUOUS' for r in targets):raise SystemExit('Expected eight unambiguously mapped exceptions')
version=subprocess.check_output(['miniprot','--version'],text=True).strip()
out=base/'step25kz17.29'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
def table(name,rows,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
seqs=fasta(src/'domain_checkpoint/02_candidates_and_references.faa');mapping=read(src/'domain_checkpoint/01_sequence_map.tsv')
# Use all 40 unflagged candidate proteins plus three reference controls; retain query identities.
flagids={r['protein_id'] for r in targets};queries=[r for r in mapping if r['protein_id'] not in flagids]
if len(queries)!=43:raise SystemExit('Expected 43 reference/support queries')
qfile=out/'01_queries.faa';qfile.write_text(''.join('>'+r['query_id']+'\n'+seqs[r['query_id']]+'\n' for r in queries));table('02_query_map.tsv',queries,list(queries[0]));qmap={r['query_id']:r for r in queries}
intervals=read(src/'GENOMIC_INTERVALS.tsv');consistency=[];alignments=[];status=[]
for i,r in enumerate(targets,1):
 sp,gene=r['species'],r['gene_id'];d=out/('L'+str(i).zfill(2));d.mkdir();source=src/sp
 vv=[x for x in intervals if x['species']==sp and x['gene_id']==gene]
 if len(vv)!=1 or vv[0]['status']!='EXTRACTED':raise SystemExit('Missing window '+gene)
 v=vv[0];offset=int(v['start_0based']);wins=fasta(source/'04_genomic_windows.fna');win=[s for h,s in wins.items() if h.split('|')[0]==gene]
 if len(win)!=1 or len(win[0])!=int(v['extracted_bp']):raise SystemExit('Window mismatch '+gene)
 genome=win[0];window=d/'window.fna';window.write_text('>WINDOW\n'+genome+'\n')
 fs=list(feats(source/'02_locus_annotations.gff3'));ms=[m for m in read(source/'01_isoform_mapping.tsv') if m['gene_id']==gene]
 proteins=fasta(source/'03_protein_isoforms.fasta');cdsseq=fasta(source/'03_cds_isoforms.fasta');targetcds=[]
 for m in ms:
  parts=[c for c,a in fs if c[2]=='CDS' and m['transcript_id'] in a.get('Parent','').split(',')];parts.sort(key=lambda c:int(c[3]),reverse=m['strand']=='-')
  if not parts:raise SystemExit('Missing CDS features '+gene)
  chunks=[]
  for c in parts:
   a,b=int(c[3]),int(c[4]);targetcds.append((a,b))
   if a<offset+1 or b>offset+len(genome):raise SystemExit('CDS outside window')
   s=genome[a-offset-1:b-offset];chunks.append(rc(s) if m['strand']=='-' else s)
  dna=''.join(chunks);phase=int(parts[0][7]) if parts[0][7] in ['0','1','2'] else 0;aa=translate(dna[phase:]);provided=cdsseq.get(m['cds_id']);protein=proteins.get(m['protein_id'])
  consistency.append(dict(protein_id=r['protein_id'],transcript_id=m['transcript_id'],cds_segments=len(parts),initial_phase=parts[0][7],genomic_CDS_bp=len(dna),provided_CDS_bp=len(provided) if provided is not None else 'MISSING',CDS_relation='EXACT' if dna==provided else ('TERMINAL_STOP_OMITTED' if dna[-3:] in ['TAA','TAG','TGA'] and dna[:-3]==provided else 'MISMATCH_OR_MISSING'),protein_matches=protein is not None and aa.removesuffix('*')==protein.removesuffix('*'),internal_stops=aa.removesuffix('*').count('*'),ambiguous_bases=sum(x not in 'ACGT' for x in dna),trailing_incomplete_codon=(len(dna)-phase)%3))
 # Annotated CDS of every candidate gene visible in this window, including the paired PO fragment.
 gene_transcripts={a['ID']:set(a.get('Parent','').split(',')) for c,a in fs if c[2] in ['mRNA','transcript'] and 'ID' in a}
 by_gene={}
 for c,a in fs:
  if c[2]=='CDS':
   for parent in a.get('Parent','').split(','):
    for g in gene_transcripts.get(parent,set()):by_gene.setdefault(g,[]).append((int(c[3]),int(c[4])))
 raw=d/'03_miniprot.local.gff3';cmd=['miniprot','-t',os.environ.get('SLURM_CPUS_PER_TASK','4'),'-G','20000','--gff','--aln','--outs=0.5',str(window),str(qfile)]
 with raw.open('w') as f,(d/'03_miniprot.log').open('w') as log:subprocess.run(cmd,stdout=f,stderr=log,check=True)
 pred=list(feats(raw));hits=0
 with (d/'04_miniprot.genomic.gff3').open('w') as f:
  f.write('##gff-version 3\n')
  for c,a in pred:
   c=c.copy();c[0]=v['seqid'];c[3]=str(int(c[3])+offset);c[4]=str(int(c[4])+offset);f.write('\t'.join(c)+'\n')
 for c,a in pred:
  if c[2]!='mRNA':continue
  hits+=1;hit=a['ID'];q,qs,qe=a['Target'].split()[:3];parts=[(int(cc[3])+offset,int(cc[4])+offset) for cc,aa in pred if cc[2]=='CDS' and hit in aa.get('Parent','').split(',')]
  overlapped=sorted(g for g,cd in by_gene.items() if any(overlap(x,y,u,w)>0 for x,y in parts for u,w in cd))
  shared=sum(overlap(x,y,u,w) for x,y in parts for u,w in set(targetcds))
  alignments.append(dict(window_id=d.name,species=sp,target_gene=gene,hit_id=hit,query_id=q,query_protein=qmap[q]['protein_id'],query_role=qmap[q]['role'],query_coverage=round((int(qe)-int(qs)+1)/len(seqs[q]),4),identity=a.get('Identity','NA'),strand=c[6],same_strand=c[6]==r['strand'],start=int(c[3])+offset,end=int(c[4])+offset,target_CDS_overlap_bp=shared,overlapping_gene_CDS=';'.join(overlapped),spans_PO_pair=sp=='PO' and {'evm.TU.LG07.1805','evm.TU.LG07.1806'}.issubset(overlapped),frameshifts=a.get('Frameshift','NA'),internal_stops=a.get('StopCodon','NA'),window_edge_distance_bp=min(int(c[3])-1,len(genome)-int(c[4])),decision='REVIEW_ONLY'))
 status.append(dict(window_id=d.name,protein_id=r['protein_id'],seqid=v['seqid'],start_0based=offset,end_exclusive=int(v['end_exclusive_actual']),reported_alignments=hits))
 shutil.copy2(source/'02_locus_annotations.gff3',d/'source_annotations.gff3')
table('05_CDS_protein_consistency.tsv',consistency,list(consistency[0]));table('06_alignment_review.tsv',alignments,['window_id','species','target_gene','hit_id','query_id','query_protein','query_role','query_coverage','identity','strand','same_strand','start','end','target_CDS_overlap_bp','overlapping_gene_CDS','spans_PO_pair','frameshifts','internal_stops','window_edge_distance_bp','decision']);table('07_window_map.tsv',status,list(status[0]))
shutil.copy2(src/'01_VPPASE_locus_ledger.tsv',out/'08_candidate_ledger_unchanged.tsv')
(out/'QC.json').write_text(json.dumps({'status':'BOUNDED_REVIEW_COMPLETE_NOT_FINAL','loci':8,'queries':43,'alignments':len(alignments),'consistency_issues':[r for r in consistency if not r['protein_matches'] or r['CDS_relation']=='MISMATCH_OR_MISSING' or r['internal_stops'] or r['ambiguous_bases'] or r['trailing_incomplete_codon']],'final_counts':False},indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps({'source_release':os.path.relpath(src,root),'source_manifest_sha256':sha(src/'SHA256SUMS.txt'),'miniprot_version':version,'alignment_options':'-G 20000 --gff --aln --outs=0.5'},indent=2)+'\n')
(out/'README.txt').write_text('Bounded evidence review only: eight existing windows, 40 unflagged candidates and three references. No annotation modification or final counts.\nGenomic GFF uses original 1-based coordinates; local GFF uses WINDOW coordinates. Both strands are retained; use same_strand and CDS overlap before interpreting hits.\nPO pair tested in two overlapping windows; these are not independent biological evidence. A bridging alignment is not proof that models should be merged.\nMissing frameshift/stop attributes are NA, not explicitly reported zero. Query coverage is alignment span.\nGenome CDS is spliced in transcript order, with initial phase applied once for translation. Terminal stop inclusion differences are recorded.\nUnflagged proteins are supporting queries, not experimentally validated controls.\nNo read evidence or new genomic windows are obtained. Retain unresolved cases after review and move to the next family.\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+str(p.relative_to(out))+'\n')
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='VPPASE_bounded_review')
print('BUNDLE='+str(bundle))
PY
