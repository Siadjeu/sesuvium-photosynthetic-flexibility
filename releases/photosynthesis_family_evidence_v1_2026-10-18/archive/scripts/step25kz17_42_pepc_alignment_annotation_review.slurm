#!/usr/bin/env bash
#SBATCH --job-name=pepc_align_review
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=pepc_align_review_%j.out
#SBATCH --error=pepc_align_review_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v mafft
python3 - <<'PY'
import csv,gzip,hashlib,json,os,re,shutil,subprocess,tarfile
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True);base=root/'08_focused_families/validation'
def need(v,msg):
 if not v:raise SystemExit(msg)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def verify(p):
 for line in (p/'SHA256SUMS.txt').read_text().splitlines():
  h,n=line.split(maxsplit=1);f=(p/n.lstrip('*')).resolve()
  need(f.is_relative_to(p) and f.is_file() and sha(f)==h,'Checksum mismatch: '+n)
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta(p):
 d={};headers={};key=None
 with p.open() as f:
  for line in f:
   if line.startswith('>'):
    key=line[1:].split()[0];need(key not in d,'Duplicate FASTA ID '+key);d[key]='';headers[key]=line[1:].strip()
   elif line.strip():
    need(key is not None,'Invalid FASTA');d[key]+=''.join(line.split()).upper()
 return d,headers
def norm(s):return s[:-1] if s.endswith('*') else s
def valid(s):return bool(s) and set(s)<=set('ACDEFGHIKLMNPQRSTVWYBXZJUO') and not set(s)<=set('ACGTN')
pp=[root/os.environ['PEPC_DOMAIN_RELEASE']] if os.environ.get('PEPC_DOMAIN_RELEASE') else [p for p in (base/'step25kz17.41').glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
need(len(pp)==1,'Set PEPC_DOMAIN_RELEASE to the intended step25kz17.41 release')
src=pp[0].resolve(strict=True);verify(src)
need(json.loads((src/'QC.json').read_text())['all_reference_controls_have_PF00311'],'Reference control failure in source')
proteins,_=fasta(src/'01_PEPC_corrected_53_candidates.faa');refs,_=fasta(src/'02_PPC_reference_controls.faa')
need(len(proteins)==53 and len(refs)==4,'Expected 53 proteins and four controls')
queries={**proteins,**refs};need(len(queries)==57 and all(valid(s) for s in queries.values()),'Invalid sequence inputs')
meta={r['protein_id']:r for r in read(src/'10_prior_provenance.tsv')}
hits=defaultdict(list)
for r in read(src/'08_all_Pfam_domain_hits.tsv'):
 if r['sequence_id'] in queries and r['accession'].split('.')[0]=='PF00311':
  for f in ['hmm_from','hmm_to','ali_from','ali_to','model_length']:r[f]=int(r[f])
  hits[r['sequence_id']].append(r)
out=base/'step25kz17.42'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'annotation').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+s+'\n' for k,s in d.items()))
reconciled=[]
for k in queries:
 hh=sorted(hits[k],key=lambda x:x['ali_from']);need(hh,'Missing domain evidence '+k)
 need(len({x['model_length'] for x in hh})==1,'Model lengths differ')
 covered=set()
 for h in hh:covered.update(range(h['hmm_from'],h['hmm_to']+1))
 cov=len(covered)/hh[0]['model_length']
 ordered=all(b['ali_from']>a['ali_to'] and b['hmm_from']>a['hmm_from'] and b['hmm_to']>a['hmm_to'] for a,b in zip(hh,hh[1:]))
 overlap=sum(max(0,min(a['hmm_to'],b['hmm_to'])-max(a['hmm_from'],b['hmm_from'])+1) for a,b in zip(hh,hh[1:]))
 category='SINGLE_HIT_BROAD_COVERAGE' if len(hh)==1 and cov>=.9 else 'ORDERED_SPLIT_BROAD_COVERAGE' if len(hh)>1 and ordered and overlap<=20 and cov>=.9 else 'PRIORITY_COMPLETENESS_REVIEW'
 reconciled.append(dict(protein_id=k,role='CANDIDATE' if k in proteins else 'REFERENCE',hit_count=len(hh),union_model_coverage=round(cov,6),ordered_nonoverlapping_protein_segments=ordered,adjacent_model_overlap_residues=overlap,protein_intersegment_gaps=';'.join(str(b['ali_from']-a['ali_to']-1) for a,b in zip(hh,hh[1:])),model_intersegment_gaps=';'.join(str(b['hmm_from']-a['hmm_to']-1) for a,b in zip(hh,hh[1:])),evidence_category=category,formal_decision='PENDING'))
table('01_reconciled_domain_evidence.tsv',reconciled,list(reconciled[0]))
priority={r['protein_id'] for r in reconciled if r['role']=='CANDIDATE' and r['evidence_category']=='PRIORITY_COMPLETENESS_REVIEW'}
need(priority=={'SP_Hex__gene_Spor0168910.1','SP_Hex__gene_Spor0380090.1'},'Unexpected priority set; inspect reconciled table before expanding annotation scope')
writefa('02_53_candidates_four_controls.faa',queries)
with (out/'03_PEPC_review_alignment.faa').open('w') as a,(out/'03_mafft.log').open('w') as e:
 subprocess.run(['mafft','--amino','--thread',os.environ.get('SLURM_CPUS_PER_TASK','4'),'--localpair','--maxiterate','1000',str(out/'02_53_candidates_four_controls.faa')],stdout=a,stderr=e,check=True)
aln,_=fasta(out/'03_PEPC_review_alignment.faa')
need(set(aln)==set(queries) and len({len(s) for s in aln.values()})==1,'Invalid alignment dimensions or IDs')
need(all(aln[k].replace('-','')==s for k,s in queries.items()),'Alignment changed a sequence')
metrics=[];gaprows=[]
for k in proteins:
 for ref in refs:
  a=aln[k];b=aln[ref];cols=[i for i,c in enumerate(b) if c!='-'];present=[j for j,i in enumerate(cols) if a[i]!='-']
  paired=len(present);matches=sum(a[i]==b[i] for i in cols if a[i]!='-')
  nmissing=present[0] if present else len(cols);cmissing=len(cols)-present[-1]-1 if present else len(cols)
  internal=sum(a[cols[j]]=='-' for j in range(present[0],present[-1]+1)) if present else 0
  insertions=sum(x!='-' and y=='-' for x,y in zip(a,b))
  metrics.append(dict(protein_id=k,reference_id=ref,reference_length=len(cols),paired_residues=paired,reference_occupancy=round(paired/len(cols),6),identity_over_paired=round(matches/paired,6) if paired else 0,N_terminal_reference_gaps=nmissing,C_terminal_reference_gaps=cmissing,internal_reference_gaps=internal,candidate_residues_in_reference_gap_columns=insertions))
  if k in priority:
   start=None
   for j in range(len(cols)+1):
    missing=j<len(cols) and a[cols[j]]=='-'
    if missing and start is None:start=j
    if not missing and start is not None:
     gaprows.append(dict(protein_id=k,reference_id=ref,reference_start_1based=start+1,reference_end_1based=j,gap_length=j-start,alignment_start_1based=cols[start]+1,alignment_end_1based=cols[j-1]+1));start=None
table('04_all_candidate_reference_alignment_metrics.tsv',metrics,list(metrics[0]))
table('05_priority_reference_gap_intervals.tsv',gaprows,['protein_id','reference_id','reference_start_1based','reference_end_1based','gap_length','alignment_start_1based','alignment_end_1based'])
# Bounded annotation collection: exact mapping aliases, two genes, +/-5 kb GFF windows.
sp='SP_Hex';mp=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv')
manifest=root/'00_admin/manifests/input_sets.tsv'
need(mp.is_file() and manifest.is_file(),'Missing mapping or input manifest')
genes={meta[k]['gene_id'] for k in priority};mapping=[r for r in read(mp) if r['gene_id'] in genes]
need({r['gene_id'] for r in mapping}==genes,'Exception gene absent from mapping')
configs=[r for r in read(manifest) if r['input_set']==sp];need(len(configs)==1,'Ambiguous SP_Hex input manifest')
conf=configs[0];aliases=defaultdict(set);windows=[];loc=[]
for gene in sorted(genes):
 mm=[r for r in mapping if r['gene_id']==gene];need(len({(r['seqid'],r['strand']) for r in mm})==1,'Ambiguous gene coordinates')
 start=min(int(r['start']) for r in mm);end=max(int(r['end']) for r in mm)
 windows.append((mm[0]['seqid'],max(1,start-5000),end+5000))
 loc.append(dict(gene_id=gene,seqid=mm[0]['seqid'],strand=mm[0]['strand'],start=start,end=end,isoform_mapping_rows=len(mm)))
 for r in mm:
  for f in ['gene_id','transcript_id','protein_id','cds_id']:
   v=r.get(f,'')
   if v and v not in {'.','NA'}:
    aliases[v].add(gene);aliases[sp+'__'+v].add(gene)
for k in priority:need(meta[k]['gene_id'] in aliases.get(k,set()),'Candidate ID not confirmed by mapping: '+k)
table('annotation/01_all_isoform_mapping.tsv',mapping,list(mapping[0]));table('annotation/02_locus_windows.tsv',loc,list(loc[0]))
def op(p):return gzip.open(p,'rt') if p.suffix=='.gz' else p.open()
def selected_sequences(p):
 d={};key=None;parts=[]
 def commit():
  if key in aliases:
   need(key not in d,'Duplicate annotation FASTA identifier');d[key]=''.join(parts).upper()
 with op(p) as f:
  for line in f:
   if line.startswith('>'):commit();key=line[1:].split()[0];parts=[]
   elif key in aliases:parts.append(''.join(line.split()))
  commit()
 return d
status=[];provenance=[];extracted={};matches=[]
for kind,col in [('gff','gff3'),('protein','protein_all'),('cds','cds_all')]:
 p=root/conf[col];source_hash=sha(p) if p.is_file() else 'MISSING'
 provenance.append(dict(kind=kind,source=str(p.relative_to(root)) if p.is_relative_to(root) else '<EXTERNAL>/'+p.name,sha256=source_hash))
 if not p.is_file():status.append(dict(kind=kind,status='MISSING',records=0));continue
 n=0
 if kind=='gff':
  with op(p) as f,(out/'annotation/03_context_windows.gff3').open('w') as w:
   w.write('##gff-version 3\n')
   for line in f:
    if line.startswith('##FASTA'):break
    if line.startswith('#'):continue
    c=line.rstrip('\n').split('\t')
    if len(c)==9 and any(c[0]==chrom and int(c[3])<=b and int(c[4])>=a for chrom,a,b in windows):w.write(line.rstrip('\n')+'\n');n+=1
 else:
  d=selected_sequences(p);extracted[kind]=d;writefa('annotation/04_'+kind+'_isoforms.fasta',d);n=len(d)
  for key,s in d.items():matches.append(dict(kind=kind,source_id=key,mapped_gene_ids=';'.join(sorted(aliases[key])),length=len(s),sequence_sha256=hashlib.sha256(s.encode()).hexdigest()))
 status.append(dict(kind=kind,status='COLLECTED' if n else 'NO_EXACT_MATCH',records=n))
 need(sha(p)==source_hash,'Annotation input changed during read')
table('annotation/05_collection_status.tsv',status,['kind','status','records']);table('annotation/06_input_provenance.tsv',provenance,['kind','source','sha256']);table('annotation/07_exact_sequence_aliases.tsv',matches,['kind','source_id','mapped_gene_ids','length','sequence_sha256'])
consistency=[]
for k in sorted(priority):
 gene=meta[k]['gene_id'];available=[(a,s) for a,s in extracted.get('protein',{}).items() if gene in aliases[a]]
 same=[a for a,s in available if norm(s)==proteins[k]]
 consistency.append(dict(protein_id=k,gene_id=gene,exact_sequence_matching_annotation_ids=';'.join(same),annotation_protein_status='EXACT_SEQUENCE_CONFIRMED' if same else 'NO_MATCH_REQUIRES_PROVENANCE_REVIEW',alignment_review='PENDING',annotation_review='PENDING',family_membership='PENDING',model_limitation=''))
table('06_priority_review_queue.tsv',consistency,list(consistency[0]))
for name in ['05_CDS_exception_ledger.tsv','11_integrated_locus_evidence_UNCHANGED.tsv','12_final_Wang_counts_PENDING.tsv','HISTORICAL_TREE_STATUS.txt']:
 shutil.copyfile(src/name,out/('prior_'+name))
qc=dict(status='BOUNDED_ALIGNMENT_AND_ANNOTATION_REVIEW_READY',candidate_proteins=53,reference_controls=4,alignment_columns=len(next(iter(aln.values()))),domain_categories=dict(Counter(r['evidence_category'] for r in reconciled if r['role']=='CANDIDATE')),priority_annotation_loci=len(priority),annotation_collection_issues=sum(r['status']!='COLLECTED' for r in status),priority_protein_provenance=dict(Counter(r['annotation_protein_status'] for r in consistency)),final_counts=False,new_phylogeny=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(prior_manifest_sha256=sha(src/'SHA256SUMS.txt'),mapping_sha256=sha(mp),input_manifest_sha256=sha(manifest),alignment_method='MAFFT --amino --localpair --maxiterate 1000',split_review_rule='Ordered, nonoverlapping protein segments; increasing model endpoints; summed adjacent model overlap <=20; union model coverage >=0.9. Review heuristic only.'),indent=2)+'\n')
(out/'README.txt').write_text('Domain reconciliation retains all profile segments. Ordered split matches with broad combined model coverage are separated from reduced-coverage cases, not automatically accepted as intact genes or assigned a subtype. All 53 candidates plus four controls are aligned de novo; no diagnostic CDS-derived fragments are included. Metrics are reported against each reference, not only the best reference, to expose reference-dependent gaps. Occupancy is not proof of structural completeness, function, or correct annotation. Gap coordinates are 1-based inclusive. GFF windows include neighbors, which are context and not automatically part of target loci. Exact mapping aliases only; no transcript suffix stripping. Missing or unmatched source sequences remain explicit. No genome translation, read validation, model repair, phylogeny, or final count promotion. Two CDS-derived historical records remain separately documented.\n')
p=out/'03_mafft.log';p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='PEPC_alignment_annotation_review')
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
