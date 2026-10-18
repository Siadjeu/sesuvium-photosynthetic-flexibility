#!/usr/bin/env bash
#SBATCH --job-name=vatp_cds_qc
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=vatp_cds_qc_%j.out
#SBATCH --error=vatp_cds_qc_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
for tool in python3; do command -v "$tool"; done
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
from urllib.parse import unquote
src=release('step25kz17.34','VATPASE_REVIEW_RELEASE')
ex=read(src/'02_exception_queue.tsv');assert len(ex)==23
# Standard nuclear genetic code; ambiguous codons are explicitly translated as X.
bases='TCAG';aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codons=[a+b+c for a in bases for b in bases for c in bases];code=dict(zip(codons,aa))
def translate(s):
 s=s.upper().replace('U','T');return ''.join(code.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
def normalize(s):
 s=s.upper();return s[:-1] if s.endswith('*') else s
assert translate('ATGGCTTAA')=='MA*' and normalize('MA**')=='MA*'
def records(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:
  h=None;seq=[]
  for l in f:
   if l.startswith('>'):
    if h is not None:yield h,''.join(seq)
    h=l[1:].strip();seq=[]
   elif l.strip():seq.append(l.strip())
  if h is not None:yield h,''.join(seq)
def attributes(s):
 d={}
 for item in s.split(';'):
  if '=' in item:
   k,v=item.split('=',1);d[k]=[unquote(x) for x in v.split(',')]
 return d
out=base/'step25kz17.35'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False)
def table(name,rr,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
# Reference proteins are the exact sequences used in the preceding subunit alignments.
original={}
for p in (src/'alignments').glob('*_input.faa'):
 for k,v in fasta(p).items():
  if k in original:assert original[k]==v
  original[k]=v
source_prov=[];recovery=[];checks=[];proteins={};cdsout={};translated={}
manifest=root/'00_admin/manifests/input_sets.tsv'
for sp in sorted({r['species'] for r in ex}):
 folder=src/'annotation'/sp;maps=read(folder/'01_all_isoform_mapping.tsv')
 pf=fasta(folder/'03_protein_isoforms.fasta');cf=fasta(folder/'03_cds_isoforms.fasta')
 for r in [r for r in ex if r['species']==sp]:
  pid=r['protein_id'];mm=[m for m in maps if m['gene_id']==r['gene_id']];assert len(mm)==1,'Expected one mapped isoform per exception'
  m=mm[0];aliases={m.get(k,'') for k in ['gene_id','transcript_id','protein_id','cds_id']}-{'','NA','.'}
  aliases|={sp+'__'+x for x in list(aliases)}
  ph=[(k,s) for k,s in pf.items() if k in aliases];ch=[(k,s) for k,s in cf.items() if k in aliases]
  if len(ph)!=1:raise SystemExit('Expected one source protein for '+pid+'; found '+str(len(ph)))
  source_protein=ph[0][1];prot=original[pid]
  if normalize(source_protein)!=normalize(prot):
   raise SystemExit('Protein sequence mismatch after case/one-terminal-stop normalization: '+pid+'; source_id='+ph[0][0]+'; source_length='+str(len(source_protein))+'; screened_length='+str(len(prot)))
  protein_relation='EXACT_BYTES' if source_protein==prot else 'MATCH_AFTER_CASE_AND_ONE_TERMINAL_STOP_NORMALIZATION'
  proteins[pid]=prot;method='COLLECTED_EXACT_ALIAS';hits=ch
  if not hits and sp=='AT':
   # Derive additional identifiers only from the mapped target transcript/CDS features.
   gff=folder/'02_annotation_windows.gff3';documented=set(aliases)
   for line in gff.read_text().splitlines():
    if line.startswith('#'):continue
    c=line.split('\t')
    if len(c)!=9 or c[2] not in ['CDS','mRNA','transcript']:continue
    a=attributes(c[8]);links=set(a.get('Parent',[])+a.get('ID',[])+a.get('protein_id',[])+a.get('transcript_id',[]))
    if links & aliases:
     for key in ['ID','Name','protein_id','transcript_id']:documented.update(a.get(key,[]))
   configs=[x for x in read(manifest) if x['input_set']=='AT'];assert len(configs)==1
   source=root/configs[0]['cds_all']
   expected=[x for x in read(src/'09_annotation_provenance.tsv') if x['species']=='AT' and x['kind']=='cds'];assert len(expected)==1
   if source.is_file():
    digest=sha(source);assert digest==expected[0]['sha256'],'AT CDS source changed'
    source_prov.append(dict(species=sp,kind='AT_original_CDS',sha256=digest,source=configs[0]['cds_all']))
    hits=[]
    for header,seq in records(source):
     first=header.split()[0];fields={unquote(v) for k,v in re.findall(r'\[(protein_id|transcript_id|ID|Name)=([^\]]+)\]',header)}
     # NCBI bracket metadata are matched exactly; no substring/suffix guessing.
     matched=({first}|fields)&documented
     if matched:
      hits.append((first,seq));recovery.append(dict(protein_id=pid,source_header=header,matched_identifiers=';'.join(sorted(matched)),method='DOCUMENTED_GFF_AND_HEADER_IDENTIFIERS'))
    method='DOCUMENTED_GFF_AND_HEADER_IDENTIFIERS'
   else:method='AT_SOURCE_MISSING'
  base_row=dict(protein_id=pid,species=sp,gene_id=r['gene_id'],subunit=r['provisional_subunit'],protein_length=len(prot),source_protein_id=ph[0][0],source_protein_length=len(source_protein),source_protein_terminal_stop='YES' if source_protein.endswith('*') else 'NO',source_vs_screened_protein=protein_relation,cds_records=len(hits),cds_source_ids=';'.join(k for k,s in hits),recovery_method=method)
  if len(hits)!=1:
   checks.append(dict(**base_row,cds_length='',length_mod3='',start_codon='',terminal_stop='',internal_stops='',ambiguous_codons='',translated_length='',mismatch_positions='',status='CDS_MISSING' if not hits else 'AMBIGUOUS_CDS_MATCHES'));continue
  cid,dna=hits[0];dna=dna.upper().replace('U','T');tr=translate(dna);pn=normalize(prot);tn=normalize(tr)
  invalid=sorted(set(dna)-set('ACGTRYSWKMBDHVN'))
  stops=tn.count('*');unknown=tr.count('X');mismatches=[str(i+1) for i in range(max(len(tn),len(pn))) if (tn[i] if i<len(tn) else '')!=(pn[i] if i<len(pn) else '')]
  flags=[]
  if invalid:flags.append('INVALID_NUCLEOTIDE_CHARACTERS')
  if len(dna)%3:flags.append('CDS_LENGTH_NOT_MULTIPLE_OF_THREE')
  if stops:flags.append('INTERNAL_STOP')
  if unknown:flags.append('AMBIGUOUS_CODONS')
  if tn!=pn:flags.append('TRANSLATION_MISMATCH')
  status=';'.join(flags) if flags else ('EXACT_MATCH_TERMINAL_STOP_PRESENT' if tr.endswith('*') else 'EXACT_MATCH_NO_TERMINAL_STOP')
  checks.append(dict(**base_row,cds_length=len(dna),length_mod3=len(dna)%3,start_codon=dna[:3],terminal_stop='YES' if tr.endswith('*') else 'NO',internal_stops=stops,ambiguous_codons=unknown,translated_length=len(tn),mismatch_positions=','.join(mismatches),status=status))
  cdsout[pid]=dna;translated[pid]=tr
assert len(checks)==23
fields=['protein_id','species','gene_id','subunit','protein_length','source_protein_id','source_protein_length','source_protein_terminal_stop','source_vs_screened_protein','cds_records','cds_source_ids','recovery_method','cds_length','length_mod3','start_codon','terminal_stop','internal_stops','ambiguous_codons','translated_length','mismatch_positions','status']
table('01_protein_CDS_consistency.tsv',checks,fields)
table('02_AT_CDS_recovery_evidence.tsv',recovery,['protein_id','source_header','matched_identifiers','method'])
writefa(out/'03_exception_proteins.faa',proteins);writefa(out/'04_recovered_CDS.fna',cdsout);writefa(out/'05_CDS_translations.faa',translated)
table('06_additional_source_provenance.tsv',source_prov,['species','kind','sha256','source'])
for name in ['01_reconciled_evidence_ledger.tsv','02_exception_queue.tsv','05_exception_alignment_completeness.tsv','06_exception_and_reference_domains.tsv','07_exception_locus_summary.tsv']:
 shutil.copy2(src/name,out/('prior_'+name))
joined=[]
for r in ex:
 c=next(c for c in checks if c['protein_id']==r['protein_id'])
 joined.append(dict(protein_id=r['protein_id'],species=r['species'],gene_id=r['gene_id'],subunit=r['provisional_subunit'],original_flags=r['flags'],protein_CDS_status=c['status'],family_decision='PENDING_EVIDENCE_REVIEW',model_limitation='ORIGINAL_ALIGNMENT_DOMAIN_CAVEATS_RETAINED',independent_read_support='NOT_ASSESSED'))
table('07_exception_review_ledger.tsv',joined,list(joined[0]))
(out/'QC.json').write_text(json.dumps(dict(status='BOUNDED_PROTEIN_CDS_CHECK_COMPLETE',exceptions=23,status_counts=dict(Counter(r['status'] for r in checks)),AT_recovery_records=len(recovery),genome_translation_checked=False,models_modified=False,final_counts=False),indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),input_manifest_sha256=sha(manifest) if manifest.is_file() else 'NOT_AVAILABLE'),indent=2)+'\n')
(out/'README.txt').write_text('Bounded consistency check of 23 V-ATPase exceptions. Standard nuclear genetic code, supplied CDS orientation and reading frame; no frame search or sequence repair. Source proteins are compared to screened proteins after case normalization and removal of at most one terminal stop; the source length and normalization status are recorded. Screened proteins are preserved unchanged in the output. One terminal stop is allowed on either protein or translation; internal stops and incomplete codons remain flagged. Ambiguous codons are X and explicitly flagged even if matching protein X.\nAT CDS recovery uses exact IDs from mapped transcript/CDS GFF attributes and original FASTA record IDs or bracketed protein_id/transcript_id/ID/Name fields. Multiple matching records remain unresolved, even if their sequences are identical. Full matched headers and original source checksum are retained.\nAn exact protein-CDS match shows internal consistency of annotation products, not independent support for a correct gene model, fusion, completeness, enzyme activity or C4/CAM function. No genome translation, new alignment, pseudogene call or gene-model repair is performed. Existing domain/alignment caveats and distinct loci are preserved.\n')
files=sorted(p for p in out.iterdir() if p.is_file())
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+p.name+'\n' for p in files))
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='VATPASE_protein_CDS_check')
print('RELEASE='+str(out));print('BUNDLE='+str(bundle));print(json.dumps(Counter(r['status'] for r in checks)))
PY
