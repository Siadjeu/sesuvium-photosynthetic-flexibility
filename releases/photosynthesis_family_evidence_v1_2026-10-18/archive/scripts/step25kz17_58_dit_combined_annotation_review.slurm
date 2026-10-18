#!/usr/bin/env bash
#SBATCH --job-name=dit_annotation
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=dit_annotation_%j.out
#SBATCH --error=dit_annotation_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
python3 - <<'PY'
import urllib.request,time,csv,gzip,hashlib,json,os,re,shutil,subprocess,tarfile
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True);base=root/'08_focused_families/validation'
def need(v,msg):
 if not v:raise SystemExit(msg)
hash_cache={}
def sha(p):
 p=p.resolve()
 if p in hash_cache:return hash_cache[p]
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 hash_cache[p]=h.hexdigest();return hash_cache[p]
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
def release(step,env):
 pp=[root/os.environ[env]] if os.environ.get(env) else [p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
 need(len(pp)==1,'Set '+env+' to the intended release');p=pp[0].resolve(strict=True);verify(p);return p



src=release('step25kz17.57','DIT_DOMAIN_RELEASE');rec=release('step25kz17.56','DIT_RECONCILIATION_RELEASE');screen=release('step25kz17.46','TRANSPORTER_RECOVERY_RELEASE')
need(json.loads((src/'PROVENANCE.json').read_text())['source_manifest_sha256']==sha(rec/'SHA256SUMS.txt'),'Domain/reconciliation provenance mismatch');need(json.loads((rec/'PROVENANCE.json').read_text())['step46_manifest_sha256']==sha(screen/'SHA256SUMS.txt'),'Screen provenance mismatch')
ledger=read(src/'08_DIT_completeness_review.tsv');queue=read(src/'09_exception_queue.tsv');exceptions={r['protein_id'] for r in queue};seqs,_=fasta(src/'01_domain_queries.faa');need(len(ledger)==38 and len(exceptions)==9,'Unexpected review scope')
mapping_hashes=json.loads((screen/'PROVENANCE.json').read_text())['mapping_hashes']
out=base/'step25kz17.58'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'annotation').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
table('01_all_candidate_domain_review.tsv',ledger,list(ledger[0]));table('02_nine_exception_queue.tsv',queue,list(queue[0]));writefa('03_exception_proteins.faa',{k:seqs[k] for k in sorted(exceptions)})
manifest=root/'00_admin/manifests/input_sets.tsv';configs=read(manifest);sources=[];status=[];mappingrows=[];checks=[]
bases='TCAG';letters='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';codons={a+b+c:v for (a,b,c),v in zip(((a,b,c) for a in bases for b in bases for c in bases),letters)}
def op(p):return gzip.open(p,'rt') if p.suffix=='.gz' else p.open()
def collect(p,aliases):
 found={};key=None;parts=[]
 def commit():
  if key in aliases:need(key not in found,'Duplicate source ID');found[key]=''.join(parts).upper()
 with op(p) as f:
  for line in f:
   if line.startswith('>'):commit();key=line[1:].split()[0];parts=[]
   elif key in aliases:parts.append(''.join(line.split()))
  commit()
 return found
for idx,k in enumerate(sorted(exceptions),1):
 r=next(x for x in ledger if x['protein_id']==k);sp=r['species'];gene=r['gene_id'];folder='annotation/E%02d'%idx;(out/folder).mkdir();mp=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv')
 need(mp.is_file(),'Missing mapping for '+sp);need(sha(mp)==mapping_hashes[sp],'Mapping provenance mismatch for '+sp);mm=[m for m in read(mp) if m['gene_id']==gene];need(mm,'No exact gene mapping '+k);need(len({(m['seqid'],m['strand']) for m in mm})==1,'Conflicting mapping coordinates');table(folder+'/mapping.tsv',mm,list(mm[0]));sources.append(dict(protein_id=k,kind='mapping',file=str(mp.relative_to(root)),sha256=sha(mp)));aliases={v for m in mm for f in ['gene_id','transcript_id','protein_id','cds_id'] if (v:=m.get(f,'')) and v not in {'.','NA'}};aliases|={sp+'__'+v for v in list(aliases)};need(k in aliases,'Candidate alias mismatch')
 conf=[c for c in configs if c['input_set']==sp];need(len(conf)==1,'Ambiguous manifest species');conf=conf[0];found={}
 for kind,col in [('protein','protein_all'),('cds','cds_all'),('gff','gff3')]:
  value=conf.get(col,'');p=root/value if value else root/'__MISSING_INPUT__';exists=p.is_file();sources.append(dict(protein_id=k,kind=kind,file=str(p.relative_to(root)) if p.is_relative_to(root) else '<EXTERNAL>/'+p.name,sha256=sha(p) if exists else 'MISSING'))
  if not exists:status.append(dict(protein_id=k,kind=kind,status='MISSING',records=0));continue
  if kind!='gff':found[kind]=collect(p,aliases);writefa(folder+'/'+kind+'_isoforms.fasta',found[kind]);n=len(found[kind])
  else:
   chrom=mm[0]['seqid'];lo=max(1,min(int(m['start']) for m in mm)-5000);hi=max(int(m['end']) for m in mm)+5000;n=0
   with op(p) as f,(out/folder/'context.gff3').open('w') as w:
    w.write('##gff-version 3\n')
    for line in f:
     if line.startswith('##FASTA'):break
     if line.startswith('#'):continue
     c=line.rstrip('\n').split('\t')
     if len(c)==9 and c[0]==chrom and int(c[3])<=hi and int(c[4])>=lo:w.write(line.rstrip('\n')+'\n');n+=1
  status.append(dict(protein_id=k,kind=kind,status='COLLECTED' if n else 'NO_EXACT_ID_MATCH',records=n))
 for m in mm:
  def lookup(kind,field):
   ident=m.get(field,'');d=found.get(kind,{});return d.get(ident,d.get(sp+'__'+ident))
  aa=lookup('protein','protein_id');dna=lookup('cds','cds_id');translation=''.join(codons.get(dna[i:i+3],'X') for i in range(0,len(dna)-len(dna)%3,3)) if dna else None
  checks.append(dict(candidate_id=k,mapped_transcript=m.get('transcript_id',''),protein_id=m.get('protein_id',''),cds_id=m.get('cds_id',''),source_protein_found=aa is not None,source_CDS_found=dna is not None,protein_matches_candidate=norm(aa)==seqs[k] if aa else 'NA',cds_nt=len(dna) if dna else 'NA',frame_remainder=len(dna)%3 if dna else 'NA',translation_matches_source_protein=norm(translation)==norm(aa) if translation and aa else 'NA',translation_matches_candidate=norm(translation)==seqs[k] if translation else 'NA',internal_stops=norm(translation).count('*') if translation else 'NA',ambiguous_amino_acids=translation.count('X') if translation else 'NA',start_codon=dna[:3] if dna else 'NA',terminal_codon=dna[-3:] if dna else 'NA',interpretation='ANNOTATION_CONSISTENCY_NOT_GENOME_OR_READ_VALIDATION'))
 mappingrows.append(dict(protein_id=k,species=sp,gene_id=gene,folder=folder))
table('06_annotation_file_map.tsv',mappingrows,list(mappingrows[0]));table('07_annotation_collection_status.tsv',status,list(status[0]));table('08_protein_CDS_consistency.tsv',checks,list(checks[0]));table('09_annotation_source_provenance.tsv',sources,list(sources[0]))
summary=[]
for r in queue:
 k=r['protein_id'];cc=[c for c in checks if c['candidate_id']==k];exact=[c for c in cc if c['protein_matches_candidate'] is True];paired=[c for c in exact if c['translation_matches_candidate'] is True and c['frame_remainder']==0 and c['internal_stops']==0 and c['ambiguous_amino_acids']==0]
 category='SOURCE_PROTEIN_UNRESOLVED' if not exact else 'CDS_CONSISTENCY_UNRESOLVED'
 if paired:category='ANNOTATION_CONSISTENT_MODEL_CAVEAT_REMAINS'
 if any(c['source_CDS_found'] is True and c['translation_matches_source_protein'] is False for c in exact):category='PROTEIN_CDS_DISCORDANCE_REVIEW'
 summary.append(dict(protein_id=k,species=r['species'],gene_id=r['gene_id'],tree_group=r['tree_group'],length_aa=r['length_aa'],domain_completeness_flags=r['flags'],mapped_isoforms=len(cc),exact_source_protein_isoforms=len(exact),consistent_candidate_CDS_isoforms=len(paired),category=category,alternative_isoforms=';'.join(c['protein_id'] for c in cc if c['source_protein_found'] is True and c['protein_matches_candidate'] is False),decision='REVIEW_ONLY',limitation='CDS_CONSISTENCY_DOES_NOT_RESOLVE_TRUNCATION_OR_FUSION',final_count='NA'))
table('10_combined_exception_summary.tsv',summary,list(summary[0]))
for n in ['source_01_existing_tree_DIT_group_audit.tsv','source_02_DIT_candidate_reconciliation.tsv','source_07_refreshed_family_coverage.tsv','source_prior_01_integrated_locus_evidence.tsv','source_prior_07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/n)
shutil.copyfile(src/'06_all_Pfam_domain_hits.tsv',out/'11_existing_domain_hits.tsv');shutil.copyfile(src/'07_domain_model_coverage.tsv',out/'12_existing_model_coverage.tsv')
(out/'README.txt').write_text('Combined bounded review of nine DIT exceptions: verifies release chains and original mapping hashes; collects exact mapped gene isoforms, proteins, CDS and GFF context +/-5kb; translates annotated CDS in supplied coding orientation with standard nuclear code; compares every mapped protein/CDS pair and candidate; reports missing inputs, disagreements, alternative isoforms and retained architecture caveats. Terminal stop stripped only for protein comparison; incomplete codons, internal stops, ambiguous residues and start/terminal codons explicitly reported. No frame hunting or sequence repair. GFF includes neighbors for context only, with no inferred merging. CDS consistency is not independent biological validation and does not establish intact transporters. Original domain/tree evidence and ledger retained. No searches, alignments, phylogenies, automatic count promotion or open-ended genome rescue. Begin review with10_combined_exception_summary.tsv and08_protein_CDS_consistency.tsv.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step57_manifest_sha256=sha(src/'SHA256SUMS.txt'),step56_manifest_sha256=sha(rec/'SHA256SUMS.txt'),step46_manifest_sha256=sha(screen/'SHA256SUMS.txt'),input_manifest_sha256=sha(manifest)),indent=2)+'\n')
qc=dict(status='DIT_COMBINED_ANNOTATION_REVIEW_READY',candidate_records=38,exception_loci=len(exceptions),summary_categories=dict(Counter(r['category'] for r in summary)),collection_status_counts=dict(Counter(r['status'] for r in status)),final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out);bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='DIT_combined_annotation_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
