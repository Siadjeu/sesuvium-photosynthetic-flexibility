#!/usr/bin/env bash
#SBATCH --job-name=carrier_review
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=carrier_review_%j.out
#SBATCH --error=carrier_review_%j.err
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
def release(step,env):
 pp=[root/os.environ[env]] if os.environ.get(env) else [p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
 need(len(pp)==1,'Set '+env+' to the intended release');p=pp[0].resolve(strict=True);verify(p);return p



src=release('step25kz17.53','CARRIER_TREE_RELEASE');align=release('step25kz17.52','CARRIER_ALIGNMENT_RELEASE');need(json.loads((src/'PROVENANCE.json').read_text())['source_manifest_sha256']==sha(align/'SHA256SUMS.txt'),'Source chain mismatch')
meta=read(src/'01_extraction_and_tree_id_map.tsv');byid={r['tree_id']:r for r in meta};seqs,_=fasta(align/'01_all_sequences.faa');refids={r['tree_id'] for r in meta if r['role']=='REFERENCE'}
# Read edge splits from unrooted Newick, including both sides of every edge.
def edges(path):
 tokens=re.findall(r'[^\s(),:;]+|[(),:;]',path.read_text());pos=0;ee=[]
 def node():
  nonlocal pos
  label=''
  if tokens[pos]=='(':
   pos+=1;leaves=set()
   while True:
    leaves|=node()
    if tokens[pos]==',':pos+=1;continue
    need(tokens[pos]==')','Malformed Newick');pos+=1;break
   if tokens[pos] not in {':',',',')',';'}:label=tokens[pos];pos+=1
  else:leaves={tokens[pos]};pos+=1
  if tokens[pos]==':':pos+=2
  ee.append((set(leaves),label));return leaves
 universe=node();need(tokens[pos]==';','Invalid tree');return universe,[(z,l) for a,l in ee for z in (a,universe-a) if 1<len(z)<len(universe)-1]
fams={'DTC':{'REF__Q9C5M0'},'DIC':{'REF__Q9SJY5','REF__Q9SB52','REF__Q9FY68'}};groups={};splitrows=[]
for mode in ['main','sensitivity_MC_partial','sensitivity_SS_partial']:
 universe,es=edges(src/mode/'tree.treefile')
 for fam,anchors in fams.items():
  options=[(z,l) for z,l in es if anchors<=z and z&refids==anchors];need(options,'No reference-exclusive group for '+fam);z,label=max(options,key=lambda x:len(x[0]));support=[float(x) for x in label.split('/')];need(len(support)==2 and support[0]>=80 and support[1]>=95,'Group support below review gate');groups[mode,fam]=z-anchors
  splitrows.append(dict(analysis=mode,family=fam,candidate_count=len(z-anchors),support_SHaLRT_UFBoot=label,tree_ids=';'.join(sorted(z-anchors))))
need(len(groups['main','DTC'])==42 and len(groups['main','DIC'])==17,'Unexpected main group sizes')
for fam in fams:
 for mode in ['sensitivity_MC_partial','sensitivity_SS_partial']:need(groups['main',fam]<=groups[mode,fam],'Main membership unstable')
exceptions={'SP_Hex__gene_Spor0293270.1','PO__evm.TU.LG05.272','PO__evm.TU.LG25.1042','SS__g18155.t1'}
out=base/'step25kz17.54'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'annotation').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
table('01_phylogenetic_group_support.tsv',splitrows,list(splitrows[0]));ledger=[]
for r in meta:
 if r['role']!='CANDIDATE':continue
 k=r['tree_id'];fam=next((f for f in fams if k in groups['main',f]),None);status='STABLE_PHYLOGENETIC_ASSIGNMENT' if fam else 'OUTSIDE_REFERENCE_DEFINED_DTC_DIC_GROUPS'
 if r['protein_id'] in exceptions:status='ARCHITECTURE_REVIEW_PENDING' if fam else 'PARTIAL_DTC_ASSOCIATION_REVIEW_ONLY'
 if r['protein_id']=='SS__g18155.t1':need(k in groups['sensitivity_SS_partial','DTC'],'SS partial placement changed');fam='DTC'
 ledger.append(dict(protein_id=r['protein_id'],species=r['species'],gene_id=r['gene_id'],tree_id=k,phylogenetic_family=fam or 'UNASSIGNED',evidence_category=status,main_tree=k in groups['main','DTC']|groups['main','DIC'],limitation='BOUNDED_PANEL;FUNCTION_AND_EXHAUSTIVE_COUNTS_NOT_ESTABLISHED',final_count='NA'))
table('02_candidate_phylogenetic_evidence.tsv',ledger,list(ledger[0]));need(len(ledger)==85,'Expected 85 records')
hits=[r for r in read(align/'03_domains_in_alignment.tsv') if r['sequence_id'] in exceptions or r['sequence_id'] in set.union(*fams.values())];table('03_exception_and_anchor_domain_hits.tsv',hits,list(hits[0]));writefa('04_exception_and_anchor_proteins.faa',{k:seqs[k] for k in sorted(exceptions|set.union(*fams.values()))})
# Quantify repeated domain overlap rather than interpreting raw hit count as repeat count.
intervals=[]
for k in sorted(exceptions|set.union(*fams.values())):
 hh=sorted([h for h in hits if h['sequence_id']==k and h['accession'].split('.')[0]=='PF00153'],key=lambda h:int(h['ali_from']));covered=set()
 for h in hh:covered.update(range(int(h['ali_from']),int(h['ali_to'])+1))
 overlaps=[]
 for i,h in enumerate(hh):
  for g in hh[i+1:]:
   n=max(0,min(int(h['ali_to']),int(g['ali_to']))-max(int(h['ali_from']),int(g['ali_from']))+1)
   if n:overlaps.append(str(n))
 intervals.append(dict(protein_id=k,length_aa=len(seqs[k]),hit_count=len(hh),intervals=';'.join(h['ali_from']+'-'+h['ali_to'] for h in hh),union_residues=len(covered),pairwise_overlap_lengths=';'.join(overlaps) or 'NONE',interpretation='REVIEW_HIT_LAYOUT_NOT_AUTOMATIC_MODEL_DEFECT'))
table('05_repeat_overlap_review.tsv',intervals,list(intervals[0]))
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
 need(mp.is_file(),'Missing mapping for '+sp);mm=[m for m in read(mp) if m['gene_id']==gene];need(mm,'No exact gene mapping '+k);need(len({(m['seqid'],m['strand']) for m in mm})==1,'Conflicting mapping coordinates');table(folder+'/mapping.tsv',mm,list(mm[0]));sources.append(dict(protein_id=k,kind='mapping',file=str(mp.relative_to(root)),sha256=sha(mp)));aliases={v for m in mm for f in ['gene_id','transcript_id','protein_id','cds_id'] if (v:=m.get(f,'')) and v not in {'.','NA'}};aliases|={sp+'__'+v for v in list(aliases)};need(k in aliases,'Candidate alias mismatch')
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
for n in ['source_prior_01_integrated_locus_evidence.tsv','source_prior_06_refreshed_family_coverage.tsv','source_prior_07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/n)
(out/'README.txt').write_text('Stable unrooted reference-exclusive edge groups recorded from three existing trees; minimum SH-aLRT80/UFBoot95. Main memberships verified preserved. 42 DTC/17 DIC phylogenetic placements; architecture review pending for three main candidates, SS partial association separate. Other prioritized candidates remain outside these reference-defined groups, not automatically excluded from every carrier family. Four exception annotations collected through exact gene/isoform mappings; all isoforms retained and paired by documented protein/CDS identifiers. GFF +/-5kb includes neighboring features for context only; no locus merging. Missing files/IDs explicit. CDS translation checks do not validate biological gene models or prove functional transport. Domain hits/overlaps retained without automatic repeat-count decisions. No searches, trees, repairs, integrated-ledger updates or final counts.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step53_manifest_sha256=sha(src/'SHA256SUMS.txt'),step52_manifest_sha256=sha(align/'SHA256SUMS.txt'),input_manifest_sha256=sha(manifest)),indent=2)+'\n');qc=dict(status='CARRIER_ASSIGNMENTS_EXCEPTION_REVIEW_READY',candidate_records=len(ledger),main_DTC=42,main_DIC=17,exception_loci=4,annotation_status=status,final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out);bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Carrier_assignments_exception_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
