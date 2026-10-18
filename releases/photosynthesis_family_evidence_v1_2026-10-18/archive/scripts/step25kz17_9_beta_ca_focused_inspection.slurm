#!/usr/bin/env bash
#SBATCH --job-name=beta_ca_focus
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=beta_ca_focus_%j.out
#SBATCH --error=beta_ca_focus_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
: "${CA_RELEASE:?Export the Step25kz17.7 release directory}"
command -v python3 >/dev/null
command -v mafft >/dev/null
python3 - <<'PY'
import os,csv,json,gzip,hashlib,tarfile,subprocess,re
from pathlib import Path
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
source=Path(os.environ['CA_RELEASE']).resolve(strict=True)
targets=['SP_Hex__gene_Spor0140320.1','SP_Hex__gene_Spor0527280.1',
         'SLP_Hap1__SP3G00376','PO__evm.TU.LG15.1221','AT__gene-LOC130820426']
threads=int(os.environ.get('SLURM_CPUS_PER_TASK','4'))
inputs={n:source/n for n in ['05_all_screened_candidates.faa','08_beta_only_review_alignment.faa','03_domain_evidence.tsv','02_CA_family_specific_review.tsv']}
for p in inputs.values():
 if not p.is_file() or not p.stat().st_size:raise SystemExit('Required source file missing: '+p.name)
def opener(p):return gzip.open(p,'rt') if p.suffix=='.gz' else p.open()
def fasta_records(p):
 with opener(p) as f:
  key=None;parts=[]
  for line in f:
   if line.startswith('>'):
    if key is not None:yield key,''.join(parts).upper()
    key=line[1:].split()[0];parts=[]
   elif line.strip():
    if key is None:raise SystemExit('FASTA data precedes header')
    parts.append(''.join(line.split()))
  if key is not None:yield key,''.join(parts).upper()
def fasta(p):
 result={}
 for k,s in fasta_records(p):
  if k in result or not s:raise SystemExit('Duplicate or empty FASTA record')
  result[k]=s
 return result
def write_fasta(p,records):
 with p.open('w') as f:
  for k,s in records:
   f.write('>'+k+'\n')
   for i in range(0,len(s),80):f.write(s[i:i+80]+'\n')
def read_tsv(p):
 with p.open(encoding='utf-8-sig',newline='') as f:
  r=csv.DictReader(f,delimiter='\t');return r.fieldnames,list(r)
def table(p,fields,rows):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def safe(s):
 s=str(s).replace(str(root),'<PROJECT_ROOT>')
 home=str(Path.home())
 if home!='/':s=s.replace(home,'<HOME>')
 s=re.sub(r'(?<![\w.:/])/(?!/)[^\s\t\n\"\'<>|;,]+','<ABS_PATH>',s)
 user=os.environ.get('USER','')
 if user:s=re.sub(r'(?<!\w)'+re.escape(user)+r'(?!\w)','<USER>',s)
 return s
# Standard nuclear genetic code; frame is the supplied CDS frame, never selected by best match.
bases='TCAG';amino='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
code=dict(zip((a+b+c for a in bases for b in bases for c in bases),amino))
def translate(seq):return ''.join(code.get(seq[i:i+3],'X') for i in range(0,len(seq)-2,3))
proteins=fasta(inputs['05_all_screened_candidates.faa'])
alignment=fasta(inputs['08_beta_only_review_alignment.faa'])
refs={k:s.replace('-','') for k,s in alignment.items() if k.startswith('REF__')}
if len(refs)!=6:raise SystemExit('Expected six beta references')
_,domains=read_tsv(inputs['03_domain_evidence.tsv']);_,screen=read_tsv(inputs['02_CA_family_specific_review.tsv'])
for pid in targets:
 if pid not in proteins or pid not in alignment:raise SystemExit('Missing focused candidate: '+pid)
 if alignment[pid].replace('-','')!=proteins[pid]:raise SystemExit('Source sequence disagreement: '+pid)
stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
out=root/'08_focused_families/validation/step25kz17.9'/('release_'+os.environ.get('SLURM_JOB_ID','manual')+'_'+stamp)
out.mkdir(parents=True,exist_ok=False)
summary=[];mapping_status=[];context_status=[];cds_checks=[];domain_fragments=[];site_checks=[]
for pid in targets:
 species=pid.split('__',1)[0];d=out/pid;d.mkdir()
 seq=proteins[pid];write_fasta(d/'01_candidate_protein.faa',[(pid,seq)])
 ds=[x for x in domains if x['protein_id']==pid and x['family']=='BETA_CA']
 table(d/'02_beta_domain_evidence.tsv',list(domains[0]),ds)
 significant=sorted([x for x in ds if float(x['independent_evalue'])<=1e-5],key=lambda x:int(x['alignment_from']))
 motifs=[]
 for pattern,label in [('C.D.R','C_x_D_x_R'),('H..C','H_x_x_C')]:
  for m in re.finditer(pattern,seq):motifs.append({'motif_pattern':label,'start_aa_1based':m.start()+1,'end_aa_1based':m.end(),'sequence':m.group()})
 table(d/'03_motif_occurrences.tsv',['motif_pattern','start_aa_1based','end_aa_1based','sequence'],motifs)
 # Whole-candidate local alignment; domain copies are also aligned separately below.
 analysis_units=[('whole',seq)]
 for i,x in enumerate(significant,1):
  start,end=int(x['alignment_from']),int(x['alignment_to'])
  if not 1<=start<=end<=len(seq):raise SystemExit('Bad domain coordinates')
  unit='domain_'+str(i);analysis_units.append((unit,seq[start-1:end]))
  domain_fragments.append({'protein_id':pid,'unit':unit,'start_aa_1based':start,'end_aa_1based':end,
                           'length_aa':end-start+1,'independent_evalue':x['independent_evalue'],
                           'interpretation':'DIAGNOSTIC_FRAGMENT_NOT_A_GENE_OR_FORMAL_TREE_TIP'})
 for unit,useq in analysis_units:
  uid=pid+'__'+unit
  inf=d/(unit+'_with_references.faa');af=d/(unit+'_alignment.faa')
  write_fasta(inf,list(refs.items())+[(uid,useq)])
  with af.open('w') as f,(d/(unit+'_mafft.log')).open('w') as log:
   subprocess.run(['mafft','--thread',str(threads),'--localpair','--maxiterate','1000',str(inf)],stdout=f,stderr=log,check=True)
  aa=fasta(af)
  if set(aa)!=(set(refs)|{uid}) or len(set(map(len,aa.values())))!=1:raise SystemExit('Alignment QC failed')
  for k,s in aa.items():
   if s.replace('-','')!=(useq if k==uid else refs[k]):raise SystemExit('Alignment changed sequence')
  anchor=next(k for k in refs if '|P27140|' in k);n=0;mp={}
  for j,c in enumerate(aa[anchor]):
   if c!='-':n+=1;mp[n]=j
  pos=[167,169,171,227,230];cols=[mp[x] for x in pos]
  consensus_ok=all(''.join(aa[k][j] for j in cols)=='CDRHC' for k in refs)
  observed=''.join(aa[uid][j] for j in cols)
  site_checks.append({'protein_id':pid,'unit':unit,'length_aa':len(useq),'reference_mapping_consistent':str(consensus_ok),
                     'aligned_residues':observed,'status':'REFERENCE_MAPPING_REVIEW' if not consensus_ok else 'MATCH_ALL_FIVE' if observed=='CDRHC' else 'RESIDUE_REVIEW_REQUIRED'})
 # Exact identifier matching in harmonized mapping tables; no substring/fuzzy match.
 hdir=root/'03_harmonized/primary_v1'/species
 aliases={pid,pid.split('__',1)[1]}
 mapping_tables=[]
 for name in [species+'.representative_gene_map.tsv',species+'.all_isoform_mapping.tsv']:
  p=hdir/name
  if not p.is_file():mapping_status.append({'protein_id':pid,'file':name,'status':'MISSING_FILE','matched_rows':0});continue
  fields,rr=read_tsv(p)
  matched=[r for r in rr if any(v in aliases for v in r.values() if isinstance(v,str))]
  mapping_tables.append((name,fields,rr,matched))
  # Expand only identifier-labelled fields from an exactly matched row.
  for r in matched:
   for k,v in r.items():
    if v and (k.lower().endswith('_id') or k.lower() in {'gene','transcript','protein','gene_id','transcript_id','protein_id'}) and not any(z in k.lower() for z in ['chrom','scaffold','contig','species']):aliases.add(v)
 for name,fields,rr,matched in mapping_tables:
  matched=[r for r in rr if any(v in aliases for v in r.values() if isinstance(v,str))]
  table(d/name,fields,[{k:safe(v) for k,v in r.items()} for r in matched])
  mapping_status.append({'protein_id':pid,'file':name,'status':'EXACT_MATCH' if matched else 'NO_EXACT_MATCH','matched_rows':len(matched)})
 table(d/'04_identifier_aliases.tsv',['identifier'],[{'identifier':v} for v in sorted(aliases)])
 # Extract candidate CDS and representative protein records without loading full genomes.
 for suffix,label in [('representative.cds.fa.gz','cds'),('representative.proteins.fa.gz','harmonized_protein')]:
  p=hdir/(species+'.'+suffix)
  matches=[]
  if p.is_file():
   for k,s in fasta_records(p):
    if k in aliases:matches.append((k,s))
  write_fasta(d/(label+'_exact_matches.fa'),matches)
  context_status.append({'protein_id':pid,'evidence':label,'status':'MISSING_FILE' if not p.is_file() else 'NO_EXACT_MATCH' if not matches else 'MATCHED' if len(matches)==1 else 'MULTIPLE_MATCHES_REVIEW','records':len(matches)})
  for k,s in matches:
   if label=='cds':
    translated=translate(s);trimmed=translated[:-1] if translated.endswith('*') else translated
    cds_checks.append({'protein_id':pid,'cds_id':k,'cds_nt':len(s),'length_modulo_3':len(s)%3,'start_codon':s[:3],
                       'terminal_stop':str(translated.endswith('*')),'internal_stop_count':trimmed.count('*'),
                       'translation_matches_screened_protein':str(trimmed==seq.rstrip('*')),
                       'ambiguous_codon_count':translated.count('X'),
                       'interpretation':'CDS_PROTEIN_CONSISTENCY_ONLY_NOT_GENOMIC_VALIDATION'})
    write_fasta(d/('translation_'+str(len(cds_checks))+'.faa'),[(k,translated)])
   else:
    context_status.append({'protein_id':pid,'evidence':'harmonized_protein_matches_screened_sequence',
                           'status':'IDENTICAL' if s.rstrip('*')==seq.rstrip('*') else 'DIFFERENT_INPUT_VERSION_OR_MODEL_REVIEW','records':1})
 # BED provides genomic placement and neighboring gene intervals, not exon structure.
 bed=hdir/(species+'.representative.genes.bed');bedrows=[]
 if bed.is_file():
  for line in bed.read_text().splitlines():
   if not line or line.startswith(('#','track','browser')):continue
   t=line.split('\t')
   if len(t)<4:continue
   try:start,end=int(t[1]),int(t[2])
   except ValueError:continue
   bedrows.append((t[0],start,end,t))
 hits=[i for i,b in enumerate(bedrows) if b[3][3] in aliases]
 selected=[]
 for i in hits:
  b=bedrows[i];chrom=sorted([v for v in bedrows if v[0]==b[0]],key=lambda v:(v[1],v[2],v[3][3]))
  idx=chrom.index(b)
  selected.extend(chrom[max(0,idx-2):idx+3])
 with (d/'05_genomic_neighborhood.bed').open('w') as f:
  for t in sorted({tuple(x[3]) for x in selected},key=lambda v:(v[0],int(v[1]),v[3])):f.write('\t'.join(t)+'\n')
 context_status.append({'protein_id':pid,'evidence':'BED_gene_context','status':'MISSING_FILE' if not bed.is_file() else 'MATCHED' if hits else 'NO_EXACT_MATCH','records':len(hits)})
 summary.append({'protein_id':pid,'length_aa':len(seq),'significant_beta_domains':len(significant),
                 'C_x_D_x_R_occurrences':sum(x['motif_pattern']=='C_x_D_x_R' for x in motifs),
                 'H_x_x_C_occurrences':sum(x['motif_pattern']=='H_x_x_C' for x in motifs),
                 'manual_decision':'PENDING','interpretation':'FOCUSED_DIAGNOSTICS_ONLY'})
table(out/'01_focused_summary.tsv',list(summary[0]),summary)
table(out/'02_local_alignment_site_checks.tsv',list(site_checks[0]),site_checks)
table(out/'03_domain_fragment_coordinates.tsv',['protein_id','unit','start_aa_1based','end_aa_1based','length_aa','independent_evalue','interpretation'],domain_fragments)
table(out/'04_CDS_translation_checks.tsv',['protein_id','cds_id','cds_nt','length_modulo_3','start_codon','terminal_stop','internal_stop_count','translation_matches_screened_protein','ambiguous_codon_count','interpretation'],cds_checks)
table(out/'05_mapping_status.tsv',['protein_id','file','status','matched_rows'],mapping_status)
table(out/'06_context_status.tsv',['protein_id','evidence','status','records'],context_status)
notes={'status':'COMPLETE_FOCUSED_INSPECTION_REVIEW_ONLY','target_count':len(targets),
 'input_sha256':{k:sha(p) for k,p in inputs.items()},
 'notes':['Local alignments are diagnostic sensitivity checks. They do not overwrite the accepted input alignment.',
 'Domain fragments are diagnostic units, not new genes or automatic phylogenetic tips. Overlapping domain boundaries are preserved and recorded.',
 'CDS translation uses the standard nuclear code and the supplied frame; partial trailing codons are reported through length_modulo_3.',
 'An exact CDS/protein match validates their consistency, not the biological correctness of the annotation.',
 'BED coordinates are zero-based, half-open and do not provide exon structure or read support.',
 'Missing/ambiguous mapping records are reported explicitly. No fuzzy matching or gene-model repair is performed.',
 'No automatic acceptance, exclusion, splitting of gene models or tree construction is performed.',
 'Genomic sequence, exon-level GFF and read support may still be required to resolve fusion or truncation.']}
(out/'FOCUSED_QC.json').write_text(json.dumps(notes,indent=2)+'\n')
(out/'COMPLETE.txt').write_text('Focused inspection complete; review required.\n')
# Sanitize paths in diagnostic tool logs before sharing.
for p in out.rglob('*.log'):p.write_text(safe(p.read_text()))
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+p.relative_to(out).as_posix()+'\n')
archive=out/'BETA_CA_focused_inspection_bundle.tar.gz'
with tarfile.open(archive,'w:gz') as t:
 for p in sorted(out.iterdir()):
  if p!=archive:t.add(p,arcname='BETA_CA_focused_inspection/'+p.name)
print('Focused inspection completed for',len(targets),'candidates.')
print('Output: <PROJECT_ROOT>/'+out.relative_to(root).as_posix())
print('Upload BETA_CA_focused_inspection_bundle.tar.gz')
PY
