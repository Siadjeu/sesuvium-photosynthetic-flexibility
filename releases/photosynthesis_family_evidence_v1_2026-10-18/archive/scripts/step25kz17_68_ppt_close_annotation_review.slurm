#!/usr/bin/env bash
#SBATCH --job-name=ppt_close_review
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=ppt_close_review_%j.out
#SBATCH --error=ppt_close_review_%j.err
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







def edges(path):
 t=re.findall(r'[^\s(),:;]+|[(),:;]',path.read_text());pos=0;ee=[]
 def node():
  nonlocal pos
  label=''
  if t[pos]=='(':
   pos+=1;z=set()
   while True:
    z|=node()
    if t[pos]==',':pos+=1;continue
    need(t[pos]==')','Malformed tree');pos+=1;break
   if t[pos] not in {':',',',')',';'}:label=t[pos];pos+=1
  else:z={t[pos]};pos+=1
  if t[pos]==':':pos+=2
  ee.append((z,label));return z
 universe=node();need(t[pos]==';','Invalid tree ending');return universe,[(z,l) for a,l in ee for z in (a,universe-a) if 1<len(z)<len(universe)-1]

src=release('step25kz17.67','PPT_TREE_RELEASE');screen=release('step25kz17.66','PPT_SCREEN_RELEASE')
prov=json.loads((src/'PROVENANCE.json').read_text());need(prov['source_manifest_sha256']==sha(screen/'SHA256SUMS.txt'),'Source chain mismatch')
rows=read(screen/'09_combined_PPT_specificity_review.tsv');old=read(src/'source_01_integrated_locus_evidence.tsv');seqs,_=fasta(screen/'05_current_candidate_proteins.faa')
need(len(rows)==451 and len(old)==1109 and not any(r['wang_family']=='PPT' for r in old),'Unexpected ledger scope')
meta={r['tree_id']:r for r in read(src/'02_tree_tip_mapping.tsv')};byprotein={r['protein_id']:r for r in meta.values()};refs={r['tree_id'] for r in meta.values() if r['role']=='REFERENCE'};anchors={'Q8RXN3','Q8H0T6'}
u,ee=edges(src/'main/tree.treefile');v,ff=edges(src/'sensitivity_all_partials/tree.treefile')
def group(edges_):
 options=sorted([(len(z),z,l) for z,l in edges_ if z&refs==anchors],key=lambda x:x[0]);need(options,'No exclusive PPT-reference group');return options[0][1:]
a,sa=group(ee);b,sb=group(ff);need(a==b&u and len(a-refs)==26 and len(b-refs)==29,'PPT association changed; inspect trees')
core={meta[k]['protein_id'] for k in a-refs};partials={meta[k]['protein_id'] for k in b-u};long='PO__evm.TU.LG01.2472'
need(partials=={'SP_Hex__gene_Spor0301920.1','SP_Tet__g13815','SP_Tet__g13816'} and long in core,'Exception scope changed')
exceptions=partials|{long};qmeta={r['protein_id']:r for r in rows};ledger=rows
out=base/'step25kz17.68'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'annotation').mkdir()
def table(n,rr,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in sorted(d.items())))
for sp in {qmeta[k]['species'] for k in exceptions}:
 mp=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv');need(sha(mp)==prov['source_provenance']['mapping_hashes'][sp],'Mapping provenance mismatch '+sp)
print('Reviewing four mapped annotation exceptions; no alignment or tree rerun',flush=True)
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
table('20_annotation_file_map.tsv',mappingrows,list(mappingrows[0]));table('21_annotation_collection_status.tsv',status,list(status[0]));table('22_protein_CDS_consistency.tsv',checks,list(checks[0]));table('23_annotation_source_provenance.tsv',sources,list(sources[0]))

# Translation agreement is a consistency check, not evidence that an incomplete model is complete.
annotation_summary=[]
for k in sorted(exceptions):
 cc=[x for x in checks if x['candidate_id']==k];matched=[x for x in cc if x['protein_matches_candidate'] is True]
 consistent=[x for x in matched if x['translation_matches_candidate'] is True and x['frame_remainder']==0 and x['internal_stops']==0 and x['ambiguous_amino_acids']==0]
 state='PROTEIN_CDS_CONSISTENT_CAVEAT_RETAINED' if consistent else 'UNRESOLVED_PROTEIN_CDS_CONSISTENCY'
 annotation_summary.append(dict(protein_id=k,mapped_isoforms=len(cc),matching_source_proteins=len(matched),consistent_CDS_records=len(consistent),annotation_status=state,remaining_limitation='PARTIAL_SEQUENCE' if k in partials else 'LONG_FULL_PROTEIN',closure='BOUNDED_REVIEW_COMPLETE_NO_READ_VALIDATION'))
table('24_exception_review_outcomes.tsv',annotation_summary,list(annotation_summary[0]));astatus={r['protein_id']:r['annotation_status'] for r in annotation_summary}
association=[];new=[]
for r in rows:
 k=r['protein_id'];tid=byprotein.get(k,{}).get('tree_id','');lim='BOUNDED_SIX_REFERENCE_PANEL;MAIN_SH_UF='+sa+';SENSITIVITY_SH_UF='+sb+';NO_FUNCTION_OR_EXHAUSTIVE_COUNT_CLAIM'
 if k in core:
  category='PPT_REFERENCE_GROUP_STABLE_MEMBERSHIP_SUPPORT_SENSITIVE';tier='LIMITED_PANEL_CAVEAT'
  if k==long:category='PPT_REFERENCE_GROUP_LONG_MODEL';tier='MODEL_CAVEAT';lim+=';LONG_FULL_PROTEIN;'+astatus[k]
 elif k in partials:category='PPT_REFERENCE_GROUP_PARTIAL_SENSITIVITY_ONLY';tier='PARTIAL';lim+=';JOINT_FRAGMENT_PLACEMENT_UNCERTAIN;'+astatus[k]
 elif r['historical_OG0000369']=='True':category='HISTORICAL_PPT_LABEL_UNSUPPORTED';tier='EXCLUDED';lim+=';NO_PANEL_HIT_OR_PF03151;NOT_A_GENOMIC_ABSENCE_CLAIM'
 else:category='OUTSIDE_BOUNDED_PPT_REFERENCE_GROUP';tier='OUTSIDE_FAMILY_SET';lim+=';NOT_PROVEN_NON_PPT;RETAIN_SCREEN_AND_DOMAIN_EVIDENCE'
 if k in exceptions and astatus[k]=='UNRESOLVED_PROTEIN_CDS_CONSISTENCY':lim+=';ANNOTATION_EVIDENCE_UNRESOLVED'
 association.append(dict(protein_id=k,tree_id=tid,main_PPT_group=k in core,sensitivity_PPT_group=k in core|partials,main_support=sa,sensitivity_support=sb,evidence_category=category,count_tier=tier))
 new.append(dict(wang_family='PPT',species=r['species'],gene_id=r['gene_id'],protein_id=k,subunit='',evidence_category=category,count_tier=tier,source_file='step25kz17.68/02_PPT_phylogenetic_associations.tsv;step25kz17.68/24_exception_review_outcomes.tsv',model_or_scope_limitation=lim,exhaustive_count='NOT_ESTABLISHED',C4_CAM_specific_function='NOT_ESTABLISHED'))
need(len({(r['species'],r['gene_id']) for r in new})==451,'Duplicate PPT locus');combined=old+new
fields=list(old[0]);table('01_integrated_locus_evidence.tsv',combined,fields);table('02_PPT_phylogenetic_associations.tsv',association,list(association[0]));table('03_PPT_closed_evidence_ledger.tsv',new,fields)
summary=[]
for sp in sorted({r['species'] for r in new}):
 rr=[r for r in new if r['species']==sp];summary.append(dict(species=sp,main_PPT_associated=sum(r['protein_id'] in core for r in rr),partial_PPT_associated=sum(r['protein_id'] in partials for r in rr),long_model_caveats=sum(r['protein_id']==long for r in rr),unsupported_historical_label=sum(r['count_tier']=='EXCLUDED' for r in rr),other_screened_records=sum(r['count_tier']=='OUTSIDE_FAMILY_SET' for r in rr),final_exhaustive_count='NA'))
table('04_PPT_species_evidence_summary.tsv',summary,list(summary[0]))
coverage=read(src/'09_refreshed_family_coverage.tsv')
for r in coverage:
 if r['wang_family']=='PPT':r.update(integrated_records=str(len(new)),current_status='BOUNDED_REVIEW_CLOSED_WITH_SUPPORT_AND_MODEL_CAVEATS',next_action='NO_ROUTINE_RERUN;REOPEN_ONLY_WITH_NEW_REFERENCE_OR_MODEL_EVIDENCE',priority='4',search_completeness='NOT_ESTABLISHED')
table('06_refreshed_family_coverage.tsv',coverage,list(coverage[0]));remaining=[r for r in coverage if int(r['integrated_records'])==0];table('08_families_without_integrated_ledger.tsv',remaining,list(coverage[0]))
shutil.copyfile(src/'source_07_final_Wang_counts_PENDING.tsv',out/'07_final_Wang_counts_PENDING.tsv')
loci=defaultdict(list)
for r in combined:loci[r['species'],r['gene_id']].append(r)
overlaps=[dict(species=sp,gene_id=g,families=';'.join(sorted({r['wang_family'] for r in rr})),evidence_categories=';'.join(sorted({r['evidence_category'] for r in rr}))) for (sp,g),rr in loci.items() if len({r['wang_family'] for r in rr})>1]
table('09_cross_family_overlap_review.tsv',overlaps,['species','gene_id','families','evidence_categories'])
for mode in ['main','sensitivity_all_partials']:
 (out/mode).mkdir()
 for n in ['tree.treefile','itol_labels.txt']:shutil.copyfile(src/mode/n,out/mode/n)
for n in ['00_OG0000369_label_correction.tsv','01_recalibrated_domain_completeness.tsv','02_tree_tip_mapping.tsv','08_tree_run_summary.tsv']:shutil.copyfile(src/n,out/('source_'+n))
(out/'VALIDATION_STOPPING_RULES.txt').write_text('One combined validation pass per family: exact input and locus provenance; panel specificity/domain evidence; one main and one sensitivity tree when necessary; one targeted annotation/CDS review of exceptions; then integrate evidence tiers and close. Do not iterate trees until support passes a threshold. Reopen only if new independent evidence could change a count or family assignment. Translation consistency cannot establish gene-model completeness, biological function, independent transcript support or genome-wide absence. Never merge partial/model/support caveats into formal counts. Remaining unintegrated families remain NA, not zero.\n')
(out/'README.txt').write_text('PPT bounded review closed with limitations. Existing trees reused without inference rerun.26 main-associated candidates have stable membership after adding fragments; support97.5/95 main versus97.2/90 sensitivity remains explicit.25 categorized LIMITED_PANEL_CAVEAT,one long PO protein MODEL_CAVEAT;3 sensitivity-only partials remain PARTIAL.37 OG0000369 labels unsupported and excluded from bounded PPT set;385 other screened loci remain OUTSIDE_FAMILY_SET, not proven non-PPT. Four exceptions receive exact annotation protein/CDS and5kb GFF-context collection; missing evidence remains unresolved. No model repair or independent read validation. New451 evidence records appended once to prior1109, not451 additional PPT genes. Final count matrix remains pending. Species summary main count includes the long-model subset and must not be summed with it.\n')
qc=dict(status='PPT_BOUNDED_REVIEW_CLOSED_WITH_CAVEATS',integrated_records=len(combined),integrated_families=len({r['wang_family'] for r in combined}),PPT_tiers=dict(Counter(r['count_tier'] for r in new)),main_associated=len(core),partial_associated=len(partials),annotation_outcomes=annotation_summary,cross_family_overlap_flags=len(overlaps),families_without_integrated_ledger=len(remaining),new_tree_runs=0,final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');(out/'PROVENANCE.json').write_text(json.dumps(dict(step67_manifest_sha256=sha(src/'SHA256SUMS.txt'),step66_manifest_sha256=sha(screen/'SHA256SUMS.txt'),input_manifest_sha256=sha(manifest),main_tree_sha256=sha(src/'main/tree.treefile'),sensitivity_tree_sha256=sha(src/'sensitivity_all_partials/tree.treefile')),indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='PPT_closed_annotation_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle))

PY
