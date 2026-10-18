#!/usr/bin/env bash
#SBATCH --job-name=tdt_dic_recover
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=tdt_dic_recover_%j.out
#SBATCH --error=tdt_dic_recover_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v diamond
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
src=release('step25kz17.45','TRANSPORTER_SCREEN_RELEASE');cor=release('step25kz17.21','CORRECTED_INPUT_RELEASE')
proteome=cor/'full_outputs/corrected_protein_candidates.faa';expected=json.loads((src/'PROVENANCE.json').read_text())['corrected_proteome_sha256']
need(sha(proteome)==expected,'Corrected protein provenance mismatch')
refs,_=fasta(src/'02_target_and_comparator_references.faa');refmeta=read(src/'01_reference_definitions.tsv');rd={r['reference_id']:r for r in refmeta}
need(len(refs)==len(rd)==8 and set(refs)==set(rd),'Expected eight verified references')
for k,s in refs.items():need(hashlib.sha256(s.encode()).hexdigest()==rd[k]['sequence_sha256'],'Reference sequence hash mismatch')
prior=read(src/'09_family_specificity_REVIEW_ONLY.tsv');oldids={r['protein_id'] for r in prior};need(len(prior)==35,'Expected 35 prior pooled candidates')
out=base/'step25kz17.46'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir();(out/'family_queues').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+s+'\n' for k,s in d.items()))
shutil.copyfile(src/'01_reference_definitions.tsv',out/'01_reference_definitions.tsv');writefa('02_reference_panel.faa',refs)
for r in prior:r.update(queue_status='PROVISIONAL_DIT_RELATED_REVIEW_ONLY',count_status='NOT_INCLUDED_IN_TDT_DTC_DIC_COUNTS')
table('03_prior_DIT_related_queue.tsv',prior,list(prior[0]))
commands=[];cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
def run(args,log):
 commands.append([str(x).replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for x in args])
 with (out/log).open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
run(['diamond','version'],'diamond_version.txt')
run(['diamond','makedb','--in',str(proteome),'--db',str(out/'work/corrected_proteome'),'--threads',cpu],'04_makedb.log')
fmt=['qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore']
run(['diamond','blastp','--query',str(out/'02_reference_panel.faa'),'--db',str(out/'work/corrected_proteome'),'--out',str(out/'04_reference_to_proteome.tsv'),'--threads',cpu,'--very-sensitive','--evalue','1e-3','--max-target-seqs','0','--max-hsps','1','--outfmt','6']+fmt,'04_discovery.log')
def hits(p):
 rr=[]
 for line in p.read_text().splitlines():
  if not line:continue
  x=line.split('\t');need(len(x)==12,'Unexpected DIAMOND columns')
  rr.append(dict(query=x[0],subject=x[1],pident=float(x[2]),qlen=int(x[4]),slen=int(x[5]),query_span=(int(x[7])-int(x[6])+1)/int(x[4]),subject_span=(int(x[9])-int(x[8])+1)/int(x[5]),evalue=float(x[10]),bitscore=float(x[11])))
 return rr
forward=hits(out/'04_reference_to_proteome.tsv');need(all(r['query'] in refs for r in forward),'Unexpected discovery reference')
byprotein=defaultdict(list)
for r in forward:byprotein[r['subject']].append(r)
wanted=set(byprotein)|oldids;seqs={};key=None;parts=[];seen=set();species_counts=Counter()
def commit():
 if key is not None:
  need(key not in seen,'Duplicate input protein ID');seen.add(key);species_counts[key.split('__',1)[0]]+=1
  if key in wanted:
   s=''.join(parts).upper();need(valid(norm(s)),'Invalid candidate protein '+key);seqs[key]=norm(s)
with proteome.open() as f:
 for line in f:
  if line.startswith('>'):commit();key=line[1:].split()[0];parts=[]
  elif key in wanted:parts.append(''.join(line.split()))
 commit()
need(set(seqs)==wanted,'Discovery/prior identifier absent from corrected proteome')
need(seqs,'No candidates available');writefa('05_all_discovered_and_prior_candidates.faa',seqs)
run(['diamond','makedb','--in',str(out/'02_reference_panel.faa'),'--db',str(out/'work/panel'),'--threads',cpu],'06_panel_makedb.log')
run(['diamond','blastp','--query',str(out/'05_all_discovered_and_prior_candidates.faa'),'--db',str(out/'work/panel'),'--out',str(out/'06_candidate_to_panel.tsv'),'--threads',cpu,'--very-sensitive','--evalue','1e-3','--max-target-seqs','0','--max-hsps','1','--outfmt','6']+fmt,'06_panel_comparison.log')
reverse=defaultdict(list)
for r in hits(out/'06_candidate_to_panel.tsv'):
 need(r['query'] in seqs and r['subject'] in rd,'Unexpected comparison identifier');reverse[r['query']].append(r)
# Exact locus mapping; evidence is retained even if mapping is unresolved.
maphash={};maps={}
for sp in {k.split('__',1)[0] for k in seqs}:
 mp=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv');maphash[sp]=sha(mp);lookup=defaultdict(set)
 for r in read(mp):
  for f in ['gene_id','transcript_id','protein_id','cds_id']:
   if r.get(f) and r[f] not in {'.','NA'}:lookup[r[f]].add(r['gene_id'])
 maps[sp]=lookup
review=[];evidence=[];targets={'TDT','DTC','DIC'}
for k,s in sorted(seqs.items()):
 sp,local=k.split('__',1);genes=maps[sp].get(local,set());hh=sorted(reverse[k],key=lambda r:(-r['bitscore'],r['subject']));best=hh[0] if hh else None
 family=rd[best['subject']]['family'] if best else '';other=next((r for r in hh if rd[r['subject']]['family']!=family),None) if best else None
 margin=best['bitscore']-other['bitscore'] if other else None
 flags=['LIMITED_REFERENCE_PANEL']
 if len(genes)!=1:flags.append('LOCUS_MAPPING_UNRESOLVED')
 if not best:flags.append('NO_REVERSE_PANEL_HIT')
 else:
  if best['query_span']<.7 or best['subject_span']<.7:flags.append('LIMITED_SINGLE_HSP_SPAN')
  if best['evalue']>1e-5:flags.append('WEAK_EVALUE')
  if margin is None:flags.append('NO_OTHER_FAMILY_HIT')
  elif margin<20:flags.append('SMALL_FAMILY_MARGIN')
 status='TARGET_FAVOURED_REVIEW' if family in targets else 'COMPARATOR_FAVOURED_REVIEW' if family else 'UNRESOLVED'
 review.append(dict(protein_id=k,species=sp,gene_id=next(iter(genes)) if len(genes)==1 else '',mapping_status='UNAMBIGUOUS' if len(genes)==1 else 'UNRESOLVED',length_aa=len(s),prior_pooled_candidate=k in oldids,discovery_families=';'.join(sorted({rd[r['query']]['family'] for r in byprotein[k]})),best_panel_family=family,best_reference=best['subject'] if best else '',bitscore=best['bitscore'] if best else '',evalue=best['evalue'] if best else '',query_span=round(best['query_span'],5) if best else '',reference_span=round(best['subject_span'],5) if best else '',next_family=rd[other['subject']]['family'] if other else '',between_family_margin=margin if margin is not None else '',review_category=status,flags=';'.join(flags),formal_family_assignment='PENDING',final_gene_count='NA'))
 for fam in sorted(targets):
  fh=[r for r in byprotein[k] if rd[r['query']]['family']==fam];rh=[r for r in hh if rd[r['subject']]['family']==fam]
  if fh or rh:evidence.append(dict(protein_id=k,family=fam,forward_reference_hits=len(fh),reverse_reference_hits=len(rh),best_forward_bitscore=max((r['bitscore'] for r in fh),default=''),best_reverse_bitscore=max((r['bitscore'] for r in rh),default=''),best_overall_family=family,decision='REVIEW_ONLY_NOT_AN_ADDITIONAL_GENE'))
table('07_all_candidate_evidence_REVIEW_ONLY.tsv',review,list(review[0]));table('08_target_family_hit_evidence.tsv',evidence,['protein_id','family','forward_reference_hits','reverse_reference_hits','best_forward_bitscore','best_reverse_bitscore','best_overall_family','decision'])
# Family queues retain every target-linked hit, including comparator-favoured ones.
# A protein may appear in several queues; this is not multiple locus counting.
for fam in sorted(targets):
 keys={r['protein_id'] for r in evidence if r['family']==fam}
 table('family_queues/'+fam+'_all_linked_hits.tsv',[r for r in review if r['protein_id'] in keys],list(review[0]))
 writefa('family_queues/'+fam+'_best_panel_candidates.faa',{r['protein_id']:seqs[r['protein_id']] for r in review if r['best_panel_family']==fam})
summary=[]
for sp in sorted(species_counts):
 for fam in ['TDT','DTC','DIC','DIT1','DIT2','UCP1']:
  rr=[r for r in review if r['species']==sp and r['best_panel_family']==fam]
  summary.append(dict(species=sp,best_panel_family=fam,screened_proteome_records=species_counts[sp],candidate_proteins=len(rr),unambiguous_loci=len({r['gene_id'] for r in rr if r['gene_id']}),formal_count='NA'))
table('09_species_screening_summary.tsv',summary,list(summary[0]))
for n in ['prior_01_integrated_locus_evidence.tsv','prior_06_refreshed_family_coverage.tsv','prior_07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/n)
need(sha(proteome)==expected,'Corrected input changed during run')
qc=dict(status='FULL_CORRECTED_PROTEOME_SIMILARITY_SCREEN_REVIEW_REQUIRED',proteome_records=len(seen),reference_controls=len(refs),reference_proteome_alignments=len(forward),discovered_proteins=len(byprotein),prior_candidates_retained=len(oldids),total_review_proteins=len(seqs),best_panel_families=dict(Counter(r['best_panel_family'] for r in review)),mapping_status=dict(Counter(r['mapping_status'] for r in review)),domain_validation='NEXT_STEP',biological_completeness='NOT_ESTABLISHED',final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),corrected_proteome_sha256=expected,mapping_hashes=maphash,commands=commands),indent=2)+'\n')
(out/'README.txt').write_text('Reference-to-full-corrected-proteome discovery followed by candidate-to-eight-reference-panel comparison. Reuses verified references, no downloads. DIAMOND very-sensitive, E<=1e-3, all reported targets, one HSP per pair. These are permissive discovery settings; 70% spans, E<=1e-5 and 20-bit margin are descriptive review flags, not formal acceptance. Search direction changes can affect scores; this is not reciprocal-best-hit orthology. All discovered and 35 prior pooled proteins retained. Separate target-family tables include all linked hits, even when comparator favoured; overlapping queue membership is not multiple gene counting. Empty family FASTAs mean no best-panel hit under this screen, not biological absence. No domain search, alignment, phylogeny, localization or function validation in this step. DIC/DTC/UCP discrimination needs broader carrier comparisons and family-specific validation; TDT/DIT share domain features. Whole annotated proteome searched, but completeness beyond sequence similarity and genome annotation is not established. Existing 875-row integrated ledger unchanged. Historical pooled DIT-related candidates remain provisional; no automatic transfer into Wang DIT2 counts.\n')
for p in out.glob('*.log'):p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Transporter_full_proteome_screen/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
