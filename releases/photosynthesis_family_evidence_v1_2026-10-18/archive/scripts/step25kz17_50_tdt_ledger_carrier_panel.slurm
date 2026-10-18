#!/usr/bin/env bash
#SBATCH --job-name=tdt_close_panel
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=tdt_close_panel_%j.out
#SBATCH --error=tdt_close_panel_%j.err
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


import urllib.parse
src=release('step25kz17.49','TD029_REVIEW_RELEASE');screen=release('step25kz17.46','TRANSPORTER_RECOVERY_RELEASE');tree=release('step25kz17.48','TDT_TREE_RELEASE')
need(json.loads((src/'PROVENANCE.json').read_text())['source_manifest_sha256']==sha(tree/'SHA256SUMS.txt'),'Step49/48 provenance mismatch')
need(json.loads((tree/'PROVENANCE.json').read_text())['step46_manifest_sha256']==sha(screen/'SHA256SUMS.txt'),'Step48/46 provenance mismatch')
splits=read(src/'05_TDT_split_comparison.tsv');need(len(splits)==4 and all(r['exact_TDT_bipartition']=='True' and r['support_SHaLRT_UFBoot']=='100/100' for r in splits),'Expected reviewed TDT split evidence changed')
meta=read(src/'07_tree_id_map.tsv');tdt=[r for r in meta if r['role']=='CANDIDATE' and r['screen_group']=='TDT'];need(len(tdt)==11,'Expected 11 TDT loci')
short='SP_Hex__gene_Spor0444110.1';qc=read(src/'08_short_candidate_CDS_check.tsv');need(len(qc)==1 and qc[0]['protein_id']==short and qc[0]['translation_matches']=='True' and qc[0]['candidate_matches']=='True','Short model consistency mismatch')
hits=read(src/'03_domain_hits.tsv');models={r['accession'].split('.')[0] for r in hits if r['sequence_id']=='TD029'};need(models=={'PF02902'},'TD029 evidence changed')
oldledger=read(src/'source_prior_01_integrated_locus_evidence.tsv');need(len(oldledger)==875 and not any(r['wang_family']=='TDT' for r in oldledger),'Unexpected prior integrated ledger')
out=base/'step25kz17.50'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'references').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
ff=list(oldledger[0]);new=[]
for r in tdt:
 caveat=r['protein_id']==short
 new.append(dict(wang_family='TDT',species=r['species'],gene_id=r['gene_id'],protein_id=r['protein_id'],subunit='',evidence_category='TDT_MODEL_COMPLETENESS_CAVEAT' if caveat else 'TDT_PHYLOGENETIC_DOMAIN_SUPPORTED',count_tier='MODEL_CAVEAT' if caveat else 'SUPPORTED',source_file='step25kz17.49/05_TDT_split_comparison.tsv;step25kz17.47/08_TDT_completeness_review.tsv',model_or_scope_limitation=('N_TERMINAL_REFERENCE_GAP_92AA;CDS_CONSISTENT;' if caveat else '')+'LIMITED_REFERENCE_PANEL;SEARCH_COMPLETENESS_NOT_ESTABLISHED;TRANSPORT_FUNCTION_NOT_TESTED',exhaustive_count='NOT_ESTABLISHED',C4_CAM_specific_function='NOT_ESTABLISHED'))
ledger=oldledger+new;need(len({(r['wang_family'],r['species'],r['gene_id']) for r in ledger})==len(ledger),'Duplicate family/locus in integrated ledger')
table('01_integrated_locus_evidence.tsv',ledger,ff);table('02_TDT_locus_evidence.tsv',new,ff)
excluded=next(r for r in meta if r['tree_id']=='TD029');ex=dict(excluded,decision='EXCLUDED_FROM_SUPPORTED_TDT_DIT_SET',reason='PF02902_Peptidase_C48;NO_PF00939_OR_PF03600;LIMITED_ORIGINAL_SIMILARITY_SPAN',precise_function='UNASSIGNED');table('03_TD029_exclusion_record.tsv',[ex],list(ex))
ss=[]
for sp in sorted({r['species'] for r in oldledger}):
 rr=[r for r in new if r['species']==sp];ss.append(dict(species=sp,supported_loci=sum(r['count_tier']=='SUPPORTED' for r in rr),model_caveat_loci=sum(r['count_tier']=='MODEL_CAVEAT' for r in rr),reviewed_TDT_loci=len(rr),final_exhaustive_count='NA'))
table('04_TDT_species_evidence_summary.tsv',ss,list(ss[0]))
coverage=read(src/'source_prior_06_refreshed_family_coverage.tsv')
for r in coverage:
 if r['wang_family']=='TDT':r.update(current_status='BOUNDED_REVIEW_CLOSED_LEDGER_INTEGRATED',integrated_records='11',search_completeness='NOT_ESTABLISHED',priority='4',next_action='RETAIN_10_SUPPORTED_1_MODEL_CAVEAT;NO_FUNCTIONAL_OR_EXHAUSTIVE_COUNT_CLAIM')
 if r['wang_family'] in {'DTC','DIC'}:r.update(current_status='SCREEN_QUEUE_AVAILABLE_REFERENCE_PANEL_PREPARATION',priority='1',next_action='REVIEW_BROAD_CARRIER_REFERENCE_PANEL_THEN_SPECIFICITY_VALIDATION',search_completeness='NOT_ESTABLISHED')
table('06_refreshed_family_coverage.tsv',coverage,list(coverage[0]));shutil.copyfile(src/'source_prior_07_final_Wang_counts_PENDING.tsv',out/'07_final_Wang_counts_PENDING.tsv')
queue=read(src/'source_10_mitochondrial_carrier_queue_UNCHANGED.tsv');need(len(queue)==724,'Unexpected carrier queue');seqs,_=fasta(screen/'05_all_discovered_and_prior_candidates.faa');need(len({r['protein_id'] for r in queue})==724,'Duplicate queued protein')
table('08_mitochondrial_carrier_queue_REVIEW_ONLY.tsv',queue,list(queue[0]));writefa('09_mitochondrial_carrier_candidates.faa',{r['protein_id']:seqs[r['protein_id']] for r in queue})
# Existing curated anchors remain immutable. New reviewed carrier-domain entries are comparators for review.
refs,_=fasta(screen/'02_reference_panel.faa');anchor=[r for r in read(screen/'01_reference_definitions.tsv') if r['family'] in {'DTC','DIC','UCP1'}];need(len(anchor)==5,'Expected five existing carrier anchors');anchorby={r['accession']:r for r in anchor};panel={r['reference_id']:refs[r['reference_id']] for r in anchor}
for r in anchor:need(hashlib.sha256(refs[r['reference_id']].encode()).hexdigest()==r['sequence_sha256'],'Anchor sequence hash mismatch')
table('references/01_existing_anchor_definitions.tsv',anchor,list(anchor[0]));writefa('references/01_existing_anchors.faa',panel)
query='(organism_id:3702) AND (reviewed:true) AND (xref:pfam-PF00153)';url='https://rest.uniprot.org/uniprotkb/stream?'+urllib.parse.urlencode({'query':query,'format':'json'});panel_status='PENDING';error='';added=[]
try:
 if os.environ.get('CARRIER_PANEL_JSON'):raw=(root/os.environ['CARRIER_PANEL_JSON']).read_bytes()
 else:
  raw=None
  for attempt in range(3):
   try:
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Sesuvium-carrier-reference-audit/1.0'}),timeout=45) as response:raw=response.read()
    break
   except Exception:
    if attempt==2:raise
    time.sleep(2)
 (out/'references/02_downloaded_reviewed_carriers.json').write_bytes(raw);entries=json.loads(raw)['results'];need(len(entries)>=10,'Too few carrier entries; inspect panel query')
 seen=set();addseq={}
 for e in entries:
  acc=e['primaryAccession'];need(acc not in seen,'Duplicate accession');seen.add(acc);need(e['organism']['taxonId']==3702 and e['entryType']=='UniProtKB reviewed (Swiss-Prot)','Unexpected species or review status')
  pfams={x['id'] for x in e.get('uniProtKBCrossReferences',[]) if x['database']=='Pfam'};need('PF00153' in pfams,'Missing carrier-domain annotation')
  seq=e['sequence']['value'];need(valid(seq),'Invalid reference protein');genes=e.get('genes',[]);loci=[x['value'] for g in genes for x in g.get('orderedLocusNames',[])];names=[g['geneName']['value'] for g in genes if 'geneName' in g];protein=e.get('proteinDescription',{}).get('recommendedName',{}).get('fullName',{}).get('value','')
  if acc in anchorby:need(seq==refs[anchorby[acc]['reference_id']],'Existing anchor sequence changed; review required')
  else:addseq['REF__'+acc]=seq
  added.append(dict(reference_id='REF__'+acc,accession=acc,loci=';'.join(loci),gene_names=';'.join(names),protein_name=protein,role='EXISTING_ANCHOR' if acc in anchorby else 'BROAD_CARRIER_COMPARATOR_REVIEW_REQUIRED',assigned_family=anchorby[acc]['family'] if acc in anchorby else 'UNASSIGNED',length_aa=len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),pfam_ids=';'.join(sorted(pfams)),source_url='https://www.uniprot.org/uniprotkb/'+acc,curation_status='SOURCE_ANNOTATION_RETAINED_FUNCTIONAL_SPECIFICITY_NOT_INFERRED'))
 need(len(addseq)>=5,'Insufficient expansion beyond existing anchors');panel.update(addseq);panel_status='BROAD_PANEL_READY_FOR_DEFINITION_REVIEW'
except (Exception,SystemExit) as exc:
 error=str(exc).replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>');panel_status='PANEL_PREPARATION_INCOMPLETE';added=[]
if added:table('references/03_broad_panel_annotations_REVIEW_REQUIRED.tsv',added,list(added[0]))
writefa('references/04_carrier_panel.faa',panel)
(out/'references/05_panel_status.json').write_text(json.dumps(dict(status=panel_status,query=query,url=url,retrieval_utc=datetime.now(timezone.utc).isoformat(),panel_sequences=len(panel),existing_anchors=5,error=error,scope='Reviewed Arabidopsis PF00153-containing proteins; includes organellar carrier comparators, not all assumed mitochondrial or experimentally characterized. Canonical entries only; not exhaustive.'),indent=2)+'\n')
for n in ['05_TDT_split_comparison.tsv','03_domain_hits.tsv','08_short_candidate_CDS_check.tsv','source_08_TDT_completeness_review.tsv']:shutil.copyfile(src/n,out/('evidence_'+n))
(out/'PROVENANCE.json').write_text(json.dumps(dict(step49_manifest_sha256=sha(src/'SHA256SUMS.txt'),step48_manifest_sha256=sha(tree/'SHA256SUMS.txt'),step46_manifest_sha256=sha(screen/'SHA256SUMS.txt'),prior_ledger_sha256=sha(src/'source_prior_01_integrated_locus_evidence.tsv')),indent=2)+'\n')
(out/'README.txt').write_text('TDT bounded review closed: 10 supported loci, one N-terminal completeness caveat with CDS consistency; stable 100/100 TDT split in four trees. TD029 excluded from supported TDT/DIT set, precise function unassigned. Integrated ledger 886 records in eight families. Final Wang matrix remains NA: evidence counts are not exhaustive or C4/CAM functional counts. 724 carrier candidates preserved with exact sequences. Expanded panel preparation uses reviewed Arabidopsis PF00153 entries, retaining source annotations and raw JSON; broad comparators remain unassigned until definition review. Reviewed database status alone is not experimental transport evidence. Panel download/validation failure is reported separately and does not discard completed ledger consolidation. Inspect QC and references/05_panel_status.json before proceeding. No similarity searches, alignments, trees or gene-model repairs performed.\n')
qc=dict(status='TDT_LEDGER_CONSOLIDATED',integrated_records=len(ledger),integrated_families=len({r['wang_family'] for r in ledger}),TDT_supported=10,TDT_model_caveat=1,TD029_excluded=True,carrier_candidates=len(queue),carrier_panel_status=panel_status,carrier_panel_sequences=len(panel),final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='TDT_ledger_carrier_panel/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
