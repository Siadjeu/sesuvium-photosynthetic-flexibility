#!/usr/bin/env bash
#SBATCH --job-name=carrier_close
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=carrier_close_%j.out
#SBATCH --error=carrier_close_%j.err
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




src=release('step25kz17.54','CARRIER_REVIEW_RELEASE');screen=release('step25kz17.51','CARRIER_SCREEN_RELEASE')
old=read(src/'source_prior_01_integrated_locus_evidence.tsv');need(len(old)==886 and not any(r['wang_family'] in {'DTC','DIC'} for r in old),'Unexpected prior ledger or already integrated families');need(sha(src/'source_prior_01_integrated_locus_evidence.tsv')==sha(screen/'prior_01_integrated_locus_evidence.tsv'),'Prior ledger provenance mismatch')
evidence=read(src/'02_candidate_phylogenetic_evidence.tsv');need(len(evidence)==85,'Expected 85 reviewed candidates');checks=read(src/'08_protein_CDS_consistency.tsv');expected={'PO__evm.TU.LG05.272','PO__evm.TU.LG25.1042','SP_Hex__gene_Spor0293270.1','SS__g18155.t1'}
need({r['candidate_id'] for r in checks}==expected,'Exception evidence mismatch')
for k in expected:need(any(r['candidate_id']==k and r['protein_matches_candidate']=='True' and r['translation_matches_candidate']=='True' and r['frame_remainder']=='0' and r['internal_stops']=='0' and r['ambiguous_amino_acids']=='0' for r in checks),'Exception consistency changed '+k)
spans={r['protein_id']:r for r in read(src/'05_repeat_overlap_review.tsv')};need(all(spans[k]['hit_count']=='4' for k in ['PO__evm.TU.LG05.272','PO__evm.TU.LG25.1042','REF__Q9FY68']),'DIC domain evidence changed');need(spans['SP_Hex__gene_Spor0293270.1']['hit_count']=='2' and spans['SS__g18155.t1']['hit_count']=='1','Partial domain evidence changed')
queue=read(screen/'09_all_candidates_specificity_REVIEW_ONLY.tsv');need(len(queue)==724,'Expected full 724-candidate screen');qmap={r['protein_id']:r for r in queue};need(len(qmap)==724,'Duplicate candidate');seqs,_=fasta(screen/'01_candidates.faa');need(set(seqs)==set(qmap),'Candidate sequence mismatch')
for r in evidence:need(r['protein_id'] in qmap and all(r[f]==qmap[r['protein_id']][f] for f in ['species','gene_id']),'Review/screen locus mismatch')
new=[];ff=list(old[0]);selected=set()
for r in evidence:
 fam=r['phylogenetic_family']
 if fam not in {'DTC','DIC'}:continue
 k=r['protein_id'];tier='SUPPORTED';category='PHYLOGENETIC_DOMAIN_SUPPORTED';lim='BOUNDED_REFERENCE_PANEL;SEARCH_COMPLETENESS_AND_TRANSPORT_FUNCTION_NOT_ESTABLISHED'
 if k=='SP_Hex__gene_Spor0293270.1':tier='MODEL_CAVEAT';category='DTC_COMPLETENESS_CAVEAT';lim+=';TWO_DETECTED_CARRIER_HITS;CDS_CONSISTENT'
 elif k=='SS__g18155.t1':tier='PARTIAL';category='DTC_ASSOCIATED_PARTIAL';lim+=';121AA_SINGLE_REPEAT;SENSITIVITY_PLACEMENT_ONLY;SUPPLIED_CDS_HAS_NO_TERMINAL_STOP;CDS_CONSISTENT'
 elif k in expected:lim+=';FOUR_HIT_LAYOUT_SIMILAR_TO_DIC_REFERENCE_Q9FY68;CDS_CONSISTENT'
 new.append(dict(wang_family=fam,species=r['species'],gene_id=r['gene_id'],protein_id=k,subunit='',evidence_category=category,count_tier=tier,source_file='step25kz17.54/02_candidate_phylogenetic_evidence.tsv;step25kz17.54/05_repeat_overlap_review.tsv;step25kz17.54/08_protein_CDS_consistency.tsv',model_or_scope_limitation=lim,exhaustive_count='NOT_ESTABLISHED',C4_CAM_specific_function='NOT_ESTABLISHED'));selected.add(k)
counts=Counter((r['wang_family'],r['count_tier']) for r in new);need(counts=={('DTC','SUPPORTED'):41,('DTC','MODEL_CAVEAT'):1,('DTC','PARTIAL'):1,('DIC','SUPPORTED'):17},'Unexpected consolidation totals')
ledger=old+new;need(len({(r['wang_family'],r['species'],r['gene_id']) for r in ledger})==len(ledger),'Duplicate family/locus')
out=base/'step25kz17.55'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
def table(n,rr,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
table('01_integrated_locus_evidence.tsv',ledger,ff);table('02_DTC_DIC_closed_evidence_ledger.tsv',new,ff)
summary=[]
for fam in ['DTC','DIC']:
 for sp in sorted({r['species'] for r in old}):
  rr=[r for r in new if r['wang_family']==fam and r['species']==sp];summary.append(dict(wang_family=fam,species=sp,supported_loci=sum(r['count_tier']=='SUPPORTED' for r in rr),model_caveat_loci=sum(r['count_tier']=='MODEL_CAVEAT' for r in rr),partial_loci=sum(r['count_tier']=='PARTIAL' for r in rr),total_evidence_loci=len(rr),final_exhaustive_count='NA'))
table('03_DTC_DIC_species_evidence_summary.tsv',summary,list(summary[0]))
# Retain all 724 screen records, with explicit distinction between reviewed outsiders and unreviewed remainder.
emap={r['protein_id']:r for r in evidence};audit=[]
for r in queue:
 k=r['protein_id'];review=emap.get(k);status='INTEGRATED_WITH_EVIDENCE_CATEGORY' if k in selected else ('OUTSIDE_REFERENCE_DEFINED_GROUPS_RETAINED' if review else 'SCREENED_NOT_IN_BOUNDED_PHYLOGENY_RETAINED')
 audit.append(dict(r,bounded_review_status=status,integrated_family=review['phylogenetic_family'] if k in selected else 'UNASSIGNED',automatic_exclusion='NO'))
table('04_all_carrier_candidates_disposition.tsv',audit,list(audit[0]));remaining=[r for r in audit if r['protein_id'] not in selected];table('05_remaining_carrier_candidates.tsv',remaining,list(audit[0]));writefa('05_remaining_carrier_candidates.faa',{r['protein_id']:seqs[r['protein_id']] for r in remaining})
coverage=read(src/'source_prior_06_refreshed_family_coverage.tsv');need(len(coverage)==23,'Expected 23 Wang family coverage rows')
for r in coverage:
 fam=r['wang_family'];r['integrated_records']=str(sum(x['wang_family']==fam for x in ledger))
 if fam in {'DTC','DIC'}:r.update(current_status='BOUNDED_REVIEW_CLOSED_LEDGER_INTEGRATED',search_completeness='NOT_ESTABLISHED',priority='4',next_action='RETAIN_EVIDENCE_TIERS_AND_UNRESOLVED_CARRIER_QUEUE;EXHAUSTIVE_COUNTS_PENDING')
table('06_refreshed_family_coverage.tsv',coverage,list(coverage[0]));shutil.copyfile(src/'source_prior_07_final_Wang_counts_PENDING.tsv',out/'07_final_Wang_counts_PENDING.tsv')
unresolved=[r for r in coverage if int(r['integrated_records'])==0];table('08_families_without_integrated_ledger.tsv',unresolved,list(coverage[0]))
# Cross-family overlap flags are for inspection, never automatic subtraction.
bylocus=defaultdict(list)
for r in ledger:bylocus[r['species'],r['gene_id']].append(r)
conflicts=[dict(species=sp,gene_id=g,families=';'.join(sorted({r['wang_family'] for r in rr})),count_tiers=';'.join(sorted({r['count_tier'] for r in rr}))) for (sp,g),rr in bylocus.items() if len({r['wang_family'] for r in rr})>1];table('09_cross_family_overlap_review.tsv',conflicts,['species','gene_id','families','count_tiers'])
for n in ['01_phylogenetic_group_support.tsv','02_candidate_phylogenetic_evidence.tsv','05_repeat_overlap_review.tsv','08_protein_CDS_consistency.tsv','09_annotation_source_provenance.tsv']:shutil.copyfile(src/n,out/('evidence_'+n))
(out/'README.txt').write_text('Bounded DTC/DIC review closed. DTC: 41 supported, one completeness caveat, one separately retained partial association. DIC:17 supported; two PO four-hit architectures retained with reference-consistent layout note. CDS agreement is annotation consistency, not genome/read validation. No repairs or biological-function claims. Integrated evidence ledger:946 records across10 families. Remaining664 carrier candidates preserved:25 reviewed outside reference-defined groups and639 not in bounded phylogeny. None automatically excluded from all carrier families. Final Wang matrix remains unchanged/pending. Coverage distinguishes absent integrated ledgers from zero genes; previously integrated families may still have unresolved categories or unknown search completeness. No searches or trees rerun.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step54_manifest_sha256=sha(src/'SHA256SUMS.txt'),step51_manifest_sha256=sha(screen/'SHA256SUMS.txt'),prior_ledger_sha256=sha(src/'source_prior_01_integrated_locus_evidence.tsv')),indent=2)+'\n');qc=dict(status='DTC_DIC_BOUNDED_REVIEW_CLOSED',integrated_records=len(ledger),integrated_families=len({r['wang_family'] for r in ledger}),new_family_tiers={f+'__'+t:n for (f,t),n in counts.items()},remaining_carrier_candidates=len(remaining),remaining_status_counts=dict(Counter(r['bounded_review_status'] for r in remaining)),families_without_integrated_ledger=len(unresolved),cross_family_overlap_flags=len(conflicts),final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out);bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='DTC_DIC_closed_Wang_coverage/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
