#!/usr/bin/env bash
#SBATCH --job-name=beta_at_isoform
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=beta_at_isoform_%j.out
#SBATCH --error=beta_at_isoform_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
command -v python3 >/dev/null
command -v mafft >/dev/null
command -v hmmsearch >/dev/null
python3 - <<'PY'
import csv,hashlib,json,os,subprocess,tarfile
from pathlib import Path
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
def release(env,step,required):
 if os.environ.get(env):choices=[Path(os.environ[env]).resolve(strict=True)]
 else:choices=[p for p in (base/step).glob('release*') if p.is_dir() and all((p/f).is_file() and (p/f).stat().st_size for f in required)]
 if len(choices)!=1:raise SystemExit('Set '+env+' to one explicit completed '+step+' release; found '+str(len(choices)))
 p=choices[0]
 for f in required:
  if not (p/f).is_file() or not (p/f).stat().st_size:raise SystemExit('Missing '+env+' input: '+f)
 return p
ca=release('CA_RELEASE','step25kz17.7',['BETA_CA/profile.hmm','ALPHA_CA/profile.hmm','BETA_CA/references.faa','02_CA_family_specific_review.tsv'])
model=release('MODEL_RELEASE','step25kz17.10',['AT/03_protein_isoforms.fasta','COMPLETE.txt'])
res=release('RESIDUE_RELEASE','step25kz17.8',['01_BETA_CA_residue_completeness_review.tsv','COMPLETE.txt'])
out=base/'step25kz17.11'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
out.mkdir(parents=True,exist_ok=False)
threads=str(int(os.environ.get('SLURM_CPUS_PER_TASK','4')))
def fasta(p):
 result={};k=None;ss=[]
 for l in p.read_text().splitlines():
  if l.startswith('>'):
   if k is not None:result[k]=''.join(ss).upper()
   k=l[1:].split()[0]
   if k in result:raise SystemExit('Duplicate FASTA ID')
   ss=[]
  elif l.strip():ss.append(l.strip())
 if k is not None:result[k]=''.join(ss).upper()
 return result
def writefa(p,seqs):p.write_text(''.join('>'+k+'\n'+s+'\n' for k,s in seqs.items()))
def rows(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,rr,fields):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def run(args,name):
 with (out/(name+'.txt')).open('w') as stdout,(out/(name+'.err')).open('w') as stderr:subprocess.run(args,stdout=stdout,stderr=stderr,check=True,cwd=out)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
proteins=fasta(model/'AT/03_protein_isoforms.fasta')
expected={'XP_057541773.1':501,'XP_057541774.1':501,'XP_057541775.1':278}
if set(proteins)!=set(expected) or any(len(proteins[k])!=n for k,n in expected.items()):raise SystemExit('Unexpected AT isoform IDs or lengths; inspect source')
if any('*' in s for s in proteins.values()):raise SystemExit('Unexpected stop symbol')
writefa(out/'01_AT_isoforms.faa',proteins)
refs=fasta(ca/'BETA_CA/references.faa')
allow={'P27140','P42737','Q94CE3','Q94CE4','Q9C6F5','Q9ZUC2'}
if len(refs)!=6 or {k.split('|')[1] for k in refs}!=allow:raise SystemExit('Unexpected beta reference set')
domains=[];full={}
for family in ['BETA_CA','ALPHA_CA']:
 # Relative copied input names keep report headers portable.
 (out/(family+'.hmm')).write_bytes((ca/family/'profile.hmm').read_bytes())
 run(['hmmsearch','--cpu',threads,'--noali','-E','10','--domE','10','--tblout',family+'.tbl','--domtblout',family+'.domtbl',family+'.hmm','01_AT_isoforms.faa'],family+'_search')
 for l in (out/(family+'.tbl')).read_text().splitlines():
  if not l or l.startswith('#'):continue
  x=l.split();full[(x[0],family)]=float(x[4])
 for l in (out/(family+'.domtbl')).read_text().splitlines():
  if not l or l.startswith('#'):continue
  x=l.split();domains.append(dict(protein_id=x[0],family=family,full_evalue=float(x[6]),independent_evalue=float(x[12]),model_coverage=(int(x[16])-int(x[15])+1)/int(x[5]),hmm_from=int(x[15]),hmm_to=int(x[16]),sequence_from=int(x[17]),sequence_to=int(x[18]),significant=float(x[6])<=1e-5 and float(x[12])<=1e-5))
site=[]
for pid,s in proteins.items():
 inp=out/(pid+'_with_references.faa');writefa(inp,{**refs,pid:s})
 with (out/(pid+'_alignment.faa')).open('w') as f,(out/(pid+'_mafft.log')).open('w') as e:
  subprocess.run(['mafft','--thread',threads,'--localpair','--maxiterate','1000',inp.name],cwd=out,stdout=f,stderr=e,check=True)
 aln=fasta(out/(pid+'_alignment.faa'));anchor=next(k for k in refs if '|P27140|' in k)
 cols={};pos=0
 for i,c in enumerate(aln[anchor]):
  if c!='-':pos+=1
  if c!='-' and pos in [167,169,171,227,230]:cols[pos]=i
 if len(cols)!=5:raise SystemExit('Anchor mapping incomplete')
 consistent=all(''.join(aln[k][cols[p]] for p in [167,169,171,227,230])=='CDRHC' for k in refs)
 residues=''.join(aln[pid][cols[p]] for p in [167,169,171,227,230])
 b=[d for d in domains if d['protein_id']==pid and d['family']=='BETA_CA' and d['significant']]
 a=[d for d in domains if d['protein_id']==pid and d['family']=='ALPHA_CA' and d['significant']]
 status='REVIEW_REQUIRED'
 if consistent and residues=='CDRHC' and len(b)==1 and b[0]['model_coverage']>=.5 and not a:status='SUPPORTED_SINGLE_DOMAIN_ISOFORM_REVIEW'
 site.append(dict(protein_id=pid,locus='AT__gene-LOC130820426',length_aa=len(s),beta_full_evalue=full.get((pid,'BETA_CA'),'NA'),alpha_full_evalue=full.get((pid,'ALPHA_CA'),'NA'),significant_beta_domains=len(b),significant_alpha_domains=len(a),best_beta_domain_coverage=max([d['model_coverage'] for d in b],default=0),reference_mapping_consistent=consistent,aligned_residues=residues,diagnostic_status=status,representative_replacement='NOT_PERFORMED'))
table(out/'02_AT_domain_evidence.tsv',domains,['protein_id','family','full_evalue','independent_evalue','model_coverage','hmm_from','hmm_to','sequence_from','sequence_to','significant'])
table(out/'03_AT_isoform_review.tsv',site,list(site[0]))
# Keep every screened candidate, including alpha-only and partial hits.
prior=rows(ca/'02_CA_family_specific_review.tsv');residue={r['protein_id']:r for r in rows(res/'01_BETA_CA_residue_completeness_review.tsv')}
if len({r['protein_id'] for r in prior})!=len(prior):raise SystemExit('Duplicate candidate IDs')
short=next(r for r in site if r['protein_id']=='XP_057541775.1')
for r in prior:
 pid=r['protein_id'];rr=residue.get(pid,{})
 r['residue_review_anchor']=rr.get('anchor_residues','NOT_ASSESSED');r['residue_review_recommendation']=rr.get('review_recommendation','NOT_ASSESSED')
 proposal='REVIEW_EXISTING_EVIDENCE'
 if r['screen_class'].startswith('ALPHA_'):proposal='KEEP_OUT_OF_BETA_SET_ALPHA_EVIDENCE'
 elif r['screen_class'].startswith('BETA_PARTIAL'):proposal='HOLD_PARTIAL_BETA_CANDIDATE'
 elif r['screen_class'].startswith('BETA_SUPPORTED'):proposal='REVIEW_FOR_BETA_INCLUSION'
 if pid in ['SP_Hex__gene_Spor0140320.1','SP_Hex__gene_Spor0527280.1']:proposal='HOLD_UNRESOLVED_RESIDUES_MAIN_TREE'
 if pid=='AT__gene-LOC130820426':proposal='REVIEW_ONE_LOCUS_TREE_REPRESENTATION;SHORT_ISOFORM='+short['diagnostic_status']
 if pid=='SLP_Hap1__SP3G00376':proposal='RETAIN_CANDIDATE_REVIEW_INTERNAL_GAP'
 if pid=='PO__evm.TU.LG15.1221':proposal='RETAIN_CANDIDATE_SHORT_LENGTH_NOT_SUFFICIENT_TO_EXCLUDE'
 r['proposed_disposition']=proposal;r['final_reviewer_decision']='PENDING';r['final_reviewer_reason']='';r['selected_tree_sequence']='PENDING';r['formal_status']='REVIEW_ONLY_NOT_EXPORTED'
table(out/'04_consolidated_CA_review.tsv',prior,list(prior[0]))
inputs=[ca/'BETA_CA/profile.hmm',ca/'ALPHA_CA/profile.hmm',ca/'BETA_CA/references.faa',ca/'02_CA_family_specific_review.tsv',model/'AT/03_protein_isoforms.fasta',res/'01_BETA_CA_residue_completeness_review.tsv']
(out/'QC.json').write_text(json.dumps({'status':'COMPLETE_REVIEW_ONLY','isoforms':len(site),'candidate_rows':len(prior),'long_isoforms_identical':proteins['XP_057541773.1']==proteins['XP_057541774.1'],'input_hashes':{str(p.relative_to(root)) if p.is_relative_to(root) else p.name:sha(p) for p in inputs},'threshold_note':'Full and independent domain E <=1e-5; coverage >=0.5 for single-domain diagnostic support. E-values use the three-isoform search space and are not directly comparable to full-proteome E-values.','limitations':['No final acceptance, representative replacement or tree export.','Motif and HMM evidence do not demonstrate enzyme activity.','Long tandem-domain proteins have alignment ambiguity; use prior domain-specific checks.']},indent=2)+'\n')
(out/'COMPLETE.txt').write_text('REVIEW_ONLY\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.iterdir()):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+p.name+'\n')
bundle=out/'BETA_CA_AT_isoform_review_bundle.tar.gz'
with tarfile.open(bundle,'w:gz') as t:
 for p in sorted(out.iterdir()):
  if p.is_file() and p!=bundle:t.add(p,arcname='BETA_CA_AT_isoform_review/'+p.name)
print('Bundle relative to PROJECT_ROOT:',bundle.relative_to(root))
PY
