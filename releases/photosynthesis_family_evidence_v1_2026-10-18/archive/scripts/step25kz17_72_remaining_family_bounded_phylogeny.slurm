#!/usr/bin/env bash
#SBATCH --job-name=remaining_batch_trees
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=remaining_batch_trees_%j.out
#SBATCH --error=remaining_batch_trees_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v mafft
command -v iqtree2 || command -v iqtree
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











import signal
src=release('step25kz17.71','REMAINING_DOMAIN_RELEASE')
proteins,_=fasta(src/'03_candidate_proteins.faa');refs,_=fasta(src/'03_reference_panel.faa');roles=read(src/'03_reference_roles.tsv');reviews=read(src/'09_combined_domain_review.tsv')
seqmap=read(src/'03_domain_sequence_map.tsv');domains=defaultdict(list)
for x in read(src/'07_Pfam_domain_hits.tsv'):
 for f in ['model_length','hmm_from','hmm_to','ali_from','ali_to']:x[f]=int(x[f])
 x['score']=float(x['score']);domains[x['sequence_id']].append(x)
lookup={}
for x in seqmap:
 kind,k=x['kind'],x['original_id'];s=(proteins if kind=='CANDIDATE' else refs)[k]
 need(hashlib.sha256(s.encode()).hexdigest()==x['sequence_sha256'],'Sequence map mismatch '+k);lookup[(kind,k)]=x['domain_sequence_id']
out=base/'step25kz17.72'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
def table(path,rows,fields):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def writefa(p,d):p.write_text(''.join('>'+k+'\n'+s+'\n' for k,s in sorted(d.items())))
cpu=os.environ.get('SLURM_CPUS_PER_TASK','4');iq=shutil.which('iqtree2') or shutil.which('iqtree');need(iq,'IQ-TREE unavailable')
commands=[];started=time.monotonic();deadline=started+6000
# Six alignment groups; four specificity trees, two limited-panel alignments.
specs=[('PEPCK',['PEPCK'],['PF01293']),('BASS',['BASS2'],['PF01758']),('Malic_enzyme',['NAD-ME','NADP-ME'],['PF00390','PF03949']),('MDH_LDH',['MDH'],['PF00056','PF02866']),('Aminotransferase',['ASPAT','AlaAT'],['PF00155']),('MPC',['MPC'],['PF03650'])]
allmembership=[];summaries=[];exception_rows=[]
def run(cmd,stdout_path,stderr_path,limit):
 remaining=deadline-time.monotonic()
 if remaining<30:raise TimeoutError('TOTAL_TIME_BUDGET_EXHAUSTED')
 commands.append(cmd)
 with stdout_path.open('w') as a,stderr_path.open('w') as b:
  proc=subprocess.Popen(cmd,stdout=a,stderr=b,start_new_session=True)
  try:rc=proc.wait(timeout=min(limit,remaining))
  except subprocess.TimeoutExpired:
   os.killpg(proc.pid,signal.SIGTERM)
   try:proc.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
   raise TimeoutError('COMMAND_TIME_BUDGET_EXHAUSTED')
  if rc:raise RuntimeError('COMMAND_EXIT_'+str(rc))
for group,fams,models in specs:
 d=out/group;d.mkdir();rr=[r for r in reviews if r['wang_family'] in fams];cand=sorted({r['protein_id'] for r in rr});ref=sorted({r['reference_id'] for r in roles if r['wang_family'] in fams});rows=[];blocks={m:{} for m in models};partial={};tiprows=[];mainids=[];refids=[]
 for kind,ids,seqs,prefix in [('CANDIDATE',cand,proteins,'C'),('REFERENCE',ref,refs,'R')]:
  for i,k in enumerate(ids,1):
   tip=prefix+str(i).zfill(5);hh=domains[lookup[(kind,k)]];chosen=[];flags=[];coverage=[]
   for m in models:
    hits=[x for x in hh if x['accession']==m];best=max(hits,key=lambda x:x['score']) if hits else None
    coverage.append((best['hmm_to']-best['hmm_from']+1)/best['model_length'] if best else 0)
    chosen.append(best)
    if len(hits)>1:flags.append('MULTIPLE_HITS_'+m)
   ordered=all(chosen[j]['ali_to']<chosen[j+1]['ali_from'] for j in range(len(chosen)-1)) if all(chosen) else False
   if len(models)==1:ordered=bool(chosen[0])
   if not ordered:flags.append('MISSING_OR_NONCANONICAL_CORE_ORDER')
   status='MAIN_CORE' if min(coverage)>=.8 and ordered and not flags else 'PARTIAL_OR_ARCHITECTURE' if max(coverage)>=.3 else 'OUTSIDE_CORE_ALIGNMENT'
   if kind=='REFERENCE' and status!='MAIN_CORE':flags.append('REFERENCE_NOT_MAIN_COMPATIBLE')
   related=[x for x in rr if x['protein_id']==k] if kind=='CANDIDATE' else []
   sp=related[0]['species'] if related else '';gene=related[0]['gene_id'] if related else ''
   if status=='MAIN_CORE':
    for m,h in zip(models,chosen):blocks[m][tip]=seqs[k][h['ali_from']-1:h['ali_to']]
    mainids.append(tip)
    if kind=='REFERENCE':refids.append(tip)
   elif status=='PARTIAL_OR_ARCHITECTURE':partial[tip]=seqs[k]
   flags.extend(sorted({f for x in related for f in x['review_flags'].split(';') if f in {'AMBIGUOUS_RESIDUES','NONOVERLAPPING_ADDITIONAL_MODEL_HEURISTIC','LENGTH_OUTLIER_VS_BEST_POSITIVE_HEURISTIC'}}))
   row=dict(group=group,tip_id=tip,kind=kind,original_id=k,species=sp,gene_id=gene,alignment_category=status,core_coverage=';'.join(m+':'+str(round(v,5)) for m,v in zip(models,coverage)),architecture_flags=';'.join(sorted(set(flags))) or 'NONE',family_evidence=';'.join(x['wang_family']+':'+x['similarity_category'] for x in related),reference_roles=';'.join(x['wang_family']+':'+x['role']+':'+x['subgroups'] for x in roles if kind=='REFERENCE' and x['reference_id']==k))
   rows.append(row);allmembership.append(row)
   tiprows.append((tip,(k+' | '+sp+' | '+row['reference_roles']).strip(' |')))
   if kind=='CANDIDATE':exception_rows.append(dict(group=group,protein_id=k,species=sp,gene_id=gene,alignment_category=status,architecture_flags=row['architecture_flags'],specificity_status='PENDING_REFERENCE_GUIDED_REVIEW' if group!='MPC' else 'LIMITED_PANEL_NO_NON_MPC_COMPARATOR',count_status='REVIEW_ONLY'))
 table(d/'01_membership.tsv',rows,list(rows[0]));writefa(d/'02_partial_or_architecture_full_proteins.faa',partial)
 (d/'itol_labels.txt').write_text('LABELS\nSEPARATOR TAB\nDATA\n'+''.join(k+'\t'+v+'\n' for k,v in tiprows))
 q=dict(group=group,candidate_records=len(cand),main_candidates=sum(x['kind']=='CANDIDATE' and x['alignment_category']=='MAIN_CORE' for x in rows),partial_or_architecture_candidates=sum(x['kind']=='CANDIDATE' and x['alignment_category']=='PARTIAL_OR_ARCHITECTURE' for x in rows),outside_core_candidates=sum(x['kind']=='CANDIDATE' and x['alignment_category']=='OUTSIDE_CORE_ALIGNMENT' for x in rows),main_references=len(refids),core_columns=0,status='PREPARED',error='')
 print('Preparing '+group+' '+str(q),flush=True)
 try:
  need(len(mainids)<=450,'Alignment group exceeds450 tips')
  positive_refs={x['reference_id'] for x in roles if x['wang_family'] in fams and x['role']=='POSITIVE'}
  mainrefs={x['original_id'] for x in rows if x['kind']=='REFERENCE' and x['alignment_category']=='MAIN_CORE'}
  if not positive_refs<=mainrefs:raise ValueError('POSITIVE_REFERENCE_CORE_INCOMPLETE_REVIEW_REQUIRED')
  if group not in {'MPC','PEPCK'} and not any(x['reference_id'] in mainrefs and x['role']=='COMPARATOR' for x in roles if x['wang_family'] in fams):raise ValueError('NO_DOMAIN_COMPATIBLE_COMPARATOR')
  if len(mainids)<4:raise ValueError('INSUFFICIENT_MAIN_TIPS')
  joined={k:'' for k in mainids};columnrows=[]
  for m in models:
   writefa(d/(m+'_input.faa'),blocks[m])
   run(['mafft','--thread',cpu,'--threadit','0','--auto',str(d/(m+'_input.faa'))],d/(m+'_aligned.faa'),d/(m+'_mafft.log'),180)
   aln,_=fasta(d/(m+'_aligned.faa'));need(set(aln)==set(mainids) and len({len(v) for v in aln.values()})==1,'Invalid alignment')
   for k in aln:need(aln[k].replace('-','')==blocks[m][k],'MAFFT sequence mismatch')
   keep=[j for j in range(len(next(iter(aln.values())))) if sum(aln[k][j]!='-' for k in refids)>=max(1,(len(refids)+1)//2) and sum(aln[k][j]!='-' for k in mainids)>=len(mainids)*.5]
   for j in keep:columnrows.append(dict(core_column=len(columnrows)+1,model=m,model_alignment_column=j+1))
   for k in mainids:joined[k]+=''.join(aln[k][j] for j in keep)
  q['core_columns']=len(columnrows);writefa(d/'03_main_core.faa',joined);table(d/'03_core_columns.tsv',columnrows,['core_column','model','model_alignment_column'])
  if len(columnrows)<60:raise ValueError('CORE_TOO_SHORT')
  occupancy=[dict(tip_id=k,non_gap_fraction=round(sum(c!='-' for c in v)/len(v),6)) for k,v in joined.items()];table(d/'03_core_occupancy.tsv',occupancy,['tip_id','non_gap_fraction'])
  if min(x['non_gap_fraction'] for x in occupancy)<.5:raise ValueError('LOW_POST_TRIM_OCCUPANCY_REVIEW_REQUIRED')
  if group in {'MPC','PEPCK'}:q['status']='ALIGNED_LIMITED_PANEL_NO_SPECIFICITY_TREE'
  else:
   run([iq,'-s',str(d/'03_main_core.faa'),'-st','AA','-m','LG+G4','-nt',cpu,'-seed','17272','-bb','1000','-alrt','1000','-bnni','-keep-ident','-pre',str(d/'main')],d/'04_iqtree_stdout.log',d/'04_iqtree_stderr.log',900)
   need((d/'main.treefile').is_file() and (d/'main.iqtree').is_file(),'Missing completed tree outputs');tree_tips=re.findall(r'(?:\(|,)\s*([CR][0-9]{5})(?=[:),])',(d/'main.treefile').read_text());need(len(tree_tips)==len(mainids) and set(tree_tips)==set(mainids),'Tree tip set mismatch');q['status']='TREE_COMPLETE_REVIEW_ONLY'
 except (TimeoutError,RuntimeError,ValueError,SystemExit) as e:q['status']='REVIEW_BLOCKED_OR_TIME_LIMIT';q['error']=str(e)
 (d/'QC.json').write_text(json.dumps(q,indent=2)+'\n');summaries.append(q)
 table(out/'04_group_summary.tsv',summaries,list(q))
# All group statuses are explicit; timed-out trees are not completed evidence.
table(out/'01_all_membership.tsv',allmembership,list(allmembership[0]));table(out/'02_separated_review_queue.tsv',exception_rows,list(exception_rows[0]))
for n in ['source_01_integrated_locus_evidence.tsv','source_07_final_Wang_counts_PENDING.tsv','03_reference_roles.tsv','09_combined_domain_review.tsv','08_reference_domain_calibration.tsv']:shutil.copyfile(src/n,out/n)
coverage=read(src/'11_refreshed_family_coverage.tsv')
for r in coverage:
 matching=[q for q in summaries for name,fams,mm in specs if q['group']==name and r['wang_family'] in fams]
 if matching:r.update(current_status=matching[0]['status'],next_action='REVIEW_BATCH_PHYLOGENY_AND_SEPARATED_EXCEPTIONS',search_completeness='NOT_ESTABLISHED')
table(out/'11_refreshed_family_coverage.tsv',coverage,list(coverage[0]))
# Versions are inexpensive; failures remain visible.
for cmd,name in [(['mafft','--version'],'mafft_version.txt'),([iq,'--version'],'iqtree_version.txt')]:
 with (out/name).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=15,check=False)
def clean(s):return s.replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>')
for p in out.rglob('*'):
 if p.is_file() and p.suffix in {'.log','.iqtree','.txt'}:p.write_text(clean(p.read_text(errors='replace')))
(out/'PROVENANCE.json').write_text(json.dumps(dict(step71_manifest_sha256=sha(src/'SHA256SUMS.txt'),upstream=json.loads((src/'PROVENANCE.json').read_text()),commands=[[clean(x) for x in c] for c in commands],model='LG+G4',seed=17272),indent=2)+'\n')
(out/'README.txt').write_text('Bounded combined alignment/phylogeny triage. No new proteome/domain searches. Shared ME and aminotransferase pools processed once each. Main cores require single hits covering at least80% of each designated HMM and canonical nonoverlapping order. This is an alignment eligibility heuristic, not full gene completeness or family acceptance. All candidate IDs retained, including outside-core and partial/architecture categories. Domain blocks aligned separately, then concatenated after retaining columns with >=50% reference and all-tip occupancy. Low post-trim occupancy blocks tree; no adaptive threshold relaxation. Extra domains/length flags retained even for core-compatible sequences. ASPAT PF01041 and MDH PF01073 not imposed universally; BASS uses common PF01758, avoiding overlapping model double-counting. Aminotransferase PF00155 tree covers only this domain-compatible subset; other folds remain outside this tree, not excluded from family evidence. Domain-incompatible comparator references retained in membership but omitted from alignment; all positive references must be main-core compatible. PEPCK ATP-core alignment only, with GTP comparators documented outside core. MPC alignment only because panel lacks independent comparator. Unrooted IQ-TREE LG+G4 fixed-model first-pass trees,1000UFBoot+1000SH-aLRT+bnni; fixed model saves model-search time and is a limitation, not selected best fit. No automatic functional assignments or counts; tree support alone does not establish C4/CAM role. Partial sequences saved separately, not placed by best hit. No automatic sensitivity-tree reruns. Four potential trees with900s per-tree cap,180s per-domain alignment and6000s total external-work budget. Timed-out outputs are explicitly incomplete. Check group statuses even when SLURM exits0. Existing1560-record ledger unchanged. Sources: https://iqtree.github.io/doc/Command-Reference ; https://mafft.cbrc.jp/alignment/software/multithreading.html .\n')
qc=dict(status='COMBINED_TRIAGE_COMPLETE_WITH_EXPLICIT_GROUP_STATUSES',groups=summaries,completed_trees=sum(q['status']=='TREE_COMPLETE_REVIEW_ONLY' for q in summaries),unique_candidates=len({x['protein_id'] for x in exception_rows}),ledger_unchanged=sha(src/'source_01_integrated_locus_evidence.tsv')==sha(out/'source_01_integrated_locus_evidence.tsv'),final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Remaining_family_bounded_phylogeny/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle))
PY
