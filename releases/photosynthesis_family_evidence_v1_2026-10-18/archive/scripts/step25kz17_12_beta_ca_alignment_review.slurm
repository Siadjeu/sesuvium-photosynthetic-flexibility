#!/usr/bin/env bash
#SBATCH --job-name=beta_ca_align
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=beta_ca_align_%j.out
#SBATCH --error=beta_ca_align_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT from the project directory}"
command -v python3 >/dev/null
command -v mafft >/dev/null
python3 - <<'PY'
import csv,hashlib,json,os,subprocess,tarfile,collections
from pathlib import Path
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
def release(env,step,files):
 choices=[Path(os.environ[env]).resolve(strict=True)] if os.environ.get(env) else [p for p in (base/step).glob('release*') if p.is_dir() and all((p/f).is_file() and (p/f).stat().st_size for f in files)]
 if len(choices)!=1:raise SystemExit('Set '+env+' to one explicit '+step+' release; found '+str(len(choices)))
 p=choices[0]
 for f in files:
  if not (p/f).is_file() or not (p/f).stat().st_size:raise SystemExit('Missing '+env+' input: '+f)
 return p
ca_files=['05_all_screened_candidates.faa','BETA_CA/references.faa','03_domain_evidence.tsv']
iso_files=['01_AT_isoforms.faa','03_AT_isoform_review.tsv','04_consolidated_CA_review.tsv','QC.json','COMPLETE.txt']
ca=release('CA_RELEASE','step25kz17.7',ca_files)
iso=release('ISOFORM_RELEASE','step25kz17.11',iso_files)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
# Confirm the selected CA source agrees with the source used in Step17.11.
q=json.loads((iso/'QC.json').read_text())
for suffix in ['BETA_CA/references.faa','02_CA_family_specific_review.tsv']:
 values=[h for p,h in q['input_hashes'].items() if p.endswith('/'+suffix)]
 if len(values)!=1 or sha(ca/suffix)!=values[0]:raise SystemExit('CA_RELEASE provenance mismatch: '+suffix)
def fasta(p):
 d={};key=None;ss=[]
 for l in p.read_text().splitlines():
  if l.startswith('>'):
   if key is not None:d[key]=''.join(ss).upper()
   key=l[1:].split()[0]
   if key in d:raise SystemExit('Duplicate FASTA identifier')
   ss=[]
  elif l.strip():ss.append(l.strip())
 if key is not None:d[key]=''.join(ss).upper()
 if any(not s for s in d.values()):raise SystemExit('Empty FASTA sequence')
 return d
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,rr,fields):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(p,d):p.write_text(''.join('>'+k+'\n'+s+'\n' for k,s in d.items()))
allseq=fasta(ca/ca_files[0]);refs=fasta(ca/ca_files[1]);isoseq=fasta(iso/'01_AT_isoforms.faa');review=read(iso/'04_consolidated_CA_review.tsv')
allow={'P27140','P42737','Q94CE3','Q94CE4','Q9C6F5','Q9ZUC2'}
if len(refs)!=6 or {k.split('|')[1] for k in refs}!=allow:raise SystemExit('Unexpected reference set')
if len(review)!=110 or len({r['protein_id'] for r in review})!=110:raise SystemExit('Unexpected candidate table')
# This is a proposed alignment selection, not final biological acceptance.
if any(r['final_reviewer_decision']!='PENDING' for r in review):raise SystemExit('Reviewer decisions have changed: reconcile selection policy before running')
hold={'SP_Hex__gene_Spor0140320.1','SP_Hex__gene_Spor0527280.1'}
at='AT__gene-LOC130820426';short='XP_057541775.1';at_tip=at+'__isoform_'+short
check=[r for r in read(iso/'03_AT_isoform_review.tsv') if r['protein_id']==short]
if len(check)!=1 or check[0]['diagnostic_status']!='SUPPORTED_SINGLE_DOMAIN_ISOFORM_REVIEW':raise SystemExit('Short AT isoform did not pass preceding diagnostic review')
selected=[r for r in review if r['screen_class']=='BETA_SUPPORTED_REVIEW' and r['protein_id'] not in hold]
if len(selected)!=45 or at not in {r['protein_id'] for r in selected}:raise SystemExit('Unexpected proposed locus set')
out=base/'step25kz17.12'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
out.mkdir(parents=True,exist_ok=False)
main={};mapping=[];selection=[]
selected_ids={r['protein_id'] for r in selected}
for r in review:
 pid=r['protein_id'];include=pid in selected_ids
 reason='SUPPORTED_BETA_PROPOSED_ALIGNMENT' if include else 'OUTSIDE_MAIN_ALIGNMENT_'+r['screen_class']
 if pid in hold:reason='HOLD_UNRESOLVED_RESIDUES'
 if pid==at:reason='DOCUMENTED_SHORT_ISOFORM_FOR_SINGLE_DOMAIN_ALIGNMENT'
 tip=at_tip if pid==at else pid
 if include:
  seq=isoseq[short] if pid==at else allseq[pid]
  if any(c not in 'ACDEFGHIKLMNPQRSTVWYBXZJUO' for c in seq):raise SystemExit('Unexpected sequence symbol: '+pid)
  main[tip]=seq;mapping.append(dict(tip_id=tip,locus_id=pid,source_sequence_id=short if pid==at else pid,role='CANDIDATE_LOCUS',length_aa=len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()))
 selection.append(dict(locus_id=pid,screen_class=r['screen_class'],alignment_selection='PROPOSED_INCLUDE' if include else 'HOLD_OUTSIDE_MAIN',reason=reason,tip_id=tip if include else '',final_biological_decision='PENDING'))
for k,s in refs.items():
 if k in main:raise SystemExit('Reference ID collision')
 main[k]=s;mapping.append(dict(tip_id=k,locus_id='',source_sequence_id=k,role='REFERENCE',length_aa=len(s),sequence_sha256=hashlib.sha256(s.encode()).hexdigest()))
assert len(main)==51
sensitivity={k:s for k,s in main.items() if k!=at_tip};sensmap=[dict(r) for r in mapping if r['tip_id']!=at_tip]
dom=sorted([r for r in read(ca/'03_domain_evidence.tsv') if r['protein_id']==at and r['family']=='BETA_CA' and float(r['independent_evalue'])<=1e-5],key=lambda r:int(r['alignment_from']))
if len(dom)!=2:raise SystemExit('Expected two AT domains')
fragments=[]
for i,r in enumerate(dom,1):
 a,b=int(r['alignment_from']),int(r['alignment_to']);tip=at+'__domain_'+str(i)
 if not 1<=a<=b<=len(allseq[at]):raise SystemExit('Invalid domain coordinates')
 seq=allseq[at][a-1:b];sensitivity[tip]=seq
 sensmap.append(dict(tip_id=tip,locus_id=at,source_sequence_id=at,role='DOMAIN_FRAGMENT_NOT_SEPARATE_GENE',length_aa=len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()))
 fragments.append(dict(tip_id=tip,locus_id=at,start_aa_1based=a,end_aa_inclusive=b,interpretation='DIAGNOSTIC_DOMAIN_NOT_NEW_LOCUS'))
assert len(sensitivity)==52
table(out/'01_selection_manifest.tsv',selection,list(selection[0]))
table(out/'02_main_tip_map.tsv',mapping,list(mapping[0]))
table(out/'03_sensitivity_tip_map.tsv',sensmap,list(sensmap[0]))
table(out/'04_AT_domain_coordinates.tsv',fragments,list(fragments[0]))
positions=[167,169,171,227,230];anchor=next(k for k in refs if '|P27140|' in k)
summary={};threads=str(int(os.environ.get('SLURM_CPUS_PER_TASK','4')))
for label,seqs,mp in [('main',main,mapping),('domain_sensitivity',sensitivity,sensmap)]:
 d=out/label;d.mkdir();writefa(d/'input.faa',seqs)
 with (d/'alignment.faa').open('w') as f,(d/'mafft.log').open('w') as e:
  subprocess.run(['mafft','--thread',threads,'--localpair','--maxiterate','1000','input.faa'],cwd=d,stdout=f,stderr=e,check=True)
 aln=fasta(d/'alignment.faa')
 if set(aln)!=set(seqs) or len({len(s) for s in aln.values()})!=1:raise SystemExit('Invalid alignment identifiers or lengths')
 if any(aln[k].replace('-','')!=s for k,s in seqs.items()):raise SystemExit('Alignment changed residues')
 length=len(aln[anchor]);cols={};pos=0
 for j,c in enumerate(aln[anchor]):
  if c!='-':
   pos+=1
   if pos in positions:cols[pos]=j
 if set(cols)!=set(positions):raise SystemExit('Incomplete reference coordinate map')
 refcheck=all(''.join(aln[k][cols[p]] for p in positions)=='CDRHC' for k in refs)
 core=[j for j in range(length) if sum(aln[k][j]!='-' for k in refs)>=5]
 if not core:raise SystemExit('Empty diagnostic core')
 metrics=[];gaps=[];motif_windows=[];bytip={r['tip_id']:r for r in mp}
 for k,s in aln.items():
  cc=''.join(s[j] for j in core);nongap=[i for i,c in enumerate(cc) if c!='-'];start=nongap[0] if nongap else len(core);end=nongap[-1] if nongap else -1
  runs=[];run=[]
  for i,c in enumerate(cc):
   if c=='-':run.append(i)
   elif run:runs.append(run);run=[]
  if run:runs.append(run)
  internal=[run for run in runs if run[0]>start and run[-1]<end]
  for run in internal:
   if len(run)>=10:gaps.append(dict(tip_id=k,core_gap_length=len(run),alignment_start_1based=core[run[0]]+1,alignment_end_1based=core[run[-1]]+1))
  residues=''.join(s[cols[p]] for p in positions)
  metrics.append(dict(tip_id=k,locus_id=bytip[k]['locus_id'],role=bytip[k]['role'],length_aa=len(seqs[k]),anchor_residues=residues,reference_mapping_consistent=refcheck,core_columns=len(core),covered_core_columns=len(nongap),core_coverage_fraction=round(len(nongap)/len(core),6),missing_core_N=start,missing_core_C=len(core)-1-end if nongap else len(core),longest_internal_core_gap=max([len(run) for run in internal],default=0),review_status='INSPECT' if not refcheck or residues!='CDRHC' or any(len(run)>=10 for run in internal) else 'RESIDUES_PASS_REVIEW_ALIGNMENT'))
  for p in positions:
   j=cols[p];motif_windows.append(dict(tip_id=k,reference_position=p,alignment_column_1based=j+1,residue=s[j],alignment_window_5_columns_each_side=s[max(0,j-5):min(length,j+6)]))
 table(d/'sequence_QC.tsv',metrics,list(metrics[0]))
 table(d/'internal_core_gaps.tsv',gaps,['tip_id','core_gap_length','alignment_start_1based','alignment_end_1based'])
 table(d/'conserved_site_windows.tsv',motif_windows,list(motif_windows[0]))
 table(d/'column_occupancy.tsv',[dict(alignment_column_1based=j+1,reference_non_gap=sum(aln[k][j]!='-' for k in refs),all_non_gap=sum(s[j]!='-' for s in aln.values()),diagnostic_core=j in set(core)) for j in range(length)],['alignment_column_1based','reference_non_gap','all_non_gap','diagnostic_core'])
 # Diagnostic core is exported for inspection only; it is not a validated tree trim.
 writefa(d/'diagnostic_reference_core.faa',{k:''.join(s[j] for j in core) for k,s in aln.items()})
 groups=collections.defaultdict(list)
 for k,s in seqs.items():groups[s].append(k)
 identical=[dict(group=i,tip_ids=','.join(ids),size=len(ids)) for i,ids in enumerate(groups.values(),1) if len(ids)>1]
 table(d/'identical_sequences.tsv',identical,['group','tip_ids','size'])
 summary[label]=dict(tips=len(seqs),alignment_columns=length,diagnostic_core_columns=len(core),reference_sites_consistent=refcheck,inspection_flags=sum(r['review_status']=='INSPECT' for r in metrics))
 table(d/'SLP_gap_focus.tsv',[r for r in metrics if r['tip_id']=='SLP_Hap1__SP3G00376'],list(metrics[0]))
(out/'README.txt').write_text('PROPOSED ALIGNMENTS FOR REVIEW; no tree inference or final family count.\nMain: 45 candidate loci + 6 references. AT uses XP_057541775.1 only here; primary representative is unchanged.\nSensitivity: 44 other loci + two AT domain fragments + 6 references. The two AT fragments remain ONE locus; boundaries may overlap.\nIdentical sequences retain their locus IDs. Do not collapse biological counts automatically.\nDiagnostic core keeps columns occupied by at least 5/6 references; this is not an approved phylogenetic trimming rule. Terminal gaps alone do not establish truncation.\nInspect sequence_QC, conserved_site_windows and internal_core_gaps, especially SLP_Hap1__SP3G00376.\nSelection manifests distinguish phylogenetic inclusion from biological gene acceptance.\n')
inputs=[ca/f for f in ca_files]+[ca/'02_CA_family_specific_review.tsv']+[iso/f for f in iso_files]
(out/'QC.json').write_text(json.dumps({'status':'COMPLETE_ALIGNMENT_REVIEW_ONLY','proposed_loci':45,'species_locus_counts':dict(collections.Counter(r['species'] for r in selected)),'alignments':summary,'inputs_sha256':{str(p.relative_to(root)) if p.is_relative_to(root) else p.name:sha(p) for p in inputs}},indent=2)+'\n')
(out/'COMPLETE.txt').write_text('ALIGNMENT_REVIEW_ONLY\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+str(p.relative_to(out))+'\n')
bundle=out/'BETA_CA_alignment_review_bundle.tar.gz'
with tarfile.open(bundle,'w:gz') as t:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p!=bundle:t.add(p,arcname='BETA_CA_alignment_review/'+str(p.relative_to(out)))
print('Bundle relative to PROJECT_ROOT:',bundle.relative_to(root))
PY
