#!/usr/bin/env bash
#SBATCH --job-name=carrier_align
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=carrier_align_%j.out
#SBATCH --error=carrier_align_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v mafft
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

src=release('step25kz17.51','CARRIER_SCREEN_RELEASE');queue=read(src/'11_phylogeny_priority_queue.tsv');proteins,_=fasta(src/'11_phylogeny_priority_candidates.faa');refs,_=fasta(src/'02_reference_panel.faa');need(len(queue)==len(proteins)==85 and len(refs)==43,'Unexpected review set');need(set(proteins)=={r['protein_id'] for r in queue},'Queue identity mismatch')
meta=read(src/'02_reference_annotations.tsv');need(set(refs)=={r['reference_id'] for r in meta},'Reference metadata mismatch')
for r in meta:need(hashlib.sha256(refs[r['reference_id']].encode()).hexdigest()==r['sequence_sha256'],'Reference hash mismatch')
allseq={**proteins,**refs};need(len(allseq)==128 and all(valid(s) for s in allseq.values()),'Invalid review sequences')
out=base/'step25kz17.52'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
writefa('01_all_sequences.faa',allseq);commands=[];cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
def run(args,log,stdout=None):
 commands.append([str(x).replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for x in args])
 with (out/log).open('w') as f:
  if stdout:
   with (out/stdout).open('w') as a:subprocess.run(args,stdout=a,stderr=f,check=True)
  else:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
run(['mafft','--amino','--thread',cpu,'--localpair','--maxiterate','1000',str(out/'01_all_sequences.faa')],'02_mafft.log','02_full_alignment.faa')
aln,_=fasta(out/'02_full_alignment.faa');need(set(aln)==set(allseq) and len({len(v) for v in aln.values()})==1,'Invalid alignment');need(all(aln[k].replace('-','')==v for k,v in allseq.items()),'Alignment changed sequences')
hits=[r for r in read(src/'06_all_domain_hits.tsv') if r['sequence_id'] in allseq];dh=defaultdict(list)
for r in hits:dh[r['sequence_id']].append(r)
# Map source Pfam residue intervals into alignment coordinates, retaining every hit.
rescols={k:[i for i,c in enumerate(v) if c!='-'] for k,v in aln.items()};mapped=[];carriercols={};arch=[]
for k,v in allseq.items():
 car=[r for r in dh[k] if r['accession'].split('.')[0]=='PF00153'];covered=set();intervals=[]
 for h in dh[k]:
  a,b=int(h['ali_from']),int(h['ali_to']);need(1<=a<=b<=len(v),'Invalid Pfam coordinates');mapped.append(dict(h,alignment_from=rescols[k][a-1]+1,alignment_to=rescols[k][b-1]+1))
  if h in car:covered.update(range(a,b+1));intervals.append((a,b))
 carriercols[k]={rescols[k][j-1] for j in covered};ordered=sorted(intervals);overlaps=sum(max(0,min(b,d)-max(a,c)+1) for i,(a,b) in enumerate(ordered) for c,d in ordered[i+1:]);extra=sorted({h['accession'] for h in dh[k] if h not in car})
 arch.append(dict(sequence_id=k,role='CANDIDATE' if k in proteins else 'REFERENCE',length_aa=len(v),carrier_hits=len(car),carrier_union_residues=len(covered),carrier_union_fraction=round(len(covered)/len(v),6),carrier_intervals=';'.join(f'{a}-{b}' for a,b in ordered),pairwise_overlap_residues=overlaps,N_terminal_outside_carrier=min(covered)-1 if covered else 'NA',C_terminal_outside_carrier=len(v)-max(covered) if covered else 'NA',additional_Pfam_models=';'.join(extra),interpretation='DESCRIPTIVE_DOMAIN_ARCHITECTURE_NOT_FUNCTION_OR_COMPLETENESS_PROOF'))
table('03_domains_in_alignment.tsv',mapped,list(hits[0])+['alignment_from','alignment_to']);table('04_domain_architecture_review.tsv',arch,list(arch[0]))
targets=['REF__Q9C5M0','REF__Q9SJY5','REF__Q9SB52','REF__Q9FY68'];need(all(carriercols[k] for k in targets),'Missing target carrier domains')
metrics=[]
for k in proteins:
 for ref in targets:
  cc=sorted(carriercols[ref]);present=[i for i in cc if aln[k][i]!='-'];identity=sum(aln[k][i]==aln[ref][i] for i in present);metrics.append(dict(protein_id=k,reference_id=ref,reference_carrier_residues=len(cc),paired_residues=len(present),reference_carrier_occupancy=round(len(present)/len(cc),6),identity_over_paired=round(identity/len(present),6) if present else 0,candidate_domain_supported_pairs=sum(i in carriercols[k] for i in present),interpretation='ALIGNMENT_REVIEW_METRIC_NOT_ACCEPTANCE_THRESHOLD'))
table('05_target_carrier_alignment_metrics.tsv',metrics,list(metrics[0]))
# Proposed tree mask only. Retain the complete alignment and explicit column decisions.
mask=[];keep=[]
for i in range(len(next(iter(aln.values())))):
 target_n=sum(i in carriercols[k] for k in targets);ref_n=sum(i in carriercols[k] for k in refs);occ=sum(aln[k][i]!='-' for k in proteins)/len(proteins);retained=target_n>=2 and occ>=.5
 if retained:keep.append(i)
 mask.append(dict(alignment_column=i+1,target_reference_domain_support=target_n,all_reference_domain_support=ref_n,candidate_occupancy=round(occ,6),proposed_core=int(retained)))
table('06_proposed_core_mask.tsv',mask,list(mask[0]));writefa('07_proposed_carrier_core_REVIEW_ONLY.faa',{k:''.join(v[i] for i in keep) for k,v in aln.items()})
coreqc=[]
for r in arch:
 k=r['sequence_id'];n=sum(aln[k][i]!='-' for i in keep);flags=[]
 if not keep or n/len(keep)<.7:flags.append('LOW_PROPOSED_CORE_OCCUPANCY')
 if r['carrier_hits']!=3:flags.append('CARRIER_HIT_COUNT_NOT_THREE')
 if r['additional_Pfam_models']:flags.append('ADDITIONAL_PFAM_MODELS')
 if r['length_aa']>600:flags.append('LONG_SEQUENCE')
 if r['length_aa']<200:flags.append('SHORT_SEQUENCE')
 coreqc.append(dict(sequence_id=k,role=r['role'],length_aa=r['length_aa'],core_residues=n,core_columns=len(keep),core_occupancy=round(n/len(keep),6) if keep else 0,flags=';'.join(flags) or 'NONE',decision='REVIEW_ONLY'))
table('08_core_occupancy_and_flags.tsv',coreqc,list(coreqc[0]));focus={'MC__gene_14767','MC__gene_5144','SS__g18155.t1'};need(focus<=set(proteins),'Expected priority exceptions missing');writefa('09_focus_exception_sequences.faa',{k:proteins[k] for k in sorted(focus)})
table('09_focus_exception_architecture.tsv',[r for r in arch if r['sequence_id'] in focus],list(arch[0]))
for n in ['02_reference_annotations.tsv','09_all_candidates_specificity_REVIEW_ONLY.tsv','10_target_competitor_comparisons.tsv','11_phylogeny_priority_queue.tsv','prior_01_integrated_locus_evidence.tsv','prior_06_refreshed_family_coverage.tsv','prior_07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/('source_'+n))
run(['mafft','--version'],'mafft_version.txt')
(out/'README.txt').write_text('Bounded review of 85 prioritized candidates plus all 43 carrier references. Full proteins aligned by MAFFT L-INS-i; unchanged sequence identities verified. Existing Pfam hits mapped from protein residues to alignment columns without rerunning domain searches. Every domain retained; overlap and repeat counts are descriptive. Candidate occupancy/identity measured separately against carrier-domain residues of each of four DTC/DIC anchors. Proposed tree core retains columns inside PF00153 alignments of >=2 target anchors and >=50% candidate residue occupancy. This is an explicit review mask, not a validated phylogenetic dataset. All 128 sequences retained even if low occupancy; inspect core QC before tree construction. No trimming of source proteins, annotation repair, functional assignment, trees or count updates. Long proteins and short fragments explicitly retained, including MC gene_14767, MC gene_5144 and SS g18155.t1. Alignment of divergent/repeated domains requires visual review; shared domain support does not establish specificity.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),commands=commands),indent=2)+'\n')
qc=dict(status='CARRIER_ALIGNMENT_ARCHITECTURE_REVIEW_READY',candidates=85,references=43,alignment_columns=len(next(iter(aln.values()))),proposed_core_columns=len(keep),core_requires_review=True,flagged_candidates=sum(r['role']=='CANDIDATE' and r['flags']!='NONE' for r in coreqc),low_core_occupancy_candidates=sum(r['role']=='CANDIDATE' and 'LOW_PROPOSED_CORE_OCCUPANCY' in r['flags'] for r in coreqc),final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
for p in out.glob('*.log'):p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Carrier_alignment_architecture_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
