#!/usr/bin/env bash
#SBATCH --job-name=carrier_tree
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=carrier_tree_%j.out
#SBATCH --error=carrier_tree_%j.err
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


src=release('step25kz17.52','CARRIER_ALIGNMENT_RELEASE');seqs,_=fasta(src/'01_all_sequences.faa');queue=read(src/'source_11_phylogeny_priority_queue.tsv');refs=read(src/'source_02_reference_annotations.tsv');need(len(queue)==85 and len(refs)==43 and len(seqs)==128,'Unexpected inputs')
refids={r['reference_id'] for r in refs};qmeta={r['protein_id']:r for r in queue};need(set(seqs)==set(qmeta)|refids,'Sequence ID mismatch')
for r in refs:need(hashlib.sha256(seqs[r['reference_id']].encode()).hexdigest()==r['sequence_sha256'],'Reference hash mismatch')
iq=next((shutil.which(x) for x in ['iqtree3','iqtree2','iqtree'] if shutil.which(x)),None);need(iq,'IQ-TREE executable not found')
partial=['MC__gene_14767','SS__g18155.t1'];need(set(partial)<=set(qmeta),'Expected partial candidates missing')
out=base/'step25kz17.53'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
# Extract one contiguous carrier span, retaining native inter-repeat linkers.
hits=read(src/'03_domains_in_alignment.tsv');mapping=[];data={};idmap={};bounds={}
for j,k in enumerate(sorted(seqs),1):
 hh=[h for h in hits if h['sequence_id']==k and h['accession'].split('.')[0]=='PF00153'];need(hh,'Missing carrier evidence '+k);a=min(int(h['ali_from']) for h in hh);b=max(int(h['ali_to']) for h in hh);need(1<=a<=b<=len(seqs[k]),'Invalid extraction')
 # No flanks for the two partial single-repeat sequences; +/-10 aa elsewhere.
 lo=a if k in partial else max(1,a-10);hi=b if k in partial else min(len(seqs[k]),b+10)
 other=[h for h in hits if h['sequence_id']==k and h['accession'].split('.')[0]!='PF00153' and int(h['ali_from'])<=hi and int(h['ali_to'])>=lo]
 need(not other,'Non-carrier domain overlaps proposed extraction: '+k+'; review before tree construction')
 tid=k if k in refids else 'CAR%03d'%j;idmap[k]=tid;data[tid]=seqs[k][lo-1:hi];bounds[tid]=(lo,hi,hh)
 mapping.append(dict(tree_id=tid,protein_id=k,role='REFERENCE' if k in refids else 'CANDIDATE',species='Arabidopsis_thaliana_REFERENCE' if k in refids else qmeta[k]['species'],gene_id=next(r['loci'] for r in refs if r['reference_id']==k) if k in refids else qmeta[k]['gene_id'],screen_group=next(r['assigned_family'] for r in refs if r['reference_id']==k) if k in refids else qmeta[k]['best_group'],original_length=len(seqs[k]),extraction_start=lo,extraction_end=hi,extracted_length=hi-lo+1,carrier_hit_count=len(hh),caveat='PARTIAL_SENSITIVITY_ONLY' if k in partial else ('LONG_MODEL_CARRIER_REGION_ONLY' if k=='MC__gene_5144' else 'SCREEN_ASSIGNMENT_REVIEW_ONLY')))
table('01_extraction_and_tree_id_map.tsv',mapping,list(mapping[0]));writefa('02_all_extracted_carrier_regions.faa',data)
main={k:v for k,v in data.items() if k not in {idmap[x] for x in partial}};writefa('03_main_alignment_input.faa',main);writefa('03_partial_fragments.faa',{idmap[k]:data[idmap[k]] for k in partial})
cpu=os.environ.get('SLURM_CPUS_PER_TASK','4');commands=[]
def run(args,log,stdout=None):
 commands.append([str(x).replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for x in args])
 with (out/log).open('w') as f:
  if stdout:
   with (out/stdout).open('w') as a:subprocess.run(args,stdout=a,stderr=f,check=True)
  else:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
run(['mafft','--amino','--thread',cpu,'--localpair','--maxiterate','1000',str(out/'03_main_alignment_input.faa')],'04_mafft_main.log','04_main_alignment.faa')
aln,_=fasta(out/'04_main_alignment.faa');need(set(aln)==set(main) and len({len(v) for v in aln.values()})==1,'Invalid main alignment');need(all(aln[k].replace('-','')==v for k,v in main.items()),'Changed main sequence')
# Preserve the main alignment and its columns while adding fragments.
run(['mafft','--amino','--thread',cpu,'--addfragments',str(out/'03_partial_fragments.faa'),'--keeplength',str(out/'04_main_alignment.faa')],'05_mafft_fragments.log','05_with_fragments_alignment.faa')
added,_=fasta(out/'05_with_fragments_alignment.faa');ncol=len(next(iter(aln.values())));need(set(added)==set(data) and all(len(v)==ncol for v in added.values()),'Invalid fragment alignment');need(all(added[k]==v for k,v in aln.items()),'Fragment addition changed backbone')
fragmentqc=[]
for k in partial:
 tid=idmap[k];retained=added[tid].replace('-','');it=iter(data[tid]);need(all(any(x==c for x in it) for c in retained),'Fragment alignment altered residue order')
 fragmentqc.append(dict(protein_id=k,tree_id=tid,input_region_length=len(data[tid]),retained_backbone_residues=len(retained),residues_omitted_by_keeplength=len(data[tid])-len(retained),caveat='PARTIAL_REPEAT_PLACEMENT_UNCERTAIN'))
table('05_fragment_retention_QC.tsv',fragmentqc,list(fragmentqc[0]))
targets={'REF__Q9C5M0','REF__Q9SJY5','REF__Q9SB52','REF__Q9FY68'};domaincols={}
for k in targets:
 lo,hi,hh=bounds[k];coords=[i for i,c in enumerate(aln[k]) if c!='-'];covered=set()
 for h in hh:
  for pos in range(int(h['ali_from']),int(h['ali_to'])+1):covered.add(coords[pos-lo])
 domaincols[k]=covered
cand=set(main)-refids;mask=[];keep=[]
for i in range(ncol):
 support=sum(i in domaincols[k] for k in targets);occ=sum(aln[k][i]!='-' for k in cand)/len(cand);ok=support>=2 and occ>=.5
 if ok:keep.append(i)
 mask.append(dict(alignment_column=i+1,target_domain_support=support,main_candidate_occupancy=round(occ,6),retained=int(ok)))
need(len(keep)>=100,'Too few carrier-core columns');table('06_fixed_core_mask.tsv',mask,list(mask[0]));core={k:''.join(v[i] for i in keep) for k,v in added.items()};need(all(v.replace('-','') for v in core.values()),'Empty core sequence')
coreqc=[dict(tree_id=k,core_residues=len(v.replace('-','')),core_columns=len(keep),occupancy=round(len(v.replace('-',''))/len(keep),6)) for k,v in core.items()];table('07_core_occupancy.tsv',coreqc,list(coreqc[0]))
runs=[];model=None
for mode,extras,seed in [('main',[],175301),('sensitivity_MC_partial',[idmap[partial[0]]],175302),('sensitivity_SS_partial',[idmap[partial[1]]],175303)]:
 (out/mode).mkdir();ids=set(main)|set(extras);dd={k:v for k,v in core.items() if k in ids};writefa(mode+'/alignment.faa',dd);prefix=out/mode/'tree'
 run([iq,'-s',str(out/mode/'alignment.faa'),'-st','AA','-m',model or 'MFP','-B','1000','-bnni','-alrt','1000','-T',cpu,'-seed',str(seed),'-keep-ident','--prefix',str(prefix)],mode+'/iqtree.stdout.log')
 tree=Path(str(prefix)+'.treefile');need(tree.is_file(),'Missing tree');tips=re.findall(r'(?:\(|,)\s*([^():,;\s]+)\s*:',tree.read_text());need(len(tips)==len(dd) and set(tips)==set(dd),'Tree tip mismatch')
 if model is None:
  report=Path(str(prefix)+'.iqtree').read_text();m=re.search(r'Best-fit model according to BIC:\s*(\S+)',report);need(m,'ModelFinder model missing');model=m.group(1)
 labels=['LABELS','SEPARATOR TAB','DATA']+[r['tree_id']+'\t'+r['protein_id']+' | '+r['screen_group']+' | '+r['caveat'] for r in mapping if r['tree_id'] in ids];(out/mode/'itol_labels.txt').write_text('\n'.join(labels)+'\n')
 runs.append(dict(analysis=mode,sequences=len(dd),columns=len(keep),model=model,seed=seed,support_order='SH-aLRT/UFBoot',rooting='UNROOTED'))
table('08_tree_run_summary.tsv',runs,list(runs[0]))
for n in ['source_02_reference_annotations.tsv','source_11_phylogeny_priority_queue.tsv','source_10_target_competitor_comparisons.tsv','04_domain_architecture_review.tsv','source_prior_01_integrated_locus_evidence.tsv','source_prior_06_refreshed_family_coverage.tsv','source_prior_07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/n)
run(['mafft','--version'],'mafft_version.txt');run([iq,'--version'],'iqtree_version.txt')
(out/'README.txt').write_text('Domain-guided carrier phylogenetic review. Contiguous span from first to last PF00153 hit plus 10-aa flanks; native inter-repeat linkers retained, source proteins unchanged. Two single-repeat partial candidates use exact hit spans only and are excluded from main backbone. MC gene_5144 retains its carrier region with architecture caveat. Main MAFFT L-INS-i: 83 candidates plus 43 references. Partial MC gene_14767 and SS g18155.t1 added with MAFFT addfragments/keeplength; backbone verified unchanged, omitted fragment residues quantified. Fixed core: >=2 target-domain support and >=50% main-candidate occupancy. Three trees: main and each partial added separately. Main model selected by ModelFinder, reused in sensitivities with parameters re-estimated. 1000 SH-aLRT/1000 UFBoot with bnni; unrooted; retain identical tips. Partial single repeats can align to alternative repeat positions: sensitivity placements alone cannot establish family membership. Repeat count flags do not automatically exclude other candidates. Limited panel and bounded candidate selection remain caveats. No ledger/count updates or gene-model repairs. Upload tree.treefile then itol_labels.txt for each tree to iTOL.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),commands=commands),indent=2)+'\n');qc=dict(status='DOMAIN_GUIDED_CARRIER_PHYLOGENIES_REVIEW_READY',extracted_sequences=len(data),main_sequences=len(main),backbone_columns=ncol,core_columns=len(keep),trees=runs,fragment_QC=fragmentqc,final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
for p in out.rglob('*'):
 if p.is_file() and p.suffix in {'.log','.iqtree','.txt'}:p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out);bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Domain_guided_carrier_phylogeny/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
