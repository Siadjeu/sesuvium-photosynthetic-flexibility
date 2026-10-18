#!/usr/bin/env bash
#SBATCH --job-name=ppt_tree_review
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=ppt_tree_review_%j.out
#SBATCH --error=ppt_tree_review_%j.err
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







src=release('step25kz17.66','PPT_SCREEN_RELEASE');iq=shutil.which('iqtree2') or shutil.which('iqtree');need(iq,'IQ-TREE executable missing')
rows=read(src/'09_combined_PPT_specificity_review.tsv');need(len(rows)==451,'Unexpected candidate scope')
meta={r['protein_id']:r for r in rows};seqs,_=fasta(src/'05_candidates_and_references.faa');panel=read(src/'reference_panel_metadata.tsv');refs={r['reference_id']:r for r in panel};refids=set(refs)
need(set(seqs)==set(meta)|refids and len(refids)==6,'Sequence/metadata mismatch')
dom=defaultdict(list)
for h in read(src/'07_Pfam_domain_hits.tsv'):dom[h['sequence_id']].append(h)
def hh(k):return [h for h in dom[k] if h['accession']=='PF03151']
def span(k):return (min(int(h['ali_from']) for h in hh(k)),max(int(h['ali_to']) for h in hh(k)))
def coverage(k):
 covered=set()
 for h in hh(k):covered.update(range(int(h['hmm_from']),int(h['hmm_to'])+1))
 return len(covered)/int(hh(k)[0]['model_length']) if hh(k) else 0.0
need(all(hh(k) for k in refs),'Missing reference PF03151')
max_ref_len=max(len(seqs[k]) for k in refs);max_ref_region=max(span(k)[1]-span(k)[0]+1 for k in refs)
out=base/'step25kz17.67'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
def table(n,rr,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in sorted(d.items())))
og=[r for r in rows if r['historical_OG0000369']=='True'];need(len(og)==37 and all(r['similarity_category']=='NO_CURRENT_PANEL_HIT' and not hh(r['protein_id']) for r in og),'OG0000369 conclusions changed')
correction=[dict(protein_id=r['protein_id'],species=r['species'],gene_id=r['gene_id'],historical_label='PPT_OG0000369',review_decision='HISTORICAL_PPT_LABEL_UNSUPPORTED',domain_evidence=r['Pfam_models'] or 'NO_GA_DOMAIN',scope='NO_SIX_REFERENCE_HIT;NO_PF03151;NO_FUNCTION_OR_ABSENCE_CLAIM') for r in og]
table('00_OG0000369_label_correction.tsv',correction,list(correction[0]))
recal=[];mapping=[];omitted=[];data={};partial=[];idmap={};bounds={}
for j,k in enumerate(sorted(seqs),1):
 isref=k in refs;target=hh(k);flags=[];overlap=[];outside=[]
 if not target:
  omitted.append(dict(protein_id=k,reason='NO_GA_SUPPORTED_PF03151_REGION',historical_OG0000369=meta.get(k,{}).get('historical_OG0000369','False')))
  recal.append(dict(protein_id=k,PF03151_model_coverage=0,region_length=0,full_length=len(seqs[k]),region_class='OUTSIDE_DOMAIN_TREE',overlapping_other_models='',nonoverlapping_other_models='',flags='NO_PF03151;RETAIN_OUTSIDE_TREE'))
  continue
 lo,hi=span(k);cov=coverage(k);need(1<=lo<=hi<=len(seqs[k]),'Invalid domain coordinates')
 for h in dom[k]:
  if h['accession']=='PF03151':continue
  a,b=int(h['ali_from']),int(h['ali_to']);ov=max(0,min(hi,b)-max(lo,a)+1)
  (overlap if ov>0 else outside).append(h['accession'])
 if cov<.8:flags.append('PARTIAL_MODEL_COVERAGE')
 if len(seqs[k])>1.5*max_ref_len:flags.append('LONG_FULL_PROTEIN_REVIEW')
 if hi-lo+1>1.5*max_ref_region:flags.append('LONG_DOMAIN_SPAN_REVIEW')
 if len(target)>1:flags.append('MULTIPLE_PF03151_HITS_REVIEW')
 if outside:flags.append('NONOVERLAPPING_OTHER_DOMAIN_REVIEW')
 if 'X' in seqs[k][lo-1:hi]:flags.append('AMBIGUOUS_DOMAIN_RESIDUES')
 # Long spans remain sensitivity-only to protect the main alignment.
 sensitivity=cov<.8 or hi-lo+1>1.5*max_ref_region
 group='PARTIAL_OR_LONG_SPAN_SENSITIVITY' if sensitivity else 'MAIN_DOMAIN_REGION'
 tid=k if isref else 'PPT%03d'%j;idmap[k]=tid;data[tid]=seqs[k][lo-1:hi];bounds[tid]=(lo,hi,target)
 if sensitivity:partial.append(k)
 mapping.append(dict(tree_id=tid,protein_id=k,role='REFERENCE' if isref else 'CANDIDATE',species='Arabidopsis_thaliana_REFERENCE' if isref else meta[k]['species'],gene_id='' if isref else meta[k]['gene_id'],screen_group=refs[k]['gene'] if isref else meta[k]['similarity_category'],extraction_start=lo,extraction_end=hi,model_coverage=round(cov,6),caveat=';'.join(flags) or 'NO_RECALIBRATED_HEURISTIC_FLAG',region_class=group))
 recal.append(dict(protein_id=k,PF03151_model_coverage=round(cov,6),region_length=hi-lo+1,full_length=len(seqs[k]),region_class=group,overlapping_other_models=';'.join(sorted(set(overlap))),nonoverlapping_other_models=';'.join(sorted(set(outside))),flags=';'.join(flags) or 'NO_RECALIBRATED_HEURISTIC_FLAG'))
need(not(refids&set(partial)),'Reference moved to sensitivity; inspect scope')
table('01_recalibrated_domain_completeness.tsv',recal,list(recal[0]));table('02_tree_tip_mapping.tsv',mapping,list(mapping[0]));table('02_omitted_no_domain.tsv',omitted,['protein_id','reason','historical_OG0000369']);writefa('02_all_extracted_domain_regions.faa',data)
main={k:v for k,v in data.items() if k not in {idmap[x] for x in partial}};need(10<len(main)<=500 and partial,'Unexpected tree scope')
writefa('03_main_alignment_input.faa',main);writefa('03_partial_fragments.faa',{idmap[k]:data[idmap[k]] for k in partial})
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
 fragmentqc.append(dict(protein_id=k,tree_id=tid,input_region_length=len(data[tid]),retained_backbone_residues=len(retained),residues_omitted_by_keeplength=len(data[tid])-len(retained),caveat='PARTIAL_REGION_PLACEMENT_UNCERTAIN'))
table('05_fragment_retention_QC.tsv',fragmentqc,list(fragmentqc[0]))
targets=refids & set(main);domaincols={}
for k in targets:
 lo,hi,hh=bounds[k];coords=[i for i,c in enumerate(aln[k]) if c!='-'];covered=set()
 for h in hh:
  for pos in range(int(h['ali_from']),int(h['ali_to'])+1):covered.add(coords[pos-lo])
 domaincols[k]=covered
cand=set(main)-refids;mask=[];keep=[]
for i in range(ncol):
 support=sum(i in domaincols[k] for k in targets);occ=sum(aln[k][i]!='-' for k in cand)/len(cand);ok=support>=3 and occ>=.5
 if ok:keep.append(i)
 mask.append(dict(alignment_column=i+1,target_domain_support=support,main_candidate_occupancy=round(occ,6),retained=int(ok)))
need(len(keep)>=100,'Too few carrier-core columns');table('06_fixed_core_mask.tsv',mask,list(mask[0]));core={k:''.join(v[i] for i in keep) for k,v in added.items()};need(all(core[k].replace('-','') for k in main),'Empty main core sequence')
coreqc=[dict(tree_id=k,core_residues=len(v.replace('-','')),core_columns=len(keep),occupancy=round(len(v.replace('-',''))/len(keep),6)) for k,v in core.items()];table('07_core_occupancy.tsv',coreqc,list(coreqc[0]));table('07_fragments_below_30_core_residues.tsv',[r for r in coreqc if r['tree_id'] not in main and r['core_residues']<30],list(coreqc[0]))
runs=[];model=None
for mode,extras,seed in [('main',[],176701),('sensitivity_all_partials',[idmap[k] for k in partial if len(core[idmap[k]].replace('-',''))>=30],176702)]:
 (out/mode).mkdir();ids=set(main)|set(extras);dd={k:v for k,v in core.items() if k in ids};writefa(mode+'/alignment.faa',dd);prefix=out/mode/'tree'
 run([iq,'-s',str(out/mode/'alignment.faa'),'-st','AA','-m',model or 'MFP','-B','1000','-bnni','-alrt','1000','-T',cpu,'-seed',str(seed),'-keep-ident','--prefix',str(prefix)],mode+'/iqtree.stdout.log')
 tree=Path(str(prefix)+'.treefile');need(tree.is_file(),'Missing tree');tips=re.findall(r'(?:\(|,)\s*([^():,;\s]+)\s*:',tree.read_text());need(len(tips)==len(dd) and set(tips)==set(dd),'Tree tip mismatch')
 if model is None:
  report=Path(str(prefix)+'.iqtree').read_text();m=re.search(r'Best-fit model according to BIC:\s*(\S+)',report);need(m,'ModelFinder model missing');model=m.group(1)
 labels=['LABELS','SEPARATOR TAB','DATA']+[r['tree_id']+'\t'+r['protein_id']+' | '+r['screen_group']+' | '+r['caveat'] for r in mapping if r['tree_id'] in ids];(out/mode/'itol_labels.txt').write_text('\n'.join(labels)+'\n')
 runs.append(dict(analysis=mode,sequences=len(dd),columns=len(keep),model=model,seed=seed,support_order='SH-aLRT/UFBoot',rooting='UNROOTED'))
table('08_tree_run_summary.tsv',runs,list(runs[0]))

# Preserve prior evidence and count state. Domain similarity is not substrate specificity.
for n in ['source_01_integrated_locus_evidence.tsv','source_07_final_Wang_counts_PENDING.tsv','09_combined_PPT_specificity_review.tsv','07_Pfam_domain_hits.tsv','reference_panel_metadata.tsv']:
 shutil.copyfile(src/n,out/n if n.startswith('source_') else 'source_'+n)
coverage_rows=read(src/'11_refreshed_family_coverage.tsv')
for r in coverage_rows:
 if r['wang_family']=='PPT':r.update(current_status='DOMAIN_PHYLOGENY_REVIEW_PENDING',next_action='REVIEW_REFERENCE_ASSOCIATION_AND_MAIN_SENSITIVITY_STABILITY;RETAIN_PARTIAL_AND_ARCHITECTURE_CAVEATS',search_completeness='BOUNDED_PANEL_SCREEN;EXHAUSTIVENESS_NOT_ESTABLISHED')
table('09_refreshed_family_coverage.tsv',coverage_rows,list(coverage_rows[0]))
run(['mafft','--version'],'mafft_version.txt');run([iq,'--version'],'iqtree_version.txt')
(out/'README.txt').write_text('PPT bounded domain phylogeny review. OG0000369 historical PPT label unsupported for37 loci, retained in correction table; not relabeled as proven enzymes and no absence claim. PF03151 is the alignment target, not proof of PPT function. PF00892 and other overlapping model hits are recorded, not required as independent domains. Full protein>1.5x maximum reference length, nonoverlapping domains, multiple target hits, and long target spans remain explicit review flags. Main includes PF03151 union HMM coverage>=0.8 and extracted span<=1.5x maximum reference span. These are alignment selection heuristics, not formal acceptance thresholds. All PF03151-bearing candidates included regardless of similarity winner, with six Arabidopsis references; no-domain candidates retained outside trees. Contiguous first-to-last PF03151 hit region extracted preserving interhit residues. MAFFT LINSI main; fragments added with keeplength; residue losses reported. Fixed core supported by>=3 references and>=50% main candidates. Fragments with<30 retained core residues omitted from sensitivity tree and explicitly listed. Main ModelFinder model reused for sensitivity;1000 SH-aLRT and1000 UFBoot bnni; identical tips retained; unrooted. Joint fragment sensitivity may be influenced by multiple partial sequences. No branch group automatically named PPT, no functional subtype assignments or final counts. Main/sensitivity trees, tip map and iTOL labels supplied. Existing ledger unchanged.\n')
qc=dict(status='PPT_DOMAIN_PHYLOGENY_REVIEW_READY',source_candidates=len(rows),historical_label_corrections=len(correction),main_candidates=len(main)-len(refids),reference_sequences=len(refids),sensitivity_only_candidates=len(partial),outside_domain_tree=len(omitted),core_columns=len(keep),tree_runs=runs,flags=dict(Counter(r['flags'] for r in recal)),final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),source_provenance=json.loads((src/'PROVENANCE.json').read_text()),commands=commands),indent=2)+'\n')
for p in out.rglob('*'):
 if p.is_file() and p.suffix in {'.log','.iqtree','.txt'}:p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='PPT_domain_phylogeny_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle))

PY
