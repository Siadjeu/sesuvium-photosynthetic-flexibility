#!/usr/bin/env bash
#SBATCH --job-name=pepc_ledger_tree
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=pepc_ledger_tree_%j.out
#SBATCH --error=pepc_ledger_tree_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v mafft
python3 - <<'PY'
import csv,gzip,hashlib,json,os,re,shutil,subprocess,tarfile
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
src=release('step25kz17.42','PEPC_REVIEW_RELEASE');dom=release('step25kz17.41','PEPC_DOMAIN_RELEASE')
need(json.loads((src/'PROVENANCE.json').read_text())['prior_manifest_sha256']==sha(dom/'SHA256SUMS.txt'),'Domain/review provenance mismatch')
iq=next((shutil.which(x) for x in ['iqtree3','iqtree2','iqtree'] if shutil.which(x)),None);need(iq,'IQ-TREE executable not found')
cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
q,_=fasta(src/'02_53_candidates_four_controls.faa');refs={k for k in q if k.startswith('REF__')};candidates=set(q)-refs
need(len(candidates)==53 and len(refs)==4,'Unexpected candidate/reference counts')
records=read(src/'01_reconciled_domain_evidence.tsv');cat={r['protein_id']:r['evidence_category'] for r in records if r['role']=='CANDIDATE'}
need(set(cat)==candidates,'Candidate evidence mismatch')
caveats={k for k in candidates if cat[k]=='PRIORITY_COMPLETENESS_REVIEW'};main=candidates-caveats
need(len(main)==51 and caveats=={'SP_Hex__gene_Spor0168910.1','SP_Hex__gene_Spor0380090.1'},'Unexpected evidence categories')
meta={r['protein_id']:r for r in read(dom/'10_prior_provenance.tsv')};exceptions=read(dom/'05_CDS_exception_ledger.tsv')
need(set(meta)==candidates|{r['protein_id'] for r in exceptions},'55-locus provenance mismatch')
out=base/'step25kz17.43'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False)
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+s+'\n' for k,s in d.items()))
# Stable safe tip IDs; the map retains original protein/locus IDs for iTOL and review.
tips={k:'PEPC_'+str(i).zfill(3) for i,k in enumerate(sorted(q),1)}
tipmap=[dict(tip_id=tips[k],original_id=k,species=meta[k]['species'] if k in meta else 'REFERENCE',gene_id=meta[k]['gene_id'] if k in meta else '',evidence_category=cat.get(k,'REFERENCE_CONTROL'),main_tree=k in main|refs,sensitivity_tree=True) for k in sorted(q)]
table('01_tree_tip_map.tsv',tipmap,list(tipmap[0]))
# Record CDS consistency for the two caveats without claiming genome/read validation.
ps,_=fasta(src/'annotation/04_protein_isoforms.fasta');cs,_=fasta(src/'annotation/04_cds_isoforms.fasta');mapping=read(src/'annotation/01_all_isoform_mapping.tsv')
code=dict(zip((a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'),'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
cons=[]
for k in sorted(caveats):
 mm=[r for r in mapping if r['gene_id']==meta[k]['gene_id']];need(len(mm)==1,'Expected one mapped isoform per exception')
 r=mm[0];aliases={r.get(f,'') for f in ['gene_id','transcript_id','protein_id','cds_id']};aliases|={meta[k]['species']+'__'+a for a in list(aliases)}
 pp=[(a,s) for a,s in ps.items() if a in aliases];cc=[(a,s) for a,s in cs.items() if a in aliases]
 need(len(pp)==len(cc)==1,'Ambiguous annotation protein/CDS pairing')
 aa=norm(''.join(code.get(cc[0][1][i:i+3],'X') for i in range(0,len(cc[0][1])-2,3)))
 ok=len(cc[0][1])%3==0 and aa==norm(pp[0][1])==q[k] and '*' not in aa and 'X' not in aa
 need(ok,'Exception CDS/protein consistency changed')
 cons.append(dict(protein_id=k,protein_source_id=pp[0][0],cds_source_id=cc[0][0],cds_length=len(cc[0][1]),translation_length=len(aa),status='EXACT_TRANSLATION_NO_INTERNAL_STOP',biological_completeness='NOT_ESTABLISHED'))
table('02_exception_CDS_consistency.tsv',cons,list(cons[0]))
notes={'SP_Hex__gene_Spor0168910.1':'Substantial near-N-terminal reference gap; 72.69% Pfam union coverage; annotation-consistent but completeness unresolved','SP_Hex__gene_Spor0380090.1':'Internal reference gaps and insertions; 79.11% Pfam union coverage; annotation-consistent but completeness unresolved'}
ledger=[]
for k,r in sorted(meta.items()):
 category='PEPC_DOMAIN_SUPPORTED_BROAD_COVERAGE' if k in main else 'PEPC_FAMILY_SUPPORTED_MODEL_CAVEAT' if k in caveats else 'EXCLUDED_HISTORICAL_CDS_RECORD'
 ledger.append(dict(wang_family='PEPC',species=r['species'],gene_id=r['gene_id'],protein_id=k,evidence_category=category,count_tier='SUPPORTED' if k in main else 'MODEL_CAVEAT' if k in caveats else 'OUTSIDE_FAMILY_SET',prior_domain_category=cat.get(k,''),limitation=notes.get(k,'Historical CDS used as protein; translated sequence lacks PF00311 and matches PF03328; precise function unassigned' if k not in candidates else 'Broad domain support does not establish function or exhaustive membership'),main_tree=k in main,sensitivity_tree=k in candidates,phylogenetic_review='PENDING' if k in candidates else 'NOT_APPLICABLE',whole_genome_completeness='NOT_ESTABLISHED',C4_CAM_function='NOT_ESTABLISHED',final_gene_count='NA'))
need(len(ledger)==55 and len({(r['species'],r['gene_id']) for r in ledger})==55,'Duplicate locus ledger')
table('03_PEPC_evidence_ledger.tsv',ledger,list(ledger[0]))
summary=[]
for sp in sorted({r['species'] for r in ledger}):
 counts=Counter(r['count_tier'] for r in ledger if r['species']==sp)
 summary.append(dict(species=sp,broad_domain_supported=counts['SUPPORTED'],model_caveat=counts['MODEL_CAVEAT'],excluded_historical_records=counts['OUTSIDE_FAMILY_SET'],final_Wang_count='NA'))
table('04_PEPC_species_evidence_summary.tsv',summary,list(summary[0]))
# Keep consolidated project ledger unchanged until tree review and cross-family integration.
shutil.copyfile(src/'prior_11_integrated_locus_evidence_UNCHANGED.tsv',out/'05_prior_integrated_ledger_UNCHANGED.tsv')
shutil.copyfile(src/'prior_12_final_Wang_counts_PENDING.tsv',out/'06_final_Wang_counts_PENDING.tsv')
commands=[]
def run(args,log,stdout=None):
 commands.append([str(x).replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for x in args])
 with (out/log).open('w') as f:
  if stdout:
   with (out/stdout).open('w') as a:subprocess.run(args,stdout=a,stderr=f,check=True)
  else:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
run([iq,'--version'],'iqtree_version.txt');run(['mafft','--version'],'mafft_version.txt')
# Main alignment excludes caveats before alignment inference. Sensitivity reuses the
# verified 57-protein corrected alignment from step17.42 (no invalid CDS inputs).
rawsens,_=fasta(src/'03_PEPC_review_alignment.faa');need(set(rawsens)==set(q),'Sensitivity alignment IDs mismatch')
need(all(rawsens[k].replace('-','')==q[k] for k in q),'Sensitivity alignment sequence mismatch')
treestats=[]
for mode,selected,seed in [('main',main|refs,174301),('sensitivity',candidates|refs,174302)]:
 folder=out/mode;folder.mkdir();selectedseq={tips[k]:q[k] for k in sorted(selected)}
 writefa(mode+'/01_input_proteins.faa',selectedseq)
 if mode=='main':
  run(['mafft','--amino','--thread',cpu,'--localpair','--maxiterate','1000',str(folder/'01_input_proteins.faa')],mode+'/02_mafft.log',mode+'/02_raw_alignment.faa')
 else:writefa(mode+'/02_raw_alignment.faa',{tips[k]:rawsens[k] for k in sorted(selected)})
 aln,_=fasta(folder/'02_raw_alignment.faa');need(set(aln)==set(selectedseq) and len({len(s) for s in aln.values()})==1,'Invalid alignment')
 need(all(aln[k].replace('-','')==selectedseq[k] for k in aln),'Alignment changed sequence')
 width=len(next(iter(aln.values())));refids={tips[k] for k in refs};mask=[];maskrows=[]
 # Reference-anchored core: retain columns occupied by any PPC reference and
 # >=50% of the selected non-reference candidates. Threshold is recorded.
 candtips=set(aln)-refids
 for i in range(width):
  occ=sum(aln[k][i]!='-' for k in candtips)/len(candtips);anchor=any(aln[k][i]!='-' for k in refids);keep=anchor and occ>=.5
  if keep:mask.append(i)
  maskrows.append(dict(raw_column_1based=i+1,reference_anchored=anchor,candidate_occupancy=round(occ,6),retained=keep))
 need(len(mask)>=500,'Unexpectedly short reference-anchored alignment; inspect before tree inference')
 core={k:''.join(s[i] for i in mask) for k,s in aln.items()};need(all(any(c!='-' for c in s) for s in core.values()),'Empty core sequence')
 writefa(mode+'/03_reference_core.faa',core);table(mode+'/03_column_mask.tsv',maskrows,list(maskrows[0]))
 identical=defaultdict(list)
 for k,s in core.items():identical[s].append(k)
 groups=[dict(group=i,tip_ids=';'.join(v),sequence_count=len(v)) for i,v in enumerate(identical.values(),1) if len(v)>1]
 table(mode+'/04_identical_sequence_groups.tsv',groups,['group','tip_ids','sequence_count'])
 prefix=folder/'PEPC_corrected'
 run([iq,'-s',str(folder/'03_reference_core.faa'),'-st','AA','-m','MFP','-B','1000','-bnni','-alrt','1000','-T',cpu,'-seed',str(seed),'-keep-ident','--prefix',str(prefix)],mode+'/05_iqtree.stdout.log')
 tree=Path(str(prefix)+'.treefile');need(tree.is_file(),'Missing completed tree')
 treeids=re.findall(r'(?<=[(,])(PEPC_\d+)(?=[:),])',tree.read_text())
 need(len(treeids)==len(selected) and set(treeids)==set(selectedseq),'Tree tip set mismatch')
 # iTOL label file maps safe tree IDs back to original identifiers.
 (folder/'iTOL_labels.txt').write_text('LABELS\nSEPARATOR TAB\nDATA\n'+''.join(tips[k]+'\t'+k+(' [MODEL_CAVEAT]' if k in caveats else '')+'\n' for k in sorted(selected)))
 treestats.append(dict(analysis=mode,candidates=len(selected)-4,references=4,raw_alignment_columns=width,core_columns=len(mask),tree_tips=len(treeids),seed=seed,rooting='UNROOTED',review_status='PENDING'))
table('07_tree_run_summary.tsv',treestats,list(treestats[0]))
(out/'COMMANDS.json').write_text(json.dumps(commands,indent=2)+'\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(review_manifest_sha256=sha(src/'SHA256SUMS.txt'),domain_manifest_sha256=sha(dom/'SHA256SUMS.txt'),sensitivity_alignment_sha256=sha(src/'03_PEPC_review_alignment.faa'),main_alignment='De novo MAFFT without two caveats',mask='Any reference residue plus >=50% candidate occupancy; independently applied per alignment'),indent=2)+'\n')
(out/'README.txt').write_text('PEPC bounded evidence ledger: 51 broadly domain-supported loci, two family-supported loci with model caveats, two excluded historical CDS records. This is not a final Wang count or whole-genome completeness claim. Two exception CDS/protein pairs are consistent; this does not validate gene models biologically. Main tree: 51 candidates + four controls, de novo alignment. Sensitivity tree: 53 candidates + controls, verified corrected step17.42 alignment. Both reference-anchored masks require >=50% candidate occupancy. Raw alignments and masks retained; between-tree differences can reflect both sample composition and alignment/mask changes. IQ-TREE ModelFinder, 1000 SH-aLRT and 1000 UFBoot replicates with BNNI; fixed seeds. Trees unrooted; no outgroup or subtype/function assignment imposed. Identical sequences retained as distinct loci; their branching order may be uninformative. Support interpretation and main/sensitivity comparison remain manual review tasks. Original 55-candidate tree superseded because of CDS contamination. iTOL treefiles and label datasets supplied for manual upload. Project-wide ledger and pending-count matrix remain unchanged.\n')
# Text logs may contain executable paths; scrub project/home/executable locations.
for p in out.rglob('*'):
 if p.is_file() and (p.suffix in {'.log','.iqtree','.txt'}):
  s=p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>').replace(str(Path(iq).parent),'<TOOL_BIN>');p.write_text(s)
qc=dict(status='PEPC_LEDGER_AND_TWO_TREES_READY_FOR_REVIEW',ledger_loci=55,categories=dict(Counter(r['evidence_category'] for r in ledger)),trees=treestats,final_counts=False,phylogenetic_acceptance='PENDING_REVIEW')
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
# Keep checkpoints on the cluster, omit compressed checkpoint/model caches from bundle.
files=sorted(p for p in out.rglob('*') if p.is_file() and not p.name.endswith(('.ckp.gz','.model.gz')))
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='PEPC_ledger_corrected_phylogeny/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
