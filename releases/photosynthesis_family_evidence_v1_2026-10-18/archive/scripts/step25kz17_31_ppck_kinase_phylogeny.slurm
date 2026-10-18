#!/usr/bin/env bash
#SBATCH --job-name=ppck_tree
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=ppck_tree_%j.out
#SBATCH --error=ppck_tree_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
command -v python3
command -v hmmsearch
command -v mafft
python3 - <<'PY'
import csv,gzip,hashlib,json,os,re,shutil,subprocess,tarfile,urllib.request
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True);base=root/'08_focused_families/validation'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def release(step,env):
 pp=[root/os.environ[env]] if os.environ.get(env) else [p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
 if len(pp)!=1:raise SystemExit('Set '+env+' explicitly')
 p=pp[0].resolve(strict=True)
 for l in (p/'SHA256SUMS.txt').read_text().splitlines():
  h,n=l.split(maxsplit=1);f=(p/n).resolve()
  if not f.is_relative_to(p) or sha(f)!=h:raise SystemExit('Checksum failure '+n)
 return p
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta(p):
 d={};h=None;s=[]
 for l in p.read_text().splitlines():
  if l.startswith('>'):
   if h is not None:d[h]=''.join(s)
   h=l[1:].split()[0];s=[]
   if h in d:raise SystemExit('Duplicate FASTA ID '+h)
  elif l.strip():s.append(l.strip())
 if h is not None:d[h]=''.join(s)
 return d
def writefa(p,d):p.write_text(''.join('>'+h+'\n'+s+'\n' for h,s in d.items()))
src=release('step25kz17.30','PPCK_RELEASE');screen=release('step25kz17.22','SCREEN_RELEASE')
iq=next((shutil.which(x) for x in ['iqtree3','iqtree2','iqtree'] if shutil.which(x)),None)
if not iq:raise SystemExit('IQ-TREE executable not found')
# Resolve the same installed Pfam source used in the previous review.
if os.environ.get('PFAM_HMM'):db=root/os.environ['PFAM_HMM']
else:
 homes=[root/'databases',root/'00_databases',root.parent/'my_interproscan']
 if os.environ.get('INTERPROSCAN_HOME'):homes.append(Path(os.environ['INTERPROSCAN_HOME']))
 exe=shutil.which('interproscan.sh')
 if exe:homes.append(Path(exe).resolve().parent)
 paths=set()
 for home in homes:
  if home.is_dir():
   for folder,dirs,files in os.walk(home):
    dirs[:]=[d for d in dirs if not d.startswith('.') and d not in ['temp','tmp','work']]
    for name in files:
     if name.lower() in ['pfam-a.hmm','pfam_a.hmm','pfam-a.hmm.gz','pfam_a.hmm.gz']:paths.add((Path(folder)/name).resolve())
 if len(paths)!=1:raise SystemExit('Set PFAM_HMM explicitly')
 db=next(iter(paths))
if sha(db)!=json.loads((src/'DATABASE_PROVENANCE.json').read_text())['compressed_or_source_sha256']:raise SystemExit('Pfam source differs from step17.30; use the same database')
out=base/'step25kz17.31'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
def table(name,rows,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
commands=[]
def run(args,log,stdout=None):
 commands.append(args)
 with (out/log).open('w') as f:
  if stdout:
   with (out/stdout).open('w') as g:subprocess.run(args,stdout=g,stderr=f,check=True)
  else:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
run([iq,'--version'],'iqtree_version.txt');run(['hmmsearch','-h'],'hmmer_version.txt')
# Keep all target-favoured and no-competitor cases from the 296-sequence batch,
# plus the five highest-scoring competitor-favoured examples per species.
queue=read(src/'10_PPCK_priority_batch.tsv');selected={r['protein_id'] for r in queue if r['limited_panel_status']!='COMPETITOR_TIED_OR_HIGHER'}
for sp in {r['species'] for r in queue}:
 rr=sorted([r for r in queue if r['species']==sp and r['limited_panel_status']=='COMPETITOR_TIED_OR_HIGHER'],key=lambda r:(-float(r['bitscore']),r['protein_id']))
 selected.update(r['protein_id'] for r in rr[:5])
high={r['protein_id'] for r in queue if float(r['bitscore'])>=247}
if len(high)!=27 or not high.issubset(selected):raise SystemExit('Expected observed 27-member high-score group; inspect input')
seqs=fasta(src/'02_candidates_and_references.faa');mapping=read(src/'01_sequence_map.tsv');tipmap=[];proteins={}
for r in mapping:
 if r['role']=='CANDIDATE' and r['protein_id'] in selected:
  q=r['query_id'];proteins[q]=seqs[q];tipmap.append(dict(tip=q,protein_id=r['protein_id'],role='HIGH_PRIORITY' if r['protein_id'] in high else 'BORDERLINE_COMPARISON',species=r['protein_id'].split('__',1)[0],full_length=len(seqs[q])))
refs=fasta(screen/'references/target_and_background.faa')
for acc,role in [('Q93VK0','PPCK_REFERENCE'),('Q9SPK4','PPCK_REFERENCE'),('Q8W490','PEPKR2_CONTROL'),('Q944Q0','WNK8_CONTROL')]:
 if acc not in refs:raise SystemExit('Missing local reference '+acc)
 q='REF_'+acc;proteins[q]=refs[acc];tipmap.append(dict(tip=q,protein_id=acc,role=role,species='Arabidopsis_thaliana',full_length=len(refs[acc])))
# Curated CPK controls downloaded as full JSON for identity and sequence provenance.
# CONTROL_DIR may contain Q06850.json and Q38872.json for offline compute nodes.
for acc,gene in [('Q06850','CPK1'),('Q38872','CPK6')]:
 path=out/(acc+'.json');cache=Path(os.environ['CONTROL_DIR'])/(acc+'.json') if os.environ.get('CONTROL_DIR') else None
 if cache and cache.is_file():shutil.copy2(cache,path)
 else:
  try:
   with urllib.request.urlopen('https://rest.uniprot.org/uniprotkb/'+acc+'.json',timeout=60) as f:path.write_bytes(f.read())
  except Exception as e:raise SystemExit('Download '+acc+'.json on a network-enabled node and set CONTROL_DIR. '+str(e))
 j=json.loads(path.read_text());names={g.get('geneName',{}).get('value') for g in j.get('genes',[])}
 if j.get('primaryAccession')!=acc or j.get('organism',{}).get('taxonId')!=3702 or gene not in names or 'reviewed' not in j.get('entryType','').lower() or 'unreviewed' in j.get('entryType','').lower():raise SystemExit('Control identity validation failed '+acc)
 q='REF_'+acc;s=j['sequence']['value'];proteins[q]=s;tipmap.append(dict(tip=q,protein_id=acc,role=gene+'_CONTROL',species='Arabidopsis_thaliana',full_length=len(s)))
writefa(out/'01_full_proteins.faa',proteins);table('02_tip_map.tsv',tipmap,list(tipmap[0]))
# Extract only PF00069. All sequences, including controls, use the same model.
blocks=[];block=[];op=gzip.open if db.suffix=='.gz' else open
with op(db,'rt') as f:
 for l in f:
  block.append(l)
  if l.strip()=='//':
   if any(x.startswith('ACC ') and x.split()[1].split('.')[0]=='PF00069' for x in block):blocks.append(''.join(block))
   block=[]
if len(blocks)!=1:raise SystemExit('Expected one PF00069 profile')
(out/'Pkinase.hmm').write_text(blocks[0]);cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'03_kinase.domtbl'),str(out/'Pkinase.hmm'),str(out/'01_full_proteins.faa')],'03_hmmsearch.log')
dom=defaultdict(list)
for l in (out/'03_kinase.domtbl').read_text().splitlines():
 if l and not l.startswith('#'):
  x=l.split();dom[x[0]].append((float(x[13]),int(x[17]),int(x[18]),int(x[15]),int(x[16]),int(x[5])))
crops={};cropmap=[]
for q,s in proteins.items():
 hits=sorted(dom[q],reverse=True)
 if not hits:
  if q.startswith('REF_'):raise SystemExit('Control lacks kinase hit '+q)
  cropmap.append(dict(tip=q,status='NO_GA_DOMAIN_NOT_EXCLUDED',domain_hits=0,start='',end='',model_coverage=''));continue
 score,a,b,ha,hb,hl=hits[0];crops[q]=s[a-1:b];cropmap.append(dict(tip=q,status='HIGHEST_SCORE_HIT_REVIEW' if len(hits)>1 else 'DOMAIN_EXTRACTED',domain_hits=len(hits),start=a,end=b,model_coverage=round((hb-ha+1)/hl,4)))
table('04_domain_crops.tsv',cropmap,['tip','status','domain_hits','start','end','model_coverage']);writefa(out/'05_kinase_domains.faa',crops)
run(['mafft','--localpair','--maxiterate','1000','--thread',cpu,str(out/'05_kinase_domains.faa')],'06_mafft.log','06_domain_alignment.faa')
aln=fasta(out/'06_domain_alignment.faa')
if set(aln)!=set(crops) or len(set(map(len,aln.values())))!=1:raise SystemExit('Invalid alignment')
for q in aln:
 if aln[q].replace('-','').upper()!=crops[q].upper():raise SystemExit('Alignment changed residues')
keep=[i for i in range(len(next(iter(aln.values())))) if aln['REF_Q93VK0'][i]!='-' or aln['REF_Q9SPK4'][i]!='-']
core={q:''.join(s[i] for i in keep) for q,s in aln.items()};writefa(out/'07_reference_core.faa',core)
(out/'07_retained_columns.txt').write_text('\n'.join(str(i+1) for i in keep)+'\n')
# Exact locus reconciliation, no transcript-suffix guessing or sequence deduplication.
loci=[];hashes={}
for sp in sorted({r['species'] for r in tipmap if not r['tip'].startswith('REF_')}):
 path=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv');rr=read(path) if path.is_file() else [];hashes[sp]=sha(path) if path.is_file() else 'MISSING'
 for t in [t for t in tipmap if t['species']==sp]:
  local=t['protein_id'].split('__',1)[1];gs={r['gene_id'] for r in rr if local in {r.get(k,'') for k in ['gene_id','transcript_id','protein_id','cds_id']}}
  loci.append(dict(tip=t['tip'],protein_id=t['protein_id'],species=sp,gene_id=';'.join(sorted(gs)),mapping_status='UNAMBIGUOUS' if len(gs)==1 else 'UNRESOLVED_OR_AMBIGUOUS',decision='PENDING_PHYLOGENETIC_REVIEW'))
table('08_locus_reconciliation.tsv',loci,list(loci[0]));table('09_core_occupancy.tsv',[dict(tip=q,core_columns=len(s),occupied=sum(c!='-' for c in s)) for q,s in core.items()],['tip','core_columns','occupied'])
# Unrooted exploratory tree; controls are not imposed as an outgroup.
run([iq,'-s',str(out/'07_reference_core.faa'),'-st','AA','-m','MFP','-B','1000','-alrt','1000','-T',cpu,'-seed','1731','-keep-ident','--prefix',str(out/'10_PPCK_kinase')],'10_iqtree.stdout.log')
if not (out/'10_PPCK_kinase.treefile').is_file():raise SystemExit('Tree output missing')
(out/'README.txt').write_text('Exploratory unrooted PPCK/related-kinase tree. No automated family acceptance or count.\nIncludes 27 observed high-score candidates, all target-favoured/no-competitor cases from the bounded batch, and five competitor-favoured cases per species. Deferred screen hits remain unclassified.\nControls: local PPCK1/2, PEPKR2 and WNK8; reviewed UniProt CPK1 and CPK6. Panel is not exhaustive.\nHighest-scoring PF00069 hit cropped per protein; inspect multi-hit and incomplete domains. Alignment mask retains columns occupied by either PPCK reference; raw alignment retained for sensitivity review.\nIdentical sequences retained; their branching order may be arbitrary. Gene counts require locus reconciliation and supported placement, not tip counting.\nReview SH-aLRT/UFBoot support, reference placement, domain completeness and competing kinase placement together. No rooting imposed.\n')
(out/'PROVENANCE.json').write_text(json.dumps({'review_manifest_sha256':sha(src/'SHA256SUMS.txt'),'screen_manifest_sha256':sha(screen/'SHA256SUMS.txt'),'pfam_sha256':sha(db),'mapping_hashes':hashes,'control_json_hashes':{a:sha(out/(a+'.json')) for a in ['Q06850','Q38872']}},indent=2)+'\n')
(out/'QC.json').write_text(json.dumps({'status':'TREE_COMPLETE_REVIEW_REQUIRED','candidate_tips_selected':len(selected),'high_priority':27,'controls':6,'tree_tips':len(core),'core_columns':len(keep),'final_counts':False},indent=2)+'\n')
# Store commands with output/project paths anonymized.
(out/'COMMANDS.json').write_text(json.dumps([[x.replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for x in cmd] for cmd in commands],indent=2)+'\n')
with (out/'SHA256SUMS.txt').open('w') as f:
 for p in sorted(out.iterdir()):
  if p.is_file() and p.name!='SHA256SUMS.txt':f.write(sha(p)+'  '+p.name+'\n')
bundle=out.parent/(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='PPCK_kinase_phylogeny')
print('BUNDLE='+str(bundle))
PY
