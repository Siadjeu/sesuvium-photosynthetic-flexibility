#!/usr/bin/env bash
#SBATCH --job-name=almt_groups_iso
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=almt_groups_iso_%j.out
#SBATCH --error=almt_groups_iso_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
command -v diamond
command -v hmmsearch
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
def edges(path):
 t=re.findall(r'[^\s(),:;]+|[(),:;]',path.read_text());pos=0;ee=[]
 def node():
  nonlocal pos
  label=''
  if t[pos]=='(':
   pos+=1;z=set()
   while True:
    z|=node()
    if t[pos]==',':pos+=1;continue
    need(t[pos]==')','Malformed tree');pos+=1;break
   if t[pos] not in {':',',',')',';'}:label=t[pos];pos+=1
  else:z={t[pos]};pos+=1
  if t[pos]==':':pos+=2
  ee.append((z,label));return z
 universe=node();need(t[pos]==';','Invalid tree ending');return universe,[(z,l) for a,l in ee for z in (a,universe-a) if 1<len(z)<len(universe)-1]
src=release('step25kz17.62','ALMT_TREE_REVIEW_RELEASE');previous=release('step25kz17.61','ALMT_RECOVERY_RELEASE');need(json.loads((src/'PROVENANCE.json').read_text())['source_manifest_sha256']==sha(previous/'SHA256SUMS.txt'),'Release chain mismatch');prov=json.loads((previous/'PROVENANCE.json').read_text())
out=base/'step25kz17.63'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
# Only documented alternative isoforms, from checksum-verified collected source files.
checks=read(src/'22_protein_CDS_consistency.tsv');filemap=read(src/'20_annotation_file_map.tsv');recovered={};isoaudit=[]
for candidate in ['SS__g16509.t2','SS__g17774.t2']:
 folder=next(r['folder'] for r in filemap if r['protein_id']==candidate);aa,_=fasta(src/folder/'protein_isoforms.fasta');cc=[r for r in checks if r['candidate_id']==candidate and r['protein_matches_candidate']=='False'];need(len(cc)==1,'Expected one alternative isoform per SS locus')
 for c in cc:
  pid=c['protein_id'];need(pid in aa and c['translation_matches_source_protein']=='True','Alternative protein/CDS mismatch');k='SS__'+pid;need(valid(norm(aa[pid])),'Invalid alternative protein');recovered[k]=norm(aa[pid]);isoaudit.append(dict(candidate_id=candidate,alternative_protein_id=k,length_aa=len(recovered[k]),translation_consistent=c['translation_matches_source_protein'],internal_stops=c['internal_stops'],terminal_codon=c['terminal_codon'],source_folder=folder))
need(len(recovered)==2,'Expected two alternatives');table('01_alternative_isoform_provenance.tsv',isoaudit,list(isoaudit[0]));writefa('03_recovered_SS_proteins.faa',recovered);shutil.copyfile(previous/'02_ALMT_positive_references.faa',out/'02_ALMT_positive_references.faa')
if os.environ.get('PFAM_HMM'):db=(root/os.environ['PFAM_HMM']).resolve(strict=True)
else:
 paths=set();homes=[root/'databases',root/'00_databases',root.parent/'my_interproscan']
 if os.environ.get('INTERPROSCAN_HOME'):homes.append(Path(os.environ['INTERPROSCAN_HOME']))
 exe=shutil.which('interproscan.sh')
 if exe:homes.append(Path(exe).resolve().parent)
 for home in homes:
  if home.is_dir():
   for folder,dirs,files in os.walk(home):
    dirs[:]=[d for d in dirs if not d.startswith('.') and d not in ['tmp','temp','work']]
    for n in files:
     if n.lower() in {'pfam-a.hmm','pfam_a.hmm','pfam-a.hmm.gz','pfam_a.hmm.gz'}:paths.add((Path(folder)/n).resolve())
 need(len(paths)==1,'Set PFAM_HMM to the intended installed Pfam-A .hmm or .hmm.gz');db=next(iter(paths))
need(sha(db)==prov['pfam_sha256'],'Pfam database differs from reused results')
dbuse=db
if db.suffix=='.gz':
 dbuse=out/'work/Pfam-A.hmm'
 with gzip.open(db,'rb') as a,dbuse.open('wb') as b:shutil.copyfileobj(a,b)
models=ga=0
with dbuse.open() as f:
 for line in f:
  if line.startswith('HMMER3/'):models+=1
  elif line.startswith('GA '):ga+=1
need(models>1000 and models==ga,'Expected complete Pfam with GA thresholds')
cpu=os.environ.get('SLURM_CPUS_PER_TASK','4');commands=[]
def run(args,log,stdout=None):
 commands.append([str(x).replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for x in args])
 with (out/log).open('w') as f:
  if stdout:
   with (out/stdout).open('w') as a:subprocess.run(args,stdout=a,stderr=f,check=True)
  else:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
if recovered:run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'04_pfam.domtbl'),str(dbuse),str(out/'03_recovered_SS_proteins.faa')],'04_hmmsearch.log')
dh=[];domains=defaultdict(list)
for line in ((out/'04_pfam.domtbl').read_text().splitlines() if recovered else []):
 if not line or line.startswith('#'):continue
 x=line.split(maxsplit=22);r=dict(sequence_id=x[0],model=x[3],accession=x[4],model_length=int(x[5]),i_evalue=float(x[12]),score=float(x[13]),hmm_from=int(x[15]),hmm_to=int(x[16]),ali_from=int(x[17]),ali_to=int(x[18]));dh.append(r);domains[x[0]].append(r)
table('06_all_Pfam_domain_hits.tsv',dh,['sequence_id','model','accession','model_length','i_evalue','score','hmm_from','hmm_to','ali_from','ali_to'])
if recovered:
 run(['diamond','makedb','--threads',cpu,'--in',str(out/'02_ALMT_positive_references.faa'),'--db',str(out/'work/ALMT')],'05_diamond_makedb.log')
 run(['diamond','blastp','--threads',cpu,'--query',str(out/'03_recovered_SS_proteins.faa'),'--db',str(out/'work/ALMT'),'--out',str(out/'05_positive_reference_hits.tsv'),'--very-sensitive','--evalue','1e-3','--max-target-seqs','0','--max-hsps','1','--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore'],'05_diamond.log')
u,main_edges=edges(src/'main/tree.treefile');v,sens_edges=edges(src/'sensitivity_all_partials/tree.treefile');need(u<=v,'Sensitivity lost main tips');refs={k for k in u if k.startswith('REF__')}
def support(label):
 try:return tuple(float(x) for x in label.split('/'))
 except ValueError:return ()
def canonical(z):return frozenset(min((z,u-z),key=lambda a:(len(a),sorted(a))))
ma={canonical(z):label for z,label in main_edges};sb=defaultdict(list)
for z,label in sens_edges:
 z=z&u
 if min(len(z),len(u-z))>=2:sb[canonical(z)].append(label)
audit=[];groups=[]
for i,(z,label) in enumerate(sorted(ma.items(),key=lambda x:(len(x[0]),sorted(x[0]))),1):
 labels=sb.get(z,[]);a=support(label);good=[l for l in labels if len(support(l))==2 and support(l)[0]>=80 and support(l)[1]>=95];stable=bool(labels);strong=stable and len(a)==2 and a[0]>=80 and a[1]>=95 and bool(good)
 audit.append(dict(split_id='S%03d'%i,main_support=label,sensitivity_supports=';'.join(labels),backbone_split_retained=stable,strong_in_both=strong,main_side_tips=';'.join(sorted(z))))
 if strong:
  for side in [set(z),u-set(z)]:
   anchors=side&refs
   if anchors and anchors!=refs and side-refs:groups.append(dict(group_id='G%03d'%(len(groups)+1),split_id='S%03d'%i,reference_ids=';'.join(sorted(anchors)),candidate_tree_ids=';'.join(sorted(side-refs)),main_support=label,sensitivity_supports=';'.join(good),interpretation='NESTED_UNROOTED_REFERENCE_ASSOCIATION_NOT_FUNCTIONAL_SUBTYPE'))
table('07_backbone_split_stability.tsv',audit,list(audit[0]));table('08_stable_reference_groups.tsv',groups,['group_id','split_id','reference_ids','candidate_tree_ids','main_support','sensitivity_supports','interpretation'])
meta=read(src/'01_extraction_and_tree_id_map.tsv');byprotein={r['protein_id']:r for r in meta};rows=read(src/'00_corrected_review_flags.tsv');evidence=[]
for r in rows:
 k=r['protein_id'];m=byprotein.get(k);tid=m['tree_id'] if m else '';gg=[g for g in groups if tid in g['candidate_tree_ids'].split(';')];gg.sort(key=lambda g:(len(g['candidate_tree_ids'].split(';')),g['group_id']));best=gg[0] if gg else None
 category='NO_PF11744_REGION_RETAINED_OUTSIDE_TREE' if not m else ('PARTIAL_SENSITIVITY_ONLY' if tid not in u else 'MAIN_TREE_DOMAIN_EVIDENCE')
 architecture='ADDITIONAL_TLD_REGION' if k in {'MC__gene_14995','PO__evm.TU.LG16.459'} else 'NO_TARGETED_TLD_ARCHITECTURE_CAVEAT'
 evidence.append(dict(protein_id=k,species=r['species'],gene_id=r['gene_id'],tree_id=tid,evidence_category=category,architecture_category=architecture,smallest_stable_reference_group=best['group_id'] if best else '',reference_ids=best['reference_ids'] if best else '',retained_completeness_flags=r['flags'],decision='REVIEW_ONLY',final_count='NA'))
table('09_candidate_evidence_categories.tsv',evidence,list(evidence[0]))
besthits={}
for l in (out/'05_positive_reference_hits.tsv').read_text().splitlines():
 x=l.split('\t')
 if x[0] not in besthits or float(x[11])>float(besthits[x[0]][11]):besthits[x[0]]=x
results=[]
for r in isoaudit:
 k=r['alternative_protein_id'];models={h['accession'].split('.')[0] for h in domains[k]};x=besthits.get(k);status='NO_ALMT_SUPPORT_IN_ALTERNATIVE' if 'PF11744' not in models and not x else 'ALTERNATIVE_ALMT_EVIDENCE_REQUIRES_REVIEW'
 results.append(dict(r,models=';'.join(sorted(models)),best_reference=x[1] if x else '',best_evalue=x[10] if x else '',review_outcome=status,limitation='BOUNDED_TEST_OF_ANNOTATED_ISOFORMS;NO_GENOME_REPAIR_OR_ABSENCE_CLAIM'))
table('10_SS_alternative_review_outcomes.tsv',results,list(results[0]))
for n in ['source_prior_01_integrated_locus_evidence.tsv','source_prior_06_refreshed_family_coverage.tsv','source_prior_07_final_Wang_counts_PENDING.tsv','22_protein_CDS_consistency.tsv','source_06_combined_domain_hits.tsv']:shutil.copyfile(src/n,out/n)
run(['hmmsearch','-h'],'hmmer_version.txt');run(['diamond','version'],'diamond_version.txt')
(out/'README.txt').write_text('Records all main unrooted splits and whether retained after pruning sensitivity tree to main tips. Strong reference-containing sides require SHaLRT>=80/UFBoot>=95 in both trees. Groups can overlap/nest; smallest group is descriptive association, not orthology or functional subtype. Partial sensitivity-only records are separate and not assigned a stable main group. Architecture flags and completeness flags retained independently. Exactly two documented alternative SS isoforms screened with full Pfam GA and14 positive references; no-result evidence does not prove genome-wide absence. No new phylogeny, model repair, count integration or automatic locus exclusion.\n')
qc=dict(status='ALMT_GROUPS_ALTERNATIVE_ISOFORMS_REVIEW_READY',candidate_records=len(evidence),main_splits=len(audit),retained_backbone_splits=sum(r['backbone_split_retained'] for r in audit),stable_reference_groups=len(groups),alternative_isoforms=len(recovered),alternative_outcomes=dict(Counter(r['review_outcome'] for r in results)),final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');(out/'PROVENANCE.json').write_text(json.dumps(dict(step62_manifest_sha256=sha(src/'SHA256SUMS.txt'),step61_manifest_sha256=sha(previous/'SHA256SUMS.txt'),pfam_sha256=sha(db),commands=commands),indent=2)+'\n')
for p in out.glob('*'):
 if p.suffix in {'.log','.domtbl'}:p.write_text(p.read_text().replace(str(dbuse),'<PFAM_HMM>').replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts);(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='ALMT_groups_isoform_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
