#!/usr/bin/env bash
#SBATCH --job-name=almt_recover
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=almt_recover_%j.out
#SBATCH --error=almt_recover_%j.err
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
src=release('step25kz17.60','ALMT_SCREEN_RELEASE');cor=release('step25kz17.21','CORRECTED_INPUT_RELEASE');prov=json.loads((src/'PROVENANCE.json').read_text());need(sha(cor/'SHA256SUMS.txt')==prov['corrected_manifest_sha256'],'Corrected release changed')
proteome=cor/'full_outputs/corrected_protein_candidates.faa';need(sha(proteome)==prov['corrected_proteome_sha256'],'Corrected proteome changed')
rows=read(src/'08_combined_ALMT_review.tsv');unresolved=[r for r in rows if r['input_status']!='READY_FOR_SCREEN'];need(len(rows)==142 and len(unresolved)==21 and all(r['species']=='SS' for r in unresolved),'Unexpected recovery scope')
mp=root/'03_harmonized/primary_v1/SS/SS.all_isoform_mapping.tsv';need(sha(mp)==prov['mapping_hashes']['SS'],'SS mapping changed');maps=read(mp);wanted=set();aliases={};selectedmaps=[]
for r in unresolved:
 mm=[m for m in maps if m['gene_id']==r['gene_id']];need(mm,'Missing exact gene mapping');selectedmaps.extend(mm);aa={}
 for m in mm:
  for field in ['protein_id','transcript_id']:
   v=m.get(field,'')
   if v and v not in {'.','NA'}:
    for a in {v,'SS__'+v}:
     aa[a]=m['gene_id'];wanted.add(a)
 aliases[r['protein_id']]=aa
current={};key=None;parts=[]
def commit():
 if key in wanted:need(key not in current,'Duplicate corrected FASTA identifier');current[key]=''.join(parts).upper()
with proteome.open() as f:
 for l in f:
  if l.startswith('>'):commit();key=l[1:].split()[0];parts=[]
  elif key in wanted:parts.append(''.join(l.split()))
 commit()
out=base/'step25kz17.61'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'work').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
recovered={};audit=[];idmap={}
for r in unresolved:
 k=r['protein_id'];matches=sorted(a for a in aliases[k] if a in current);status='UNRESOLVED_NO_MAPPED_PROTEIN' if not matches else 'UNRESOLVED_MULTIPLE_MAPPED_PROTEINS';chosen=''
 if len(matches)==1:
  chosen=matches[0];status='RECOVERED_EXACT_DOCUMENTED_ALIAS' if valid(norm(current[chosen])) else 'UNRESOLVED_INVALID_PROTEIN'
  if status.startswith('RECOVERED'):need(chosen not in recovered,'Same recovered protein assigned to multiple genes');recovered[chosen]=norm(current[chosen]);idmap[k]=chosen
 audit.append(dict(historical_id=k,gene_id=r['gene_id'],matched_corrected_ids=';'.join(matches),resolved_protein_id=chosen,status=status))
table('01_SS_identifier_recovery.tsv',audit,list(audit[0]));table('02_SS_documented_mapping.tsv',selectedmaps,list(selectedmaps[0]));writefa('03_recovered_SS_proteins.faa',recovered)
oldseq,_=fasta(src/'03_reconciled_candidates.faa');positive,_=fasta(src/'02_ALMT_positive_references.faa');need(not(set(oldseq)&set(recovered)),'Recovered sequence already screened');writefa('03_all_resolved_candidates.faa',{**oldseq,**recovered});shutil.copyfile(src/'02_ALMT_positive_references.faa',out/'02_ALMT_positive_references.faa')
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
allhits=read(src/'06_all_Pfam_domain_hits.tsv')+dh
for h in allhits:
 for f in ['model_length','hmm_from','hmm_to','ali_from','ali_to']:h[f]=int(h[f])
domains=defaultdict(list)
for h in allhits:domains[h['sequence_id']].append(h)
table('06_combined_domain_hits.tsv',allhits,list(allhits[0]))
lines=(src/'05_positive_reference_hits.tsv').read_text().splitlines()+((out/'05_positive_reference_hits.tsv').read_text().splitlines() if recovered else []);(out/'05_combined_positive_reference_hits.tsv').write_text('\n'.join(lines)+'\n');best={}
for l in lines:
 x=l.split('\t')
 if x[0] not in best or float(x[11])>float(best[x[0]][11]):best[x[0]]=x
coverage={};layouts={}
for k,hh in domains.items():
 for model in {h['accession'].split('.')[0] for h in hh}:
  hits=[h for h in hh if h['accession'].split('.')[0]==model];covered=set()
  for h in hits:covered.update(range(h['hmm_from'],h['hmm_to']+1))
  coverage[k,model]=len(covered)/hits[0]['model_length'];layouts[k,model]=';'.join(str(h['ali_from'])+'-'+str(h['ali_to']) for h in hits)
refmodels={m for k,m in coverage if k in positive};calibration=[]
for k in positive:
 for m in sorted(refmodels):calibration.append(dict(reference_id=k,model=m,coverage=round(coverage.get((k,m),0),6),protein_intervals=layouts.get((k,m),'')))
table('07_reference_domain_calibration.tsv',calibration,list(calibration[0]));review=[]
for r in rows:
 old=r['protein_id'];k=idmap.get(old,old);available=k in oldseq or k in recovered;x=best.get(k);mm={h['accession'].split('.')[0] for h in domains[k]};flags=[];ref=x[1] if x else '';comparisons=[]
 if not available:flags.append('UNRESOLVED_INPUT')
 else:
  if not x:flags.append('NO_POSITIVE_REFERENCE_SIMILARITY')
  else:
   if (int(x[9])-int(x[8])+1)/int(x[5])<.8:flags.append('REFERENCE_SPAN_LT_80_PERCENT')
   for m in sorted({h['accession'].split('.')[0] for h in domains[ref]}):
    rc=coverage.get((ref,m),0);cc=coverage.get((k,m),0);comparisons.append(m+':'+str(round(cc/rc,4)) if rc else m+':NA')
    if m not in mm:flags.append('MISSING_BEST_REFERENCE_MODEL_'+m)
    elif rc and cc/rc<.8:flags.append('COVERAGE_RELATIVE_TO_REFERENCE_LT_80_PERCENT_'+m)
  if mm-refmodels:flags.append('MODEL_NOT_SEEN_IN_POSITIVE_PANEL')
  if 'PF11744' not in mm:flags.append('NO_PF11744_DETECTED')
 if r['historical_status']!='PASS':flags.append('HISTORICAL_FAILURE_RETAINED')
 review.append(dict(historical_id=old,protein_id=k if available else '',species=r['species'],gene_id=r['gene_id'],input_status='RESOLVED' if available else 'UNRESOLVED',length_aa=len((recovered.get(k) or oldseq.get(k,''))),historical_status=r['historical_status'],best_reference=ref,models=';'.join(sorted(mm)),model_coverage_ratios_to_best_reference=';'.join(comparisons),flags=';'.join(flags) or 'NO_RECALIBRATED_HEURISTIC_EXCEPTION',decision='REVIEW_ONLY',final_count='NA'))
table('08_recalibrated_all_142_review.tsv',review,list(review[0]));table('09_remaining_review_queue.tsv',[r for r in review if r['flags']!='NO_RECALIBRATED_HEURISTIC_EXCEPTION'],list(review[0]))
for n in ['prior_01_integrated_locus_evidence.tsv','prior_06_refreshed_family_coverage.tsv','prior_07_final_Wang_counts_PENDING.tsv','01_reference_role_audit.tsv','02_other_reference_controls.faa']:shutil.copyfile(src/n,out/n)
run(['hmmsearch','-h'],'hmmer_version.txt');run(['diamond','version'],'diamond_version.txt')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step60_manifest_sha256=sha(src/'SHA256SUMS.txt'),corrected_proteome_sha256=sha(proteome),SS_mapping_sha256=sha(mp),pfam_sha256=sha(db),commands=commands),indent=2)+'\n')
qc=dict(status='ALMT_SS_RECOVERY_RECALIBRATION_READY',historical_records=142,recovered_SS=len(recovered),unresolved_SS=21-len(recovered),resolved_candidates=len(oldseq)+len(recovered),remaining_flagged_records=sum(r['flags']!='NO_RECALIBRATED_HEURISTIC_EXCEPTION' for r in review),final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'README.txt').write_text('Exact SS gene-to-documented-protein/transcript alias recovery; no guessed suffixes or longest-isoform selection. Multiple matches remain unresolved. Reuses121 prior candidate results and reference controls; only recovered SS proteins searched with identical Pfam database. Coverage compared to best DIAMOND positive reference per model;80% ratios and HSP span are descriptive review heuristics, not acceptance thresholds. PF04632 occurrence anywhere in positive panel no longer automatically abnormal; missing model relative to best reference retained explicitly. Overlapping models are not independent domains; full intervals/raw hits preserved. All142 historical records and failures retained. No subtype/function/final-count claims.\n')
for p in out.glob('*'):
 if p.suffix in {'.log','.domtbl'}:p.write_text(p.read_text().replace(str(dbuse),'<PFAM_HMM>').replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts);(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='ALMT_SS_recovery_recalibration/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
