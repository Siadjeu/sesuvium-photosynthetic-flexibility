#!/usr/bin/env bash
#SBATCH --job-name=tdt_dic_screen
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=tdt_dic_screen_%j.out
#SBATCH --error=tdt_dic_screen_%j.err
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
src=release('step25kz17.44','INTEGRATED_RELEASE');cor=release('step25kz17.21','CORRECTED_INPUT_RELEASE')
proteome=cor/'full_outputs/corrected_protein_candidates.faa'
exp=[r['sha256'] for r in read(cor/'08_full_output_checksums.tsv') if r['file']=='full_outputs/corrected_protein_candidates.faa']
need(len(exp)==1 and sha(proteome)==exp[0],'Corrected protein checksum mismatch')
manifest=read(src/'09_pooled_source_manifest.tsv')
def snapshot(suffix):
 rr=[r for r in manifest if r['source'].endswith(suffix)];need(len(rr)==1,'Missing/ambiguous pooled source '+suffix);return src/rr[0]['snapshot']
registry=read(snapshot('/candidate_validation.tsv'));old,_=fasta(snapshot('/TDT_DTC_DIC.accepted.proteins.fa'))
need(len(registry)==35 and len(old)==34,'Expected 35 registry records and 34 historical proteins')
out=base/'step25kz17.45'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False);(out/'references').mkdir();(out/'work').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+s+'\n' for k,s in d.items()))
# Family definitions anchored to exact Arabidopsis loci, not ambiguous gene names.
# DIC1-3 have historical PUMP/UCP synonyms. Arabidopsis controls are not project AT.
panel=[('Q8LG88','AT5G47560','TDT','TARGET'),('Q9C5M0','AT5G19760','DTC','TARGET'),('Q9SJY5','AT2G22500','DIC','TARGET'),('Q9SB52','AT4G24570','DIC','TARGET'),('Q9FY68','AT5G09470','DIC','TARGET'),('Q9LXV3','AT5G12860','DIT1','COMPARATOR'),('Q9FMF7','AT5G64290','DIT2','COMPARATOR'),('O81845','AT3G54110','UCP1','COMPARATOR')]
# Optional offline cache: accession.json files from UniProt REST, identity checked below.
refs={};refmeta=[]
for accession,locus,family,role in panel:
 url='https://rest.uniprot.org/uniprotkb/'+accession+'.json';raw=None
 if os.environ.get('REFERENCE_JSON_DIR'):
  raw=(root/os.environ['REFERENCE_JSON_DIR']/(accession+'.json')).read_bytes()
 else:
  for attempt in range(3):
   try:
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'bounded-family-review/1.0'}),timeout=45) as response:raw=response.read()
    break
   except Exception:
    if attempt==2:raise SystemExit('Reference retrieval failed for '+accession+'; set REFERENCE_JSON_DIR to an offline UniProt JSON cache, or retry with network access')
    time.sleep(2)
  
 data=json.loads(raw);need(data.get('primaryAccession')==accession and data.get('organism',{}).get('taxonId')==3702,'Reference accession/species mismatch '+accession)
 loci={x['value'].upper() for g in data.get('genes',[]) for x in g.get('orderedLocusNames',[])}
 need(locus in loci,'Reference locus mismatch '+accession)
 s=data['sequence']['value'].upper();need(valid(s),'Invalid reference sequence '+accession)
 (out/'references'/(accession+'.json')).write_bytes(raw);rid='REF__'+accession;refs[rid]=s
 refmeta.append(dict(reference_id=rid,accession=accession,locus=locus,family=family,role=role,length_aa=len(s),entry_type=data.get('entryType',''),source_url=url,json_sha256=hashlib.sha256(raw).hexdigest(),sequence_sha256=hashlib.sha256(s.encode()).hexdigest(),definition_status='CURATED_ANCHOR_NOT_CANDIDATE_FUNCTION_PROOF'))
table('01_reference_definitions.tsv',refmeta,list(refmeta[0]));writefa('02_target_and_comparator_references.faa',refs)
# Exact identifier reconciliation. No suffix stripping or sequence-only locus merging.
mappings={};maphash={};aliases={};mappinginfo={};wanted=set()
for r in registry:
 sp,sep,k=r['sequence_id'].partition('|');need(sep and k.startswith(sp+'__'),'Unexpected pooled ID format')
 if sp not in mappings:
  mp=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv');mappings[sp]=read(mp);maphash[sp]=sha(mp)
 local=k.split('__',1)[1];mm=[m for m in mappings[sp] if local in {m.get(f,'') for f in ['gene_id','transcript_id','protein_id','cds_id']}]
 genes={m['gene_id'] for m in mm};mappinginfo[k]=(sp,genes)
 aa={k}
 if len(genes)==1:
  for m in mappings[sp]:
   if m['gene_id'] in genes:
    aa|={sp+'__'+m[f] for f in ['gene_id','transcript_id','protein_id','cds_id'] if m.get(f) and m[f] not in {'.','NA'}}
 aliases[k]=aa;wanted|=aa
# Stream corrected proteome retaining only requested exact aliases.
selected={};key=None;parts=[]
def commit():
 if key in wanted:
  need(key not in selected,'Duplicate corrected protein ID');selected[key]=''.join(parts).upper()
with proteome.open() as f:
 for line in f:
  if line.startswith('>'):commit();key=line[1:].split()[0];parts=[]
  elif key in wanted:parts.append(''.join(line.split()))
 commit()
audit=[];queries={};canonical_to_source=defaultdict(list)
for r in registry:
 sp,_,k=r['sequence_id'].partition('|');genes=mappinginfo[k][1];matched=sorted(aliases[k]&set(selected));historical=old.get(r['sequence_id'],'')
 status='UNAMBIGUOUS_CORRECTED_PROTEIN' if len(matched)==1 and len(genes)==1 else 'MAPPING_OR_SEQUENCE_UNRESOLVED'
 ck=matched[0] if len(matched)==1 else '';s=selected.get(ck,'');same=bool(historical and s and norm(historical)==norm(s))
 if status=='UNAMBIGUOUS_CORRECTED_PROTEIN':
  need(valid(norm(s)),'Invalid corrected protein');queries[ck]=norm(s);canonical_to_source[ck].append(k)
 audit.append(dict(historical_id=k,species=sp,gene_ids=';'.join(sorted(genes)),registry_status=r['status'],corrected_id=ck,corrected_length=len(s) if s else '',historical_length=len(historical) if historical else '',historical_sequence_matches_corrected=same if historical and s else 'NOT_COMPARABLE',status=status,family_assignment='PENDING'))
table('03_candidate_provenance.tsv',audit,list(audit[0]));need(queries,'No valid corrected candidates')
writefa('04_corrected_pooled_candidates.faa',queries)
# Installed complete Pfam, same locations as prior family steps.
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
def run(args,name):
 commands.append([str(x).replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for x in args])
 with (out/name).open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
writefa('05_candidates_and_controls.faa',{**queries,**refs})
run(['hmmsearch','--cpu',cpu,'--cut_ga','--noali','--domtblout',str(out/'06_pfam.domtbl'),str(dbuse),str(out/'05_candidates_and_controls.faa')],'06_hmmsearch.log')
run(['diamond','makedb','--in',str(out/'02_target_and_comparator_references.faa'),'--db',str(out/'work/reference_panel'),'--threads',cpu],'07_diamond_makedb.log')
run(['diamond','blastp','--query',str(out/'04_corrected_pooled_candidates.faa'),'--db',str(out/'work/reference_panel'),'--out',str(out/'07_similarity.tsv'),'--threads',cpu,'--very-sensitive','--evalue','1e-3','--max-target-seqs','0','--max-hsps','1','--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore'],'07_diamond.log')
domains=defaultdict(list);dh=[]
for line in (out/'06_pfam.domtbl').read_text().splitlines():
 if not line or line.startswith('#'):continue
 x=line.split(maxsplit=22);r=dict(sequence_id=x[0],model=x[3],accession=x[4],model_length=int(x[5]),i_evalue=float(x[12]),score=float(x[13]),hmm_from=int(x[15]),hmm_to=int(x[16]),ali_from=int(x[17]),ali_to=int(x[18]));dh.append(r);domains[x[0]].append(r)
table('08_all_domain_hits.tsv',dh,['sequence_id','model','accession','model_length','i_evalue','score','hmm_from','hmm_to','ali_from','ali_to'])
rd={r['reference_id']:r for r in refmeta};sim=defaultdict(list)
for line in (out/'07_similarity.tsv').read_text().splitlines():
 x=line.split('\t');need(len(x)==12,'Unexpected DIAMOND output')
 sim[x[0]].append(dict(reference=x[1],family=rd[x[1]]['family'],role=rd[x[1]]['role'],pident=float(x[2]),qcov=(int(x[7])-int(x[6])+1)/int(x[4]),scov=(int(x[9])-int(x[8])+1)/int(x[5]),evalue=float(x[10]),bitscore=float(x[11])))
review=[]
for k,s in queries.items():
 hh=sorted(sim[k],key=lambda x:(-x['bitscore'],x['reference']));best=hh[0] if hh else None;other=next((x for x in hh if x['family']!=best['family']),None) if best else None
 refmodels={d['accession'].split('.')[0] for d in domains[best['reference']]} if best else set();candidate_models={d['accession'].split('.')[0] for d in domains[k]};shared=refmodels&candidate_models
 margin=best['bitscore']-other['bitscore'] if other else ''
 flags=[]
 if not best:flags.append('NO_PANEL_HIT')
 else:
  if best['role']=='COMPARATOR':flags.append('COMPARATOR_FAVOURED')
  if best['qcov']<.7 or best['scov']<.7:flags.append('LIMITED_SINGLE_HSP_COVERAGE')
  if not shared:flags.append('NO_SHARED_REFERENCE_PFAM_MODEL')
  if not other:flags.append('NO_OTHER_FAMILY_HIT_IN_LIMITED_PANEL')
  elif margin<20:flags.append('SMALL_BETWEEN_FAMILY_SCORE_MARGIN')
 flags.append('LIMITED_REFERENCE_PANEL_REQUIRES_REVIEW')
 review.append(dict(protein_id=k,historical_ids=';'.join(canonical_to_source[k]),length_aa=len(s),best_panel_family=best['family'] if best else '',best_reference=best['reference'] if best else '',best_bitscore=best['bitscore'] if best else '',best_evalue=best['evalue'] if best else '',query_span=round(best['qcov'],5) if best else '',reference_span=round(best['scov'],5) if best else '',next_family=other['family'] if other else '',between_family_margin=margin,shared_reference_Pfam_models=';'.join(sorted(shared)),all_candidate_models=';'.join(sorted(candidate_models)),review_flags=';'.join(flags),formal_family_assignment='PENDING',final_gene_count='NA'))
table('09_family_specificity_REVIEW_ONLY.tsv',review,list(review[0]))
controls=[dict(reference_id=k,family=rd[k]['family'],Pfam_models=';'.join(sorted({h['accession']+':'+h['model'] for h in domains[k]})),domain_hits=len(domains[k])) for k in refs]
table('10_reference_domain_controls.tsv',controls,list(controls[0]))
for n in ['01_integrated_locus_evidence.tsv','06_refreshed_family_coverage.tsv','07_final_Wang_counts_PENDING.tsv']:shutil.copyfile(src/n,out/('prior_'+n))
for tool,args in [('diamond',['version']),('hmmsearch',['-h'])]:run([tool]+args,tool+'_version.txt')
(out/'PROVENANCE.json').write_text(json.dumps(dict(prior_manifest_sha256=sha(src/'SHA256SUMS.txt'),corrected_proteome_sha256=exp[0],mapping_hashes=maphash,pfam_sha256=sha(db),pfam_models=models,reference_retrieved_UTC=datetime.now(timezone.utc).isoformat(),commands=commands),indent=2)+'\n')
(out/'README.txt').write_text('TDT: tonoplast dicarboxylate transporter anchored to At5g47560. DTC: mitochondrial dicarboxylate/tricarboxylate carrier At5g19760. DIC: mitochondrial dicarboxylate carrier group At2g22500/At4g24570/At5g09470, also called PUMP/UCP isoforms in UniProt. DIT1/DIT2 chloroplast transporters and UCP1 are comparators, not targets for these three counts. Arabidopsis controls are not project species AT. Identity checked by accession, taxon 3702 and ordered locus. Sources: https://www.ncbi.nlm.nih.gov/gene/834806 ; https://www.uniprot.org/uniprotkb/Q9C5M0/entry ; https://pubmed.ncbi.nlm.nih.gov/18039180/ . Panel is limited, not exhaustive. Shared mitochondrial-carrier domains cannot distinguish DIC from DTC or other carriers. Family labels are best-panel-hit suggestions only, not accepted assignments. 70% single-HSP spans and 20-bit margins are review flags only. All 35 historical registry entries retained, including prior failure; only unambiguously mapped corrected proteins screened. Missing/changed sequences and mapping ambiguity remain explicit, never silently rescued. No full-proteome completeness search, tree, localization or function validation. No pooled candidate counted three times. Existing integrated ledger unchanged; final counts NA.\n')
qc=dict(status='TDT_DTC_DIC_BOUNDED_SCREEN_REVIEW_REQUIRED',historical_registry_records=35,historical_accepted_records=34,corrected_screened_proteins=len(queries),provenance_status=dict(Counter(r['status'] for r in audit)),best_panel_families=dict(Counter(r['best_panel_family'] for r in review)),reference_controls=len(refs),reference_controls_without_Pfam_hits=sum(not domains[k] for k in refs),whole_genome_completeness='NOT_TESTED',final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
for p in out.glob('*'):
 if p.suffix in {'.log','.domtbl'}:p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(dbuse),'<PFAM_HMM>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file() and 'work' not in p.relative_to(out).parts)
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='TDT_DTC_DIC_definition_screen/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
