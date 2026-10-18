#!/usr/bin/env bash
#SBATCH --job-name=almt_tree_review
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=almt_tree_review_%j.out
#SBATCH --error=almt_tree_review_%j.err
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
src=release('step25kz17.61','ALMT_RECOVERY_RELEASE');iq=shutil.which('iqtree2') or shutil.which('iqtree');need(iq,'IQ-TREE executable missing')
ledger=read(src/'08_recalibrated_all_142_review.tsv');need(len(ledger)==142 and all(r['input_status']=='RESOLVED' for r in ledger),'Expected142 resolved candidates');seqs,_=fasta(src/'03_all_resolved_candidates.faa');refs,_=fasta(src/'02_ALMT_positive_references.faa');seqs.update(refs);refids=set(refs);qmeta={r['protein_id']:r for r in ledger};hits=read(src/'06_combined_domain_hits.tsv')
out=base/'step25kz17.62'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False);(out/'annotation').mkdir()
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(n,d):(out/n).write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
for r in ledger:
 flags=[f for f in r['flags'].split(';') if f not in {'MISSING_BEST_REFERENCE_MODEL_PF04632','NO_RECALIBRATED_HEURISTIC_EXCEPTION'}];r['flags']=';'.join(flags) or 'NO_RECALIBRATED_HEURISTIC_EXCEPTION'
table('00_corrected_review_flags.tsv',ledger,list(ledger[0]))
mapping=[];data={};bounds={};idmap={};partial=[];omitted=[]
for j,k in enumerate(sorted(seqs),1):
 hh=[h for h in hits if h['sequence_id']==k and h['accession'].split('.')[0]=='PF11744']
 if not hh:omitted.append(dict(protein_id=k,reason='NO_PF11744_REGION_FOR_HOMOLOGOUS_ALIGNMENT'));continue
 covered=set()
 for h in hh:covered.update(range(int(h['hmm_from']),int(h['hmm_to'])+1))
 cov=len(covered)/int(hh[0]['model_length']);lo=min(int(h['ali_from']) for h in hh);hi=max(int(h['ali_to']) for h in hh);need(1<=lo<=hi<=len(seqs[k]),'Invalid extraction');tid=k if k in refs else 'ALM%03d'%j
 # Extract only detected ALMT region; no unrelated appended region in tree.
 need(not any(h['sequence_id']==k and h['accession'].split('.')[0]=='PF07534' and int(h['ali_from'])<=hi and int(h['ali_to'])>=lo for h in hits),'TLD overlaps extracted ALMT region; review required')
 idmap[k]=tid;data[tid]=seqs[k][lo-1:hi];bounds[tid]=(lo,hi,hh)
 if cov<.65:partial.append(k)
 mapping.append(dict(tree_id=tid,protein_id=k,role='REFERENCE' if k in refs else 'CANDIDATE',species='Arabidopsis_thaliana_REFERENCE' if k in refs else qmeta[k]['species'],gene_id='' if k in refs else qmeta[k]['gene_id'],screen_group='ALMT_REFERENCE' if k in refs else 'ALMT_REVIEW',extraction_start=lo,extraction_end=hi,model_coverage=round(cov,6),caveat='PARTIAL_SENSITIVITY_ONLY' if cov<.65 else 'DOMAIN_REGION_ONLY_REVIEW'))
table('01_extraction_and_tree_id_map.tsv',mapping,list(mapping[0]));table('01_not_in_domain_tree.tsv',omitted,['protein_id','reason']);writefa('02_all_extracted_carrier_regions.faa',data)
need(partial and len(data)-len(partial)>10,'Unexpected main/partial scope')
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
need(len(keep)>=100,'Too few carrier-core columns');table('06_fixed_core_mask.tsv',mask,list(mask[0]));core={k:''.join(v[i] for i in keep) for k,v in added.items()};need(all(v.replace('-','') for v in core.values()),'Empty core sequence')
coreqc=[dict(tree_id=k,core_residues=len(v.replace('-','')),core_columns=len(keep),occupancy=round(len(v.replace('-',''))/len(keep),6)) for k,v in core.items()];table('07_core_occupancy.tsv',coreqc,list(coreqc[0]))
runs=[];model=None
for mode,extras,seed in [('main',[],176201),('sensitivity_all_partials',[idmap[k] for k in partial],176202)]:
 (out/mode).mkdir();ids=set(main)|set(extras);dd={k:v for k,v in core.items() if k in ids};writefa(mode+'/alignment.faa',dd);prefix=out/mode/'tree'
 run([iq,'-s',str(out/mode/'alignment.faa'),'-st','AA','-m',model or 'MFP','-B','1000','-bnni','-alrt','1000','-T',cpu,'-seed',str(seed),'-keep-ident','--prefix',str(prefix)],mode+'/iqtree.stdout.log')
 tree=Path(str(prefix)+'.treefile');need(tree.is_file(),'Missing tree');tips=re.findall(r'(?:\(|,)\s*([^():,;\s]+)\s*:',tree.read_text());need(len(tips)==len(dd) and set(tips)==set(dd),'Tree tip mismatch')
 if model is None:
  report=Path(str(prefix)+'.iqtree').read_text();m=re.search(r'Best-fit model according to BIC:\s*(\S+)',report);need(m,'ModelFinder model missing');model=m.group(1)
 labels=['LABELS','SEPARATOR TAB','DATA']+[r['tree_id']+'\t'+r['protein_id']+' | '+r['screen_group']+' | '+r['caveat'] for r in mapping if r['tree_id'] in ids];(out/mode/'itol_labels.txt').write_text('\n'.join(labels)+'\n')
 runs.append(dict(analysis=mode,sequences=len(dd),columns=len(keep),model=model,seed=seed,support_order='SH-aLRT/UFBoot',rooting='UNROOTED'))
table('08_tree_run_summary.tsv',runs,list(runs[0]))
exceptions={'SLP_Hap1__SP4G00453','SS__g16509.t2','SS__g17774.t2','MC__gene_14995','PO__evm.TU.LG16.459'}
need(exceptions<=set(qmeta),'Expected architecture exceptions missing')
manifest=root/'00_admin/manifests/input_sets.tsv';configs=read(manifest);sources=[];status=[];mappingrows=[];checks=[]
bases='TCAG';letters='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';codons={a+b+c:v for (a,b,c),v in zip(((a,b,c) for a in bases for b in bases for c in bases),letters)}
def op(p):return gzip.open(p,'rt') if p.suffix=='.gz' else p.open()
def collect(p,aliases):
 found={};key=None;parts=[]
 def commit():
  if key in aliases:need(key not in found,'Duplicate source ID');found[key]=''.join(parts).upper()
 with op(p) as f:
  for line in f:
   if line.startswith('>'):commit();key=line[1:].split()[0];parts=[]
   elif key in aliases:parts.append(''.join(line.split()))
  commit()
 return found
for idx,k in enumerate(sorted(exceptions),1):
 r=next(x for x in ledger if x['protein_id']==k);sp=r['species'];gene=r['gene_id'];folder='annotation/E%02d'%idx;(out/folder).mkdir();mp=root/'03_harmonized/primary_v1'/sp/(sp+'.all_isoform_mapping.tsv')
 need(mp.is_file(),'Missing mapping for '+sp);mm=[m for m in read(mp) if m['gene_id']==gene];need(mm,'No exact gene mapping '+k);need(len({(m['seqid'],m['strand']) for m in mm})==1,'Conflicting mapping coordinates');table(folder+'/mapping.tsv',mm,list(mm[0]));sources.append(dict(protein_id=k,kind='mapping',file=str(mp.relative_to(root)),sha256=sha(mp)));aliases={v for m in mm for f in ['gene_id','transcript_id','protein_id','cds_id'] if (v:=m.get(f,'')) and v not in {'.','NA'}};aliases|={sp+'__'+v for v in list(aliases)};need(k in aliases,'Candidate alias mismatch')
 conf=[c for c in configs if c['input_set']==sp];need(len(conf)==1,'Ambiguous manifest species');conf=conf[0];found={}
 for kind,col in [('protein','protein_all'),('cds','cds_all'),('gff','gff3')]:
  value=conf.get(col,'');p=root/value if value else root/'__MISSING_INPUT__';exists=p.is_file();sources.append(dict(protein_id=k,kind=kind,file=str(p.relative_to(root)) if p.is_relative_to(root) else '<EXTERNAL>/'+p.name,sha256=sha(p) if exists else 'MISSING'))
  if not exists:status.append(dict(protein_id=k,kind=kind,status='MISSING',records=0));continue
  if kind!='gff':found[kind]=collect(p,aliases);writefa(folder+'/'+kind+'_isoforms.fasta',found[kind]);n=len(found[kind])
  else:
   chrom=mm[0]['seqid'];lo=max(1,min(int(m['start']) for m in mm)-5000);hi=max(int(m['end']) for m in mm)+5000;n=0
   with op(p) as f,(out/folder/'context.gff3').open('w') as w:
    w.write('##gff-version 3\n')
    for line in f:
     if line.startswith('##FASTA'):break
     if line.startswith('#'):continue
     c=line.rstrip('\n').split('\t')
     if len(c)==9 and c[0]==chrom and int(c[3])<=hi and int(c[4])>=lo:w.write(line.rstrip('\n')+'\n');n+=1
  status.append(dict(protein_id=k,kind=kind,status='COLLECTED' if n else 'NO_EXACT_ID_MATCH',records=n))
 for m in mm:
  def lookup(kind,field):
   ident=m.get(field,'');d=found.get(kind,{});return d.get(ident,d.get(sp+'__'+ident))
  aa=lookup('protein','protein_id');dna=lookup('cds','cds_id');translation=''.join(codons.get(dna[i:i+3],'X') for i in range(0,len(dna)-len(dna)%3,3)) if dna else None
  checks.append(dict(candidate_id=k,mapped_transcript=m.get('transcript_id',''),protein_id=m.get('protein_id',''),cds_id=m.get('cds_id',''),source_protein_found=aa is not None,source_CDS_found=dna is not None,protein_matches_candidate=norm(aa)==seqs[k] if aa else 'NA',cds_nt=len(dna) if dna else 'NA',frame_remainder=len(dna)%3 if dna else 'NA',translation_matches_source_protein=norm(translation)==norm(aa) if translation and aa else 'NA',translation_matches_candidate=norm(translation)==seqs[k] if translation else 'NA',internal_stops=norm(translation).count('*') if translation else 'NA',ambiguous_amino_acids=translation.count('X') if translation else 'NA',start_codon=dna[:3] if dna else 'NA',terminal_codon=dna[-3:] if dna else 'NA',interpretation='ANNOTATION_CONSISTENCY_NOT_GENOME_OR_READ_VALIDATION'))
 mappingrows.append(dict(protein_id=k,species=sp,gene_id=gene,folder=folder))
table('20_annotation_file_map.tsv',mappingrows,list(mappingrows[0]));table('21_annotation_collection_status.tsv',status,list(status[0]));table('22_protein_CDS_consistency.tsv',checks,list(checks[0]));table('23_annotation_source_provenance.tsv',sources,list(sources[0]))
for n in ['prior_01_integrated_locus_evidence.tsv','prior_06_refreshed_family_coverage.tsv','prior_07_final_Wang_counts_PENDING.tsv','06_combined_domain_hits.tsv']:shutil.copyfile(src/n,out/('source_'+n))
run(['mafft','--version'],'mafft_version.txt');run([iq,'--version'],'iqtree_version.txt')
(out/'README.txt').write_text('Combined ALMT review: PF04632-only absence flags removed. Extract contiguous first-to-last PF11744 hit region, preserve intervening sequence; no appended TLD region included. Main includes models with>=65% HMM union coverage; lower coverage candidates/reference retained in joint partial sensitivity only. No-PF11744 candidates explicitly omitted, not declared absent genes. MAFFT main LINSI; addfragments keeplength for partials. Fixed columns supported by>=3 main positive reference domain spans and>=50% main candidate occupancy. Main ModelFinder model reused in sensitivity;1000SHaLRT/1000UFBoot bnni, identical tips retained, unrooted. Short reference Q3E9Z9 treated by same partial rule. Placement of partials uncertain; joint sensitivity can be influenced by multiple fragments. Limited Arabidopsis panel does not establish subtype or transport function. Five exceptions collected through exact mapped isoforms, source protein/CDS comparison and local GFF context; CDS consistency does not validate gene models. No counts or model repairs.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(source_manifest_sha256=sha(src/'SHA256SUMS.txt'),input_manifest_sha256=sha(manifest),commands=commands),indent=2)+'\n');qc=dict(status='ALMT_PHYLOGENY_ANNOTATION_REVIEW_READY',historical_candidates=142,main_sequences=len(main),partial_sequences=len(partial),omitted_no_domain=len(omitted),core_columns=len(keep),tree_runs=runs,annotation_status=status,final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
for p in out.rglob('*'):
 if p.is_file() and p.suffix in {'.log','.iqtree','.txt'}:p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out);bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='ALMT_phylogeny_annotation_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
