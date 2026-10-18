#!/usr/bin/env bash
#SBATCH --job-name=tdt_dit_tree
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=tdt_dit_tree_%j.out
#SBATCH --error=tdt_dit_tree_%j.err
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

src=release('step25kz17.47','TDT_REVIEW_RELEASE');screen=release('step25kz17.46','TRANSPORTER_RECOVERY_RELEASE')
need(json.loads((src/'PROVENANCE.json').read_text())['source_manifest_sha256']==sha(screen/'SHA256SUMS.txt'),'Source release chain mismatch')
iq=next((shutil.which(x) for x in ['iqtree3','iqtree2','iqtree'] if shutil.which(x)),None);need(iq,'IQ-TREE executable not found')
rows=read(screen/'07_all_candidate_evidence_REVIEW_ONLY.tsv');selected=[r for r in rows if r['best_panel_family'] in {'TDT','DIT1','DIT2'}]
need(Counter(r['best_panel_family'] for r in selected)=={'TDT':11,'DIT1':18,'DIT2':21},'Unexpected candidate panel')
need(all(r['mapping_status']=='UNAMBIGUOUS' for r in selected),'Ambiguous locus mapping')
need(len({(r['species'],r['gene_id']) for r in selected})==len(selected),'Duplicate locus')
seqs,_=fasta(screen/'05_all_discovered_and_prior_candidates.faa');refs,_=fasta(screen/'02_reference_panel.faa');meta=[r for r in read(screen/'01_reference_definitions.tsv') if r['family'] in {'TDT','DIT1','DIT2'}]
need(len(meta)==3,'Expected three curated anchors')
for r in meta:need(hashlib.sha256(refs[r['reference_id']].encode()).hexdigest()==r['sequence_sha256'],'Reference hash mismatch')
tdtseq,_=fasta(src/'01_TDT_candidates.faa');need(all(seqs[k]==v for k,v in tdtseq.items()),'TDT sequence provenance mismatch')
short='SP_Hex__gene_Spor0444110.1';need(short in tdtseq,'Short candidate missing')
out=base/'step25kz17.48'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
def table(n,rr,ff):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,ff,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def writefa(p,d):p.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in d.items()))
mapping=[];data={}
for i,r in enumerate(sorted(selected,key=lambda x:x['protein_id']),1):
 k='TD%03d'%i;data[k]=seqs[r['protein_id']];mapping.append(dict(tree_id=k,protein_id=r['protein_id'],species=r['species'],gene_id=r['gene_id'],screen_group=r['best_panel_family'],role='CANDIDATE',caveat='SHORT_N_TERMINUS_CDS_CONSISTENT' if r['protein_id']==short else 'SCREEN_ASSIGNMENT_NOT_VALIDATED'))
for r in meta:
 k=r['reference_id'];data[k]=refs[k];mapping.append(dict(tree_id=k,protein_id=k,species='Arabidopsis_thaliana_REFERENCE',gene_id=r['locus'],screen_group=r['family'],role='CURATED_REFERENCE',caveat='LIMITED_REFERENCE_PANEL'))
need(all(valid(s) for s in data.values()),'Invalid protein sequence')
shortid=next(r['tree_id'] for r in mapping if r['protein_id']==short);refids={r['reference_id'] for r in meta};candids=set(data)-refids
writefa(out/'01_all_sequences.faa',data);table('02_tree_id_map.tsv',mapping,list(mapping[0]));table('03_candidate_screen_evidence.tsv',selected,list(selected[0]));table('04_reference_definitions.tsv',meta,list(meta[0]))
commands=[];cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
def run(args,log,stdout=None):
 commands.append([str(x).replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>') for x in args])
 with (out/log).open('w') as f:
  if stdout:
   with (out/stdout).open('w') as a:subprocess.run(args,stdout=a,stderr=f,check=True)
  else:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
run(['mafft','--version'],'mafft_version.txt');run([iq,'--version'],'iqtree_version.txt')
run(['mafft','--amino','--thread',cpu,'--localpair','--maxiterate','1000',str(out/'01_all_sequences.faa')],'05_mafft.log','05_all_alignment.faa')
aln,_=fasta(out/'05_all_alignment.faa');need(set(aln)==set(data) and len({len(s) for s in aln.values()})==1,'Invalid alignment')
need(all(aln[k].replace('-','')==s for k,s in data.items()),'Alignment changed sequence identities')
# One fixed mask for both trees: avoids confounding taxon removal with realignment.
mask=[];kept=[]
for i in range(len(next(iter(aln.values())))):
 ref_n=sum(aln[k][i]!='-' for k in refids);occ=sum(aln[k][i]!='-' for k in candids)/len(candids);keep=ref_n>0 and occ>=.5
 if keep:kept.append(i)
 mask.append(dict(alignment_column=i+1,reference_residues=ref_n,candidate_occupancy=round(occ,6),retained=int(keep)))
need(len(kept)>=100,'Insufficient shared alignment core')
table('06_alignment_mask.tsv',mask,list(mask[0]));core={k:''.join(v[i] for i in kept) for k,v in aln.items()};need(all(any(c!='-' for c in s) for s in core.values()),'Empty core sequence')
summary=[]
for mode,omit,seed in [('main_without_short',True,174801),('sensitivity_with_short',False,174802)]:
 folder=out/mode;folder.mkdir();dd={k:v for k,v in core.items() if not(omit and k==shortid)};writefa(folder/'alignment.faa',dd)
 prefix=folder/'tree';run([iq,'-s',str(folder/'alignment.faa'),'-st','AA','-m','MFP','-B','1000','-bnni','-alrt','1000','-T',cpu,'-seed',str(seed),'-keep-ident','--prefix',str(prefix)],mode+'/iqtree.stdout.log')
 tree=Path(str(prefix)+'.treefile');need(tree.is_file(),'Missing tree');tips=re.findall(r'(?:\(|,)\s*([^():,;\s]+)\s*:',tree.read_text());need(len(tips)==len(dd) and set(tips)==set(dd),'Tree tip mismatch')
 labels=['LABELS','SEPARATOR TAB','DATA']+[r['tree_id']+'\t'+r['protein_id']+' | '+r['screen_group']+' | '+r['role']+(' | SHORT_MODEL' if r['tree_id']==shortid else '') for r in mapping if r['tree_id'] in dd];(folder/'itol_labels.txt').write_text('\n'.join(labels)+'\n')
 summary.append(dict(analysis=mode,sequences=len(dd),core_columns=len(kept),short_candidate_included=not omit,seed=seed,tree_rooting='UNROOTED',support_order='SH-aLRT/UFBoot',formal_counts='NA'))
table('07_tree_run_summary.tsv',summary,list(summary[0]))
# Preserve the bounded CDS check as reproducible evidence.
cds,_=fasta(src/'annotation/03_cds_isoforms.fasta');prots,_=fasta(src/'annotation/03_protein_isoforms.fasta');need(len(cds)==len(prots)==1,'Unexpected annotation isoforms')
dna=next(iter(cds.values()));aa=next(iter(prots.values()));bases='TCAG';letters='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';codons={a+b+c:v for (a,b,c),v in zip(((a,b,c) for a in bases for b in bases for c in bases),letters)};translation=''.join(codons.get(dna[i:i+3],'X') for i in range(0,len(dna),3))
check=dict(protein_id=short,cds_length=len(dna),protein_length=len(aa),frame_remainder=len(dna)%3,start_codon=dna[:3],stop_codon=dna[-3:],internal_stops=norm(translation).count('*'),ambiguous_residues=translation.count('X'),translation_matches=norm(translation)==norm(aa),candidate_matches=norm(aa)==tdtseq[short],interpretation='ANNOTATION_CONSISTENCY_ONLY_NOT_GENOME_OR_READ_VALIDATION');table('08_short_candidate_CDS_check.tsv',[check],list(check))
for n in ['prior_01_integrated_locus_evidence.tsv','prior_06_refreshed_family_coverage.tsv','prior_07_final_Wang_counts_PENDING.tsv','08_TDT_completeness_review.tsv','07_domain_model_coverage.tsv','10_mitochondrial_carrier_queue_UNCHANGED.tsv']:shutil.copyfile(src/n,out/('source_'+n))
(out/'README.txt').write_text('Bounded TDT-versus-DIT phylogenetic review. Includes all 50 TDT/DIT best-panel candidates from corrected-proteome screening and three curated Arabidopsis anchors. DIT candidate labels are provisional comparison groups, not validated functions. Main tree: 49 candidates plus three anchors, excluding shorter SP_Hex model. Sensitivity tree: all 50 plus three anchors. Same MAFFT alignment and fixed reference-supported, >=50% candidate occupancy mask used for both trees. ModelFinder selects models independently; seeds recorded. Trees are unrooted; internal supports SH-aLRT/UFBoot. Review separation of TDT-reference-containing group from both DIT anchors, support and placement stability before ledger decisions. Limited references, shared Pfam models, and a phylogenetic grouping do not prove transport specificity or C4/CAM function. Short-model CDS consistency is not genome/read validation. No counts updated. Upload each tree.treefile to iTOL followed by its itol_labels.txt.\n')
(out/'PROVENANCE.json').write_text(json.dumps(dict(step47_manifest_sha256=sha(src/'SHA256SUMS.txt'),step46_manifest_sha256=sha(screen/'SHA256SUMS.txt'),commands=commands),indent=2)+'\n')
qc=dict(status='TDT_DIT_PHYLOGENY_REVIEW_READY',candidate_groups=dict(Counter(r['best_panel_family'] for r in selected)),references=len(meta),raw_columns=len(next(iter(aln.values()))),core_columns=len(kept),trees=summary,final_counts=False);(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
for p in out.rglob('*'):
 if p.is_file() and p.suffix in {'.log','.iqtree','.txt'}:p.write_text(p.read_text().replace(str(out),'<RELEASE>').replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>'))
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='TDT_DIT_phylogenetic_review/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle.relative_to(root)))
PY
