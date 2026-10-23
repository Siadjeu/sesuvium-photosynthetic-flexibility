#!/usr/bin/env python3
"""Exon-aware correction of historical PEPC references and reproducible trees."""
import itertools,os,sys,re,json,math,shutil,subprocess,tarfile,collections,urllib.request,urllib.parse,xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime,timezone
import tree_utils as U
from subtype_utils import evidence
BASE='08_focused_families/downstream_supported/pepc_reference_review'
TREE=BASE+'/tree_itol_20835576_20261022T211258984458Z'
AUDIT=BASE+'/e1_reference_audit_20837851_20261022T232648722417Z'
RAW=['08_focused_families/validation/step25kz12/references/pepc_36_genes.fasta','08_focused_families/validation/step25kz12/references/pepc_wang_recovered_mixed.fasta']
CODE=dict(zip((a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'),'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
IUPAC=dict(A='A',C='C',G='G',T='T',R='AG',Y='CT',S='GC',W='AT',K='GT',M='AC',B='CGT',D='AGT',H='ACT',V='ACG',N='ACGT')
def translate(s):
 result=[]
 for i in range(0,len(s)-2,3):
  codon=s[i:i+3];values={CODE[''.join(c)] for c in itertools.product(*(IUPAC[b] for b in codon))}
  result.append(next(iter(values)) if len(values)==1 else 'X')
 return ''.join(result)
def split_args(s):
 parts=[];start=depth=0
 for i,c in enumerate(s):
  depth+=(c=='(')-(c==')')
  if c==',' and depth==0:parts.append(s[start:i]);start=i+1
 return parts+[s[start:]]
def extract(location,nt):
 """Deliberately reject remote, order(), fuzzy ranges beyond endpoint < >."""
 s=location.replace(' ','')
 if s.startswith('complement(') and s.endswith(')'):
  seq=extract(s[11:-1],nt);return seq.translate(str.maketrans('ACGTRYMKBDHVN','TGCAYRKMVHDBN'))[::-1]
 if s.startswith('join(') and s.endswith(')'):return ''.join(extract(x,nt) for x in split_args(s[5:-1]))
 m=re.fullmatch(r'[<>]?(\d+)\.\.[<>]?(\d+)',s)
 if m:
  a,b=map(int,m.groups());U.need(1<=a<=b<=len(nt),'CDS interval outside sequence');return nt[a-1:b]
 if re.fullmatch(r'\d+',s):
  a=int(s);U.need(1<=a<=len(nt),'CDS point outside sequence');return nt[a-1:a]
 raise ValueError('Unsupported CDS location: '+s)
def records(path):return {x.findtext('GBSeq_accession-version'):x for x in ET.parse(path).findall('.//GBSeq')}
def qualifiers(f):return {q.findtext('GBQualifier_name'):q.findtext('GBQualifier_value','') for q in f.findall('./GBFeature_quals/GBQualifier')}
def fetch_worker(ids,path):
 url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode(dict(db='nuccore',id=ids,rettype='gb',retmode='xml',tool='pepc_cds_rebuild'))
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'PEPC-CDS-rebuild/1.0'}),timeout=25) as f:data=f.read(16*1024*1024+1)
 U.need(len(data)<=16*1024*1024,'Response exceeds size bound');Path(path).write_bytes(data);records(path)
def raw_fasta(p):
 d={};key=None
 for l in p.read_text().splitlines():
  if l.startswith('>'):key=l[1:];U.need(key not in d,'Duplicate raw header');d[key]=''
  elif l.strip():U.need(key is not None,'Invalid FASTA');d[key]+=''.join(l.split()).upper()
 return d
def prepare(root,out):
 src=root/TREE;audit=root/AUDIT
 names=['tip_metadata.tsv','representative_proteins.faa','representative_to_locus_mapping.tsv','historical_reference_metadata.tsv']
 prov=U.checked(src,names);prov.update(U.checked(audit,['public_records/NCBI_GenBank.xml']))
 old=U.fasta(src/'representative_proteins.faa');oldmeta=U.readtsv(src/'tip_metadata.tsv');hist=U.readtsv(src/'historical_reference_metadata.tsv')
 proteins={};meta={};oldids={}
 for r in oldmeta:
  k=r['accession'] if r['role']=='REFERENCE' else r['tip_id'];oldids[r['tip_id']]=k
  rr=dict(r,tip_id=k);rr['represented_loci']=k if r['role']=='REFERENCE' else r['represented_loci']
  proteins[k]=old[r['tip_id']];meta[k]=rr
 U.need(len(proteins)==80 and sum(r['role']=='SUPPORTED' for r in meta.values())==49,'Unexpected predecessor membership')
 raw=[]
 for name in RAW:
  p=root/name;U.need(p.is_file(),'Missing raw reference input: '+name);prov[str(p)]=U.sha(p)
  raw.extend((name,k,v) for k,v in raw_fasta(p).items())
 (out/'public_records').mkdir(exist_ok=True)
 cache=out/'public_records/previous_E1.xml';shutil.copyfile(audit/'public_records/NCBI_GenBank.xml',cache);gb=records(cache)
 ids={r['sequence_id'] for r in hist if re.fullmatch(r'KJ\d+\.\d+',r['sequence_id'])};U.need(len(ids)==18,'Expected 18 historical KJ records')
 extra=out/'public_records/additional_records.xml'
 if os.environ.get('PEPC_GENBANK_XML'):
  p=Path(os.environ['PEPC_GENBANK_XML']);p=p if p.is_absolute() else root/p
  U.need(p.resolve()!=extra.resolve(),'Cache override must be external to output');shutil.copyfile(p,extra);prov[str(p.resolve())]=U.sha(p)
 if extra.exists():gb.update(records(extra))
 missing=sorted(ids-set(gb))
 if missing:
  print('FETCHING_NCBI_ACCESSIONS:',','.join(missing),flush=True)
  subprocess.run([sys.executable,str(Path(__file__).resolve()),'--fetch',','.join(missing),str(extra)],check=True,timeout=90)
  gb.update(records(extra))
 U.need(ids<=set(gb),'Missing exact accession.version records; inspect public_records')
 audits=[];corrected_cds={}
 for h in hist:
  k=h['sequence_id'];previous=proteins.get(k,'');species=h['species'].replace('_',' ')
  if k not in meta:
   meta[k]=dict(tip_id=k,display_label=species+' | '+k,role='REFERENCE',species_code='REFERENCE',species=species,accession=k,original_id='WANG_S7__'+k,represented_loci=k,historical_subtype=h['explicit_subtype'],metadata_status=h['metadata_status'],metadata_source='historical_reference_metadata.tsv',historical_species=species,note='Previously excluded; reassessed after CDS correction')
  a=dict(accession=k,old_aa=len(previous),new_aa='',status='',reason='',CDS_location='',codon_start='',gene='',public_species='',public_protein_id='',exact_old_match='',naive_genomic_translation_matches_old='',subtype_anchor_eligible=0)
  try:
   if k in ids:
    g=gb[k];nt=g.findtext('GBSeq_sequence','').upper();fs=[f for f in g.findall('./GBSeq_feature-table/GBFeature') if f.findtext('GBFeature_key')=='CDS' and qualifiers(f).get('translation')]
    U.need(len(fs)==1,'Expected exactly one translated CDS');f=fs[0];q=qualifiers(f);loc=f.findtext('GBFeature_location','');off=int(q.get('codon_start','1'))-1
    U.need(off in (0,1,2) and q.get('transl_table','1')=='1','Unsupported codon start/genetic code');U.need('transl_except' not in q,'Translation exception requires review')
    cds=extract(loc,nt)[off:];new=''.join(q['translation'].split()).upper().rstrip('*');computed=translate(cds).rstrip('*')
    if len(cds)%3 and len(new)==len(computed)+1 and ('<' in loc or '>' in loc):
     computed=translate(cds+'N'*(3-len(cds)%3)).rstrip('*');a['reason']='Partial terminal codon checked with ambiguous N completion against annotated translation'
    U.need(computed==new,'Exon-aware CDS translation disagrees with annotation');U.need(new and set(new)<=U.AA|{'X'},'Invalid annotated protein')
    corrected_cds[k]=cds;proteins[k]=new
    species=g.findtext('GBSeq_organism','');meta[k]['species']=species;meta[k]['display_label']=species+' | '+k+' | '+q.get('gene','PEPC')
    # Verify subtype identity from the gene annotation; retain historical conflict as a caveat.
    gene=q.get('gene','');annotated='PEPC'+re.sub(r'^ppc[-_]?','',gene,flags=re.I) if re.match(r'^ppc[-_]?',gene,re.I) else ''
    eligible=h['metadata_status']=='S7_EXPLICIT_SUBTYPE' and annotated.lower()==h['explicit_subtype'].lower()
    meta[k]['metadata_status']='CDS_AND_SUBTYPE_VERIFIED' if eligible else 'CDS_VERIFIED_SUBTYPE_REVIEW'
    meta[k]['metadata_source']='NCBI accession.version CDS; cached XML'
    meta[k]['note']+='; original metadata status: '+h['metadata_status']
    a.update(status='VERIFIED_ANNOTATED_CDS',new_aa=len(new),CDS_location=loc,codon_start=off+1,gene=gene,public_species=species,public_protein_id=q.get('protein_id',''),exact_old_match=int(previous==new) if previous else 'NOT_PREVIOUSLY_RETAINED',naive_genomic_translation_matches_old=int(''.join(CODE.get(nt[i:i+3],'X') for i in range(off,len(nt)-2,3)).rstrip('*').replace('*','X')==previous) if previous else '',subtype_anchor_eligible=int(eligible))
   else:
    matches=[(name,head,s.replace('-','').replace('.','').rstrip('*')) for name,head,s in raw if re.search(r'(?<![A-Za-z0-9])'+re.escape(k)+r'(?![A-Za-z0-9])',head)]
    good={s for _,_,s in matches if set(s)<=U.AA|{'X'} and set(s)-set('ACGTRYSWKMBDHVN')}
    U.need(len(good)==1,'No unique original protein record; nucleotide-only/unresolved source requires CDS annotation')
    new=next(iter(good));proteins[k]=new
    a.update(status='ORIGINAL_PROTEIN_RECORD',new_aa=len(new),exact_old_match=int(previous==new) if previous else 'NOT_PREVIOUSLY_RETAINED',reason='Original amino-acid FASTA; no nucleotide translation performed')
   U.need(len(proteins[k])>=30,'Reference protein too short')
  except (RuntimeError,ValueError) as e:
   proteins.pop(k,None);a.update(status='EXCLUDED_UNVERIFIED_REFERENCE',reason=str(e))
  audits.append(a)
 U.table(out/'01_reference_correction_audit.tsv',audits);U.writefa(out/'02_verified_reference_CDS.fna',corrected_cds)
 U.need(all(a['status']=='VERIFIED_ANNOTATED_CDS' for a in audits if a['accession'] in ids),'One or more KJ references failed CDS verification; inspect audit before tree inference')
 U.writefa(out/'03_all_input_proteins.faa',proteins)
 mafft=shutil.which('mafft');U.need(mafft,'MAFFT unavailable in activated environment')
 cpu=os.environ.get('SLURM_CPUS_PER_TASK','4')
 with (out/'04_all_alignment.faa').open('w') as f,(out/'MAFFT.log').open('w') as err:subprocess.run([mafft,'--auto','--thread',cpu,str(out/'03_all_input_proteins.faa')],stdout=f,stderr=err,check=True)
 aln=U.fasta(out/'04_all_alignment.faa');U.need(set(aln)==set(proteins) and len({len(s) for s in aln.values()})==1,'Invalid alignment')
 U.need(all(U.ungap(aln[k])==s for k,s in proteins.items()),'MAFFT changed protein residues')
 anchor='Q84KR7';U.need(len(proteins[anchor])==970,'Unexpected maize coordinate reference length');coords={};n=0
 for i,c in enumerate(aln[anchor]):
  if c not in '-.':n+=1;coords[n]=i
 decisions=[];wins=[];keep=[]
 for k,s in aln.items():
  covered=sum(s[coords[p]] in U.AA for p in U.POS);retained=covered>=2 or k in {'SS__g29856','SS__g7519'}
  decisions.append(dict(tip_id=k,role=meta[k]['role'],covered_positions=covered,retained=int(retained),reason='USER_SS_EXCEPTION' if covered<2 and retained else 'AT_LEAST_TWO_POSITIONS' if retained else 'INSUFFICIENT_FOCAL_COVERAGE'))
  if retained:keep.append(k)
  for p in U.POS:
   cols=[coords[j] for j in range(p-2,p+3)];w=''.join(s[c] for c in cols)
   wins.append(dict(tip_id=k,maize_position=p,window=w,center_residue=s[coords[p]],center_covered=int(s[coords[p]] in U.AA),full_window_covered=int(all(a in U.AA for a in w)),source_alignment_columns_1based=','.join(str(c+1) for c in cols)))
 U.table(out/'05_coverage_decisions.tsv',decisions);U.table(out/'06_all_five_residue_windows.tsv',wins)
 U.need(anchor in keep,'Maize anchor lost');candidates=[k for k in keep if meta[k]['role']=='SUPPORTED'];refs=[k for k in keep if meta[k]['role']=='REFERENCE']
 U.need(len(candidates)>=4 and len(refs)>=2,'Insufficient retained taxa')
 U.writefa(out/'source_alignment.faa',{k:aln[k] for k in keep});U.writefa(out/'representative_proteins.faa',{k:proteins[k] for k in keep});U.table(out/'tip_metadata.tsv',[meta[k] for k in keep]);U.table(out/'maize_windows.tsv',[r for r in wins if r['tip_id'] in keep])
 reps=U.readtsv(src/'representative_to_locus_mapping.tsv')
 for r in reps:
  r['tree_representative']=oldids[r['tree_representative']];r['represented_tip']=oldids.get(r['represented_tip'],r['represented_tip'])
 U.table(out/'representative_to_locus_mapping.tsv',reps)
 masks=[];main=[];sens=[]
 for i in range(len(aln[anchor])):
  nc=sum(aln[k][i] in U.AA for k in candidates);nr=sum(aln[k][i] in U.AA for k in refs);m=nc>=math.ceil(len(candidates)/2) and nr>=1;ss=nc+nr>0
  if m:main.append(i)
  if ss:sens.append(i)
  masks.append(dict(source_column_1based=i+1,supported_resolved=nc,reference_resolved=nr,main_keep=int(m),sensitivity_keep=int(ss)))
 U.need(len(main)>=300,'Main alignment too short; review required')
 for mode,cols in [('main',main),('sensitivity',sens)]:
  d={k:''.join(aln[k][i] for i in cols) for k in keep};U.need(all(any(a in U.AA for a in s) for s in d.values()),'Empty masked sequence');U.writefa(out/(mode+'_alignment.faa'),d)
 U.table(out/'column_mask.tsv',masks)
 U.table(out/'INPUT_PROVENANCE.tsv',[dict(project_relative_path=str(Path(p).relative_to(root)),sha256=h) for p,h in prov.items()])
 return dict(input_tips=len(keep),supported_representatives=len(candidates),supported_loci=51,references=len(refs),source_columns=len(aln[anchor]),main_columns=len(main),sensitivity_columns=len(sens),coordinate_reference=anchor,positions=list(U.POS),historical_KJ_references_verified=len(ids),excluded_unverified_references=[r['accession'] for r in audits if r['status']=='EXCLUDED_UNVERIFIED_REFERENCE'])
def classify(out):
 meta={r['tip_id']:r for r in U.readtsv(out/'tip_metadata.tsv')};audit=U.readtsv(out/'01_reference_correction_audit.tsv')
 anchors={r['accession']:meta[r['accession']]['historical_subtype'] for r in audit if r['subtype_anchor_eligible']=='1' and r['accession'] in meta}
 U.need(anchors,'No retained verified subtype anchors');per={};groups=[];evrows=[]
 for mode in ['main','sensitivity']:
  t=U.parse_newick((out/'trees'/f'{mode}.treefile').read_text());ev,gg,_=evidence(t,meta,anchors);per[mode]=ev;groups.extend(dict(tree=mode,**r) for r in gg);evrows.extend(dict(tree=mode,tip_id=k,**v) for k,v in ev.items())
 U.table(out/'07_reference_group_coherence.tsv',groups);U.table(out/'08_per_tree_subtype_evidence.tsv',evrows)
 final={};locus=[]
 for r in U.readtsv(out/'representative_to_locus_mapping.tsv'):
  k=r['tree_representative'];represented=r['represented_tip']
  if not represented.startswith(tuple(s+'__' for s in U.SPECIES)):continue
  if k not in final:
   a=per['main'].get(k,{});b=per['sensitivity'].get(k,{})
   same=bool(a.get('subtype') and a.get('subtype')==b.get('subtype'));strong=same and a.get('status')==b.get('status')=='SUPPORTED_PLACEMENT'
   final[k]=dict(tip_id=k,subtype=a['subtype'] if same else '',classification='SUPPORTED_SUBTYPE' if strong else 'PROVISIONAL_SUBTYPE' if same else 'UNRESOLVED',reason='BOTH_TREES_MEET_REFERENCE_AND_SUPPORT_CRITERIA' if strong else 'CONCORDANT_ASSOCIATION_WITH_INCOMPLETE_SUPPORT' if same else 'NOT_IN_TREE_AFTER_COVERAGE_FILTER' if k not in meta else 'INSUFFICIENT_OR_DISCORDANT_PLACEMENT',functional_C4_CAM_assignment='NOT_ESTABLISHED')
  sp,gene=represented.split('__',1);locus.append(dict(species=sp,gene_id=gene,tree_representative=k,**{kk:vv for kk,vv in final[k].items() if kk!='tip_id'}))
 U.need(len(locus)==51 and len({(r['species'],r['gene_id']) for r in locus})==51,'Locus accounting changed')
 U.table(out/'09_representative_subtypes.tsv',list(final.values()));U.table(out/'10_locus_subtypes.tsv',locus)
 counts=collections.Counter((r['species'],r['classification'],r['subtype']) for r in locus)
 U.table(out/'11_species_subtype_summary.tsv',[dict(species=k[0],classification=k[1],subtype=k[2],loci=v) for k,v in sorted(counts.items())])
 for mode in ['main','sensitivity']:
  folder=out/'itol'/mode;order=U.leaves(U.parse_newick((folder/'01_PEPC.treefile').read_text()))
  U.dataset_text(folder/'08_subtype_evidence.txt','Reference-guided subtype; not functional proof',order,{k:(r['subtype']+' | '+r['classification']).strip(' |') for k,r in final.items() if k in meta})
  for p in folder.glob('*.txt'):
   lines=p.read_text().splitlines()
   if 'DATA' in lines:
    ids=[x.split('\t')[0] for x in lines[lines.index('DATA')+1:] if x and not x.startswith('#')];U.need(len(ids)==len(set(ids)) and set(ids)<=set(meta),'Invalid iTOL IDs')
 return dict(eligible_retained_subtype_anchors=len(anchors),locus_categories=dict(collections.Counter(r['classification'] for r in locus)))
def main():
 root=Path(os.environ['PROJECT_ROOT']).resolve();pkg=Path(__file__).resolve().parent
 iq=next((shutil.which(x) for x in ['iqtree3','iqtree2','iqtree'] if shutil.which(x)),None);U.need(iq and shutil.which('mafft'),'Activate an environment with Python3, MAFFT and IQ-TREE2/3 before submission')
 base=root/BASE;resume=os.environ.get('PEPC_CDS_RUN');base.mkdir(parents=True,exist_ok=True)
 out=(root/resume).resolve() if resume else base/('cds_rebuild_'+os.environ.get('SLURM_JOB_ID','manual')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
 out.relative_to(base.resolve())
 if not resume:out.mkdir(exist_ok=False)
 U.need(out.is_dir() and not (out/'COMPLETE').exists(),'Output missing or already complete')
 print('OUTPUT_DIRECTORY:',out,flush=True)
 print('Resume an interrupted run by exporting PEPC_CDS_RUN='+str(out.relative_to(root)),flush=True)
 package_hash={p.name:U.sha(p) for p in pkg.iterdir() if p.is_file() and p.suffix in {'.py','.slurm','.txt'}}
 init=out/'PACKAGE_STATE.json'
 if init.exists():U.need(json.loads(init.read_text())==package_hash,'Resume requires unchanged package')
 else:
  init.write_text(json.dumps(package_hash,indent=2)+'\n')
  for name in package_hash:shutil.copyfile(pkg/name,out/name)
 statefile=out/'RUN_STATE.json'
 if statefile.exists():
  state=json.loads(statefile.read_text());prep=state['preparation']
  for n,h in state['hashes'].items():U.need(U.sha(out/n)==h,'Prepared input changed: '+n)
 else:
  prep=prepare(root,out);state=dict(preparation=prep,hashes={str(p.relative_to(out)):U.sha(p) for p in out.rglob('*') if p.is_file()});statefile.write_text(json.dumps(state,indent=2)+'\n')
 version=subprocess.run([iq,'--version'],capture_output=True,text=True,check=True);v=version.stdout+version.stderr;U.need(re.search(r'version\s+[23]\.',v,re.I),'IQ-TREE2/3 required')
 vp=out/'IQTREE_VERSION.txt'
 if vp.exists():U.need(vp.read_text()==v,'IQ-TREE changed on resume')
 else:vp.write_text(v)
 (out/'trees').mkdir(exist_ok=True)
 for mode,seed in [('main',519780),('sensitivity',519890)]:
  flag=out/'trees'/f'{mode}.SUCCESS';tree=out/'trees'/f'{mode}.treefile'
  if flag.exists():U.need(flag.read_text().strip()==U.sha(tree),'Completed tree changed');continue
  cmd=[iq,'-s',mode+'_alignment.faa','-st','AA','-m','MFP','-B','1000','-bnni','-alrt','1000','-T',os.environ.get('SLURM_CPUS_PER_TASK','4'),'-seed',str(seed),'-keep-ident','-pre','trees/'+mode]
  (out/'trees'/f'{mode}.command.json').write_text(json.dumps([Path(iq).name]+cmd[1:],indent=2)+'\n')
  print('RUNNING_TREE:',mode,flush=True)
  with (out/'trees'/f'{mode}.stdout.log').open('a') as log:subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
  flag.write_text(U.sha(tree)+'\n')
 qc=U.export(out,prep);qc.update(classify(out));qc.update(status='CORRECTED_CDS_TREES_AND_SUBTYPE_REVIEW_COMPLETE',subtype_assignments_to_supported_genes='RECOMPUTED_FROM_CORRECTED_REFERENCES',functional_validation=False,previous_reference_dependent_conclusions='SUPERSEDED_FOR_REASSESSMENT',reference_species_conflict_retained=True)
 (out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
 for r in U.readtsv(out/'INPUT_PROVENANCE.tsv'):U.need(U.sha(root/r['project_relative_path'])==r['sha256'],'Source changed during job')
 for n,h in state['hashes'].items():U.need(U.sha(out/n)==h,'Prepared input changed during inference')
 (out/'COMPLETE').write_text('New corrected reference release; inspect QC and subtype uncertainty before reporting.\n')
 with (out/'SHA256SUMS.txt').open('w') as f:
  for p in sorted(out.rglob('*')):
   if p.is_file() and p.name!='SHA256SUMS.txt':f.write(U.sha(p)+'  '+str(p.relative_to(out))+'\n')
 bundle=out.with_name(out.name+'_bundle.tar.gz')
 with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname=out.name)
 print(json.dumps(qc,indent=2));print('UPLOAD_BUNDLE:',bundle,flush=True)
if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='--fetch':fetch_worker(sys.argv[2],sys.argv[3])
 else:main()
