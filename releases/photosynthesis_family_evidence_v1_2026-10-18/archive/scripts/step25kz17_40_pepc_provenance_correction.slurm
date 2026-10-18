#!/usr/bin/env bash
#SBATCH --job-name=pepc_prov_fix
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=pepc_prov_fix_%j.out
#SBATCH --error=pepc_prov_fix_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
python3 - <<'PY'
import csv, gzip, hashlib, json, os, shutil, tarfile
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime, timezone
root=Path(os.environ['PROJECT_ROOT']).resolve(strict=True)
base=root/'08_focused_families/validation'
def require(ok,msg):
    if not ok: raise SystemExit(msg)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()
def verify(p):
    manifest=p/'SHA256SUMS.txt'
    require(manifest.is_file(),'Missing release checksums')
    for line in manifest.read_text().splitlines():
        h,n=line.split(maxsplit=1); n=n.lstrip('*'); f=(p/n).resolve()
        require(f.is_relative_to(p) and f.is_file() and sha(f)==h,'Release checksum mismatch: '+n)
def release(step,env):
    pp=[root/os.environ[env]] if os.environ.get(env) else [p for p in (base/step).glob('release_*') if p.is_dir() and (p/'QC.json').is_file()]
    require(len(pp)==1,'Expected one release; set '+env+' explicitly')
    p=pp[0].resolve(strict=True); verify(p); return p
def rows(p):
    with p.open() as f: return list(csv.DictReader(f,delimiter='\t'))
def label(p):
    return str(p.relative_to(root)) if p.is_relative_to(root) else '<EXTERNAL>/'+p.name
def clean(s): return str(s).replace(str(root),'<PROJECT_ROOT>').replace(str(Path.home()),'<HOME>')
src=release('step25kz17.39','PEPC_AUDIT_RELEASE')
cor=release('step25kz17.21','CORRECTED_INPUT_RELEASE')
candidates=rows(src/'03_PEPC_provenance_audit.tsv')
ids={r['protein_id']:r for r in candidates}
require(len(candidates)==len(ids)==55,'Expected 55 distinct PEPC candidates')
ledger=rows(src/'01_integrated_locus_evidence.tsv')
require(len(ledger)==820,'Unexpected integrated ledger size')
proteome=cor/'full_outputs/corrected_protein_candidates.faa'
expected=[r['sha256'] for r in rows(cor/'08_full_output_checksums.tsv') if r['file']=='full_outputs/corrected_protein_candidates.faa']
require(len(expected)==1 and proteome.is_file() and sha(proteome)==expected[0],'Corrected protein provenance mismatch')
out=base/'step25kz17.40'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
out.mkdir(parents=True,exist_ok=False)
def table(name,rr,fields):
    with (out/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n'); w.writeheader();w.writerows(rr)
issues=[]; provenance=[]
def record(role,p):
    provenance.append(dict(role=role,path=label(p),sha256=sha(p),bytes=p.stat().st_size))
def locate(env,rel): return (root/os.environ.get(env,rel)).resolve()
registry=locate('PEPC_REGISTRY','08_focused_families/validation/step25kn/results/PEPC/candidate_validation.tsv')
classification=locate('PEPC_CLASSIFICATION','08_focused_families/validation/step25kz16.7/release/02_all_species_pepc_classification.tsv')
treefa=locate('PEPC_TREE_FASTA','08_focused_families/validation/step25kz16.7b/release/work/corrected_all_species_PEPC_primary_with_refs.faa')
cds=locate('SS_RECOVERY_CDS','06_positive_selection/sequences/positive_selection_sequence_recovery_primary_v1/species_sequences/SS.candidate.cds.fa')
def normalize(v):
    if v in ids: return v,'EXACT'
    prefix,sep,key=v.partition('|')
    if sep and key in ids and prefix==ids[key]['species']: return key,'SPECIES_PREFIX_CONFIRMED'
    return None,None
reg=defaultdict(list); clas=defaultdict(list); extracted=[]
for role,p,idfield,target in [('registry',registry,'sequence_id',reg),('classification',classification,'tree_id',clas)]:
    if not p.is_file():
        issues.append(dict(source=role,issue='FILE_MISSING'));continue
    record(role,p)
    with p.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        require(idfield in (reader.fieldnames or []),'Unexpected '+role+' identifier column')
        for line,r in enumerate(reader,2):
            k,method=normalize(r.get(idfield,''))
            if k:
                target[k].append(r)
                extracted.append(dict(protein_id=k,source=role,line=line,match_method=method,original_id=r[idfield],record_json=clean(json.dumps(r,sort_keys=True))))
# Stream only designated FASTAs; retain exact requested IDs, never strip transcript suffixes.
def read_selected(p,wanted):
    result={}; key=None; parts=[]
    def commit():
        if key in wanted:
            require(key not in result,'Duplicate requested FASTA identifier: '+key)
            result[key]=''.join(parts).upper()
    op=gzip.open if p.suffix=='.gz' else open
    with op(p,'rt') as f:
        for line in f:
            if line.startswith('>'):
                commit();key=line[1:].split()[0];parts=[]
            elif key in wanted: parts.append(''.join(line.split()))
        commit()
    return result
def load(role,p,wanted,required=False):
    if not p.is_file():
        require(not required,'Missing '+role)
        issues.append(dict(source=role,issue='FILE_MISSING'));return {}
    record(role,p);return read_selected(p,wanted)
historical=load('historical_tree_input',treefa,set(ids))
current=load('corrected_proteome',proteome,set(ids),True)
# Local recovery identifiers are scoped to SS only; exact aliases are recorded.
ss_alias={}
for k,r in ids.items():
    if r['species']=='SS':
        for alias in (k,k.split('__',1)[1],'SS|'+k): ss_alias[alias]=k
recovery=load('SS_recovery_CDS',cds,set(ss_alias))
recovery_by_id=defaultdict(list)
for k,s in recovery.items(): recovery_by_id[ss_alias[k]].append((k,s))
def norm(s): return s[:-1] if s.endswith('*') else s
def seqhash(s): return hashlib.sha256(s.encode()).hexdigest() if s else ''
# Standard nuclear code, frame 1 only; no alternative ORF selection or model repair.
bases='TCAG'; aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
code=dict(zip((a+b+c for a in bases for b in bases for c in bases),aas))
def translate(s): return ''.join(code.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
comparison=[]; audit=[]; fa=[]; recovery_audit=[]
for k,r in sorted(ids.items()):
    old=historical.get(k,''); new=current.get(k,''); rg=reg[k]; cl=clas[k]
    relation='EXACT_NORMALIZED_MATCH' if old and new and norm(old)==norm(new) else 'SEQUENCE_DIFFERENCE' if old and new else 'MISSING_CORRECTED_PROTEIN' if not new else 'MISSING_HISTORICAL_SEQUENCE'
    dna_like=bool(old) and set(norm(old))<=set('ACGTN')
    exact_cds=[]
    for alias,s in recovery_by_id[k]:
        tr=translate(s); matches=bool(old) and norm(old)==norm(s)
        if matches: exact_cds.append(alias)
        recovery_audit.append(dict(protein_id=k,recovery_id=alias,cds_length=len(s),length_mod3=len(s)%3,cds_sha256=seqhash(s),historical_equals_CDS=matches,translation_length=len(norm(tr)),translation_internal_stops=norm(tr).count('*'),translation_equals_corrected_protein=bool(new) and norm(tr)==norm(new),translation_equals_historical=bool(old) and norm(tr)==norm(old)))
        fa.append(('RECOVERY_CDS|'+alias,s))
    if exact_cds: relation='HISTORICAL_SEQUENCE_EQUALS_RECOVERY_CDS'
    comparison.append(dict(protein_id=k,species=r['species'],historical_length=len(old) if old else '',corrected_protein_length=len(new) if new else '',historical_sha256=seqhash(old),corrected_sha256=seqhash(new),historical_ACGTN_only=dna_like,historical_internal_stops=norm(old).count('*') if old else '',corrected_internal_stops=norm(new).count('*') if new else '',recovery_CDS_exact_ids=';'.join(exact_cds),sequence_relation=relation))
    flags=[]
    if len(rg)!=1: flags.append('REGISTRY_MISSING' if not rg else 'MULTIPLE_REGISTRY_RECORDS')
    if rg and any(x.get('status')!='PASS' for x in rg): flags.append('REGISTRY_NON_PASS')
    if len(cl)!=1: flags.append('CLASSIFICATION_MISSING' if not cl else 'MULTIPLE_CLASSIFICATION_RECORDS')
    if relation!='EXACT_NORMALIZED_MATCH': flags.append(relation)
    if dna_like: flags.append('HISTORICAL_NUCLEOTIDE_ALPHABET_REQUIRES_REVIEW')
    if new and norm(new).count('*'):flags.append('CORRECTED_PROTEIN_INTERNAL_STOP')
    audit.append(dict(protein_id=k,species=r['species'],gene_id=r['gene_id'],registry_record_count=len(rg),registry_statuses=';'.join(sorted({x.get('status','') for x in rg})),registry_evidence=';'.join(sorted({x.get('evidence','') for x in rg})),classification_record_count=len(cl),historical_relaxed_similarity=';'.join(sorted({x.get('relaxed_similarity','') for x in cl})),historical_registry_flag=';'.join(sorted({x.get('validated_registry','') for x in cl})),provenance_category='PROVENANCE_EXCEPTION' if flags else 'TRACEABLE_SEQUENCE_AND_REGISTRY',review_flags=';'.join(flags),membership_decision='PENDING_PEPC_SPECIFIC_VALIDATION',final_count='NA'))
    if old:fa.append(('HISTORICAL|'+k,old))
    if new:fa.append(('CORRECTED_PROTEIN|'+k,new))
table('01_PEPC_corrected_provenance.tsv',audit,list(audit[0]))
table('02_exact_source_records.tsv',extracted,['protein_id','source','line','match_method','original_id','record_json'])
table('03_sequence_comparison.tsv',comparison,list(comparison[0]))
table('04_SS_recovery_CDS_comparison.tsv',recovery_audit,['protein_id','recovery_id','cds_length','length_mod3','cds_sha256','historical_equals_CDS','translation_length','translation_internal_stops','translation_equals_corrected_protein','translation_equals_historical'])
(out/'05_collected_sequences.fasta').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in fa))
# Preserve the integrated ledger byte-for-byte: this step corrects provenance only.
shutil.copyfile(src/'01_integrated_locus_evidence.tsv',out/'06_integrated_locus_evidence_UNCHANGED.tsv')
shutil.copyfile(src/'10_final_Wang_counts_PENDING.tsv',out/'07_final_Wang_counts_PENDING.tsv')
table('08_collection_issues.tsv',issues,['source','issue'])
record('prior_audit_manifest',src/'SHA256SUMS.txt');record('corrected_input_manifest',cor/'SHA256SUMS.txt')
table('09_input_provenance.tsv',provenance,['role','path','sha256','bytes'])
qc=dict(status='BOUNDED_PEPC_PROVENANCE_CORRECTION_COMPLETE_REVIEW_REQUIRED',PEPC_loci=len(audit),registry_matched_loci=len([k for k in ids if reg[k]]),classification_matched_loci=len([k for k in ids if clas[k]]),historical_sequences=len(historical),corrected_proteins=len(current),sequence_relations=dict(Counter(x['sequence_relation'] for x in comparison)),provenance_categories=dict(Counter(x['provenance_category'] for x in audit)),collection_issues=len(issues),integrated_ledger_unchanged=sha(out/'06_integrated_locus_evidence_UNCHANGED.tsv')==sha(src/'01_integrated_locus_evidence.tsv'),analysis_rerun=False,final_counts=False)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n')
(out/'README.txt').write_text('Bounded PEPC provenance correction. No alignment, search, tree inference, model repair, or formal membership promotion. Registry PASS alone does not establish PEPC specificity. Historical registry inclusion ignored PASS/FAIL; corrected tree metadata labels do not independently validate family membership. Species prefixes are removed only with an exact species-confirmed match; transcript suffixes are never stripped. Corrected-proteome absence is retained as unresolved, not biological absence. ACGTN-only composition is a flag; exact equality to the documented CDS is reported separately. CDS translations use standard nuclear code, frame 1, with one terminal stop ignored for comparisons. Original sequence hashes and lengths retained. Missing inputs and multiple rows remain explicit. The mixed collected FASTA is for inspection only, never an alignment input. Previously integrated 820 records and pending final counts remain unchanged.\n')
# Recheck sources after reads to detect an input changing during this audit.
for p in (registry,classification,treefa,cds,proteome):
    if p.is_file():
        expected_hashes={x['sha256'] for x in provenance if x['path']==label(p)}
        require(expected_hashes=={sha(p)},'Input changed during audit: '+label(p))
files=sorted(p for p in out.rglob('*') if p.is_file())
(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files))
verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:t.add(out,arcname='PEPC_provenance_correction')
print(json.dumps(qc,indent=2));print('RELEASE='+label(out));print('BUNDLE='+label(bundle))
PY
