#!/usr/bin/env python3
"""Rebuild the frozen Wang-layout evidence tables. Python 3 standard library only."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

TIERS = ('SUPPORTED','MODEL_CAVEAT','LIMITED_PANEL_CAVEAT','PARTIAL','UNRESOLVED','OUTSIDE_FAMILY_SET','EXCLUDED')
ID_TIERS = TIERS[:3]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return list(reader), reader.fieldnames

def write(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def require(condition, message):
    if not condition:
        raise ValueError(message)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--out', type=Path, required=True, help='New directory; refuses to overwrite existing output')
    args = parser.parse_args()
    release = args.release.resolve()
    cfg = json.loads((release/'provenance/freeze.json').read_text())
    ledger = release/'inputs'/cfg['ledger_file']
    matrix = release/'inputs/matrix_SUPPORTED.tsv'
    require(digest(ledger) == cfg['ledger_sha256'], 'Frozen ledger checksum mismatch')
    require(digest(matrix) == cfg['matrix_sha256'], 'Frozen matrix checksum mismatch')
    rows, fields = read(ledger)
    reference, matrix_fields = read(matrix)
    species, families = matrix_fields[1:], [r['wang_family'] for r in reference]
    blocked = set(cfg['blocked_families'])
    require(len(set(families)) == 23 and len(species) == 8, 'Unexpected reporting scope')
    require(len(rows) == cfg['expected_records'], 'Unexpected ledger record count')
    keys = [(r['wang_family'], r['species'], r['gene_id']) for r in rows]
    require(len(keys) == len(set(keys)), 'Duplicate family/species/locus key')
    for r in rows:
        require(r['wang_family'] in families and r['wang_family'] not in blocked, 'Invalid family')
        require(r['species'] in species and r['count_tier'] in TIERS, 'Invalid species or tier')
        for field in ('gene_id','protein_id'):
            require(bool(r[field]) and not any(c.isspace() for c in r[field]), 'Missing or whitespace-containing ID')
    totals = collections.Counter(r['count_tier'] for r in rows)
    require(dict(totals) == cfg['expected_tier_totals'], 'Tier totals differ from freeze')
    counts = collections.Counter((r['wang_family'],r['species'],r['count_tier']) for r in rows)
    for r in reference:
        for sp in species:
            expected = 'NA' if r['wang_family'] in blocked else str(counts[r['wang_family'],sp,'SUPPORTED'])
            require(r[sp] == expected, 'Count/ledger mismatch')
    rows.sort(key=lambda r:(families.index(r['wang_family']),species.index(r['species']),r['gene_id']))
    args.out.mkdir(parents=True, exist_ok=False)
    write(args.out/'gene_list_all_evidence.tsv',fields,rows)
    for tier in TIERS:
        write(args.out/f'gene_list_{tier.lower()}.tsv',fields,[r for r in rows if r['count_tier']==tier])
        table=[dict(Species=sp,**{fam: 'NA' if fam in blocked else counts[fam,sp,tier] for fam in families}) for sp in species]
        write(args.out/f'counts_{tier.lower()}.tsv',['Species']+families,table)
    evidence=[]
    manifest=[]
    for fam in families:
        for sp in species:
            evidence.append(dict(wang_family=fam,species=sp,**{t:'NA' if fam in blocked else counts[fam,sp,t] for t in TIERS},formal_Wang_count='NA'))
            for tier in ID_TIERS:
                rr=[r for r in rows if (r['wang_family'],r['species'],r['count_tier'])==(fam,sp,tier)]
                item=dict(family=fam,species=sp,tier=tier,locus_count='NA' if fam in blocked else len(rr),status='DEFINITION_UNRESOLVED' if fam in blocked else 'RECORDS_AVAILABLE' if rr else 'NO_RECORDS_IN_TIER_NOT_ABSENCE')
                for field in ('gene_id','protein_id'):
                    path=Path('id_lists')/tier/fam/(sp+'.'+field+'s.txt')
                    item[field+'_file']='NA' if fam in blocked else path.as_posix()
                    if fam not in blocked:
                        values=[r[field] for r in rr]
                        require(len(values)==len(set(values)), 'Repeated identifier within ID list')
                        target=args.out/path; target.parent.mkdir(parents=True,exist_ok=True)
                        target.write_text(''.join(v+'\n' for v in values),encoding='utf-8')
                manifest.append(item)
    write(args.out/'evidence_counts.tsv',list(evidence[0]),evidence)
    write(args.out/'id_list_manifest.tsv',list(manifest[0]),manifest)
    qc={'status':'PASS','release_id':cfg['release_id'],'family_species_cells':len(families)*len(species),'ledger_records':len(rows),'tier_totals':dict(sorted(totals.items())),'unique_family_species_locus_keys':len(set(keys)),'formal_counts':False,'blocked_families':sorted(blocked),'ledger_sha256':digest(ledger)}
    (args.out/'QC.json').write_text(json.dumps(qc,indent=2,sort_keys=True)+'\n')
    paths=sorted(p for p in args.out.rglob('*') if p.is_file())
    (args.out/'SHA256SUMS.txt').write_text(''.join(digest(p)+'  '+p.relative_to(args.out).as_posix()+'\n' for p in paths))
    print(json.dumps(qc,indent=2))

if __name__=='__main__':
    main()
