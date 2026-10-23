#!/usr/bin/env python3
"""Conservative unrooted reference-guided PEPC subtype review."""
import csv,hashlib,json,os,re,shutil,tarfile,collections
from pathlib import Path
from datetime import datetime,timezone
from tree_utils import parse_newick,leaves,need
SOURCE='08_focused_families/downstream_supported/pepc_reference_review/tree_itol_20835576_20261022T211258984458Z'
SHA_MIN=80.;BOOT_MIN=95.
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,rows,fields=None):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def edges(tree):
 universe=frozenset(leaves(tree));result={}
 def walk(n):
  side=frozenset(leaves(n));other=universe-side
  if len(side)>1 and len(other)>1:
   canonical=min(tuple(sorted(side)),tuple(sorted(other)))
   try:alrt,boot=map(float,n['name'].split('/'))
   except (ValueError,KeyError):alrt=boot=None
   result[canonical]=(side,other,alrt,boot)
  for c in n['children']:walk(c)
 walk(tree);return result

def evidence(tree,meta,anchors):
 universe=set(leaves(tree));es=edges(tree);groups=collections.defaultdict(set)
 for k,v in anchors.items():groups[v].add(k)
 blockers={k:r['historical_subtype'] for k,r in meta.items() if r['role']=='REFERENCE' and r['historical_subtype'].startswith('PEPC')}
 candidates=[k for k,r in meta.items() if r['role']=='SUPPORTED'];possible=collections.defaultdict(lambda:collections.defaultdict(list));group_rows=[]
 for subtype,refids in sorted(groups.items()):
  pure=[]
  for edge_id,(a,b,alrt,boot) in enumerate(es.values(),1):
   for side in (a,b):
    hits=set(side)&set(anchors)
    if hits and {anchors[k] for k in hits}=={subtype} and not any(blockers[k]!=subtype for k in set(side)&set(blockers)):
     row=dict(subtype=subtype,side=side,anchors=hits,alrt=alrt,boot=boot,all_refs=hits==refids,distinct_species=len({meta[k]['species'] for k in hits}),edge_id=edge_id)
     pure.append(row)
     for k in set(side)&set(candidates):possible[k][subtype].append(row)
  full=[r for r in pure if r['all_refs']]
  group_rows.append(dict(subtype=subtype,eligible_reference_count=len(refids),reference_species_count=len({meta[k]['species'] for k in refids}),exclusive_full_reference_group_present=int(bool(full)),reference_ids=';'.join(sorted(refids))))
 result={}
 for k in candidates:
  choices=possible[k]
  if len(choices)!=1:
   result[k]=dict(subtype='',status='AMBIGUOUS_REFERENCE_GROUPS' if choices else 'NO_EXCLUSIVE_REFERENCE_GROUP',candidate_subtypes=';'.join(sorted(choices)),alrt='',boot='',anchor_count=0,species_count=0,all_refs=0,side='',anchor_ids='');continue
  subtype=next(iter(choices));rr=choices[subtype]
  # Prefer a full-reference group. Within this set choose the tightest enclosing side.
  full=[r for r in rr if r['all_refs']];pool=full or rr
  best=min(pool,key=lambda r:(len(r['side']),tuple(sorted(r['side']))))
  supported=bool(best['all_refs'] and len(best['anchors'])>=2 and best['distinct_species']>=2 and best['alrt'] is not None and best['alrt']>=SHA_MIN and best['boot']>=BOOT_MIN)
  result[k]=dict(subtype=subtype,status='SUPPORTED_PLACEMENT' if supported else 'PROVISIONAL_PLACEMENT',candidate_subtypes=subtype,alrt=best['alrt'] if best['alrt'] is not None else '',boot=best['boot'] if best['boot'] is not None else '',anchor_count=len(best['anchors']),species_count=best['distinct_species'],all_refs=int(best['all_refs']),side=';'.join(sorted(best['side'])),anchor_ids=';'.join(sorted(best['anchors'])))
 return result,group_rows,es

