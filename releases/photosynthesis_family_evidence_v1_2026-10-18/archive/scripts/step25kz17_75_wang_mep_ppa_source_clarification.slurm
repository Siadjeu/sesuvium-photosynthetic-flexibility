#!/usr/bin/env bash
#SBATCH --job-name=wang_source_audit
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=wang_source_audit_%j.out
#SBATCH --error=wang_source_audit_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
command -v python3
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














import urllib.parse,xml.etree.ElementTree as ET,zipfile,io,posixpath
from html.parser import HTMLParser
src=release('step25kz17.74','WANG_REPORT_RELEASE')
out=base/'step25kz17.75'/('release_'+os.environ.get('SLURM_JOB_ID','local')+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
(out/'sources').mkdir();(out/'extracted').mkdir()
def table(n,rows,fields):
 with (out/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
# Small, predefined primary-source set. No genome analyses, no recursive web crawl.
requests=[('Wang2023_article','WANG_PRIMARY','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10663116/fullTextXML','xml'),('Wang2023_supplements','WANG_PRIMARY','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10663116/supplementaryFiles','zip'),('Maize2012_MEP_article','RELATED_PRIMARY_NOT_WANG_MEMBERSHIP','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3440195/fullTextXML','xml'),('Ferrari2020_PPA_article','RELATED_PRIMARY_NOT_WANG_MEMBERSHIP','https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.16265','html')]
retrieval=[];inventory=[];passages=[];identifiers=[];references=[];parse_status=[];processed=set();totalbytes=0;MAX_BYTES=60*1024**2;MAX_FILE=20*1024**2;start=time.monotonic()
patterns={'MEP':re.compile(r'\bMEP(?:[123])?\b|proton.{0,12}pyruvate|mesophyll.{0,20}envelope',re.I),'PPA':re.compile(r'\bPPA(?:[-_][A-Za-z0-9.]+)?\b|inorganic pyrophosphatase|soluble.{0,20}pyrophosphatase',re.I)}
idpat=re.compile(r'\b(?:AT[1-5CM]G\d{5}(?:\.\d+)?|GRMZM\dG\d+(?:_T\d+)?|Zm\d+[a-z]+\d+|evm\.TU\.[A-Za-z0-9_.-]+|(?:NP|XP)_\d+(?:\.\d+)?|(?:MEP|PPA)[-_][A-Za-z0-9.]+)\b',re.I)
def emit(source,origin,location,text,context='DIRECT_TEXT'):
 text=' '.join(text.split());families=[fam for fam,pat in patterns.items() if pat.search(text)]
 if not families:return
 if len(passages)>=1500:return
 # Keep bounded contextual excerpts, not arbitrary large article/chapter copies.
 spans=[]
 for fam in families:
  for m in patterns[fam].finditer(text):spans.append((max(0,m.start()-350),min(len(text),m.end()+700)))
 spans=sorted(spans);merged=[]
 for a,b in spans:
  if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
  else:merged.append((a,b))
 excerpt=' [...] '.join(text[a:b] for a,b in merged)[:5000]
 pid='P'+str(len(passages)+1).zfill(5);found=sorted(set(idpat.findall(excerpt)))
 passages.append(dict(passage_id=pid,source_file=source,source_class=origin,location=location,family_mentions=';'.join(families),context=context,excerpt=excerpt,identifier_tokens=';'.join(found),interpretation='SOURCE_LEAD_NOT_VALIDATED_COUNTED_MEMBERSHIP'))
 for ident in found:identifiers.append(dict(passage_id=pid,source_class=origin,family_mentions=';'.join(families),identifier=ident,status='CONTEXT_TOKEN_REQUIRES_EXPLICIT_FAMILY_AND_SEQUENCE_MAPPING'))
def fetch(tag,origin,url,ext):
 global totalbytes
 p=out/'sources'/(tag+'.'+ext)
 try:
  if time.monotonic()-start>600:raise ValueError('NETWORK_TIME_BUDGET_EXCEEDED')
  req=urllib.request.Request(url,headers={'User-Agent':'Wang-source-clarification/1.0','Accept':'*/*'})
  with urllib.request.urlopen(req,timeout=45) as response:
   data=response.read(MAX_FILE+1);final=response.url;mime=response.headers.get('Content-Type','')
  if len(data)>MAX_FILE or totalbytes+len(data)>MAX_BYTES:raise ValueError('DOWNLOAD_SIZE_LIMIT')
  if ext=='zip' and not zipfile.is_zipfile(io.BytesIO(data)):raise ValueError('RESPONSE_NOT_ZIP')
  if ext=='xml':
   doc=ET.fromstring(data)
   if doc.tag.split('}')[-1]!='article':raise ValueError('RESPONSE_NOT_ARTICLE_XML')
  if ext=='html' and any(x in data[:10000].lower() for x in [b'just a moment',b'access denied',b'captcha',b'cf-chl']):raise ValueError('ACCESS_CHALLENGE_NOT_ARTICLE')
  p.write_bytes(data);totalbytes+=len(data);retrieval.append(dict(tag=tag,source_class=origin,url=url,final_url=final,status='DOWNLOADED',bytes=len(data),sha256=sha(p),error=''));return p
 except (OSError,ValueError,ET.ParseError) as e:
  retrieval.append(dict(tag=tag,source_class=origin,url=url,final_url='',status='UNAVAILABLE_OR_BLOCKED',bytes=0,sha256='',error=str(e)));return None
class TextHTML(HTMLParser):
 def __init__(self):super().__init__();self.parts=[];self.skip=0
 def handle_starttag(self,tag,attrs):
  if tag in ['script','style']:self.skip+=1
  if tag in ['p','div','br','tr','li','h1','h2','h3']:self.parts.append('\n')
 def handle_endtag(self,tag):
  if tag in ['script','style']:self.skip=max(0,self.skip-1)
 def handle_data(self,data):
  if not self.skip:self.parts.append(data)
def analyze(p,origin,provenance,depth=0):
 global totalbytes
 digest=sha(p)
 if digest in processed:return
 processed.add(digest);rel=str(p.relative_to(out));suffix=p.suffix.lower();inventory.append(dict(file=rel,source_class=origin,provenance=provenance,bytes=p.stat().st_size,sha256=digest))
 try:
  if suffix=='.zip':
   if depth>=2:raise ValueError('NESTED_ARCHIVE_DEPTH_LIMIT')
   with zipfile.ZipFile(p) as z:
    if len(z.infolist())>300:raise ValueError('ARCHIVE_MEMBER_LIMIT')
    for j,info in enumerate(z.infolist()):
     if info.is_dir():continue
     ext=Path(info.filename).suffix.lower()
     if ext not in {'.pdf','.xlsx','.docx','.xml','.txt','.tsv','.csv','.zip','.html','.htm'}:continue
     if info.file_size>MAX_FILE or totalbytes+info.file_size>MAX_BYTES:parse_status.append(dict(file=rel+'::'+info.filename,status='SIZE_LIMIT',detail=''));continue
     data=z.read(info);totalbytes+=len(data);child=out/'extracted'/(digest[:10]+'_'+str(j)+ext);child.write_bytes(data);analyze(child,origin,rel+'::'+info.filename,depth+1)
  elif suffix=='.xlsx':
   ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
   with zipfile.ZipFile(p) as z:
    need(sum(x.file_size for x in z.infolist())<MAX_FILE*3,'Workbook expansion exceeds60MB')
    ss=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',ns)] if 'xl/sharedStrings.xml' in z.namelist() else []
    rels={x.attrib['Id']:x.attrib['Target'] for x in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    for sheet in ET.fromstring(z.read('xl/workbook.xml')).findall('m:sheets/m:sheet',ns):
     target=rels[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']];name=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
     for row in ET.fromstring(z.read(name)).findall('.//m:row',ns):
      cells=[]
      for c in row:
       val=c.find('m:v',ns);v=val.text or '' if val is not None else ''.join(c.find('m:is',ns).itertext()) if c.find('m:is',ns) is not None else ''
       if c.attrib.get('t')=='s':v=ss[int(v)]
       if v:cells.append(c.attrib['r']+'='+v)
      emit(rel,origin,sheet.attrib['name']+' row'+row.attrib.get('r','?'),' | '.join(cells),'WORKSHEET_ROW_HEADER_OR_DATA')
  elif suffix=='.docx':
   with zipfile.ZipFile(p) as z:
    need(sum(x.file_size for x in z.infolist())<MAX_FILE*3,'DOCX expansion exceeds60MB')
    doc=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    for i,node in enumerate(doc.findall('.//w:p',ns),1):emit(rel,origin,'paragraph'+str(i),' '.join(x.text or '' for x in node.findall('.//w:t',ns)))
  elif suffix=='.pdf':
   exe=shutil.which('pdftotext')
   if not exe:raise ValueError('PDF_TEXT_EXTRACTION_UNAVAILABLE_INSTALL_POPPLER_OR_REVIEW_PDF_MANUALLY')
   result=subprocess.run([exe,'-layout',str(p),'-'],capture_output=True,text=True,timeout=60,check=True)
   for i,page in enumerate(result.stdout.split('\f'),1):
    for j,para in enumerate(re.split(r'\n\s*\n',page),1):emit(rel,origin,'page'+str(i)+' block'+str(j),para)
  elif suffix=='.xml':
   doc=ET.parse(p).getroot()
   if doc.tag.split('}')[-1]=='article':
    for i,node in enumerate(doc.iter(),1):
     tag=node.tag.split('}')[-1]
     if tag in ['p','caption','table-wrap','article-title','ref']:emit(rel,origin,tag+' '+node.attrib.get('id',str(i)),' '.join(node.itertext()))
     if tag=='ref':references.append(dict(source_file=rel,reference_id=node.attrib.get('id',''),citation=' '.join(' '.join(node.itertext()).split())))
    for node in doc.iter():
     if node.tag.split('}')[-1] in ['supplementary-material','media','ext-link']:
      href=node.attrib.get('{http://www.w3.org/1999/xlink}href','')
      if href:references.append(dict(source_file=rel,reference_id='LINK',citation=href))
   else:
    for i,node in enumerate(doc,1):emit(rel,origin,'xml_element'+str(i),' '.join(node.itertext()))
  elif suffix in ['.html','.htm']:
   parser=TextHTML();parser.feed(p.read_text(errors='replace'))
   for i,part in enumerate(''.join(parser.parts).splitlines(),1):emit(rel,origin,'text_block'+str(i),part)
  else:
   for i,line in enumerate(p.read_text(errors='replace').splitlines(),1):emit(rel,origin,'line'+str(i),line)
  parse_status.append(dict(file=rel,status='PARSED',detail='Keyword matches are leads, not accepted biological definitions'))
 except (OSError,ValueError,KeyError,IndexError,ET.ParseError,zipfile.BadZipFile,subprocess.SubprocessError,SystemExit) as e:parse_status.append(dict(file=rel,status='PARSE_BLOCKED',detail=str(e)))
for tag,origin,url,ext in requests:
 print('Retrieving '+tag,flush=True);p=fetch(tag,origin,url,ext)
 if p:analyze(p,origin,url)
# Optional explicitly supplied folder only; no search through all project files.
if os.environ.get('WANG_SOURCE_DIR'):
 folder=(root/os.environ['WANG_SOURCE_DIR']).resolve(strict=True);need(folder.is_dir(),'WANG_SOURCE_DIR is not a directory');paths=sorted(p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in {'.pdf','.xlsx','.docx','.xml','.txt','.csv','.tsv','.zip'})
 need(len(paths)<=30,'More than30 local source files; supply a focused folder')
 for i,p in enumerate(paths):
  if p.stat().st_size>MAX_FILE or totalbytes+p.stat().st_size>MAX_BYTES:continue
  dest=out/'sources'/('local_'+str(i)+p.suffix.lower());shutil.copyfile(p,dest);totalbytes+=p.stat().st_size;analyze(dest,'USER_SUPPLIED_SOURCE_ATTRIBUTION_REQUIRES_REVIEW',p.name)
# Carry established source findings forward; new keyword hits never overwrite them automatically.
for n in ['01_integrated_locus_evidence.tsv','07_remaining_definition_decisions.tsv','06_source_document_provenance.tsv','06_Wang_S6_published_counts.tsv','08_refreshed_family_coverage.tsv']:shutil.copyfile(src/n,out/('source_'+n))
known=[dict(family='MEP',status='LABEL_RESOLVED_MEMBERSHIP_UNRESOLVED',source='Wang Supplemental FigureS9 legend PDFpage8; source hashes in source_06_source_document_provenance.tsv',finding='Proton:pyruvate cotransporter label; original counted loci and reference sequence membership not established'),dict(family='PPA',status='WANG_SCOPE_UNRESOLVED',source='Wang TableS6; Ferrari2020 doi10.1111/nph.16265 is a contextual lead only',finding='PPA count column separate from V_Ppase; related Portulaca literature uses Ppa for inorganic pyrophosphatase, but this does not establish Wang exact inclusion rules')]
table('01_retrieval_status.tsv',retrieval,['tag','source_class','url','final_url','status','bytes','sha256','error']);table('02_source_inventory.tsv',inventory,['file','source_class','provenance','bytes','sha256']);table('03_source_passages.tsv',passages,['passage_id','source_file','source_class','location','family_mentions','context','excerpt','identifier_tokens','interpretation']);table('04_identifier_leads.tsv',identifiers,['passage_id','source_class','family_mentions','identifier','status']);table('05_article_references_and_links.tsv',references,['source_file','reference_id','citation']);table('06_parse_status.tsv',parse_status,['file','status','detail']);table('07_definition_status.tsv',known,['family','status','source','finding'])
review=[]
for fam in ['MEP','PPA']:
 direct=[p for p in passages if p['source_class']=='WANG_PRIMARY' and fam in p['family_mentions'].split(';')];related=[p for p in passages if p['source_class']!='WANG_PRIMARY' and fam in p['family_mentions'].split(';')]
 review.append(dict(family=fam,Wang_matching_passages=len(direct),related_matching_passages=len(related),Wang_identifier_tokens=len({r['identifier'] for r in identifiers if r['source_class']=='WANG_PRIMARY' and fam in r['family_mentions'].split(';')}),accepted_definition='PENDING_SOURCE_REVIEW',accepted_counted_identifiers='NONE_AUTOMATICALLY',next_action='REVIEW_PRIMARY_PASSAGES_AND_CONTEXTUAL_ID_LEADS'))
table('08_clarification_summary.tsv',review,list(review[0]))
(out/'README.txt').write_text('Source clarification only. Primary Wang full text and supplement archive retrieved via EuropePMC, with one related maize MEP paper and one related Portulaca PPA paper. Related-paper terminology is not automatically Wang membership. Downloads and archive expansion bounded60MB aggregate,20MB per file,30 explicitly supplied local files,2 archive nesting levels; no project-wide scan. Archive members use generated safe filenames. XML paragraphs, PDF page blocks, workbook rows and DOCX paragraphs retain source locations. MEP pathway/ammonium/protease uses can be false keyword matches. TableS6 numbers and header mentions do not identify counted loci. Identifier patterns find contextual tokens, not proof of family function. No automatic definition or candidate promotion, no sequence searches or trees, ledger copied unchanged. Download/parse failures explicit; a zero match after failed retrieval is not source absence. Review03_source_passages.tsv and04_identifier_leads.tsv with source files before proposing reference panels. OptionalWANG_SOURCE_DIR supplies only a small focused folder; raw source content is never executed.\n')
(out/'AUTHOR_QUERY_DRAFT.txt').write_text('Subject: Clarification of MEP and PPA gene sets in Supplemental TableS6 (doi10.1093/plphys/kiad451)\n\nDear authors,\nWe are comparing photosynthesis-related gene families in Sesuvium and related species using your TableS6 as a reference. Could you share the gene identifiers and reference/query sequences used for MEP and PPA, ideally for P. oleracea and A. thaliana? For MEP, does the counted set correspond to a particular mesophyll envelope protein family, and which reference loci define it? For PPA, which biochemical family and subfamilies were included, and how was it distinguished from the separately listed V_Ppase category? It would also help to know whether counts represent genes or transcripts and how partial models were handled.\nThank you for your help.\n\nDraft only; nothing has been sent.\n')
qc=dict(status='SOURCE_CLARIFICATION_COLLECTION_COMPLETE_REVIEW_REQUIRED',downloaded=sum(r['status']=='DOWNLOADED' for r in retrieval),retrieval_failures=sum(r['status']!='DOWNLOADED' for r in retrieval),source_files=len(inventory),parse_failures=sum(r['status']!='PARSED' for r in parse_status),passages=len(passages),identifier_leads=len(identifiers),passage_limit_reached=len(passages)>=1500,ledger_unchanged=sha(src/'01_integrated_locus_evidence.tsv')==sha(out/'source_01_integrated_locus_evidence.tsv'),automatic_definition_acceptance=False,new_searches=0,new_trees=0)
(out/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');(out/'PROVENANCE.json').write_text(json.dumps(dict(step74_manifest_sha256=sha(src/'SHA256SUMS.txt'),retrieval_time_UTC=datetime.now(timezone.utc).isoformat(),primary_DOI='10.1093/plphys/kiad451'),indent=2)+'\n')
files=sorted(p for p in out.rglob('*') if p.is_file());(out/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(out))+'\n' for p in files));verify(out)
bundle=out.with_name(out.name+'_bundle.tar.gz')
with tarfile.open(bundle,'w:gz') as t:
 for p in files+[out/'SHA256SUMS.txt']:t.add(p,arcname='Wang_MEP_PPA_source_clarification/'+str(p.relative_to(out)))
print(json.dumps(qc,indent=2));print('BUNDLE='+str(bundle))
PY
