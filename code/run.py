#!/usr/bin/env python3
"""Run the bulk sensitivity and signature evidence analyses from processed inputs."""
from pathlib import Path, PurePosixPath
import argparse,csv,hashlib,os,shutil,subprocess,sys,zipfile
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

def verify():
 with (ROOT/'MANIFEST.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 for r in rows:
  p=ROOT/r['path']
  if not p.is_file() or p.stat().st_size!=int(r['bytes']) or digest(p)!=r['sha256']:raise ValueError('Missing or altered file: '+r['path'])
 return rows

def prepare(destination):
 if destination.exists():raise ValueError('Output directory already exists: '+str(destination))
 destination.mkdir(parents=True)
 with (ROOT/'data/CONTENTS.tsv').open() as f:records={r['path']:r for r in csv.DictReader(f,delimiter='\t')}
 for archive in sorted({r['archive'] for r in records.values()}):
  subset={n:r for n,r in records.items() if r['archive']==archive}
  with zipfile.ZipFile(ROOT/'data'/archive) as z:
   if {r['member'] for r in subset.values()}!=set(z.namelist()):raise ValueError('Dataset archive and index differ')
   for name,r in subset.items():
    parts=PurePosixPath(name)
    if parts.is_absolute() or '..' in parts.parts:raise ValueError('Invalid dataset path')
    data=z.read(r['member'])
    if len(data)!=int(r['bytes']) or hashlib.sha256(data).hexdigest()!=r['sha256']:raise ValueError('Altered dataset: '+name)
    p=destination/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
 for sub in ['config','scripts','analysis']:
  shutil.copytree(ROOT/'code'/sub,destination/sub,dirs_exist_ok=True)
 sig=destination/'results/tables/bulk/frozen_multicellular_signature.tsv';sig.parent.mkdir(parents=True,exist_ok=True)
 shutil.copy2(ROOT/'data/results/tables/bulk/frozen_multicellular_signature.tsv',sig)
 return destination

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',default='results/analysis');ap.add_argument('--prepare-only',action='store_true');ap.add_argument('--cells-only',action='store_true');args=ap.parse_args()
 verify();workspace=prepare(ROOT/args.output)
 if args.prepare_only:print('Prepared scientific inputs: '+str(workspace));return
 env=dict(os.environ,AIS_DATA_ROOT=str(workspace),AIS_OUTPUT=str(workspace/'computed'),PYTHONDONTWRITEBYTECODE='1')
 for script in ([] if args.cells_only else ['bulk_sensitivity.R'])+['signature_evidence.R']:
  subprocess.run(['Rscript',str(workspace/'analysis'/script)],cwd=workspace,env=env,check=True)
 verify();print('Analysis outputs: '+str(workspace/'computed'))
if __name__=='__main__':main()
