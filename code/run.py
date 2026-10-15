#!/usr/bin/env python3
"""Verify inputs, reproduce all retained results/figures, then run independent checks."""
from __future__ import annotations
import csv,gzip,hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path,PurePosixPath
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code/scripts'))
import retained_reproduction as reproduction
from audit_release import root_git_metadata


def safe(root,name):
    path=PurePosixPath(name)
    if not name or path.is_absolute() or str(path)!=name or '..' in path.parts or '\\' in name or '.git' in path.parts:
        raise ValueError('Unsafe path: '+name)
    target=root
    for part in path.parts:
        target=target/part
        if target.is_symlink():raise ValueError('Symlink path: '+name)
    if root.resolve() not in target.resolve().parents:raise ValueError('Path escapes root')
    return target


def verify(root=ROOT):
    root=Path(root).resolve()
    with safe(root,'MANIFEST.tsv').open(newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if reader.fieldnames!=['path','runtime','bytes','sha256']:raise ValueError('Invalid manifest header')
        records=list(reader)
    seen=set();targets=set()
    for row in records:
        if None in row or any(v is None for v in row.values()):raise ValueError('Malformed manifest')
        name=row['path'];p=safe(root,name)
        if name in seen or not p.is_file() or p.stat().st_size!=int(row['bytes']) or reproduction.digest(p)!=row['sha256']:
            raise ValueError('Changed, missing or duplicate source: '+name)
        seen.add(name)
        if row['runtime']:
            safe(root,row['runtime'])
            if row['runtime'] in targets:raise ValueError('Duplicate runtime path')
            targets.add(row['runtime'])
        if name.startswith('data/') and p.suffix in {'.py','.R','.sh'}:raise ValueError('Code inside data directory')
    git_kind=root_git_metadata(root)
    for directory,dirs,files in os.walk(root,followlinks=False):
        parent=Path(directory)
        for name in list(dirs):
            p=parent/name
            if p.is_symlink():raise ValueError('Symlink directory')
            if name=='.git':
                if parent!=root or git_kind!='directory':raise ValueError('Unexpected Git metadata')
                dirs.remove(name)
            elif name=='__pycache__' or (parent==root and name=='results'):dirs.remove(name)
        for name in files:
            p=parent/name;rel=p.relative_to(root).as_posix()
            if p.is_symlink():raise ValueError('Symlink file')
            if parent==root and name=='.git' and git_kind=='file':continue
            if name.endswith('.pyc') and '__pycache__' in p.parts:continue
            if rel!='MANIFEST.tsv' and rel not in seen:raise ValueError('Uninventoried source: '+rel)
    return records


def restore_table(root,name,table):
    path=safe(root,name)
    if path.exists() or not name.endswith(('.tsv','.tsv.gz')):raise ValueError('Duplicate/invalid expected table')
    columns=table['columns'];rows=table['rows']
    if not columns or len(set(columns))!=len(columns) or any(len(r)!=len(columns) for r in rows):raise ValueError('Malformed expected table')
    if not all(isinstance(v,str) for row in [columns,*rows] for v in row):raise ValueError('Expected cells must retain exact string values')
    path.parent.mkdir(parents=True,exist_ok=True)
    with (gzip.open(path,'wt',newline='') if path.suffix=='.gz' else path.open('w',newline='')) as f:
        writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow(columns);writer.writerows(rows)


def prepare(root,reference):
    records=verify(root);reference.mkdir()
    for row in records:
        if not row['runtime']:continue
        dst=safe(reference,row['runtime']);dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(safe(root,row['path']),dst)
    with gzip.open(root/'data/expected_results.json.gz','rt') as f:expected=json.load(f)
    if expected['format']!=1:raise ValueError('Unsupported expected-result format')
    for name,table in expected['tables'].items():restore_table(reference,name,table)
    with (reference/'SHA256SUMS.tsv').open('w',newline='') as f:
        writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow(['path','bytes','sha256'])
        for p in sorted(reference.rglob('*')):
            if p.is_file() and p.name!='SHA256SUMS.tsv':writer.writerow([p.relative_to(reference).as_posix(),p.stat().st_size,reproduction.digest(p)])
    reproduction.manifest(reference)
    return reference


def main():
    verify()
    output=ROOT/'results';output.mkdir(exist_ok=True)
    run=Path(tempfile.mkdtemp(prefix='run-',dir=output))
    reference=prepare(ROOT,run/'reference')
    workspace=reproduction.stage(reference)
    # Keep one visible, predictable output location without changing baseline checks.
    destination=run/'workspace'
    shutil.move(str(workspace.parent/'initial_file_states.json'),run/'initial_file_states.json')
    shutil.move(str(workspace),destination);workspace=destination
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
    subprocess.run(['bash','scripts/lib_reproduce_revision_stages.sh'],cwd=workspace,env=env,check=True)
    comparison=reproduction.compare(reference,workspace)
    (run/'numerical_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    subprocess.run(['bash','validate_release.sh'],cwd=workspace,env=env,check=True)
    verify()
    receipt={'status':'PASS','source_manifest_sha256':reproduction.digest(ROOT/'MANIFEST.tsv'),
             'reference_files':comparison['retained_reference_files_verified'],
             'scope':'Independent retained-input result/figure reproduction, source joins, model checks and regression tests; no new raw-read alignment',
             'outputs':'workspace','numerical_comparison':'numerical_comparison.json'}
    (run/'PASS.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS. Results and verification: '+str(run))
if __name__=='__main__':main()
