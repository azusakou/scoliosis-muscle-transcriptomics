"""Stage retained re-analysis without overwriting checksum-pinned reference tables."""
import argparse,csv,gzip,hashlib,json,math,re,shutil,tempfile
from decimal import Decimal,InvalidOperation
from pathlib import Path,PurePosixPath

REL_TOL=1e-10
ABS_TOL=1e-12
# The same scientific model-refit tolerance already used by
# verify_retained_models.R; only these explicit derived model fields use it.
# This does not relax arithmetic, counts, inputs, null grids or decision gates.
MODEL_TOL=1e-7
GATE_MODEL_COLUMNS={
    'analysis/bulk_gate_observed_genes_locked_original_TMM.tsv.gz':
        {'logFC','logCPM','F','PValue','FDR'},
    'analysis/bulk_gate_observed_genes_retained_library_reTMM.tsv.gz':
        {'logFC','logCPM','F','PValue','FDR'},
    'analysis/bulk_gate_B32_logFC_locked_original_TMM.tsv.gz':
        {f'configuration_{i:03d}' for i in range(1,33)},
    'analysis/bulk_gate_B32_logFC_retained_library_reTMM.tsv.gz':
        {f'configuration_{i:03d}' for i in range(1,33)},
    'analysis/bulk_gate_ANGPT1_LEP_fixed_candidate_results.tsv': {
        'fitted_shrunk_logFC','MLE_unshrunk_logFC','QL95_lower','QL95_upper',
        'P_reproduction_absolute_error','QL_F','QL_P','genome_wide_FDR_original_17961',
        'BH_two_fixed_candidates_within_normalization',
        'BH_four_candidate_normalization_sensitivity_family'},
}
DERIVED_EXPRESSION_AUDITS={
    'analysis/audit_gse311180_all_column_expression_comparisons.tsv',
    'analysis/audit_gse311180_same_name_and_closest_expression.tsv',
}
DERIVED_GATE_GZIP_OUTPUTS={
    'analysis/bulk_gate_observed_genes_locked_original_TMM.tsv.gz',
    'analysis/bulk_gate_observed_genes_retained_library_reTMM.tsv.gz',
    'analysis/bulk_gate_B32_logFC_locked_original_TMM.tsv.gz',
    'analysis/bulk_gate_B32_logFC_retained_library_reTMM.tsv.gz',
}

IDENTIFIER_COLUMNS={'id','gene','patient_id','sample_id','barcode','cohort','side',
    'cell_type','lineage','assay','normalization','cluster','configuration'}

def probability_column(column):
    key=column.lower().replace('.','_')
    if key.startswith(('n_','p_reproduction_')):return False
    return bool(re.search(r'(^|_)(p|pvalue|pval|p_value|p_val|p_adj|padj|p_adjusted|qvalue|q_value|qval|q_val|fdr|bh)(_|$)|_tail$',key))

def effect_or_interval_column(column):
    key=column.lower()
    return ('logfc' in key or key.startswith(('ql95_','configuration_')) or
            key=='patient_paired_logcpm_mean' or
            ('ci' in key and any(t in key for t in ('lower','upper','low','high','95_'))))

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def manifest(root):
    with (root/'SHA256SUMS.tsv').open(newline='') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    seen=set()
    for r in rows:
        name=PurePosixPath(r['path']);p=root/name
        if name.is_absolute() or '..' in name.parts or '\\' in r['path'] or r['path'] in seen:
            raise ValueError('Unsafe or duplicate reference path')
        seen.add(r['path'])
        if p.is_symlink() or root.resolve() not in p.resolve().parents:
            raise ValueError('Reference escapes release root')
        if not p.is_file() or p.stat().st_size!=int(r['bytes']) or digest(p)!=r['sha256']:
            raise ValueError('Changed retained reference: '+r['path'])
    return rows

def compare_tsv(reference,reproduced,model_columns=frozenset()):
    opener=lambda p: gzip.open(p,'rt',newline='') if p.suffix=='.gz' else p.open(newline='')
    with opener(reference) as f:a=list(csv.reader(f,delimiter='\t'))
    with opener(reproduced) as f:b=list(csv.reader(f,delimiter='\t'))
    if not a or len(a)!=len(b) or a[0]!=b[0]:raise ValueError('Changed table rows/header')
    if not model_columns.issubset(a[0]):raise ValueError('Missing declared model field')
    differences=0;maximum=0.
    for x,y in zip(a[1:],b[1:]):
        if len(x)!=len(a[0]) or len(y)!=len(a[0]):raise ValueError('Changed table columns')
        for column,old,new in zip(a[0],x,y):
            if old==new:continue
            # IDs are labels, even when their spelling happens to be numeric.
            if column.lower() in IDENTIFIER_COLUMNS:
                raise ValueError('Changed table identifier/text label')
            try:u,v=Decimal(old),Decimal(new)
            except InvalidOperation:raise ValueError('Changed table text')
            if not u.is_finite() or not v.is_finite():raise ValueError('Changed nonfinite value')
            # Qualitative decisions are protected for ordinary arithmetic too,
            # not only for the explicitly declared model fields.
            if probability_column(column) and (not (0<=u<=1 and 0<=v<=1) or
                    (u<Decimal('.05'))!=(v<Decimal('.05'))):
                raise ValueError('Changed probability decision or range')
            if effect_or_interval_column(column) and ((u>0)-(u<0))!=((v>0)-(v<0)):
                raise ValueError('Changed effect or interval sign')
            if column in model_columns:
                # Model F/P fields can serialize as integer tokens (F=0, P=1),
                # unlike sample counts, which remain exact outside this list.
                if abs(u-v)>Decimal(str(MODEL_TOL))*max(Decimal(1),abs(u)):
                    raise ValueError('Changed model result outside existing model-refit tolerance')
            # Exact integer/count tokens must not pass via float rounding or relative tolerance.
            elif re.fullmatch(r'[+-]?\d+',old) or re.fullmatch(r'[+-]?\d+',new):
                if u!=v:raise ValueError('Changed integer/count value')
            elif not math.isclose(float(u),float(v),rel_tol=REL_TOL,abs_tol=ABS_TOL):
                raise ValueError('Changed numerical result outside fixed tolerance')
            differences+=1;maximum=max(maximum,float(abs(u-v)))
    return {'rows':len(a)-1,'numeric_cells_differing':differences,'max_absolute_difference':maximum,
        'model_columns_using_existing_refit_tolerance':sorted(model_columns)}

def stage(root):
    root=root.resolve();rows=manifest(root)
    log=root/'logs/reproductions';log.mkdir(parents=True,exist_ok=True)
    run=Path(tempfile.mkdtemp(prefix='run-',dir=log));workspace=run/'workspace';workspace.mkdir()
    for r in rows:
        dst=workspace/r['path'];dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(root/r['path'],dst)
    shutil.copy2(root/'SHA256SUMS.tsv',workspace/'SHA256SUMS.tsv')
    manifest(workspace)
    # This derived table is deliberately not distributed twice: its published
    # reference is the fixed parameter file.  It must be freshly generated.
    generated_normalization=workspace/'analysis/bulk_gate_normalization_factors.tsv'
    if generated_normalization.exists():
        raise ValueError('Unexpected copied normalization output; use the canonical parameter input only')
    (run/'initial_file_states.json').write_text(json.dumps({
        r['path']:{'mtime_ns':(workspace/r['path']).stat().st_mtime_ns,
                   'sha256':r['sha256']} for r in rows},indent=2)+'\n')
    return workspace

def prepare_raw(root,workspace):
    """Create a mutable raw-analysis checkout while retaining the release baseline."""
    root=root.resolve();workspace=workspace.resolve();rows=manifest(root)
    if workspace==root or root in workspace.parents or workspace in root.parents:
        raise ValueError('Raw workspace must be separate from the release root')
    if workspace.exists():
        raise ValueError('Raw workspace destination already exists; refusing overwrite')
    workspace.mkdir(parents=True)
    for r in rows:
        dst=workspace/r['path'];dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(root/r['path'],dst)
    baseline=workspace/'logs/release_baseline';baseline.mkdir(parents=True)
    shutil.copy2(root/'SHA256SUMS.tsv',baseline/'SHA256SUMS.tsv')
    # A raw run can change metadata path fields and numerical artifacts. It must
    # not present the old immutable manifest as a hash of those new outputs.
    (baseline/'SCOPE.txt').write_text('Original supplied reference manifest only; raw reconstruction outputs are mutable and are not certified by these old hashes.\n')
    return workspace

def compare(root,workspace):
    root=root.resolve();workspace=workspace.resolve();rows=manifest(root)
    changed=[];exact=0
    for r in rows:
        p=workspace/r['path']
        if not p.is_file() or p.is_symlink() or workspace not in p.resolve().parents:
            raise ValueError('Missing or unsafe reproduced artifact')
        if digest(p)==r['sha256']:exact+=1;continue
        path=PurePosixPath(r['path'])
        # Only retained revision output TSVs may vary numerically. Code, metadata,
        # gzip inputs, primary results and native RNA audits stay exact. The two
        # explicitly named expression-correlation audits are derived floats.
        allowed=(path.parent==PurePosixPath('analysis') and path.name.startswith(('bulk_','signature_'))) or (path.parent==PurePosixPath('analysis/tables') and path.name.startswith('figure_')) or r['path'] in DERIVED_EXPRESSION_AUDITS
        if not ((allowed and path.suffix=='.tsv') or r['path'] in DERIVED_GATE_GZIP_OUTPUTS):
            raise ValueError('Changed immutable code/input: '+r['path'])
        changed.append({'path':r['path'],**compare_tsv(root/r['path'],p,
            GATE_MODEL_COLUMNS.get(r['path'],frozenset()))})
    normalization=verify_generated_normalization(root,workspace)
    figures=None
    if (root/'config/replayed_figures.json').is_file():
        from verify_replayed_figures import verify as verify_figures
        figures=verify_figures(workspace)
    # This receipt distinguishes actually written outputs from immutable input
    # tables carried into the workspace; identical hashes alone cannot do so.
    initial_path=workspace.parent/'initial_file_states.json'
    initial=json.loads(initial_path.read_text()) if initial_path.is_file() else {}
    regenerated=[r['path'] for r in rows if r['path'] in initial and
                 (workspace/r['path']).stat().st_mtime_ns!=initial[r['path']]['mtime_ns']]
    return {'status':'PASS','retained_reference_files_verified':len(rows),
        'reproduced_files_byte_identical':exact,'numerically_equivalent_tables':changed,
        'fresh_normalization_output_compared_to_canonical_published_input':normalization,
        'fresh_figure_outputs':figures,
        'actually_rewritten_reference_paths':regenerated,
        'relative_tolerance':REL_TOL,'absolute_tolerance':ABS_TOL,
        'derived_expression_audit_allowlist':sorted(DERIVED_EXPRESSION_AUDITS),
        'derived_gate_gzip_output_allowlist':sorted(DERIVED_GATE_GZIP_OUTPUTS),
        'model_refit_tolerance':MODEL_TOL,
        'model_refit_field_allowlist':{k:sorted(v) for k,v in GATE_MODEL_COLUMNS.items()},
        'integer_text_structure_checks':'Exact; any changed count, text, header or row count fails',
        'scope':'Retained-table reproduction only; scientific model-refit tolerances unchanged; no raw RNA alignment or new DNA identity claim'}

def verify_generated_normalization(root,workspace):
    reference=root/'analysis/bulk_gate_inputs/published_retained_TMM_factors.tsv'
    if not reference.is_file():
        # compare() also serves small non-gate regression fixtures. The full
        # public release pins this required reference in its file manifest.
        return None
    generated=workspace/'analysis/bulk_gate_normalization_factors.tsv'
    if not generated.is_file() or generated.is_symlink():
        raise ValueError('Missing freshly generated normalization output')
    return compare_tsv(reference,generated)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['stage','compare','prepare-raw'])
    p.add_argument('--root',type=Path,required=True);p.add_argument('--workspace',type=Path)
    a=p.parse_args()
    if a.mode=='stage':print(stage(a.root))
    elif a.mode=='prepare-raw':
        if not a.workspace:p.error('--workspace required')
        print(prepare_raw(a.root,a.workspace))
    else:
        if not a.workspace:p.error('--workspace required')
        out=compare(a.root,a.workspace)
        (a.workspace.parent/'numerical_comparison.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(out,indent=2))
if __name__=='__main__':main()
