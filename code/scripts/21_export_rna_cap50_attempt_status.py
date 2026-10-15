#!/usr/bin/env python3
"""Export terminal cap50 attempt metadata, including unsuccessful alignments.

This verifies source-bound receipts and reports their actual scope. It does not
turn timeout/partial outputs into successful alignments or establish DNA identity.
Native allele-count/BAM checks for completed pairs are a separate stage18 run.
"""
import argparse,importlib,json,re
from pathlib import Path
from input_integrity import digest,validate_pairs
m=importlib.import_module('14_rna_pairing_failure_diagnostics')
ROOT=Path(__file__).resolve().parents[1]
SAMPLES=['DR118','DR144','DR145','DR148','DR73','DR75']

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',nargs='+',type=Path,required=True)
    p.add_argument('--primary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--manifest',type=Path,default=ROOT/'metadata/audit/gse311180_requant_manifest.tsv')
    a=p.parse_args();meta=m.rows(a.manifest);validate_pairs(meta,expected_pairs=7)
    patient={r['sample_id']:r['patient_id'] for r in meta}
    qc={r['sample_id']:r for r in m.rows(ROOT/'analysis/audit_gse311180_raw_RNA_mapping_QC.tsv')}
    source={r['sample_id']:r for r in m.rows(ROOT/'analysis/audit_gse311180_completed_raw_source_sha256.tsv')}
    rows=[]
    for name in SAMPLES:
        candidates=[r/name for r in a.runs if (r/name/'completed_diagnostic.json').is_file()]
        if len(candidates)!=1:raise ValueError('Require one actual terminal receipt per requested library: '+name)
        d=candidates[0];r=json.loads((d/'completed_diagnostic.json').read_text());inputs=r['inputs']
        assert inputs['sample']==name and inputs['seedPerWindowNmax']==50
        assert inputs['primary_fingerprint_sha256']==digest(a.primary/'fingerprints'/f'{name}.somalier')==source[name]['RNA_fingerprint_sha256']
        assert inputs['prior_verified_FASTQ_sha256']==[source[name]['mate1_sha256'],source[name]['mate2_sha256']]
        assert inputs['prior_source_provenance_sha256']==source[name]['FASTQ_provenance_sha256']
        assert inputs['reference_binding_sha256']==digest(a.primary/'reference_artifact_binding.json')
        alignment=d/'alignment.log';native_log=d/'Log.out'
        native_complete=(native_log.is_file() and 'ALL DONE!' in native_log.read_text() and alignment.is_file() and 'finished successfully' in alignment.read_text())
        if r['status']=='PASS':
            status='COMPLETE';assert native_complete and int(r['STAR_input_pairs'])==int(qc[name]['validated_read_pairs'])
            assert digest(d/f'{name}.somalier')==r['fingerprint_sha256']
        elif r['status']=='FAILED_PRESERVED':
            status='FAILED_TIMEOUT' if 'timed out' in r['error'] else 'FAILED_ERROR'
        else:raise ValueError('Unsupported/unfinished receipt: '+name)
        rows.append(dict(sample_id=name,patient_id=patient[name],cap50_attempt_status=status,
            controller_reported_status=r['status'],native_alignment_complete=native_complete,
            recovered_after_legacy_wrapper_timeout=bool(r.get('recovery')),
            initial_controller_status=r.get('recovery',{}).get('original_controller_status',r['status']),
            validated_FASTQ_pairs=qc[name]['validated_read_pairs'],completed_STAR_input_pairs=r.get('STAR_input_pairs',''),
            seedPerWindowNmax=50,producer_script_sha256=inputs['source_script_sha256'],
            receipt_sha256=digest(d/'completed_diagnostic.json'),alignment_log_sha256=digest(alignment),native_log_sha256=digest(native_log),
            progress_log_sha256=digest(d/'Log.progress.out'),
            failed_stage=('STAR_alignment' if '--twopassMode' in r.get('error','') else 'other_native_stage') if status!='COMPLETE' else '',
            stage_timeout_seconds=(re.search(r'timed out after ([0-9.]+) seconds',r.get('error','')).group(1) if re.search(r'timed out after ([0-9.]+) seconds',r.get('error','')) else ''),
            failure_type='stage_timeout' if status=='FAILED_TIMEOUT' else ('native_stage_error' if status=='FAILED_ERROR' else ''),
            interpretation_scope='Attempt metadata only; failed/partial output excluded from pairing statistics'))
    m.write(a.output/'audit_RNA_cap50_alignment_attempts.tsv',rows)
    proof={'status':'PASS','requested_libraries':6,'source_bound_terminal_receipts':6,
        'complete_libraries':sum(r['cap50_attempt_status']=='COMPLETE' for r in rows),
        'failed_libraries':sum(r['cap50_attempt_status']!='COMPLETE' for r in rows),
        'script_sha256':digest(Path(__file__)),'scope':'Terminal source-bound attempt metadata; not BAM/pair verification or DNA identity'}
    (a.output/'cap50_attempt_status_QA.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(proof,indent=2))
if __name__=='__main__':main()
