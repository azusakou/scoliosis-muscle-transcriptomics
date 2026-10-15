#!/usr/bin/env python3
"""Run all12 post hoc native RNA depth/balance settings and uncalibrated proxies.

All 14 primary sketches remain private. All1092 native integer comparisons are
checked against a separate format2 decoder; original gates are never changed.
"""
import argparse,importlib,itertools,json,os,subprocess
from pathlib import Path
from input_integrity import digest,validate_pairs
from project_paths import executable
m=importlib.import_module('14_rna_pairing_failure_diagnostics')
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['identity','sites','output']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--manifest',type=Path,default=ROOT/'metadata/audit/gse311180_requant_manifest.tsv')
    a=p.parse_args()
    if a.output.exists() and any(a.output.iterdir()):raise ValueError('Preserve previous diagnostic outputs; choose a new directory')
    a.output.mkdir(parents=True,exist_ok=True)
    meta=m.rows(a.manifest);validate_pairs(meta,expected_pairs=7)
    patient={r['sample_id']:r['patient_id'] for r in meta};names=sorted(patient)
    tool=executable('somalier','AIS_SOMALIER');ver=subprocess.run([tool,'--help'],capture_output=True,text=True)
    assert 'somalier version: 0.3.5' in ver.stdout+ver.stderr
    binding=json.loads((a.identity/'reference_artifact_binding.json').read_text())['binding']['references']
    assert digest(a.sites)==binding['compressed']['sites.hg38.rna.vcf.gz']
    sketches=[a.identity/'fingerprints'/f'{n}.somalier' for n in names];counts={}
    for n,file in zip(names,sketches):
        name,size,cs=m.decode(file);assert name==n and size[0]==len(cs)==16469
        original=json.loads((a.identity/f'{n}.provenance.json').read_text())
        assert digest(file)==original['outputs']['fingerprint_sha256'];counts[n]=cs
    env=dict(os.environ,SOMALIER_REPORT_ALL_PAIRS='1');proof=[]
    for depth,ab in itertools.product([3,7,15,30],[.1,.2,.3]):
        prefix=a.output/f'd{depth}_ab{ab:g}'
        with Path(str(prefix)+'.log').open('w') as log:
            subprocess.run([tool,'relate','--sites',str(a.sites),'--groups',str(a.identity/'expected_pairs.txt'),
                '--min-depth',str(depth),'--min-ab',str(ab),'--output-prefix',str(prefix),*map(str,sketches)],
                check=True,env=env,stdout=log,stderr=subprocess.STDOUT)
        native=m.rows(Path(str(prefix)+'.pairs.tsv'));assert len(native)==91
        calls={n:[m.gt(c,depth,ab) for c in counts[n]] for n in names}
        for r in native:
            v=m.compare(calls[r['sample_a']],calls[r['sample_b']])
            for k in ['n','ibs0','ibs2','shared_hets','hets_ab','shared_hom_alts']:assert v[k]==int(r[k])
            assert abs(v['relatedness']-float(r['relatedness']))<=.0005001
        proof.append({'depth':depth,'min_AB':ab,'all91_integer_statistics':'PASS',
                      'native_pairs_sha256':digest(Path(str(prefix)+'.pairs.tsv'))})
    samples=[];pairs=[]
    for depth in [15,30]:
        prefix=a.output/f'native_contamination_d{depth}'
        with Path(str(prefix)+'.log').open('w') as log:
            subprocess.run([tool,'contamination','--sites',str(a.sites),'--min-depth',str(depth),
                '--output',str(prefix),*map(str,sketches)],check=True,stdout=log,stderr=subprocess.STDOUT)
        ss=m.rows(Path(str(prefix)+'.samples.tsv'));pp=m.rows(Path(str(prefix)+'.pairs.tsv'))
        assert len(ss)==14 and len(pp)==182
        for r in ss:samples.append(dict(depth=depth,sample_id=r['sample_name'],patient_id=patient[r['sample_name']],
            n_sites_usable=r['n_sites_usable'],native_CHARR_proxy=r['contamination_charr'],RNA_contamination_calibrated=False))
        for r in pp:pairs.append(dict(depth=depth,sample_id=r['sample_name'],anchor_sample=r['anchor_sample'],
            expected_metadata_pair=patient[r['sample_name']]==patient[r['anchor_sample']],n_sites_usable=r['n_sites_usable'],
            native_paired_MLE_proxy=r['contamination_mle'],RNA_contamination_calibrated=False))
    m.write(a.output/'audit_RNA_CHARR_proxy.tsv',samples);m.write(a.output/'audit_RNA_paired_MLE_proxy.tsv',pairs)
    result={'status':'PASS','tool_sha256':digest(Path(tool)),'settings':proof,'total_integer_pair_comparisons':1092,
        'native_proxy_sample_records':28,'native_proxy_ordered_pair_records':364,'script_sha256':digest(Path(__file__)),
        'decoder_sha256':digest(Path(m.__file__)),'scope':'Post hoc; all original depth7/AB0.2 and interpretation gates unchanged; RNA proxies not calibrated contamination percentages'}
    (a.output/'native_grid_independent_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print('ALL1092_NATIVE_DIAGNOSTIC_COMPARISONS_AND_PROXY_COMMANDS_PASS')
if __name__=='__main__':main()
