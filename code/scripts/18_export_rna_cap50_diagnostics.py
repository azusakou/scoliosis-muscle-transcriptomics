#!/usr/bin/env python3
"""Validate and export only aggregate full-read cap50 sensitivity outputs.

Keep all private BAMs/counts/sketches on the compute host. This post hoc
six-library alignment comparison does not replace the uniform primary recipe.
"""
import argparse,importlib,itertools,json,os,subprocess
from pathlib import Path
from input_integrity import digest
from project_paths import executable
m=importlib.import_module('14_rna_pairing_failure_diagnostics')
t=importlib.import_module('19_rna_read_filter_transitions')
f=importlib.import_module('15_rna_pairing_cap50_sensitivity')
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--primary',type=Path,required=True)
    p.add_argument('--completed',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--manifest',type=Path,default=ROOT/'metadata/audit/gse311180_requant_manifest.tsv')
    p.add_argument('--generate-native-pairs',action='store_true',help='Generate native comparisons for an explicitly selected completed subset in a new output directory')
    p.add_argument('--sites',type=Path,help='Pinned RNA sites required with --generate-native-pairs')
    p.add_argument('--samples',nargs='+',default=['DR144','DR73','DR118','DR148','DR145','DR75'],help='Only complete metadata pairs with verified completed receipts; default all six requested libraries')
    a=p.parse_args();meta=m.rows(a.manifest);patient={r['sample_id']:r['patient_id'] for r in meta}
    groups=f.selected_pairs(meta,a.samples);names=sorted(a.samples);counts={};qc=[]
    assert set(names)<=set(['DR144','DR73','DR118','DR148','DR145','DR75'])
    comparisons=len(names)*(len(names)-1)//2
    for n in names:
        folder=a.completed/n;r=json.loads((folder/'completed_diagnostic.json').read_text())
        assert r['status']=='PASS' and r['inputs']['sample']==n and r['inputs']['seedPerWindowNmax']==50
        fp=folder/f'{n}.somalier';sample,size,cs=m.decode(fp)
        assert sample==n and size[0]==len(cs)==16469 and digest(fp)==r['fingerprint_sha256']
        assert digest(folder/'Aligned.sortedByCoord.out.bam')==r['BAM_sha256']
        counts[n]=[m.gt(c) for c in cs]
        qc.append(dict(sample_id=n,patient_id=patient[n],status=r['status'],seedPerWindowNmax=50,
            input_read_pairs=r['STAR_input_pairs'],uniquely_mapped_percent=r['STAR_uniquely_mapped_percent'],
            BAM_sha256=r['BAM_sha256'],fingerprint_sha256=r['fingerprint_sha256'],
            mate1_sha256=r['inputs']['prior_verified_FASTQ_sha256'][0],mate2_sha256=r['inputs']['prior_verified_FASTQ_sha256'][1],
            reference_binding_sha256=r['inputs']['reference_binding_sha256'],producer_script_sha256=r['inputs']['source_script_sha256']))
    native_file=a.completed/'cap50.pairs.tsv'
    if a.generate_native_pairs:
        assert a.sites and digest(a.sites)==json.loads((a.primary/'reference_artifact_binding.json').read_text())['binding']['references']['compressed']['sites.hg38.rna.vcf.gz']
        if a.output.exists() and any(a.output.iterdir()):raise ValueError('Generate comparisons in a new private output directory')
        a.output.mkdir(parents=True,exist_ok=True)
        som=executable('somalier','AIS_SOMALIER')
        version=subprocess.run([som,'--help'],capture_output=True,text=True)
        assert 'somalier version: 0.3.5' in version.stdout+version.stderr
        groups_file=a.output/'selected_complete_pairs.txt';groups_file.write_text(''.join(','.join(v)+'\n' for v in groups))
        with (a.output/'native_relate.log').open('w') as log:
            subprocess.run([som,'relate','--sites',str(a.sites),'--groups',str(groups_file),'--min-depth','7','--min-ab','0.2','--output-prefix',str(a.output/'native_cap50'),*[str(a.completed/n/f'{n}.somalier') for n in names]],check=True,env=dict(os.environ,SOMALIER_REPORT_ALL_PAIRS='1'),stdout=log,stderr=subprocess.STDOUT)
        native_file=a.output/'native_cap50.pairs.tsv'
    native=m.rows(native_file)
    assert len(native)==comparisons and {tuple(sorted([r['sample_a'],r['sample_b']])) for r in native}==set(itertools.combinations(names,2))
    out=[]
    for r in native:
        x,y=r['sample_a'],r['sample_b'];v=m.compare(counts[x],counts[y])
        for k in ['n','ibs0','ibs2','shared_hets','hets_ab','shared_hom_alts']:assert v[k]==int(r[k])
        assert abs(v['relatedness']-float(r['relatedness']))<=.0005001
        out.append(dict(sample_a=x,sample_b=y,patient_a=patient[x],patient_b=patient[y],
            expected_metadata_pair=patient[x]==patient[y],somalier_reported_relatedness=r['relatedness'],**v))
    m.write(a.output/'audit_RNA_pairing_full_read_cap50.tsv',out);m.write(a.output/'audit_RNA_cap50_source_QC.tsv',qc)
    primary={}
    for n in names:
        file=a.primary/'fingerprints'/f'{n}.somalier';sample,size,cs=m.decode(file)
        original=json.loads((a.completed/n/'completed_diagnostic.json').read_text())
        assert sample==n and size[0]==len(cs)==16469 and digest(file)==original['inputs']['primary_fingerprint_sha256']
        primary[n]=[m.gt(c) for c in cs]
    transitions=[]
    for x,y in itertools.combinations(names,2):
        if patient[x]!=patient[y]:continue
        before=[t.category(a,b) for a,b in zip(primary[x],primary[y])]
        after=[t.category(a,b) for a,b in zip(counts[x],counts[y])]
        for start in ['uncalled','concordant','opposite_homozygotes','het_hom_discordance']:
            for end in ['uncalled','concordant','opposite_homozygotes','het_hom_discordance']:
                transitions.append(dict(patient_id=patient[x],sample_a=x,sample_b=y,primary_category=start,cap50_category=end,
                    sites=sum(a==start and b==end for a,b in zip(before,after))))
    assert len(transitions)==16*len(groups);m.write(a.output/'audit_RNA_cap50_call_transitions.tsv',transitions)

    proof={'status':'PASS','samples':len(names),'complete_read_alignments':len(names),'native_and_independent_comparisons':comparisons,'sample_ids':names,'native_comparisons_generated_by_exporter':a.generate_native_pairs,
        'original_gates_changed':False,'producer_script_sha256s':sorted({r['producer_script_sha256'] for r in qc}),
        'export_script_sha256':digest(Path(__file__)),'scope':'Post hoc cap50 sensitivity of explicitly selected completed pairs; no absent/failed alignment is counted; no DNA identity or side authentication'}
    (a.output/'independent_full_read_cap50_QA.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(proof,indent=2))
if __name__=='__main__':main()
