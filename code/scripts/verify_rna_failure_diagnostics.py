#!/usr/bin/env python3
"""Independently check retained aggregate RNA-failure diagnostic arithmetic.

This checks published summaries, not private per-marker counts or DNA identity.
Raw/private-count checks are separate explicitly recorded compute-host runs.
"""
import argparse,csv,itertools,json,math
from pathlib import Path
from project_paths import resolve_analysis_input
ROOT=Path(__file__).resolve().parents[1]

def read(path):
    with resolve_analysis_input(path.relative_to(ROOT), ROOT).open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def close(a,b):
    assert math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-12),(a,b)
def arithmetic(r):
    n,ibs0,ibs2,h,hab=map(lambda k:int(r[k]),['n','ibs0','ibs2','shared_hets','hets_ab'])
    disc=int(r['heterozygote_homozygote_discordance'])
    assert 0<=ibs0+ibs2<=n and disc==n-ibs0-ibs2
    assert hab==2*h+disc and 0<=h+int(r['shared_hom_alts'])<=ibs2
    close(r['relatedness'],2*(h-2*ibs0)/max(1,hab))
    close(r['loss_from_het_hom'],disc/max(1,hab));close(r['loss_from_IBS0'],4*ibs0/max(1,hab))
    if hab:close(1-float(r['relatedness']),float(r['loss_from_het_hom'])+float(r['loss_from_IBS0']))
def matching(data,expected_names):
    assert {tuple(sorted([r['sample_a'],r['sample_b']])) for r in data}==set(itertools.combinations(sorted(expected_names),2))
    closest={}
    for name in expected_names:
        options=[(float(r['relatedness']),r['sample_b'] if r['sample_a']==name else r['sample_a']) for r in data if name in (r['sample_a'],r['sample_b'])]
        best=max(v for v,n in options);hits=[n for v,n in options if v==best]
        assert len(hits)==1
        closest[name]=hits[0]
    for r in data:
        if r['expected_metadata_pair']=='True':
            assert closest[r['sample_a']]==r['sample_b'] and closest[r['sample_b']]==r['sample_a']

def verify(root):
    rev=root/'analysis';meta=read(root/'metadata/audit/gse311180_requant_manifest.tsv')
    patient={r['sample_id']:r['patient_id'] for r in meta};names=sorted(patient)
    combos={tuple(x) for x in itertools.combinations(names,2)}
    primary={tuple(sorted([r['sample_a'],r['sample_b']])):r for r in read(rev/'audit_gse311180_rna_all_pairs.tsv')}
    loss=read(rev/'audit_RNA_pairing_loss_decomposition.tsv');assert len(loss)==7
    for r in loss:
        arithmetic(r);base=primary[tuple(sorted([r['sample_a'],r['sample_b']]))]
        for k in ['n','ibs0','ibs2','shared_hets','hets_ab','shared_hom_alts']:assert r[k]==base[k]
        close(r['relatedness'],base['relatedness'])
        assert int(r['het_a_hom_b'])+int(r['hom_a_het_b'])==int(r['heterozygote_homozygote_discordance'])
        close(r['fraction_of_relatedness_loss_from_het_hom'],float(r['loss_from_het_hom'])/(1-float(r['relatedness'])))
    coverage=read(rev/'audit_RNA_fingerprint_site_coverage.tsv');assert len(coverage)==14
    for r in coverage:
        assert int(r['all_markers'])==16469
        assert int(r['covered_markers'])+int(r['zero_depth_sites'])==16469
        assert sum(int(r[k]) for k in ['callable_sites','rejected_low_depth','rejected_other_base','rejected_intermediate_AB'])==16469
        close(r['callable_fraction'],int(r['callable_sites'])/16469)
        assert 0<float(r['top1percent_markers_fraction_of_depth'])<=1
    for fname,size in [('audit_RNA_pairing_depth_strata.tsv',28),('audit_RNA_pairing_chromosome_strata.tsv',154)]:
        data=read(rev/fname);assert len(data)==size
        for r in data:arithmetic(r)
        for base in loss:
            group=[r for r in data if r['patient_id']==base['patient_id']]
            for k in ['n','ibs0','ibs2','shared_hets','hets_ab','shared_hom_alts','heterozygote_homozygote_discordance']:
                assert sum(int(r[k]) for r in group)==int(base[k]),(fname,base['patient_id'],k)
    grid=read(rev/'audit_RNA_pairing_depth_AB_grid.tsv');assert len(grid)==1092
    for dep,ab in itertools.product([3,7,15,30],[.1,.2,.3]):
        rs=[r for r in grid if int(r['depth'])==dep and float(r['min_het_AB'])==ab]
        assert {tuple(sorted([r['sample_a'],r['sample_b']])) for r in rs}==combos and len(rs)==91
        matching(rs,names)
        for r in rs:
            arithmetic(r)
            assert (r['expected_metadata_pair']=='True')==(patient[r['sample_a']]==patient[r['sample_b']])
            if dep==7 and ab==.2:
                base=primary[tuple(sorted([r['sample_a'],r['sample_b']]))]
                for k in ['n','ibs0','ibs2','shared_hets','hets_ab','shared_hom_alts']:assert r[k]==base[k]
                close(r['relatedness'],base['relatedness'])
    direction=read(rev/'audit_RNA_pairing_directional_dropout.tsv');assert len(direction)==14
    for base in loss:
        group=[r for r in direction if r['patient_id']==base['patient_id']]
        assert len(group)==2 and sum(int(r['discordant_sites']) for r in group)==int(base['heterozygote_homozygote_discordance'])
        for r in group:
            assert sum(int(r[k]) for k in ['hom_depth_7_14','hom_depth_15_29','hom_depth_30_99','hom_depth_100_plus'])==int(r['discordant_sites'])
            close(r['hom_call_side_zero_minor_fraction'],int(r['zero_minor_allele_sites'])/int(r['discordant_sites']))
            assert int(r['hom_reference_allele_sites'])+int(r['hom_alternate_allele_sites'])==int(r['discordant_sites'])
    fragment=read(rev/'audit_RNA_discordant_fragment_coordinate_diversity.tsv');assert len(fragment)==14
    for r in fragment:
        d=next(x for x in direction if x['patient_id']==r['patient_id'] and x['homozygous_call_sample']==r['homozygous_call_sample'])
        assert int(r['original_discordant_loci'])==int(d['discordant_sites'])
        loci=int(r['original_discordant_loci'])
        assert 0<=int(r['homo_loci_with_HQ_proper_templates'])<=loci
        assert 0<=int(r['het_loci_with_HQ_proper_templates'])<=loci
        assert int(r['homo_loci_with_1_or_2_coordinates'])+int(r['homo_loci_with_7plus_coordinates'])<=int(r['homo_loci_with_HQ_proper_templates'])
        assert 0<=int(r['homo_loci_with_other_allele_coordinates'])<=int(r['homo_loci_with_HQ_proper_templates'])
        assert 0<=int(r['het_loci_with_both_alleles_coordinates'])<=int(r['het_loci_with_HQ_proper_templates'])
        for prefix in ['homo','het']:
            assert 0<=float(r[prefix+'_distinct_coordinates_median'])<=float(r[prefix+'_qname_templates_median'])<=float(r[prefix+'_HQ_read_bases_median'])
        assert 0<=int(r['homo_total_distinct_coordinates'])<=int(r['homo_total_qname_templates'])
    transitions=read(rev/'audit_RNA_read_filter_call_transitions.tsv');assert len(transitions)==112
    strict=read(rev/'audit_RNA_pairing_unique_high_quality_reads.tsv');assert len(strict)==91
    for base in loss:
        group=[r for r in transitions if r['patient_id']==base['patient_id']]
        assert len(group)==16 and sum(int(r['sites']) for r in group)==16469
        assert sum(int(r['sites']) for r in group if r['primary_category']=='het_hom_discordance')==int(base['heterozygote_homozygote_discordance'])
        q=next(r for r in strict if r['expected_metadata_pair']=='True' and r['patient_a']==base['patient_id'])
        assert sum(int(r['sites']) for r in group if r['strict_category']=='het_hom_discordance')==int(q['heterozygote_homozygote_discordance'])
        assert sum(int(r['sites']) for r in group if r['strict_category']!='uncalled')==int(q['n'])

    assert {tuple(sorted([r['sample_a'],r['sample_b']])) for r in strict}==combos
    matching(strict,names)
    for r in strict:arithmetic(r)
    retained=read(rev/'audit_RNA_pairing_unique_high_quality_reads_overlap_retained.tsv');assert len(retained)==91
    assert {tuple(sorted([r['sample_a'],r['sample_b']])) for r in retained}==combos
    matching(retained,names)
    for r in retained:arithmetic(r)
    retained_transitions=read(rev/'audit_RNA_read_filter_call_transitions_overlap_retained.tsv');assert len(retained_transitions)==112
    for base in loss:
        group=[r for r in retained_transitions if r['patient_id']==base['patient_id']]
        assert len(group)==16 and sum(int(r['sites']) for r in group)==16469
        q=next(r for r in retained if r['expected_metadata_pair']=='True' and r['patient_a']==base['patient_id'])
        assert sum(int(r['sites']) for r in group if r['strict_category']=='het_hom_discordance')==int(q['heterozygote_homozygote_discordance'])
        assert sum(int(r['sites']) for r in group if r['strict_category']!='uncalled')==int(q['n'])
    cap=read(rev/'audit_RNA_pairing_full_read_cap50.tsv')
    capqc=read(rev/'audit_RNA_cap50_source_QC.tsv')
    capnames=sorted(r['sample_id'] for r in capqc)
    assert len(capnames) in (2,4,6) and len(set(capnames))==len(capnames)
    assert set(capnames)<=set(['DR118','DR144','DR145','DR148','DR73','DR75'])
    assert all(sum(patient[n]==pid for n in capnames)==2 for pid in {patient[n] for n in capnames})
    assert len(cap)==len(capnames)*(len(capnames)-1)//2
    matching(cap,capnames)
    for r in cap:
        arithmetic(r)
        assert (r['expected_metadata_pair']=='True')==(patient[r['sample_a']]==patient[r['sample_b']])
    captrans=read(rev/'audit_RNA_cap50_call_transitions.tsv');assert len(captrans)==8*len(capnames)
    for q in [r for r in cap if r['expected_metadata_pair']=='True']:
        group=[r for r in captrans if r['patient_id']==q['patient_a']]
        assert len(group)==16 and sum(int(r['sites']) for r in group)==16469
        base=next(r for r in loss if r['patient_id']==q['patient_a'])
        assert sum(int(r['sites']) for r in group if r['primary_category']=='het_hom_discordance')==int(base['heterozygote_homozygote_discordance'])
        assert sum(int(r['sites']) for r in group if r['cap50_category']!='uncalled')==int(q['n'])
        assert sum(int(r['sites']) for r in group if r['cap50_category']=='het_hom_discordance')==int(q['heterozygote_homozygote_discordance'])
    attempts=read(rev/'audit_RNA_cap50_alignment_attempts.tsv')
    assert len(attempts)==6 and {r['sample_id'] for r in attempts}==set(['DR118','DR144','DR145','DR148','DR73','DR75'])
    assert all(r['cap50_attempt_status'] in {'COMPLETE','FAILED_TIMEOUT','FAILED_ERROR'} for r in attempts)
    completed={r['sample_id'] for r in attempts if r['cap50_attempt_status']=='COMPLETE'}
    paired_complete={n for n in completed if sum(patient[n]==patient[y] for y in completed)==2}
    assert set(capnames)==paired_complete
    for r in attempts:
        assert int(r['seedPerWindowNmax'])==50
        if r['cap50_attempt_status']=='COMPLETE':
            assert r['native_alignment_complete']=='True'
            assert int(r['validated_FASTQ_pairs'])==int(r['completed_STAR_input_pairs'])
    assert len(capqc)==len(capnames)
    assert {r['sample_id'] for r in capqc}==set(capnames)
    primaryqc={r['sample_id']:r for r in read(rev/'audit_gse311180_raw_RNA_mapping_QC.tsv')}
    for r in capqc:
        assert r['status']=='PASS' and int(r['seedPerWindowNmax'])==50
        assert len(r['BAM_sha256'])==len(r['fingerprint_sha256'])==64
        # Whole-read-count agreement is checked against the primary producer audit.
        assert int(r['input_read_pairs'])==int(primaryqc[r['sample_id']]['STAR_input_pairs'])
    for fname,n in [('audit_RNA_CHARR_proxy.tsv',28),('audit_RNA_paired_MLE_proxy.tsv',364)]:
        data=read(rev/fname);assert len(data)==n
        assert all(r['RNA_contamination_calibrated']=='False' for r in data)
    expected=read(rev/'audit_gse311180_rna_expected_pairs.tsv')
    failed={r['patient_id'] for r in expected if r['interpretation']!='consistent_RNA_profiles'}
    assert failed=={'STB-084','STB-100'}
    return {'status':'PASS','primary_gates_changed':False,'default_loss_pairs':7,'marker_coverage_samples':14,'fragment_coordinate_directions':14,
      'depth_AB_comparisons':1092,'strict_read_comparisons_per_overlap_mode':91,'strict_read_overlap_modes':2,'full_read_cap50_comparisons':len(cap),'full_read_cap50_completed_libraries':len(capnames),
      'scope':'Independent retained aggregate arithmetic/coverage/source-count checks; private count/native comparisons are separately executed; not DNA identity, anatomy, or calibrated contamination'}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--report',type=Path)
    a=p.parse_args();out=verify(a.root);text=json.dumps(out,indent=2)+'\n'
    if a.report:a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(text)
    print(text)
