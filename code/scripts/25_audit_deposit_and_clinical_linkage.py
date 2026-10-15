#!/usr/bin/env python3
"""Post hoc all-column expression comparison and source-only clinical audit.

Expression similarity and sex markers cannot identify anatomical side or
authorize clinical-ID merges. Original labels and primary results stay fixed.
"""
import argparse
import csv
import gzip
import json
import math
import statistics
import xml.etree.ElementTree as ET
from pathlib import Path
from input_integrity import digest
from verify_results import rank_unit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis'


def read(p):
    with (gzip.open(p, 'rt') if str(p).endswith('.gz') else p.open()) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(name, rows):
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def unit(x):
    mean = statistics.mean(x); centered = [v-mean for v in x]
    norm = math.sqrt(math.fsum(v*v for v in centered))
    if not norm:
        raise ValueError('Constant expression vector')
    return [v/norm for v in centered]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--clinical-bioc', type=Path)
    a = p.parse_args()
    candidates_path = OUT / 'audit_gse254300_published_clinical_candidates.tsv'
    if a.clinical_bioc:
        candidates = []
        for collection in json.loads(a.clinical_bioc.read_text()):
            for doc in collection['documents']:
                for passage in doc['passages']:
                    info = passage.get('infons', {})
                    if info.get('type') != 'table' or info.get('id') not in ('jsp21318-tbl-0001', 'jsp21318-tbl-0002'):
                        continue
                    for row in ET.fromstring(info['xml']).find('tbody'):
                        v = [''.join(c.itertext()).strip() for c in row]
                        if len(v) != 9:
                            raise ValueError('Unexpected source clinical table structure')
                        candidates.append(dict(published_case_id=v[0], diagnosis='AIS' if v[0].startswith('IS') else 'CS',
                            age_years=v[1], sex=v[2], classification=v[3], major_curve_range=v[4],
                            curve_length_vertebrae=v[5], cobb_angle_degrees=v[6], apex=v[7], application=v[8],
                            source_DOI='10.1002/jsp2.1318', source_BioC_sha256=digest(a.clinical_bioc),
                            expression_alias_merge='not authenticated; source clinical record only'))
        if len(candidates) != 30:
            raise ValueError('Require both 15-patient source tables')
        write(candidates_path.name, candidates)
    else:
        candidates = read(candidates_path)
    summaries = []
    for diagnosis in ('AIS', 'CS'):
        rows = [r for r in candidates if r['diagnosis']==diagnosis and r['application']=='RNA‐seq']
        if len(rows) != 5:
            raise ValueError('Require five source-reported RNA-seq cases per diagnosis')
        ages = [float(r['age_years']) for r in rows]; angles = [float(r['cobb_angle_degrees']) for r in rows]
        summaries.append(dict(diagnosis=diagnosis, source_reported_RNAseq_cases=5,
            female_cases=sum(r['sex']=='Female' for r in rows), male_cases=sum(r['sex']=='Male' for r in rows),
            age_mean=statistics.mean(ages), age_min=min(ages), age_max=max(ages),
            cobb_mean=statistics.mean(angles), cobb_min=min(angles), cobb_max=max(angles),
            classification_counts='; '.join(k+':'+str(sum(r['classification']==k for r in rows)) for k in sorted({r['classification'] for r in rows})),
            scope='Source-publication cohort description; not verified expression-patient covariate linkage'))
    write('audit_gse254300_source_cohort_clinical_summary.tsv', summaries)
    aliases = read(ROOT/'analysis/tables/audit_gse254300_alias_clinical_link.tsv')
    cases = {r['published_case_id']:r for r in candidates}; linkage = []
    ycounts = read(OUT/'audit_gse254300_unfiltered_sex_marker_counts.tsv')
    if len(ycounts)!=5:
        raise ValueError('Require five full-count Y marker rows')
    for row in aliases:
        alias = row['analysis_patient_id']; submitted = row['GEO_concave_title'].split(' [')[0][:-1]
        for basis, candidate in [('processed_column_name', alias), ('GEO_submitted_title', submitted)]:
            c = cases.get(candidate, {})
            left = sum(float(r[alias+'A']) for r in ycounts)
            right = sum(float(r[alias+'T']) for r in ycounts)
            linkage.append(dict(analysis_patient_id=alias, candidate_basis=basis, candidate_case_id=candidate,
                candidate_exists_in_source_table=bool(c), candidate_sex=c.get('sex','unavailable'),
                candidate_age=c.get('age_years','unavailable'), candidate_application=c.get('application','unavailable'),
                concave_five_Y_marker_count_sum=left,convex_five_Y_marker_count_sum=right,
                sex_marker_caveat='Y expression conflicts with female-labelled candidate' if left>0 and right>0 and c.get('sex')=='Female'
                    else 'Sex-marker compatibility/non-detection is not authenticated identity or clinical mapping',
                linkage_status='not authenticated; same string or sex compatibility does not establish identity'))
    write('audit_gse254300_clinical_alias_candidates.tsv', linkage)

    full = read(OUT/'bulk_followup_full_gene_TPM_counts.tsv.gz')
    deposited_path = OUT/'followup_inputs/GSE311180_deposited_TPM.txt.gz'
    if digest(deposited_path) != 'b344f80f289d7f3ca9450ed4b8da80f62fb35d24fbc6957d5e2f5f5605522fd5':
        raise ValueError('Deposited source checksum mismatch')
    deposited = read(deposited_path)
    named = [r for r in deposited if r['gene_name']]
    if len(named) != len({r['gene_name'] for r in named}):
        raise ValueError('Duplicate deposited gene names')
    raw = {r['gene']:r for r in full}; dep = {r['gene_name']:r for r in named}
    genes = sorted(set(raw)&set(dep)); names = sorted(k[:-4] for k in full[0] if k.endswith('_TPM'))
    columns = [c for c in deposited[0] if c!='gene_name']
    if len(names)!=14 or len(columns)!=24 or len(genes)<10000:
        raise ValueError('Incomplete comparison universe')
    x = {n:[float(raw[g][n+'_TPM']) for g in genes] for n in names}
    y = {n:[float(dep[g][n]) for g in genes] for n in columns}
    for v in [*x.values(),*y.values()]:
        if any(not math.isfinite(z) or z<0 for z in v):
            raise ValueError('Invalid gene TPM')
    xr = {n:rank_unit(v) for n,v in x.items()}; yr = {n:rank_unit(v) for n,v in y.items()}
    xl = {n:unit([math.log1p(z) for z in v]) for n,v in x.items()}
    yl = {n:unit([math.log1p(z) for z in v]) for n,v in y.items()}
    comparisons = []
    for n in names:
        for c in columns:
            comparisons.append(dict(raw_sample=n, deposited_column=c, same_column_name=n==c,
                shared_named_genes=len(genes), spearman_TPM=math.fsum(a*b for a,b in zip(xr[n],yr[c])),
                pearson_log1p_TPM=math.fsum(a*b for a,b in zip(xl[n],yl[c])),
                scope='Descriptive different-reference expression comparison; no inferred relabelling/side identity'))
    write('audit_gse311180_all_column_expression_comparisons.tsv', comparisons)
    best = []
    for n in names:
        rows = sorted((r for r in comparisons if r['raw_sample']==n), key=lambda r:r['spearman_TPM'], reverse=True)
        own = next(r for r in rows if r['deposited_column']==n)
        best.append(dict(raw_sample=n, same_name_spearman=own['spearman_TPM'],
            highest_spearman_deposited_column=rows[0]['deposited_column'], highest_spearman=rows[0]['spearman_TPM'],
            runner_up_column=rows[1]['deposited_column'], runner_up_spearman=rows[1]['spearman_TPM'],
            relabelled=False, scope='No expression-based patient/side authentication or replacement map'))
    write('audit_gse311180_same_name_and_closest_expression.tsv', best)
    meta = read(ROOT/'metadata/audit/gse311180_requant_manifest.tsv'); by_sample={r['sample_id']:r for r in meta}
    markers = []
    for n in names:
        for gene in ('XIST','RPS4Y1','DDX3Y','KDM5D','UTY','EIF1AY','CKM','ACAN','COL2A1','ADIPOQ','ANGPT1','LEP'):
            if gene not in raw:
                raise ValueError('Missing unfiltered nominated gene: '+gene)
            markers.append(dict(sample_id=n, patient_id=by_sample[n]['patient_id'], side=by_sample[n]['side'],
                gene=gene, raw_requant_TPM=raw[gene][n+'_TPM'], raw_requant_NumReads=raw[gene][n+'_NumReads'],
                deposited_TPM=dep.get(gene,{}).get(n,'gene_not_present'),
                scope='Unfiltered expression measurement; no DNA/anatomical authentication'))
    write('audit_gse311180_unfiltered_nominated_markers.tsv', markers)
    print(json.dumps(dict(status='PASS', all_column_comparisons=len(comparisons), common_named_genes=len(genes),
        unnamed_deposited_rows=len(deposited)-len(named), clinical_candidates=len(candidates),
        labels_changed=False),indent=2))


if __name__ == '__main__':
    main()
