#!/usr/bin/env python3
"""Independent arithmetic/rank checks for the available-data follow-up tables."""
import csv,gzip,json,math,statistics
from pathlib import Path
from project_paths import resolve_analysis_input
from verify_results import close,t_critical,t_two_sided,bh
ROOT=Path(__file__).resolve().parents[1]

def rows(name):
    return source_rows('analysis/'+name)

def source_rows(name):
    p=resolve_analysis_input(name, ROOT)
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else p.open()) as f:
        r=csv.DictReader(f,delimiter='\t')
        if not r.fieldnames or len(set(r.fieldnames))!=len(r.fieldnames):raise ValueError('Duplicate/missing header '+name)
        out=list(r)
    if not out or any(None in z or any(v is None for v in z.values()) for z in out):raise ValueError('Incomplete table '+name)
    return out

def unique_index(records,fields,context):
    out={}
    for r in records:
        key=tuple(r[k] for k in fields)
        if key in out:raise ValueError('Duplicate '+context+' identity '+str(key))
        out[key]=r
    return out

def verify_proxy_input_links(effects):
    """Reconstruct each exported patient contrast from the actual expression matrix."""
    observed=unique_index(effects,['cohort','patient_id'],'proxy patient')
    ma=source_rows('metadata/audit/bulk_samples.tsv')
    mb=source_rows('metadata/audit/gse311180_requant_manifest.tsv')
    expected={}
    recipes=[('GSE254300_AIS','results/tables/bulk/GSE254300_AIS_TMM_logCPM.tsv.gz',
              [r for r in ma if r['cohort']=='AIS'],'processed_sample_id',5),
             ('GSE311180_AIS','results/tables/bulk_gse311180_requant/GSE311180_TMM_logCPM.tsv.gz',mb,'sample_id',7),
             ('GSE254300_CS','results/tables/disease_control_gse254300_cs/GSE254300_AIS_CS_joint_TMM_logCPM.tsv.gz',
              [r for r in ma if r['cohort']=='congenital_scoliosis'],'processed_sample_id',5)]
    for cohort,path,metadata,sample_key,n in recipes:
        pairs=unique_index(metadata,['patient_id','side'],'source sample')
        patients={key[0] for key in pairs}
        if len(patients)!=n or len(pairs)!=2*n:raise ValueError('Incomplete proxy source cohort')
        with gzip.open(ROOT/path,'rt') as f:
            rd=csv.DictReader(f,delimiter='\t');gene_key=rd.fieldnames[0]
            matrix={r[gene_key]:r for r in rd if r[gene_key] in {'ADIPOQ','ANGPT1','LEP'}}
        if set(matrix)!={'ADIPOQ','ANGPT1','LEP'}:raise ValueError('Missing proxy source gene')
        for pid in patients:
            a=pairs[(pid,'concave')][sample_key];b=pairs[(pid,'convex')][sample_key]
            expected[(cohort,pid)]={gene+'_difference':float(r[a])-float(r[b]) for gene,r in matrix.items()}
    if set(observed)!=set(expected):raise ValueError('Proxy patient/cohort identity mismatch')
    for key,values in expected.items():
        for field,value in values.items():close(observed[key][field],value,'source-linked proxy '+str(key)+' '+field)

def verify_clinical_input_links(patients):
    """Check covariates and scores against their separate source tables."""
    source=unique_index([r for r in source_rows('analysis/tables/audit_patient_clinical_metadata_19_pairs.tsv')
                         if r['dataset']=='GSE311180'],['analysis_patient_id'],'clinical source')
    score=unique_index(source_rows('results/tables/bulk_gse311180_requant/frozen_axis_paired_effects.tsv'),['patient_id'],'score source')
    observed=unique_index(patients,['analysis_patient_id'],'clinical patient')
    if len(source)!=7 or set(observed)!=set(source) or set(score)!=set(source):raise ValueError('Clinical/score identity mismatch')
    for key,r in observed.items():
        for field in ['age_years','cobb_angle_degrees']:close(r[field],source[key][field],'source-linked '+field)
        close(r['score_difference'],score[key]['paired_difference'],'source-linked score')

def verify_composition_input_links(comp):
    original=source_rows('results/tables/bulk/deconvolution_two_methods.tsv')
    source=unique_index(original,['method','patient_id','cell_type','side'],'composition source')
    expected={}
    for method,pid,cell_type,side in source:
        key=(method,pid)
        if cell_type not in expected.setdefault(key,{}):
            expected[key][cell_type]=float(source[(method,pid,cell_type,'concave')]['proportion'])-float(source[(method,pid,cell_type,'convex')]['proportion'])
    observed=unique_index(comp,['method','patient_id'],'composition patient')
    if len(expected)!=10 or set(observed)!=set(expected):raise ValueError('Composition patient/method mismatch')
    for key,values in expected.items():
        if set(observed[key])-{'method','patient_id'}!=set(values):raise ValueError('Composition lineage identity mismatch')
        for field,value in values.items():close(observed[key][field],value,'source-linked composition '+str(key)+' '+field)

def eliminate(a,b=None,tol=1e-10):
    m=[list(map(float,r))+([float(v)] if b is not None else []) for r,v in zip(a,b if b is not None else [0]*len(a))]
    columns=len(a[0]);rank=0;pivots=[]
    for col in range(columns):
        if rank==len(m):break
        pivot=max(range(rank,len(m)),key=lambda i:abs(m[i][col]))
        if abs(m[pivot][col])<=tol:continue
        m[rank],m[pivot]=m[pivot],m[rank];v=m[rank][col];m[rank]=[x/v for x in m[rank]]
        for i in range(len(m)):
            if i!=rank:
                v=m[i][col];m[i]=[x-v*y for x,y in zip(m[i],m[rank])]
        pivots.append(col);rank+=1
    if b is None:return rank
    if rank!=columns or len(a)!=columns:raise ValueError('Singular linear system')
    return [m[i][-1] for i in range(columns)]

def ols(x,y):
    p=len(x[0]);n=len(y)
    xtx=[[math.fsum(r[i]*r[j] for r in x) for j in range(p)] for i in range(p)]
    beta=eliminate(xtx,[math.fsum(r[i]*v for r,v in zip(x,y)) for i in range(p)])
    inverse=[eliminate(xtx,[float(i==j) for i in range(p)]) for j in range(p)]
    rss=math.fsum((v-math.fsum(a*b for a,b in zip(r,beta)))**2 for r,v in zip(x,y))
    df=n-p;se=[math.sqrt(inverse[i][i]*rss/df) for i in range(p)]
    return beta,se,df

def main():
    candidates=rows('audit_gse254300_published_clinical_candidates.tsv')
    if len(candidates)!=30 or len({r['published_case_id'] for r in candidates})!=30:raise ValueError('Source cases')
    for s in rows('audit_gse254300_source_cohort_clinical_summary.tsv'):
        rr=[r for r in candidates if r['diagnosis']==s['diagnosis'] and r['application']=='RNA‐seq']
        if len(rr)!=5:raise ValueError('RNAseq case selection')
        for field,prefix in [('age_years','age'),('cobb_angle_degrees','cobb')]:
            v=[float(r[field]) for r in rr]
            for statistic,value in [('mean',statistics.mean(v)),('min',min(v)),('max',max(v))]:close(s[prefix+'_'+statistic],value,'source '+prefix)
        close(s['female_cases'],sum(r['sex']=='Female' for r in rr),'source sex')
    effects=rows('bulk_followup_ADIPOQ_proxy_patient_effects.tsv');models=rows('bulk_followup_ADIPOQ_proxy_models.tsv')
    verify_proxy_input_links(effects)
    if len(effects)!=17 or len(models)!=6 or len({r['cohort'] for r in effects})!=3:raise ValueError('Missing diagnosis/cohort')
    if set(unique_index(models,['cohort','target'],'proxy model'))!={(c,g) for c in {r['cohort'] for r in effects} for g in ['ANGPT1','LEP']}:raise ValueError('Missing/duplicated proxy model')
    probabilities=[]
    for r in models:
        group=[v for v in effects if v['cohort']==r['cohort']]
        x=[[1.,float(v['ADIPOQ_difference'])] for v in group];y=[float(v[r['target']+'_difference']) for v in group]
        beta,se,df=ols(x,y);p=t_two_sided(beta[0]/se[0],df);probabilities.append(p)
        xv=[z[1] for z in x];xm=statistics.mean(xv);ym=statistics.mean(y)
        correlation=math.fsum((a-xm)*(b-ym) for a,b in zip(xv,y))/math.sqrt(math.fsum((a-xm)**2 for a in xv)*math.fsum((b-ym)**2 for b in y))
        for field,value in [('n_pairs',len(group)),('unadjusted_mean_logCPM_difference',ym),('intercept_at_zero_ADIPOQ_difference',beta[0]),('SE',se[0]),('residual_df',df),('ADIPOQ_slope',beta[1]),('pearson_paired_difference_correlation',correlation),('CI95_lower',beta[0]-t_critical(df)*se[0]),('CI95_upper',beta[0]+t_critical(df)*se[0]),('nominal_t_p',p)]:close(r[field],value,'proxy '+field)
    for r,value in zip(models,bh(probabilities)):close(r['BH_six_proxy_intercept_tests'],value,'six-test BH')
    patients=rows('bulk_followup_available_clinical_patient_effects.tsv')
    verify_clinical_input_links(patients)
    gates=rows('audit_gse311180_rna_expected_pairs.tsv')
    passing={r['patient_id'] for r in gates if r['interpretation']=='consistent_RNA_profiles'}
    if len(passing)!=5 or {r['patient_id'] for r in gates}-{r['analysis_patient_id'] for r in patients}:raise ValueError('RNA-gate patient source mismatch')
    modifier_rows=rows('bulk_followup_available_clinical_modifiers.tsv')
    if set(unique_index(modifier_rows,['subset','term'],'modifier model'))!={(s,t) for s in ['all_seven_metadata_pairs','five_RNA_gate_passing_pairs'] for t in ['(Intercept)','age_centered','cobb_centered']}:raise ValueError('Missing/duplicated modifier model')
    for label in ['all_seven_metadata_pairs','five_RNA_gate_passing_pairs']:
        group=[r for r in patients if label.startswith('all') or r['analysis_patient_id'] in passing]
        age=[float(r['age_years']) for r in group];cobb=[float(r['cobb_angle_degrees']) for r in group]
        x=[[1.,a-statistics.mean(age),c-statistics.mean(cobb)] for a,c in zip(age,cobb)]
        beta,se,df=ols(x,[float(r['score_difference']) for r in group])
        for r in [v for v in modifier_rows if v['subset']==label]:
            i=['(Intercept)','age_centered','cobb_centered'].index(r['term'])
            for field,value in [('n_pairs',len(group)),('estimate',beta[i]),('SE',se[i]),('residual_df',df),('CI95_lower',beta[i]-t_critical(df)*se[i]),('CI95_upper',beta[i]+t_critical(df)*se[i]),('nominal_t_p',t_two_sided(beta[i]/se[i],df))]:close(r[field],value,'clinical '+field)
    comp=rows('bulk_followup_composition_patient_differences.tsv')
    verify_composition_input_links(comp)
    slopes=[r for r in modifier_rows if r['term']!='(Intercept)']
    if len(slopes)!=4:raise ValueError('Clinical modifier family')
    for r,value in zip(slopes,bh([float(z['nominal_t_p']) for z in slopes])):close(r['BH_four_modifier_tests'],value,'modifier BH')
    for r in rows('bulk_followup_full_composition_rank.tsv'):
        rr=[z for z in comp if z['method']==r['method']];cols=[k for k in rr[0] if k not in {'method','patient_id'}]
        cols=[k for k in cols if statistics.stdev(float(z[k]) for z in rr)>1e-8]
        nuisance_rank=eliminate([[float(z[k]) for k in cols] for z in rr],tol=1e-7)
        rank=eliminate([[1.]+[float(z[k]) for k in cols] for z in rr],tol=1e-7)
        close(r['augmented_full_composition_design_rank'],rank,'composition rank');close(r['residual_df'],len(rr)-rank,'rank df')
        close(r['variable_paired_components'],len(cols),'variable components');close(r['nuisance_design_rank'],nuisance_rank,'nuisance rank')
        if r['conditional_intercept_identifiable']!=str(rank>nuisance_rank).upper():raise ValueError('Incorrect intercept identifiability')
    print(json.dumps(dict(status='PASS',clinical_source_cases=30,source_cohorts=2,proxy_models=6,
        source_linked_proxy_patients=17,source_linked_clinical_patients=7,source_linked_composition_rows=10,
        modifier_models=2,composition_rank_checks=2,scope='Independent expression/metadata/score source joins and arithmetic/OLS/rank; no clinical/DNA authentication'),indent=2))

if __name__=='__main__':main()
