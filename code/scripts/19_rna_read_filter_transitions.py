#!/usr/bin/env python3
"""Aggregate original-marker transitions under the strict RNA read pileup.

All per-marker genotypes stay private. All seven expected pairs are controls;
changes describe calling sensitivity, not DNA identity or molecular causation.
"""
import argparse,importlib,itertools
from pathlib import Path
m=importlib.import_module('14_rna_pairing_failure_diagnostics')
ROOT=Path(__file__).resolve().parents[1]

def category(x,y):
    if x<0 or y<0:return 'uncalled'
    if x==y:return 'concordant'
    if x in (0,2) and y in (0,2):return 'opposite_homozygotes'
    return 'het_hom_discordance'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['primary','strict','output']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--manifest',type=Path,default=ROOT/'metadata/audit/gse311180_requant_manifest.tsv')
    p.add_argument('--overlap-retained',action='store_true')
    a=p.parse_args();meta=m.rows(a.manifest);patient={r['sample_id']:r['patient_id'] for r in meta}
    old={};new={}
    for n in patient:
        for output,folder in [(old,a.primary/'fingerprints'),(new,a.strict/'private_fingerprints')]:
            name,size,cs=m.decode(folder/f'{n}.somalier')
            assert name==n and len(cs)==size[0]==16469
            output[n]=[m.gt(c) for c in cs]
    rows=[]
    for x,y in itertools.combinations(sorted(patient),2):
        if patient[x]!=patient[y]:continue
        before=[category(a,b) for a,b in zip(old[x],old[y])]
        after=[category(a,b) for a,b in zip(new[x],new[y])]
        for start in ['uncalled','concordant','opposite_homozygotes','het_hom_discordance']:
            for end in ['uncalled','concordant','opposite_homozygotes','het_hom_discordance']:
                rows.append(dict(patient_id=patient[x],sample_a=x,sample_b=y,primary_category=start,
                    strict_category=end,sites=sum(a==start and b==end for a,b in zip(before,after))))
    assert len(rows)==112;m.write(a.output/('audit_RNA_read_filter_call_transitions'+('_overlap_retained' if a.overlap_retained else '')+'.tsv'),rows)
    print('ALL_SEVEN_PAIRS_ORIGINAL_MARKER_TRANSITIONS_COMPLETE')
if __name__=='__main__':main()
