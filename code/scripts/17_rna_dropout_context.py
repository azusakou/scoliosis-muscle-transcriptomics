#!/usr/bin/env python3
"""Aggregate directional, depth and chromosome context of RNA-call discordance.

This describes observed RNA allele dropout, not its biological cause or DNA
identity. Private marker counts/sketches are never exported by this command.
"""
from __future__ import annotations
import argparse, gzip, importlib, itertools, statistics
from pathlib import Path
m=importlib.import_module('14_rna_pairing_failure_diagnostics')
ROOT=Path(__file__).resolve().parents[1]

def summarize(patient,a,b,ca,cb,markers,reference_alleles=None):
    ga=[m.gt(c) for c in ca];gb=[m.gt(c) for c in cb]
    out=[]
    for hom_name,het_name,hom,het,gh,gt in [(a,b,ca,cb,ga,gb),(b,a,cb,ca,gb,ga)]:
        idx=[i for i,(x,y) in enumerate(zip(gh,gt)) if x in (0,2) and y==1]
        hd=[sum(hom[i][:2]) for i in idx];td=[sum(het[i][:2]) for i in idx]
        out.append(dict(patient_id=patient,homozygous_call_sample=hom_name,heterozygous_call_sample=het_name,
            discordant_sites=len(idx),zero_minor_allele_sites=sum(min(hom[i][:2])==0 for i in idx),
            hom_depth_7_14=sum(7<=d<15 for d in hd),hom_depth_15_29=sum(15<=d<30 for d in hd),
            hom_depth_30_99=sum(30<=d<100 for d in hd),hom_depth_100_plus=sum(d>=100 for d in hd),
            hom_depth_median=statistics.median(hd) if hd else 0,
            opposite_het_depth_median=statistics.median(td) if td else 0,
            chromosomes_with_discordance=len({markers[i][0] for i in idx}),
            hom_call_side_zero_minor_fraction=sum(min(hom[i][:2])==0 for i in idx)/len(idx) if idx else 0))
    if reference_alleles is not None:
        for r,hom,het,gh,gt in zip(out,[ca,cb],[cb,ca],[ga,gb],[gb,ga]):
            idx=[i for i,(x,y) in enumerate(zip(gh,gt)) if x in (0,2) and y==1]
            r['hom_reference_allele_sites']=sum(markers[i][2 if gh[i]==0 else 3]==reference_alleles[i] for i in idx)
            r['hom_alternate_allele_sites']=len(idx)-r['hom_reference_allele_sites']
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['identity','sites','output']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--manifest',type=Path,default=ROOT/'metadata/audit/gse311180_requant_manifest.tsv')
    a=p.parse_args();meta=m.rows(a.manifest);assert len(meta)==14
    patient={r['sample_id']:r['patient_id'] for r in meta};markers=m.sites(a.sites);counts={}
    for n in patient:
        sample,size,cs=m.decode(a.identity/'fingerprints'/f'{n}.somalier')
        assert sample==n and len(cs)==size[0]==len(markers)==16469
        counts[n]=cs
    refs={}
    with gzip.open(a.sites,'rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            fields=line.split('\t');refs[(fields[0],int(fields[1]))]=fields[3]
    reference_alleles=[refs[(c,pos)] for c,pos,x,y in markers]
    output=[]
    for x,y in itertools.combinations(sorted(patient),2):
        if patient[x]==patient[y]:output.extend(summarize(patient[x],x,y,counts[x],counts[y],markers,reference_alleles))
    assert len(output)==14
    m.write(a.output/'audit_RNA_pairing_directional_dropout.tsv',output)
    print('ALL_SEVEN_PAIRS_DIRECTIONAL_RNA_DROPOUT_CONTEXT_COMPLETE')
if __name__=='__main__':main()
