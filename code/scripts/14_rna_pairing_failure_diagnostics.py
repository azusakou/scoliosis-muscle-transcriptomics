#!/usr/bin/env python3
"""Post hoc RNA pairing diagnostics; primary consistency gates are never replaced.

Decode Somalier v0.3.5 format2 on the compute host. Publish aggregate diagnostics
only, never the sketches, genotypes or per-site allele counts. Verify all default
integer pair statistics before interpreting sensitivities or loss decomposition.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,itertools,json,math,statistics,struct
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def rows(path):
    with path.open(newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        reader.fieldnames=[n.lstrip('#$ ') for n in reader.fieldnames]
        return list(reader)


def write(path,data):
    if not data: raise ValueError('No aggregate rows: '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)


def decode(path):
    raw=path.read_bytes()
    if len(raw)<8 or raw[0]!=2:raise ValueError('Require Somalier format2')
    length=raw[1];sample=raw[2:2+length].decode('utf8');offset=2+length
    sizes=struct.unpack_from('<HHH',raw,offset);offset+=6
    if len(raw)!=offset+12*sum(sizes):raise ValueError('Truncated or trailing sketch bytes')
    counts=list(struct.iter_unpack('<III',raw[offset:]))
    return sample,sizes,counts[:sizes[0]]


def sites(path):
    result=[]
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            r=line.rstrip().split('\t')
            if r[0] in {'X','chrX','Y','chrY','NC_000023.10','NC_000023.11','NC_000024.9','NC_000024.10'}:continue
            if len(r[3])!=1 or len(r[4])!=1:raise ValueError('Non-SNP marker')
            result.append((r[0],int(r[1]),*sorted((r[3],r[4]))))
    return sorted(result,key=lambda r:(r[0],r[1]))


def other_fraction(count):
    # The version0.3.5 binary rounds this ratio to float32 before comparing
    # with float64 0.1. Validate this against all native default/grid counts.
    if not count[2]:return 0.0
    return struct.unpack("<f",struct.pack("<f",count[2]/sum(count)))[0]


def gt(count,depth=7,min_ab=.2):
    a,b,other=count;n=a+b
    if n<depth or other_fraction(count)>.1:return -1
    ab=b/n
    if ab<.01:return 0
    if ab>.99:return 2
    if min_ab<=ab<=1-min_ab:return 1
    return -1


def compare(a,b,indices=None):
    if indices is None:indices=range(len(a))
    cats=Counter((a[i],b[i]) for i in indices if a[i]>=0 and b[i]>=0)
    n=sum(cats.values());h=cats[1,1]
    ibs0=cats[0,2]+cats[2,0]
    hetab=sum(n*((x==1)+(y==1)) for (x,y),n in cats.items())
    discord=cats[0,1]+cats[1,0]+cats[1,2]+cats[2,1]
    return dict(n=n,ibs0=ibs0,ibs2=sum(n for (x,y),n in cats.items() if x==y),shared_hets=h,
                hets_ab=hetab,shared_hom_alts=cats[2,2],
                heterozygote_homozygote_discordance=discord,
                relatedness=2*(h-2*ibs0)/max(1,hetab),
                loss_from_het_hom=discord/max(1,hetab),loss_from_IBS0=4*ibs0/max(1,hetab))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--identity',type=Path,required=True);ap.add_argument('--sites',type=Path,required=True)
    ap.add_argument('--manifest',type=Path,default=ROOT/'metadata/audit/gse311180_requant_manifest.tsv')
    ap.add_argument('--default-pairs',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    manifest=rows(args.manifest);assert len(manifest)==14
    names=sorted(r['sample_id'] for r in manifest);patient={r['sample_id']:r['patient_id'] for r in manifest}
    marker=sites(args.sites);counts={};hashes={}
    for name in names:
        p=args.identity/'fingerprints'/f'{name}.somalier';sample,size,count=decode(p)
        assert sample==name and size[0]==len(marker)==16469
        counts[name]=count;hashes[name]=hashlib.sha256(p.read_bytes()).hexdigest()
    baseline={n:[gt(c) for c in counts[n]] for n in names}
    original=rows(args.default_pairs);assert len(original)==91
    for r in original:
        value=compare(baseline[r['sample_a']],baseline[r['sample_b']])
        for k in ['n','ibs0','ibs2','shared_hets','hets_ab','shared_hom_alts']:
            assert value[k]==int(r[k]),(r['sample_a'],r['sample_b'],k,value[k],r[k])
        assert abs(value['relatedness']-float(r['relatedness']))<.0005001
    expected=[(a,b) for a,b in itertools.combinations(names,2) if patient[a]==patient[b]]
    assert len(expected)==7
    loss=[];coverage=[];chrom=[];sampleq=[];sensitivity=[]
    for name in names:
        cs=counts[name];g=baseline[name];dp=[c[0]+c[1] for c in cs];positive=[x for x in dp if x]
        total=sum(dp);sorted_dp=sorted(dp,reverse=True);top=max(1,math.ceil(.01*len(dp)))
        sampleq.append(dict(sample_id=name,patient_id=patient[name],all_markers=len(dp),covered_markers=len(positive),
            callable_sites=sum(x>=0 for x in g),heterozygous_sites=g.count(1),zero_depth_sites=dp.count(0),
            median_all_marker_depth=statistics.median(dp),median_positive_marker_depth=statistics.median(positive),
            callable_fraction=sum(x>=0 for x in g)/len(g),top1percent_markers_fraction_of_depth=sum(sorted_dp[:top])/total,
            rejected_low_depth=sum(x<7 for x in dp),
            rejected_other_base=sum(c[0]+c[1]>=7 and other_fraction(c)>.1 for c in cs),
            rejected_intermediate_AB=sum(d>=7 and other_fraction(c)<=.1 and .01<=c[1]/d<=.99 and not .2<=c[1]/d<=.8 for c,d in zip(cs,dp))))
    for a,b in expected:
        v=compare(baseline[a],baseline[b]);assert abs((1-v['relatedness'])-v['loss_from_het_hom']-v['loss_from_IBS0'])<1e-12
        loss.append(dict(patient_id=patient[a],sample_a=a,sample_b=b,**v,
            fraction_of_relatedness_loss_from_het_hom=v['loss_from_het_hom']/max(1e-15,1-v['relatedness']),
            het_a_hom_b=sum(x==1 and y in (0,2) for x,y in zip(baseline[a],baseline[b])),
            hom_a_het_b=sum(y==1 and x in (0,2) for x,y in zip(baseline[a],baseline[b])),
            opposite_homozygotes_fraction=v['ibs0']/v['n']))
        for lo,hi,label in [(7,15,'7-14'),(15,30,'15-29'),(30,100,'30-99'),(100,10**15,'100+')]:
            idx=[i for i in range(len(marker)) if lo<=min(sum(counts[a][i][:2]),sum(counts[b][i][:2]))<hi]
            q=compare(baseline[a],baseline[b],idx)
            coverage.append(dict(patient_id=patient[a],stratum=label,**q))
        for c in sorted({r[0] for r in marker}):
            idx=[i for i,m in enumerate(marker) if m[0]==c]
            chrom.append(dict(patient_id=patient[a],chromosome=c,**compare(baseline[a],baseline[b],idx)))
    for depth in [3,7,15,30]:
        for ab in [.1,.2,.3]:
            calls={n:[gt(c,depth,ab) for c in counts[n]] for n in names}
            for a,b in itertools.combinations(names,2):
                sensitivity.append(dict(depth=depth,min_het_AB=ab,sample_a=a,sample_b=b,
                    expected_metadata_pair=patient[a]==patient[b],patient_a=patient[a],patient_b=patient[b],
                    **compare(calls[a],calls[b])))
    for name,data in [('audit_RNA_pairing_loss_decomposition.tsv',loss),('audit_RNA_fingerprint_site_coverage.tsv',sampleq),
        ('audit_RNA_pairing_depth_strata.tsv',coverage),('audit_RNA_pairing_chromosome_strata.tsv',chrom),
        ('audit_RNA_pairing_depth_AB_grid.tsv',sensitivity)]:write(args.output/name,data)
    summary={'status':'PASS','scope':'Independent default integer reproduction and post hoc aggregate RNA diagnostics; not DNA identity or calibrated contamination',
        'default_native_comparisons_reproduced':91,'samples':14,'autosomal_markers':len(marker),
        'primary_gates_changed':False,'sketch_sha256':hashes,'sites_sha256':hashlib.sha256(args.sites.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'native_site_count_order':'v0.3.5 common.nim: lexical chromosome then position; alphabetical allele counts',
        'source_format_reference':'https://github.com/brentp/somalier/blob/v0.3.5/src/somalierpkg/common.nim',
        'genotype_reference':'https://github.com/brentp/somalier/blob/v0.3.5/src/somalierpkg/relate.nim'}
    (args.output/'independent_baseline_allele_count_reproduction.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
