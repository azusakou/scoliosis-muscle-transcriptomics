#!/usr/bin/env python3
"""Descriptive fragment-coordinate diversity at original RNA het/hom discordances.

All seven expected pairs / all fourteen libraries are retained as controls.
MAPQ60/baseQ20 primary proper pairs only; overlap mates are collapsed by read
pair name. Coordinates are alignment-template proxies, not UMI molecules or
calibrated PCR duplicates. Per-read names/coordinates/counts stay private.
"""
import argparse,bisect,importlib,json,os,re,signal,statistics,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from collections import defaultdict
from input_integrity import digest,validate_pairs
from project_paths import executable
m=importlib.import_module('14_rna_pairing_failure_diagnostics')
ROOT=Path(__file__).resolve().parents[1]

def aligned_bases(pos,cigar,seq,qual,targets):
    tokens=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if ''.join(n+op for n,op in tokens)!=cigar or not tokens:raise ValueError('Malformed CIGAR')
    if len(seq)!=len(qual):raise ValueError('Sequence/quality length mismatch')
    ref=pos;query=0
    for value,op in tokens:
        n=int(value)
        if n<=0:raise ValueError('Zero-length CIGAR operation')
        if op in 'M=X':
            for target in targets[bisect.bisect_left(targets,ref):bisect.bisect_left(targets,ref+n)]:
                i=query+target-ref
                if ord(qual[i])-33>=20:yield target,seq[i].upper()
            ref+=n;query+=n
        elif op in 'IS':query+=n
        elif op in 'DN':ref+=n
    if query!=len(seq):raise ValueError('CIGAR/query length mismatch')

def summarize_templates(records,alleles):
    valid=[(coord,next(iter(bases))) for (qname,coord),bases in records.items() if len(bases)==1 and next(iter(bases)) in alleles]
    coords={coord for coord,base in valid}
    support={base:{coord for coord,b in valid if b==base} for base in alleles}
    return dict(qname_templates=len(valid),distinct_coordinates=len(coords),allele_coordinates=support,
        conflicting_or_other_templates=len(records)-len(valid))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['identity','sites','output']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--manifest',type=Path,default=ROOT/'metadata/audit/gse311180_requant_manifest.tsv')
    p.add_argument('--jobs',type=int,default=2)
    a=p.parse_args();assert a.jobs>0
    if a.output.exists():raise ValueError('Preserve prior diagnostics; use a new private output directory')
    a.output.mkdir();meta=m.rows(a.manifest);validate_pairs(meta,expected_pairs=7)
    patient={r['sample_id']:r['patient_id'] for r in meta};names=sorted(patient)
    binding=json.loads((a.identity/'reference_artifact_binding.json').read_text())['binding']['references']
    assert digest(a.sites)==binding['compressed']['sites.hg38.rna.vcf.gz']
    markers=m.sites(a.sites);assert len(markers)==16469
    calls={};provenance={}
    for n in names:
        provenance[n]=json.loads((a.identity/f'{n}.provenance.json').read_text())
        file=a.identity/'fingerprints'/f'{n}.somalier';name,size,cs=m.decode(file)
        assert name==n and size[0]==16469 and digest(file)==provenance[n]['outputs']['fingerprint_sha256']
        calls[n]=[m.gt(c) for c in cs]
    pairs=[sorted(r['sample_id'] for r in meta if r['patient_id']==pid) for pid in sorted(set(patient.values()))]
    target={};directions=[]
    for x,y in pairs:
        idx=[i for i,(u,v) in enumerate(zip(calls[x],calls[y])) if (u==1 and v in (0,2)) or (v==1 and u in (0,2))]
        target[x]=target[y]=idx
        for hom,het in [(x,y),(y,x)]:
            directions.append((patient[x],hom,het,[i for i in idx if calls[hom][i] in (0,2)]))
    sam=executable('samtools','AIS_SAMTOOLS');assert subprocess.check_output([sam,'--version'],text=True).splitlines()[0]=='samtools 1.24'
    def sample(n):
        idx=target[n];lookup={(markers[i][0],markers[i][1]):i for i in idx};positions=defaultdict(list)
        for c,pos in lookup:positions[c].append(pos)
        for values in positions.values():values.sort()
        bed=a.output/f'{n}.discordant_sites.bed';bed.write_text(''.join(f'{c}\t{pos-1}\t{pos}\n' for c,pos in sorted(lookup)))
        bam=a.identity/n/'Aligned.sortedByCoord.out.bam';assert digest(bam)==provenance[n]['outputs']['bam_sha256']
        templates={i:defaultdict(set) for i in idx};read_counts=defaultdict(int)
        with (a.output/f'{n}.view.log').open('w') as log:
            process=subprocess.Popen([sam,'view','-M','-L',str(bed),'-q','60','-f','2','-F','3844',str(bam)],stdout=subprocess.PIPE,stderr=log,text=True,start_new_session=True)
            try:
                for line in process.stdout:
                    row=line.rstrip('\n').split('\t');assert len(row)>=11
                    name,flag,chrom,pos,mate,tlen=row[0],int(row[1]),row[2],int(row[3]),int(row[7]),int(row[8])
                    if row[6] not in ('=',chrom) or mate<=0 or tlen==0 or row[9]=='*' or row[10]=='*':continue
                    coord=(chrom,min(pos,mate),abs(tlen))
                    for locus,base in aligned_bases(pos,row[5],row[9],row[10],positions[chrom]):
                        i=lookup[chrom,locus];templates[i][name,coord].add(base);read_counts[i]+=1
                assert process.wait()==0,'Native coordinate-view stage failed'
            except BaseException:
                try:os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                process.wait();raise
        data={i:dict(summarize_templates(templates[i],markers[i][2:]),read_bases=read_counts[i]) for i in idx}
        print(n,'DISCORDANT_FRAGMENT_COORDINATE_AUDIT_PASS',flush=True)
        return n,data
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:data=dict(pool.map(sample,names))
    out=[]
    for pid,hom,het,idx in directions:
        hs=[data[hom][i] for i in idx];ts=[data[het][i] for i in idx]
        minority=[len(data[hom][i]['allele_coordinates'][markers[i][3] if calls[hom][i]==0 else markers[i][2]]) for i in idx]
        med=lambda rs,k:statistics.median(r[k] for r in rs) if rs else ''
        out.append(dict(patient_id=pid,homozygous_call_sample=hom,heterozygous_call_sample=het,original_discordant_loci=len(idx),
            homo_loci_with_HQ_proper_templates=sum(r['qname_templates']>0 for r in hs),het_loci_with_HQ_proper_templates=sum(r['qname_templates']>0 for r in ts),
            homo_HQ_read_bases_median=med(hs,'read_bases'),het_HQ_read_bases_median=med(ts,'read_bases'),
            homo_qname_templates_median=med(hs,'qname_templates'),het_qname_templates_median=med(ts,'qname_templates'),
            homo_distinct_coordinates_median=med(hs,'distinct_coordinates'),het_distinct_coordinates_median=med(ts,'distinct_coordinates'),
            homo_loci_with_1_or_2_coordinates=sum(0<r['distinct_coordinates']<=2 for r in hs),
            homo_loci_with_7plus_coordinates=sum(r['distinct_coordinates']>=7 for r in hs),
            homo_loci_with_other_allele_coordinates=sum(v>0 for v in minority),
            het_loci_with_both_alleles_coordinates=sum(all(data[het][i]['allele_coordinates'][b] for b in markers[i][2:]) for i in idx),
            homo_total_qname_templates=sum(r['qname_templates'] for r in hs),homo_total_distinct_coordinates=sum(r['distinct_coordinates'] for r in hs),
            scope='Coordinate diversity proxy; not molecular/PCR duplicate count, genotype reclassification or DNA identity'))
    assert len(out)==14;m.write(a.output/'audit_RNA_discordant_fragment_coordinate_diversity.tsv',out)
    proof={'status':'PASS','samples':14,'expected_pairs':7,'directional_summaries':14,'original_discordant_locus_union_by_pair':sum(len(target[x]) for x,y in pairs),
        'MAPQ':60,'individual_baseQ':20,'primary_proper_pair_filter':True,'mate_overlap':'QNAME/template-coordinate collapse; ambiguous alleles excluded',
        'coordinate_key':'chromosome,min(POS,PNEXT),abs(TLEN)','molecular_or_PCR_duplicate_calibrated':False,
        'script_sha256':digest(Path(__file__)),'decoder_sha256':digest(Path(m.__file__)),'scope':'All original het/hom-discordant loci, all seven pairs; aggregate proxy, not primary calling replacement'}
    (a.output/'independent_fragment_coordinate_QA.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof,indent=2))
if __name__=='__main__':main()
