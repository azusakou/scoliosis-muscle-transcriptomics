#!/usr/bin/env python3
"""Post hoc MAPQ60/baseQ20 indexed RNA marker pileups, all libraries as controls.

Only aggregate pair tables are exported. Per-marker counts and constructed
sketches stay in the specified private compute output. This does not change
primary calls or provide DNA-level/clinical/anatomical identity confirmation.
"""
from __future__ import annotations
import argparse,csv,gzip,importlib,json,os,struct,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from input_integrity import digest,validate_pairs
from project_paths import executable
m=importlib.import_module('14_rna_pairing_failure_diagnostics')


def parse_bases(value,ref):
    result={x:0 for x in 'ACGTN'};i=0
    while i<len(value):
        c=value[i]
        if c=='^':
            if i+1>=len(value):raise ValueError('Truncated mapping-quality token')
            i+=2;continue
        if c=='$':i+=1;continue
        if c in '+-':
            i+=1;j=i
            while i<len(value) and value[i].isdigit():i+=1
            if i==j:raise ValueError('Missing pileup indel length')
            n=int(value[j:i])
            if n<=0:raise ValueError('Invalid pileup indel length')
            i+=n
            if i>len(value):raise ValueError('Truncated pileup indel')
            continue
        if c in '.,':result[ref]+=1
        elif c.upper() in result:result[c.upper()]+=1
        elif c not in '*#<>':raise ValueError('Unknown pileup token: '+c)
        i+=1
    return result


def pileup_arguments(sam,bed,genome,overlap_mode):
    if overlap_mode not in {'remove','retain'}:raise ValueError('Unknown paired-overlap mode')
    return [sam,'mpileup','-B',*(['-x'] if overlap_mode=='retain' else []),'-A','-q','60','-Q','20','-d','0','-l',str(bed),'-f',str(genome),'-']

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['identity','reference','output']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--manifest',type=Path,default=Path(__file__).resolve().parents[1]/'metadata/audit/gse311180_requant_manifest.tsv')
    p.add_argument('--paired-overlap',choices=['remove','retain'],default='remove')
    p.add_argument('--jobs',type=int,default=2);p.add_argument('--view-threads',type=int,default=2)
    a=p.parse_args()
    assert a.jobs>0 and a.view_threads>0
    if a.output.exists() and any(a.output.iterdir()):raise ValueError('Preserve existing private pileup output; choose a new directory')
    a.output.mkdir(parents=True,exist_ok=True)
    sites=a.reference/'sites.hg38.rna.vcf.gz';marker=m.sites(sites);assert len(marker)==16469
    bed=a.output/'autosomal_marker_intervals.bed'
    bed.write_text(''.join(f'{c}\t{pos-1}\t{pos}\n' for c,pos,x,y in marker))
    auto_sites=a.output/'autosomal_marker_sites.vcf.gz'
    with gzip.open(sites,'rt') as src,gzip.open(auto_sites,'wt') as out:
        for line in src:
            if line.startswith('#') or line.split('\t',1)[0] not in {'chrX','chrY','X','Y'}:out.write(line)
    manifest=m.rows(a.manifest);validate_pairs(manifest,expected_pairs=7);names=sorted(r['sample_id'] for r in manifest);patient={r['sample_id']:r['patient_id'] for r in manifest}
    sam=executable('samtools','AIS_SAMTOOLS');som=executable('somalier','AIS_SOMALIER')
    assert subprocess.check_output([sam,'--version'],text=True).splitlines()[0]=='samtools 1.24'
    version=subprocess.run([som,'--help'],capture_output=True,text=True)
    assert 'somalier version: 0.3.5' in version.stdout+version.stderr
    binding=json.loads((a.identity/'reference_artifact_binding.json').read_text())['binding']['references']
    assert digest(sites)==binding['compressed']['sites.hg38.rna.vcf.gz']
    assert digest(a.reference/'GRCh38.primary_assembly.genome.fa')==binding['uncompressed']['GRCh38.primary_assembly.genome.fa']
    fingerprints=a.output/'private_fingerprints';fingerprints.mkdir(exist_ok=True)
    def sample(name):
        count={};bam=a.identity/name/'Aligned.sortedByCoord.out.bam'
        original=json.loads((a.identity/f'{name}.provenance.json').read_text())
        assert original['inputs']['sample']==name and digest(bam)==original['outputs']['bam_sha256']
        log=a.output/f'{name}.pileup.log'
        with log.open('w') as stderr:
            view=subprocess.Popen([sam,'view','-M','-L',str(bed),'-q','60','-F','1796','-@',str(a.view_threads),'-u',str(bam)],stdout=subprocess.PIPE,stderr=stderr)
            pileup=subprocess.Popen(pileup_arguments(sam,bed,a.reference/'GRCh38.primary_assembly.genome.fa',a.paired_overlap),stdin=view.stdout,stdout=subprocess.PIPE,stderr=stderr,text=True)
            view.stdout.close()
            for line in pileup.stdout:
                row=line.rstrip().split('\t');assert len(row)==6
                key=row[0],int(row[1]);assert key not in count
                count[key]=(row[2],parse_bases(row[4],row[2]))
            assert pileup.wait()==0 and view.wait()==0,'Native pileup failed: '+name
        counts=[]
        for c,pos,x,y in marker:
            ref,d=count.get((c,pos),(x,{x:0,y:0}))
            assert ref in (x,y),(c,pos,ref,x,y)
            counts.append((d.get(x,0),d.get(y,0),sum(d.values())-d.get(x,0)-d.get(y,0)))
        data=bytes([2,len(name)])+name.encode()+struct.pack('<HHH',len(counts),0,0)+b''.join(struct.pack('<III',*r) for r in counts)
        (fingerprints/f'{name}.somalier').write_bytes(data)
        print(name,'UNIQUE_HIGH_QUALITY_MARKER_PILEUP_PASS',flush=True)
        return name,counts
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:counts=dict(pool.map(sample,names))
    env=dict(os.environ,SOMALIER_REPORT_ALL_PAIRS='1')
    with (a.output/'native_relate.log').open('w') as log:
        subprocess.run([som,'relate','--sites',str(auto_sites),'--groups',str(a.identity/'expected_pairs.txt'),'--min-depth','7','--min-ab','0.2','--output-prefix',str(a.output/'strict_unique'),*[str(fingerprints/f'{n}.somalier') for n in names]],check=True,env=env,stdout=log,stderr=subprocess.STDOUT)
    pairs=m.rows(a.output/'strict_unique.pairs.tsv');assert len(pairs)==91
    calls={n:[m.gt(c) for c in counts[n]] for n in names};aggregate=[]
    for r in pairs:
        val=m.compare(calls[r['sample_a']],calls[r['sample_b']])
        for k in ['n','ibs0','ibs2','shared_hets','hets_ab','shared_hom_alts']:assert val[k]==int(r[k])
        assert abs(val['relatedness']-float(r['relatedness']))<=.0005001
        aggregate.append(dict(sample_a=r['sample_a'],sample_b=r['sample_b'],patient_a=patient[r['sample_a']],patient_b=patient[r['sample_b']],expected_metadata_pair=patient[r['sample_a']]==patient[r['sample_b']],**val))
    name='audit_RNA_pairing_unique_high_quality_reads'+('_overlap_retained' if a.paired_overlap=='retain' else '')+'.tsv'
    m.write(a.output/name,aggregate)
    report={'status':'PASS','samples':len(names),'markers':len(marker),'native_and_independent_pair_counts_verified':91,'mapq':60,'baseq':20,'BAQ':'disabled','paired_overlap_mode':a.paired_overlap,'paired_overlap_removal':'samtools default quality adjustment' if a.paired_overlap=='remove' else 'disabled; both overlapping mate bases counted','excluded_flags':1796,'read_selection':'indexed multi-region iterator, deduplicated overlapping target intervals','somalier_binary_sha256':digest(Path(som)),'samtools_version':'1.24','script_sha256':digest(Path(__file__)),'allele_decoder_script_sha256':digest(Path(m.__file__)),'scope':'Post hoc read-quality sensitivity; original RNA gates unchanged; no calibrated contamination or DNA/anatomy claim'}
    (a.output/'independent_unique_high_quality_pileup_QA.json').write_text(json.dumps(report,indent=2)+'\n')
    print('ALL14_UNIQUE_HIGH_QUALITY_READ_DIAGNOSTICS_PASS',flush=True)
if __name__=='__main__':main()
