#!/usr/bin/env python3
"""Finish a fully completed legacy STAR output after its wrapper timed out.

Never resume an incomplete alignment. Verify full raw bytes/reference/header,
ALL DONE and input totals, then sort/extract in a separate new namespace. Keep
the original FAILED receipt and all prior files unchanged. No pair verdict here.
"""
import argparse,json,shlex,shutil,time
from pathlib import Path
import subprocess,importlib
from input_integrity import digest
from rna_reference_integrity import verify_reference
from project_paths import executable
f=importlib.import_module('15_rna_pairing_cap50_sensitivity')
ROOT=Path(__file__).resolve().parents[1]

def require_native_completion(directory):
    # STAR writes ALL DONE to Log.out, not its redirected stdout.
    native=directory/'Log.out'
    assert 'ALL DONE!' in native.read_text()
    assert 'finished successfully' in (directory/'alignment.log').read_text()
    return native

def check_header(header,sample,reference,raw,original):
    pg=[dict(x.split(':',1) for x in line.split('\t')[1:]) for line in header.splitlines() if line.startswith('@PG\t')]
    star=next(x for x in pg if x.get('PN')=='STAR');assert star['VN']=='2.7.11b'
    args=shlex.split(star['CL'])
    required={'--runThreadN':['8'],'--genomeDir':[str(reference/'STAR_2.7.11b')],
        '--readFilesIn':[str(raw/f'bulk_gse311180/fastq/{sample}_{i}.fastq') for i in [1,2]],
        '--twopassMode':['Basic'],'--outSAMtype':['BAM','Unsorted'],'--outSAMmapqUnique':['60'],
        '--seedPerWindowNmax':['50'],'--outSAMattrRGline':[f'ID:{sample}',f'SM:{sample}','PL:ILLUMINA'],
        '--outFileNamePrefix':[str(original)+'/']}
    for flag,value in required.items():
        assert args.count(flag)==1 and args[args.index(flag)+1:args.index(flag)+1+len(value)]==value,flag
    rg=[dict(x.split(':',1) for x in line.split('\t')[1:]) for line in header.splitlines() if line.startswith('@RG\t')]
    assert len(rg)==1 and rg[0]['SM']==rg[0]['ID']==sample

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['original','output','raw','reference','primary']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--sample',required=True);p.add_argument('--timeout',type=int,default=3600)
    a=p.parse_args();n=a.sample;assert n in {'DR118','DR148'}
    if a.output.exists():raise ValueError('Use a new recovery namespace; never replace legacy failure artifacts')
    r=json.loads((a.original/'completed_diagnostic.json').read_text());assert r['status']=='FAILED_PRESERVED' and 'timed out' in r['error'] and '--twopassMode' in r['error']
    inputs=r['inputs'];assert inputs['sample']==n and inputs['seedPerWindowNmax']==50
    allowed={digest(ROOT/'scripts/execution_records'/name) for name in ['rna_cap50_alignment_producer_20261013.txt','rna_cap50_guarded_alignment_producer_20261013.txt']}
    assert inputs['source_script_sha256'] in allowed
    assert inputs['reference_binding_sha256']==digest(a.primary/'reference_artifact_binding.json')
    assert verify_reference(a.reference)==json.loads((a.primary/'reference_artifact_binding.json').read_text())['binding']['references']
    assert inputs['primary_fingerprint_sha256']==digest(a.primary/'fingerprints'/f'{n}.somalier')
    fq=a.raw/f'bulk_gse311180/fastq/{n}.provenance.json';source=json.loads(fq.read_text())
    assert digest(fq)==inputs['prior_source_provenance_sha256'] and source['inputs']['sample_id']==n
    for i in [1,2]:assert digest(a.raw/f'bulk_gse311180/fastq/{n}_{i}.fastq')==inputs['prior_verified_FASTQ_sha256'][i-1]
    log=require_native_completion(a.original)
    final=a.original/'Log.final.out';values={x.split('|',1)[0].strip():x.split('|',1)[1].strip() for x in final.read_text().splitlines() if '|' in x}
    assert int(values['Number of input reads'])==source['outputs']['read_pairs']
    sam=executable('samtools','AIS_SAMTOOLS');som=executable('somalier','AIS_SOMALIER')
    assert subprocess.check_output([sam,'--version'],text=True).splitlines()[0]=='samtools 1.24'
    version=subprocess.run([som,'--help'],capture_output=True,text=True);assert 'somalier version: 0.3.5' in version.stdout+version.stderr
    unsorted=a.original/'Aligned.out.bam';subprocess.run([sam,'quickcheck','-v',str(unsorted)],check=True)
    check_header(subprocess.check_output([sam,'view','-H',str(unsorted)],text=True),n,a.reference,a.raw,a.original)
    original_bam_sha=digest(unsorted);a.output.mkdir(parents=True);start=time.time();bam=a.output/'Aligned.sortedByCoord.out.bam'
    f.run_checked_stage([sam,'sort','-@','8','-m','512M','-o',str(bam),str(unsorted)],a.output/'sort.log',a.timeout)
    assert digest(unsorted)==original_bam_sha,'Legacy output changed during recovery; do not accept it'
    subprocess.run([sam,'quickcheck','-v',str(bam)],check=True)
    f.run_checked_stage([sam,'index','-@','8',str(bam)],a.output/'index.log',a.timeout)
    f.run_checked_stage([som,'extract','--sites',str(a.reference/'sites.hg38.rna.vcf.gz'),'--fasta',str(a.reference/'GRCh38.primary_assembly.genome.fa'),'--out-dir',str(a.output),str(bam)],a.output/'extract.log',a.timeout)
    for name in ['alignment.log','Log.out','Log.final.out','Log.progress.out']:shutil.copy2(a.original/name,a.output/name)
    recovery={'original_controller_status':r['status'],'original_receipt_sha256':digest(a.original/'completed_diagnostic.json'),
        'completed_native_alignment_log_sha256':digest(log),'verified_unsorted_BAM_sha256':original_bam_sha,
        'recovery_script_sha256':digest(Path(__file__)),'whole_sources_and_reference_reverified':True,
        'scope':'Legacy native STAR fully completed after wrapper timeout; separate verified sort/extraction, not an incomplete alignment resume'}
    out={'status':'PASS','inputs':inputs,'STAR_input_pairs':int(values['Number of input reads']),
        'STAR_uniquely_mapped_percent':float(values['Uniquely mapped reads %'].rstrip('%')),
        'BAM_sha256':digest(bam),'fingerprint_sha256':digest(a.output/f'{n}.somalier'),
        'elapsed_seconds':time.time()-start,'recovery':recovery}
    (a.output/'completed_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
