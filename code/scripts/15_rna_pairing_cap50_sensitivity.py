#!/usr/bin/env python3
"""Full-read post hoc STAR cap50 sensitivity for failed pairs and one control.

Preserve the uniform cap20 primary outputs. Results are diagnostic only; no
replacement gate or sample relabelling. Requires the original source-validated
FASTQs and same verified STAR index/reference, tools and read totals.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,json,os,signal,subprocess,time,collections
from input_integrity import digest,validate_pairs
from rna_reference_integrity import verify_reference
from project_paths import executable


def run_checked_stage(argv,logfile,timeout):
 # STAR can be a shell wrapper with a native child. On failure/timeout, stop
 # the whole isolated group so an orphan cannot keep writing a failed stage.
 with logfile.open('w') as log:
  process=subprocess.Popen(argv,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  try:
   code=process.wait(timeout=timeout)
   if code:raise subprocess.CalledProcessError(code,argv)
  except BaseException:
   try:os.killpg(process.pid,signal.SIGKILL)
   except ProcessLookupError:pass
   process.wait()
   raise


def selected_pairs(manifest,samples):
 if len(set(samples))!=len(samples) or any(n not in {r['sample_id'] for r in manifest} for n in samples):
  raise ValueError('Duplicate or unknown diagnostic sample')
 selected=[r for r in manifest if r['sample_id'] in samples];validate_pairs(selected)
 pairs=collections.defaultdict(list)
 for r in selected:pairs[r['patient_id']].append(r['sample_id'])
 return [sorted(v) for k,v in sorted(pairs.items())]

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ['raw','reference','identity','output']:p.add_argument('--'+k,type=Path,required=True)
 p.add_argument('--samples',nargs='+',default=['DR144','DR73','DR118','DR148','DR145','DR75'])
 p.add_argument('--manifest',type=Path,default=Path(__file__).resolve().parents[1]/'metadata/audit/gse311180_requant_manifest.tsv')
 p.add_argument('--validate-inputs-only',action='store_true')
 p.add_argument('--threads',type=int,default=8);p.add_argument('--jobs',type=int,default=2);p.add_argument('--timeout',type=int,default=10800)
 a=p.parse_args();assert a.threads==8 and a.jobs>0 and a.timeout>0
 import csv
 with a.manifest.open(newline='') as handle:manifest=list(csv.DictReader(handle,delimiter='\t'))
 validate_pairs(manifest,expected_pairs=7);groups_selected=selected_pairs(manifest,a.samples)
 if a.output.exists() and any(a.output.iterdir()):raise ValueError('Preserve previous diagnostic stages; choose a new output directory')
 a.output.mkdir(parents=True,exist_ok=True)
 star=executable('STAR','AIS_STAR');sam=executable('samtools','AIS_SAMTOOLS');som=executable('somalier','AIS_SOMALIER')
 assert subprocess.check_output([star,'--version'],text=True).strip()=='2.7.11b'
 assert subprocess.check_output([sam,'--version'],text=True).splitlines()[0]=='samtools 1.24'
 version=subprocess.run([som,'--help'],capture_output=True,text=True)
 assert 'somalier version: 0.3.5' in version.stdout+version.stderr
 verified=verify_reference(a.reference)
 assert verified==json.loads((a.identity/'reference_artifact_binding.json').read_text())['binding']['references']
 def run(argv,logfile):
  run_checked_stage(argv,logfile,a.timeout)
 def sample(n):
  dst=a.output/n
  assert not dst.exists(),'Preserve previous diagnostic stage: '+n
  dst.mkdir();stamp=json.loads((a.identity/f'{n}.provenance.json').read_text());fq=json.loads((a.raw/f'bulk_gse311180/fastq/{n}.provenance.json').read_text())
  mates=[a.raw/f'bulk_gse311180/fastq/{n}_{i}.fastq' for i in [1,2]]
  row=next(r for r in manifest if r['sample_id']==n)
  assert fq['inputs']['sample_id']==n and fq['inputs']['sra_run']==row['sra_run']
  assert stamp['inputs']['sample']==n and stamp['inputs']['fastq_provenance_sha256']==digest(a.raw/f'bulk_gse311180/fastq/{n}.provenance.json')
  assert digest(a.identity/'fingerprints'/f'{n}.somalier')==stamp['outputs']['fingerprint_sha256']
  for key,value in {'STAR':'2.7.11b','Somalier':'0.3.5','samtools':'1.24','seedPerWindowNmax':20,'alignment_threads':8,'twopassMode':'Basic','unique_mapq':60}.items():assert stamp['inputs']['recipe'][key]==value
  inputs={'sample':n,'primary_fingerprint_sha256':digest(a.identity/'fingerprints'/f'{n}.somalier'),
     'prior_verified_FASTQ_sha256':[fq['outputs']['mate1_sha256'],fq['outputs']['mate2_sha256']],
     'prior_source_provenance_sha256':digest(a.raw/f'bulk_gse311180/fastq/{n}.provenance.json'),
     'reference_binding_sha256':digest(a.identity/'reference_artifact_binding.json'),
     'seedPerWindowNmax':50,'source_script_sha256':digest(Path(__file__))}
  start=time.time()
  try:
   # Revalidate read bytes; preserve original primary artifacts in a distinct output.
   for i,file in enumerate(mates):assert digest(file)==inputs['prior_verified_FASTQ_sha256'][i]
   if a.validate_inputs_only:
    result={'status':'PASS','inputs':inputs,'scope':'Whole source/reference/tool validation only; no new alignment'}
    (dst/'verified_inputs.json').write_text(json.dumps(result,indent=2)+'\n');return result
   run([star,'--runThreadN',str(a.threads),'--genomeDir',str(a.reference/'STAR_2.7.11b'),'--readFilesIn',*map(str,mates),'--twopassMode','Basic','--outSAMtype','BAM','Unsorted','--outSAMmapqUnique','60','--seedPerWindowNmax','50','--outSAMattrRGline',f'ID:{n}',f'SM:{n}','PL:ILLUMINA','--outFileNamePrefix',str(dst)+'/' ],dst/'alignment.log')
   lines={x.split('|',1)[0].strip():x.split('|',1)[1].strip() for x in (dst/'Log.final.out').read_text().splitlines() if '|' in x}
   assert int(lines['Number of input reads'])==fq['outputs']['read_pairs']
   bam=dst/'Aligned.sortedByCoord.out.bam'
   run([sam,'sort','-@',str(a.threads),'-m','512M','-o',str(bam),str(dst/'Aligned.out.bam')],dst/'sort.log')
   subprocess.run([sam,'quickcheck','-v',str(bam)],check=True)
   run([sam,'index','-@',str(a.threads),str(bam)],dst/'index.log')
   run([som,'extract','--sites',str(a.reference/'sites.hg38.rna.vcf.gz'),'--fasta',str(a.reference/'GRCh38.primary_assembly.genome.fa'),'--out-dir',str(dst),str(bam)],dst/'extract.log')
   (dst/'Aligned.out.bam').unlink()
   result={'status':'PASS','inputs':inputs,'STAR_input_pairs':int(lines['Number of input reads']),
       'STAR_uniquely_mapped_percent':float(lines['Uniquely mapped reads %'].rstrip('%')),
       'BAM_sha256':digest(bam),'fingerprint_sha256':digest(dst/f'{n}.somalier'),'elapsed_seconds':time.time()-start}
  except Exception as e:
   result={'status':'FAILED_PRESERVED','inputs':inputs,'error':str(e),'elapsed_seconds':time.time()-start}
  (dst/'completed_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
  print(n,result['status'],flush=True);return result
 with ThreadPoolExecutor(max_workers=a.jobs) as pool:results=list(pool.map(sample,a.samples))
 (a.output/('input_validation_summary.json' if a.validate_inputs_only else 'full_read_cap50_summary.json')).write_text(json.dumps(results,indent=2)+'\n')
 if any(r['status']!='PASS' for r in results):raise RuntimeError('Some diagnostics failed; preserve and report them rather than accepting partial output')
 if a.validate_inputs_only:
  print('ALL_DIAGNOSTIC_SOURCE_REFERENCE_TOOL_GUARDS_PASS',flush=True);return
 os.environ['SOMALIER_REPORT_ALL_PAIRS']='1'
 groups=a.output/'expected_pairs.txt';groups.write_text(''.join(','.join(v)+'\n' for v in groups_selected))
 run([som,'relate','--sites',str(a.reference/'sites.hg38.rna.vcf.gz'),'--groups',str(groups),'--min-depth','7','--min-ab','0.2','--output-prefix',str(a.output/'cap50'),*[str(a.output/n/f'{n}.somalier') for n in a.samples]],a.output/'relate.log')
 print(f'CAP50_COMPLETE_FULL_READ_{len(a.samples)}_LIBRARY_DIAGNOSTIC_PASS',flush=True)
if __name__=='__main__':main()
