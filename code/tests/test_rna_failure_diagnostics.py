import importlib
import struct
import os
import json
import subprocess
import signal
import time
import tempfile
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
d=importlib.import_module('14_rna_pairing_failure_diagnostics')
p=importlib.import_module('16_rna_unique_read_pileup_diagnostics')
c=importlib.import_module('17_rna_dropout_context')
v=importlib.import_module('verify_rna_failure_diagnostics')
f=importlib.import_module('15_rna_pairing_cap50_sensitivity')
g=importlib.import_module('23_rna_fragment_coordinate_diversity')
h=importlib.import_module('22_finish_completed_legacy_cap50_alignment')

class FailureDiagnosticsTests(unittest.TestCase):
    def test_legacy_completion_requires_native_log_and_successful_stdout(self):
        with tempfile.TemporaryDirectory() as folder:
            d=Path(folder)
            (d/'Log.out').write_text('ALL DONE!\n')
            (d/'alignment.log').write_text('finished successfully\n')
            self.assertEqual(h.require_native_completion(d),d/'Log.out')
            (d/'Log.out').write_text('still mapping\n')
            with self.assertRaises(AssertionError):h.require_native_completion(d)
            (d/'Log.out').write_text('ALL DONE!\n')
            (d/'alignment.log').write_text('still mapping\n')
            with self.assertRaises(AssertionError):h.require_native_completion(d)
    def test_legacy_complete_header_rejects_wrong_seed_pair_or_recipe(self):
        ref=Path('ref');raw=Path('raw');original=Path('old/DR118')
        command='STAR --runThreadN 8 --genomeDir ref/STAR_2.7.11b --readFilesIn raw/bulk_gse311180/fastq/DR118_1.fastq raw/bulk_gse311180/fastq/DR118_2.fastq --twopassMode Basic --outSAMtype BAM Unsorted --outSAMmapqUnique 60 --seedPerWindowNmax 50 --outSAMattrRGline ID:DR118 SM:DR118 PL:ILLUMINA --outFileNamePrefix old/DR118/'
        header='@PG\tID:STAR\tPN:STAR\tVN:2.7.11b\tCL:'+command+'\n@RG\tID:DR118\tSM:DR118\tPL:ILLUMINA\n'
        h.check_header(header,'DR118',ref,raw,original)
        for bad in [header.replace('--seedPerWindowNmax 50','--seedPerWindowNmax 20'),header.replace('--twopassMode Basic','--twopassMode None'),header.replace('SM:DR118','SM:DR148')]:
            with self.assertRaises(AssertionError):h.check_header(bad,'DR118',ref,raw,original)

    def test_fragment_proxy_cigar_splicing_clipping_indels_and_base_quality(self):
        seq='AACGTACGTAA';qual=list('I'*len(seq));qual[3]='4'
        calls=list(g.aligned_bases(10,'2S3M2N2M1I2M1D1M',seq,''.join(qual),list(range(10,25))))
        self.assertEqual(calls,[(10,'C'),(12,'T'),(15,'A'),(16,'C'),(17,'T'),(18,'A'),(20,'A')])
        for cigar in ['0M','2Z','3M']:
            with self.assertRaises(ValueError):list(g.aligned_bases(10,cigar,seq,''.join(qual),[10]))
    def test_fragment_coordinate_collisions_are_not_molecule_counts(self):
        one=('chr1',10,150);two=('chr1',20,150)
        result=g.summarize_templates({('q1',one):{'A'},('q2',one):{'A'},('q3',two):{'G'},('q4',two):{'A','G'},('q5',two):{'N'}},('A','G'))
        self.assertEqual((result['qname_templates'],result['distinct_coordinates'],result['conflicting_or_other_templates']),(3,2,2))
        self.assertEqual(result['allele_coordinates'],{'A':{one},'G':{two}})

    @unittest.skipUnless(os.name=='posix','Native STAR workflow uses POSIX process groups')
    def test_timed_out_native_stage_stops_wrapper_and_child(self):
        with tempfile.TemporaryDirectory() as folder:
            pidfile=Path(folder)/'pids.json'
            heartbeat=Path(folder)/'child_heartbeat.txt'
            childcode="import sys,time; from pathlib import Path; p=Path(sys.argv[1]);\nwhile True: p.write_text(str(time.monotonic_ns())); time.sleep(.01)"
            code="import os,json,time,subprocess,sys; from pathlib import Path; child=subprocess.Popen([sys.executable,'-c',sys.argv[3],sys.argv[2]]); Path(sys.argv[1]).write_text(json.dumps([os.getpid(),child.pid])); time.sleep(30)"
            try:
                with self.assertRaises(subprocess.TimeoutExpired):
                    f.run_checked_stage([sys.executable,'-c',code,str(pidfile),str(heartbeat),childcode],Path(folder)/'stage.log',1)
                parent,child=json.loads(pidfile.read_text())
                time.sleep(.05)
                before=heartbeat.read_bytes()
                time.sleep(.15)
                self.assertEqual(before,heartbeat.read_bytes(),'Timed-out native child still writes failed stage')
            finally:
                if pidfile.exists():
                    parent,child=json.loads(pidfile.read_text())
                    try:os.killpg(parent,signal.SIGKILL)
                    except ProcessLookupError:pass

    def test_cap50_groups_reject_partial_unknown_and_duplicate_samples(self):
        meta=[dict(sample_id='A',patient_id='P',side='concave'),dict(sample_id='B',patient_id='P',side='convex')]
        self.assertEqual(f.selected_pairs(meta,['B','A']),[['A','B']])
        for bad in [['A'],['A','X'],['A','A','B']]:
            with self.assertRaises(ValueError):f.selected_pairs(meta,bad)
    def test_quality_pileup_explicitly_distinguishes_mate_overlap_modes(self):
        removed=p.pileup_arguments('samtools','sites.bed','genome.fa','remove')
        retained=p.pileup_arguments('samtools','sites.bed','genome.fa','retain')
        self.assertNotIn('-x',removed);self.assertIn('-x',retained)
        self.assertEqual([x for x in retained if x!='-x'],removed)
        with self.assertRaises(ValueError):p.pileup_arguments('s','b','f','unspecified')
    def test_native_float32_other_base_boundary(self):
        self.assertGreater(d.other_fraction((9,0,1)),.1)
        self.assertEqual(d.gt((9,0,1)),-1)
        self.assertEqual(d.gt((10,0,1)),0)
        self.assertEqual(d.gt((4,3,0)),1)
    def test_binary_truncation_and_extra_bytes_rejected(self):
        raw=bytes([2,1])+b'S'+struct.pack('<HHH',1,0,0)+struct.pack('<III',8,7,0)
        with tempfile.TemporaryDirectory() as folder:
            file=Path(folder)/'x.somalier';file.write_bytes(raw)
            self.assertEqual(d.decode(file),('S',(1,0,0),[(8,7,0)]))
            for damaged in [raw[:-1],raw+b'!']:
                file.write_bytes(damaged)
                with self.assertRaises(ValueError):d.decode(file)
    def test_loss_decomposition_counts_both_discordance_directions(self):
        r=d.compare([1,0,2,0,1],[1,1,1,2,1])
        self.assertEqual((r['n'],r['ibs0'],r['heterozygote_homozygote_discordance']),(5,1,2))
        self.assertAlmostEqual(1-r['relatedness'],r['loss_from_het_hom']+r['loss_from_IBS0'])
    def test_retained_verifier_rejects_wrong_relatedness_and_discordance(self):
        good=d.compare([1,0,1,0],[1,1,1,2]);v.arithmetic(good)
        for key,value in [('relatedness',.9),('heterozygote_homozygote_discordance',0),('hets_ab',8)]:
            bad=dict(good);bad[key]=value
            with self.assertRaises(AssertionError):v.arithmetic(bad)
    def test_pileup_controls_indels_and_strand_not_double_counted(self):
        self.assertEqual(p.parse_bases('^].,+2AGa-3tttC$*#<>','G'),{'A':1,'C':1,'G':2,'T':0,'N':0})
        for bad in ['.+AG','.-4AC','.?','^','.-0']:
            with self.assertRaises(ValueError):p.parse_bases(bad,'A')
    def test_directional_aggregation_preserves_homozygous_side(self):
        rows=c.summarize('P','A','B',[(7,0,0),(4,4,0)],[(10,10,0),(0,8,0)],[('chr1',1,'A','G'),('chr2',1,'A','G')])
        self.assertEqual(rows[0]['homozygous_call_sample'],'A')
        self.assertEqual(rows[0]['hom_depth_median'],7)
        self.assertEqual(rows[1]['homozygous_call_sample'],'B')
        self.assertEqual(rows[1]['zero_minor_allele_sites'],1)
if __name__=='__main__':unittest.main()
