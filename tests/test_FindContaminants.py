"""
Unit tests for FindContaminants
"""
# Info
__author__ = 'Gisela Gabernet'

# Imports
import os
import shutil
import sys
import tempfile
import time
import unittest
from collections import OrderedDict
from Bio import SeqIO

# Presto imports
from presto.Annotation import parseAnnotation, parseLog
from presto.Defaults import default_out_args

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import FindContaminants


class TestFindContaminants(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        self.labels = ['A', 'B', 'C']
        self.files_1 = [os.path.join(data_path, 'contam_sample%s_R1.fastq' % s) for s in self.labels]
        self.files_2 = [os.path.join(data_path, 'contam_sample%s_R2.fastq' % s) for s in self.labels]
        self.out_dir = tempfile.mkdtemp()
        self.start = time.time()

    def tearDown(self):
        shutil.rmtree(self.out_dir)
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    def _write_fastq(self, name, records):
        path = os.path.join(self.out_dir, name)
        with open(path, 'w') as h:
            for title, seq in records:
                h.write('@%s\n%s\n+\n%s\n' % (title, seq, 'I' * len(seq)))
        return path

    def _run(self, **kwargs):
        out_dir = tempfile.mkdtemp(dir=self.out_dir)
        out_args = dict(default_out_args)
        out_args.update({'out_dir': out_dir, 'failed': True, 'gzip_output': False})
        out_args['log_file'] = os.path.join(out_dir, 'log.txt')
        FindContaminants.findContaminants(self.files_1, self.files_2, labels=self.labels,
                                          out_args=out_args, **kwargs)
        result = {}
        for s in self.labels:
            for r in ('R1', 'R2'):
                for label in ('pass', 'fail'):
                    f = os.path.join(out_dir, 'contam_sample%s_%s_contam-%s.fastq' % (s, r, label))
                    result[(s, r, label)] = OrderedDict(
                        (parseAnnotation(rec.id)['ID'], parseAnnotation(rec.id))
                        for rec in SeqIO.parse(f, 'fastq'))
        with open(out_args['log_file']) as h:
            log = [parseLog(r) for r in h.read().split('\n\n') if r.strip()]
        return result, log

    def test_matchSeq(self):
        match = FindContaminants.matchSeq
        self.assertTrue(match('ACGTACGTAC', 'ACGTACGTAC'))
        # Ns are ignored, but different lengths never match
        self.assertTrue(match('ACGTACGTAC', 'ACGNACGTAC'))
        self.assertFalse(match('ACGTACGTAC', 'ACGTACGTACGGGG'))
        self.assertFalse(match('ACGTACGTAC', 'ACGTACGTA', max_error=0.5))
        self.assertFalse(match('ACGTACGTAC', 'ACGTACGTAA'))
        self.assertTrue(match('ACGTACGTAC', 'ACGTACGTAA', max_error=0.1))
        self.assertFalse(match('ACGTACGTAC', 'ACGTACGTTT', max_error=0.1))
        # Nothing comparable never matches
        self.assertFalse(match('NNNNN', 'NNNNN'))
        self.assertFalse(match('ACGT', 'NNNN'))

    def test_findSharedReads(self):
        reads = OrderedDict([('A', [('a1', 'AAAA', 'CCCC'), ('a2', 'AAAA', 'CCCC'),
                                    ('a3', 'AAAA', 'GGGG')]),
                             ('B', [('b1', 'AAAA', 'CCCC')]),
                             ('C', [('c1', 'AAAA', 'TTTT'), ('c2', 'AAAA', 'GGGG')])])
        shared = FindContaminants.findSharedReads(reads)
        self.assertEqual({'a1': ['B'], 'a2': ['B'], 'a3': ['C']}, shared['A'])
        self.assertEqual({'b1': ['A']}, shared['B'])
        self.assertEqual({'c2': ['A']}, shared['C'])

    def test_findSharedReadsHashing(self):
        # Hash based matching must agree with scoring every pair of read pairs
        import random
        from itertools import combinations
        rng = random.Random(7)
        pool = ['ACGTACGT', 'ACGTACG', 'ACGTTCGT', 'ACGNACGT', 'TTGCA', 'TTGCAAG', 'TTGCT',
                'NNNN', '']
        for trial in range(400):
            max_error = 0.0 if trial % 2 == 0 else 0.2
            reads = OrderedDict()
            for s in 'ABCD':
                reads[s] = [('%s%i' % (s, i), rng.choice(pool), rng.choice(pool))
                            for i in range(rng.randint(0, 6))]
            expected = {s: {} for s in reads}
            for a, b in combinations(reads, 2):
                for ka, a1, a2 in reads[a]:
                    for kb, b1, b2 in reads[b]:
                        if FindContaminants.matchSeq(a1, b1, max_error) and \
                                FindContaminants.matchSeq(a2, b2, max_error):
                            expected[a].setdefault(ka, set()).add(b)
                            expected[b].setdefault(kb, set()).add(a)
            expected = {s: {k: sorted(v) for k, v in d.items()} for s, d in expected.items()}
            self.assertEqual(expected, FindContaminants.findSharedReads(reads, max_error=max_error),
                             msg=str(reads))

    def test_collectSampleReads(self):
        # Barcode on read 1 only, plus an unpaired read 1
        f1 = self._write_fastq('r1.fastq', [('S1|BARCODE=AAA', 'ACGT'), ('S2|BARCODE=CCC', 'TTTT'),
                                            ('S3|BARCODE=AAA', 'GGGG')])
        f2 = self._write_fastq('r2.fastq', [('S1', 'CCCC'), ('S2', 'AAAA')])
        groups, conflicts = FindContaminants.collectSampleReads(f1, f2)
        self.assertEqual({'AAA': [('S1', 'ACGT', 'CCCC')], 'CCC': [('S2', 'TTTT', 'AAAA')]}, groups)
        self.assertEqual(0, conflicts)
        groups, _ = FindContaminants.collectSampleReads(f1, f2, candidates={'CCC'})
        self.assertEqual({'CCC': [('S2', 'TTTT', 'AAAA')]}, groups)

        # Barcode on read 2 only, and a conflict where read 1 takes precedence
        f1 = self._write_fastq('r1b.fastq', [('S1', 'acgt'), ('S2|BARCODE=GGG', 'TTTT'), ('S3', 'GGGG')])
        f2 = self._write_fastq('r2b.fastq', [('S1|BARCODE=AAA', 'CCCC'), ('S2|BARCODE=CCC', 'AAAA'),
                                             ('S4|BARCODE=AAA', 'TTTT')])
        groups, conflicts = FindContaminants.collectSampleReads(f1, f2)
        self.assertEqual({'AAA': [('S1', 'ACGT', 'CCCC')], 'GGG': [('S2', 'TTTT', 'AAAA')]}, groups)
        self.assertEqual(1, conflicts)
        self.assertEqual({'AAA', 'CCC', 'GGG'},
                         FindContaminants.collectSampleBarcodes(f1, f2))

    def test_findContaminants(self):
        result, log = self._run()
        for r in ('R1', 'R2'):
            # Only the A read matching the B read with barcode CCCCAAAAGGGG is removed
            self.assertEqual(['A_READ01', 'A_READ02', 'A_READ03', 'A_READ08'],
                             list(result[('A', r, 'fail')]))
            self.assertEqual(['B_READ01', 'B_READ02', 'B_READ03', 'B_READ07'],
                             list(result[('B', r, 'fail')]))
            self.assertEqual(['C_READ02', 'C_READ03'], list(result[('C', r, 'fail')]))
            self.assertIn('A_READ09', result[('A', r, 'pass')])
        self.assertEqual('B', result[('A', 'R2', 'fail')]['A_READ01']['SHARED'])
        self.assertEqual('C', result[('B', 'R1', 'fail')]['B_READ02']['SHARED'])
        # UMI collision with different R2 and one R1 mismatch are not flagged
        self.assertIn('A_READ04', result[('A', 'R1', 'pass')])
        self.assertIn('A_READ05', result[('A', 'R1', 'pass')])
        self.assertIn('C_READ01', result[('C', 'R1', 'pass')])
        # R2 read without mate or barcode passes through
        self.assertIn('A_READ99', result[('A', 'R2', 'pass')])
        # Log holds one record per contaminated barcode with per-sample counts
        found = {r['BARCODE']: (r['SAMPLES'], r['SEQCOUNT'], r['GROUPCOUNT']) for r in log}
        self.assertEqual({'AAAACCCCGGGG': ('A,B', '3,1', '3,1'),
                          'GGGGTTTTAAAA': ('B,C', '2,2', '2,2'),
                          'CCCCAAAAGGGG': ('A,B', '1,1', '2,1')}, found)

    def test_findContaminantsMaxError(self):
        result, _ = self._run(max_error=0.05)
        # One mismatch over 30 bases is now allowed
        self.assertIn('A_READ05', result[('A', 'R1', 'fail')])
        self.assertIn('B_READ04', result[('B', 'R1', 'fail')])
        self.assertNotIn('B_READ04', result[('B', 'R1', 'pass')])


if __name__ == '__main__':
    unittest.main()
