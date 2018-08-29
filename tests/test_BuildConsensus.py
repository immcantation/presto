"""
Unit tests for BuildConsensus
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import os
import sys
import time
import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Sequence import calculateSetError, deleteSeqPositions, findGapPositions, \
                            frequencyConsensus, qualityConsensus

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import BuildConsensus


class TestBuildConsensus(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        # Test DNA
        seq_dna = [Seq('CGGCGTAA'),
                   Seq('CCNNGTAA'),
                   Seq('CGGC--TA'),
                   Seq('CGNN--TA'),
                   Seq('CGGC--AA')]
        seq_qual = [30, 40, 20, 20, 20, 20, 40, 40]
        self.records_dna = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='',
                                      letter_annotations={'phred_quality':seq_qual})
                            for i, s in enumerate(seq_dna, start=1)]
        self.nonmiss_dna = 8 + 6 + 6 + 4 + 6

        # Consensus sequences
        self.cons_loose = 'CGGCGTAA'
        self.cons_freq = 'CGGCGTNA'
        self.cons_qual = 'CGGCNNAA'
        self.cons_freqqual = 'CGGCNNNA'
        self.cons_gap = 'CGGCAA'

        # Set error against cons_freq and cons_qual
        self.freq_error = 1.0 - 27.0 / 36
        self.qual_error = 1.0 - 23.0 / 30

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    #@unittest.skip('-> qualityConsensus() skipped\n')
    def test_qualityConsensus(self):
        result = qualityConsensus(self.records_dna, min_qual=0)
        print('  MIN_QUAL=0> %s %s' % (result.seq, result.letter_annotations['phred_quality']))
        self.assertEqual(self.cons_loose, str(result.seq))

        result = qualityConsensus(self.records_dna, min_qual=20)
        print(' MIN_QUAL=20> %s %s' % (result.seq, result.letter_annotations['phred_quality']))
        self.assertEqual(self.cons_qual, str(result.seq))

        result = qualityConsensus(self.records_dna, min_qual=20, min_freq=0.8)
        print('MIN_FREQ=0.8> %s %s' % (result.seq, result.letter_annotations['phred_quality']))
        self.assertEqual(self.cons_freqqual, str(result.seq))

    #@unittest.skip('-> frequencyConsensus() skipped\n')
    def test_frequencyConsensus(self):
        result = frequencyConsensus(self.records_dna, min_freq=0.2)
        print('MIN_FREQ=0.2> %s' % result.seq)
        self.assertEqual(self.cons_loose, str(result.seq))

        result = frequencyConsensus(self.records_dna, min_freq=0.8)
        print('MIN_FREQ=0.8> %s' % result.seq)
        self.assertEqual(self.cons_freq, str(result.seq))

    #@unittest.skip('-> calculateSetError() skipped\n')
    def test_calculateSetError(self):
        cons = frequencyConsensus(self.records_dna)
        error = calculateSetError(self.records_dna, cons)
        print('Frequency consensus error>')
        print('  REF> %s' % cons.seq)
        for x in self.records_dna:  print(' %s> %s' %(x.id, x.seq))
        print('ERROR> %f' % error)
        self.assertAlmostEqual(self.freq_error, error, places=4)

        print('Quality consensus error>')
        cons = qualityConsensus(self.records_dna, min_qual=20)
        error = calculateSetError(self.records_dna, cons)
        print('  REF> %s' % cons.seq)
        for x in self.records_dna:  print(' %s> %s' %(x.id, x.seq))
        print('ERROR> %f' % error)
        self.assertAlmostEqual(self.qual_error, error, places=4)

    #@unittest.skip('-> findGapPositions() skipped\n')
    def test_findGapPositions(self):
        result = findGapPositions(self.records_dna, max_gap=0.4)
        print('MAX_GAP=0.4> %s' % result)
        self.assertEqual([4, 5], result)

        result = findGapPositions(self.records_dna, max_gap=0.8)
        print('MAX_GAP=0.8> %s' % result)
        self.assertEqual([], result)

    #@unittest.skip('-> deleteSeqPositions() skipped\n')
    def test_deleteSeqPositions(self):
        pos = findGapPositions(self.records_dna, max_gap=0.4)
        cons = frequencyConsensus(self.records_dna)
        result = deleteSeqPositions(cons, pos)
        print('MAX_GAP=0.4> %s' % result.seq)
        self.assertEqual(self.cons_gap, str(result.seq))

        pos = findGapPositions(self.records_dna, max_gap=0.8)
        cons = frequencyConsensus(self.records_dna)
        result = deleteSeqPositions(cons, pos)
        print('MAX_GAP=0.8> %s' % result.seq)
        self.assertEqual(self.cons_loose, str(result.seq))


if __name__ == '__main__':
    unittest.main()