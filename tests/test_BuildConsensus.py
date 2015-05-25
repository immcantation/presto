"""
Unit tests for BuildConsensus
"""
# Future
from __future__ import absolute_import, division, print_function

# Imports
import time
import unittest
import presto.Sequence
from bin import BuildConsensus as script
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Info
__author__ = 'Jason Anthony Vander Heiden'


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
        self.cons_loose = 'CGGCGTAA'
        self.cons_freq = 'CGGCGTNA'
        self.cons_qual = 'CGGCNNAA'
        self.cons_freqqual = 'CGGCNNNA'
        self.cons_gap = 'CGGCAA'

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    #@unittest.skip('-> qualityConsensus() skipped\n')
    def test_qualityConsensus(self):
        result = presto.Sequence.qualityConsensus(self.records_dna, min_qual=0)
        print('  MIN_QUAL=0> %s %s' % (result.seq, result.letter_annotations['phred_quality']))
        self.assertEqual(self.cons_loose, str(result.seq))

        result = presto.Sequence.qualityConsensus(self.records_dna, min_qual=20)
        print(' MIN_QUAL=20> %s %s' % (result.seq, result.letter_annotations['phred_quality']))
        self.assertEqual(self.cons_qual, str(result.seq))

        result = presto.Sequence.qualityConsensus(self.records_dna, min_qual=20, min_freq=0.8)
        print('MIN_FREQ=0.8> %s %s' % (result.seq, result.letter_annotations['phred_quality']))
        self.assertEqual(self.cons_freqqual, str(result.seq))

    #@unittest.skip('-> frequencyConsensus() skipped\n')
    def test_frequencyConsensus(self):
        result = presto.Sequence.frequencyConsensus(self.records_dna, min_freq=0.2)
        print('MIN_FREQ=0.2> %s' % result.seq)
        self.assertEqual(self.cons_loose, str(result.seq))

        result = presto.Sequence.frequencyConsensus(self.records_dna, min_freq=0.8)
        print('MIN_FREQ=0.8> %s' % result.seq)
        self.assertEqual(self.cons_freq, str(result.seq))

    #@unittest.skip('-> calculateSetError() skipped\n')
    def test_calculateSetError(self):
        cons = presto.Sequence.frequencyConsensus(self.records_dna)
        error = script.calculateSetError(self.records_dna, cons)
        print('  REF> %s' % cons.seq)
        print('ERROR> %f' % error)
        self.assertAlmostEqual(1 - 27 / 30, error, places=4)

        cons = presto.Sequence.qualityConsensus(self.records_dna)
        error = script.calculateSetError(self.records_dna, cons)
        print('  REF> %s' % cons.seq)
        print('ERROR> %f' % error)
        self.assertAlmostEqual(1 - 23 / 28, error, places=4)

    #@unittest.skip('-> findGapPositions() skipped\n')
    def test_findGapPositions(self):
        result = script.findGapPositions(self.records_dna, max_gap=0.4)
        print('MAX_GAP=0.4> %s' % result)
        self.assertEqual([4, 5], result)

        result = script.findGapPositions(self.records_dna, max_gap=0.8)
        print('MAX_GAP=0.8> %s' % result)
        self.assertEqual([], result)

    #@unittest.skip('-> deleteSeqPositions() skipped\n')
    def test_deleteSeqPositions(self):
        pos = script.findGapPositions(self.records_dna, max_gap=0.4)
        cons = presto.Sequence.frequencyConsensus(self.records_dna)
        result = script.deleteSeqPositions(cons, pos)
        print('MAX_GAP=0.4> %s' % result.seq)
        self.assertEqual(self.cons_gap, str(result.seq))

        pos = script.findGapPositions(self.records_dna, max_gap=0.8)
        cons = presto.Sequence.frequencyConsensus(self.records_dna)
        result = script.deleteSeqPositions(cons, pos)
        print('MAX_GAP=0.8> %s' % result.seq)
        self.assertEqual(self.cons_loose, str(result.seq))