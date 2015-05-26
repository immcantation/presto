"""
Unit tests for BuildConsensus
"""
# Future
from __future__ import absolute_import, division, print_function

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.5'
__date__      = '2015.04.20'

# Imports
import itertools, time, unittest
import IgCore as core
import BuildConsensus as mod
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
        result = core.qualityConsensus(self.records_dna, min_qual=0)
        print('  MIN_QUAL=0> %s %s' % (result.seq, result.letter_annotations['phred_quality']))
        self.assertEqual(self.cons_loose, str(result.seq))

        result = core.qualityConsensus(self.records_dna, min_qual=20)
        print(' MIN_QUAL=20> %s %s' % (result.seq, result.letter_annotations['phred_quality']))
        self.assertEqual(self.cons_qual, str(result.seq))

        result = core.qualityConsensus(self.records_dna, min_qual=20, min_freq=0.8)
        print('MIN_FREQ=0.8> %s %s' % (result.seq, result.letter_annotations['phred_quality']))
        self.assertEqual(self.cons_freqqual, str(result.seq))

    #@unittest.skip('-> frequencyConsensus() skipped\n')
    def test_frequencyConsensus(self):
        result = core.frequencyConsensus(self.records_dna, min_freq=0.2)
        print('MIN_FREQ=0.2> %s' % result.seq)
        self.assertEqual(self.cons_loose, str(result.seq))

        result = core.frequencyConsensus(self.records_dna, min_freq=0.8)
        print('MIN_FREQ=0.8> %s' % result.seq)
        self.assertEqual(self.cons_freq, str(result.seq))

    #@unittest.skip('-> calculateSetError() skipped\n')
    def test_calculateSetError(self):
        cons = core.frequencyConsensus(self.records_dna)
        error = core.calculateSetError(self.records_dna, cons)
        print('  REF> %s' % cons.seq)
        print('ERROR> %f' % error)
        self.assertAlmostEqual(1 - 27 / 30, error, places=4)

        cons = core.qualityConsensus(self.records_dna)
        error = core.calculateSetError(self.records_dna, cons)
        print('  REF> %s' % cons.seq)
        print('ERROR> %f' % error)
        self.assertAlmostEqual(1 - 23 / 28, error, places=4)

    #@unittest.skip('-> findGapPositions() skipped\n')
    def test_findGapPositions(self):
        result = core.findGapPositions(self.records_dna, max_gap=0.4)
        print('MAX_GAP=0.4> %s' % result)
        self.assertEqual([4, 5], result)

        result = core.findGapPositions(self.records_dna, max_gap=0.8)
        print('MAX_GAP=0.8> %s' % result)
        self.assertEqual([], result)

    #@unittest.skip('-> deleteSeqPositions() skipped\n')
    def test_deleteSeqPositions(self):
        pos = core.findGapPositions(self.records_dna, max_gap=0.4)
        cons = core.frequencyConsensus(self.records_dna)
        result = core.deleteSeqPositions(cons, pos)
        print('MAX_GAP=0.4> %s' % result.seq)
        self.assertEqual(self.cons_gap, str(result.seq))

        pos = core.findGapPositions(self.records_dna, max_gap=0.8)
        cons = core.frequencyConsensus(self.records_dna)
        result = core.deleteSeqPositions(cons, pos)
        print('MAX_GAP=0.8> %s' % result.seq)
        self.assertEqual(self.cons_loose, str(result.seq))