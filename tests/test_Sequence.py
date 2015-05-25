"""
Unit tests for Sequence module
"""

# Imports
import time
import itertools
import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import presto.Sequence

# Info
__author__ = 'Jason Anthony Vander Heiden'


class TestIgCore(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName
        # Test DNA
        seq_dna = [Seq('CGGCGTAA'),
                   Seq('CGNNGTAA'),
                   Seq('CGGC--AA'),
                   Seq('CGNN--AA')]
        self.records_dna = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                            for i, s in enumerate(seq_dna, start=1)]
        self.weight_dna = [8, 6, 6, 4]
        # scoreSeqPair returns (score, minimum weight, error rate)
        # TODO:  this 0.3333 case is funny (SEQ2 vs SEQ3). fix is not obvious.
        self.error_dna = [0.0, 0.0, 0.0, 0.3333, 0.0, 0.0]

        #Test Amino Acids
        seq_aa = [Seq('PQRRRWQQ'),
                  Seq('PQXRRWQX'),
                  Seq('PQR--WQQ'),
                  Seq('PQX--WQX')]
        self.records_aa = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                           for i, s in enumerate(seq_aa, start=1)]
        self.weight_aa = [8, 6, 6, 4]

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print '<- %s() %.3f' % (self._testMethodName, t)

    #@unittest.skip('-> weightDNA() skipped\n')
    def test_weightDNA(self):
        scores = [presto.Sequence.weightDNA(x) for x in self.records_dna]
        print 'DNA Weight>'
        for x, s in zip(self.records_dna, scores):
            print '  %s> %s' % (x.id, s)

        self.assertEqual(scores, self.weight_dna)

    #@unittest.skip('-> weightAA() skipped\n')
    def test_weightAA(self):
        scores = [presto.Sequence.weightAA(x) for x in self.records_aa]
        print 'AA Weight>'
        for x, s in zip(self.records_aa, scores):
            print '  %s> %s' % (x.id, s)

        self.assertEqual(scores, self.weight_aa)

    #@unittest.skip('-> scoreSeqPair() skipped\n')
    def test_scoreSeqPair(self):
        pairs = list(itertools.combinations(self.records_dna, 2))
        scores = [presto.Sequence.scoreSeqPair(x, y) for x, y in pairs]
        print 'DNA Scores>'
        for (x, y), s in zip(pairs, scores):
            print '  %s/%s> %s' % (x.id, y.id, s)

        self.assertEqual([round(s[2], 4) for s in scores],
                         [round(s, 4) for s in self.error_dna])

