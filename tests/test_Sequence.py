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
        # Test DNA sequences
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

        #Test amino acid sequences
        seq_aa = [Seq('PQRRRWQQ'),
                  Seq('PQXRRWQX'),
                  Seq('PQR--WQQ'),
                  Seq('PQX--WQX')]
        self.records_aa = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                           for i, s in enumerate(seq_aa, start=1)]
        self.weight_aa = [8, 6, 6, 4]

        # Test DNA characters
        self.pairs_dna_chars = [('A', 'A'),
                                ('A', 'C'),
                                ('A', 'N'),
                                ('N', 'A'),
                                ('A', '-'),
                                ('-', 'A')]

        # Test amino acid characters
        self.pairs_aa_chars = [('P', 'P'),
                               ('P', 'Q'),
                               ('P', 'X'),
                               ('X', 'P'),
                               ('P', '-'),
                               ('-', 'P')]

        # Character pair scores for default case
        self.pairs_scores_def = [1, 0, 1, 1, 0, 0]
        # Character pair scores for symmetric case
        self.pairs_scores_sym = [1, 0, 1, 1, 1, 1]
        # Character pair scores for asymmetric case where N/gaps in character one are a mismatch
        self.pairs_scores_asym = [1, 0, 1, 0, 1, 0]


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
        scores = [presto.Sequence.scoreDNASeqPair(x, y) for x, y in pairs]
        print 'DNA Scores>'
        for (x, y), s in zip(pairs, scores):
            print '  %s/%s> %s' % (x.id, y.id, s)

        self.assertEqual([round(s[2], 4) for s in scores],
                         [round(s, 4) for s in self.error_dna])

    #@unittest.skip('-> scoreDNA() skipped\n')
    def test_scoreDNA(self):
        scores = [presto.Sequence.scoreDNA(a, b) for a, b in self.pairs_dna_chars]
        print 'Default DNA Scores>'
        for (a, b), s in zip(self.pairs_dna_chars, scores):
            print '  %s==%s> %s' % (a, b, s)

        self.assertSequenceEqual(self.pairs_scores_def, scores)

        scores = [presto.Sequence.scoreDNA(a, b, n_score=(1, 1), gap_score=(1, 1)) \
                  for a, b in self.pairs_dna_chars]
        print 'Symmetric DNA Scores>'
        for (a, b), s in zip(self.pairs_dna_chars, scores):
            print '  %s==%s> %s' % (a, b, s)

        self.assertSequenceEqual(self.pairs_scores_sym, scores)

        scores = [presto.Sequence.scoreDNA(a, b, n_score=(0, 1), gap_score=(0, 1)) \
                  for a, b in self.pairs_dna_chars]
        print 'Asymmetric DNA Scores>'
        for (a, b), s in zip(self.pairs_dna_chars, scores):
            print '  %s==%s> %s' % (a, b, s)

        self.assertSequenceEqual(self.pairs_scores_asym, scores)