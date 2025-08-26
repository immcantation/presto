"""
Unit tests for Sequence module
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import math
import os
import sys
import time
import unittest
from itertools import combinations
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Sequence import getDNAScoreDict, scoreDNA, scoreSeqPair, weightSeq, \
                            calculateSetError, meanQuality, filterQuality
from presto.Multiprocessing import SeqData

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))


class TestSequence(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        # Test DNA sequences
        seq_dna = [Seq('CGGCGTAA'),
                   Seq('CGNNGTAG'),
                   Seq('CGGC--AA'),
                   Seq('CGNN--AG'),
                   Seq('NNNNNNNN'),
                   Seq('NNNNNNAA'),
                   Seq('--------'),
                   Seq('CG------')]
        ref_dna = [Seq('CGGCGTAA'),
                   Seq('NNNNNNNN')]
        qual_dna = [[30, 30, 20, 20, 30, 30, 40, 40],
                    [30, 30,  0,  0, 30, 30, 20, 20],
                    [ 0,  0,  0,  0,  0,  0,  0,  0]]
        self.records_dna = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                            for i, s in enumerate(seq_dna, start=1)]
        self.reference_dna = [SeqRecord(s, id='REF%i' % i, name='REF%i' % i, description='')
                              for i, s in enumerate(ref_dna, start=1)]
        self.phred_dna = [SeqRecord(seq_dna[i], id='SEQ%i' % i, name='SEQ%i' % i, description='',
                                    letter_annotations={'phred_quality':qual_dna[i]})
                          for i in range(len(qual_dna))]

        # Mean quality
        self.qual_mean = list()
        for qual in qual_dna:
            m = math.floor(sum(qual) / len(qual))
            self.qual_mean.append(m)

        # Make sequence pairs
        self.seq_pairs = list(combinations(self.records_dna, 2))

        # Weights
        self.weight_dna_mask = [8, 6, 8, 6, 0, 2, 8, 8]

        # Set error rates against reference
        self.error_ref = [1.0 - float(8 + 5 + 6 + 3 + 0 + 2 + 0 + 2) / (8 + 6 + 8 + 6 + 0 + 2 + 8 + 8),
                          1.0]

        # Pairwise rrror rates
        # scoreSeqPair returns (score, minimum weight, error rate)
        # Asymmetric case is with score_dict = getDNAScoreDict(n_score=(0, 1), gap_score=(0, 1))
        # 1vs2, 1vs3, 1vs4, 1vs5, 1vs6, 1v7, 1v8
        # 2vs3, 2vs4, 2vs5, 2vs6, 2v7, 2v8,
        # 3vs4, 3vs5, 3vs6, 3vs7, 3v8,
        # 4vs5, 4vs6, 4vs7, 4vs8,
        # 5vs6, 5vs7, 5v8,
        # 6vs7, 6vs8,
        # 7vs8
        self.error_dna_def = [1.0/8, 2.0/8, 3.0/8, 0.0/8.0, 0.0/8.0, 8.0/8, 6.0/8,
                              3.0/8, 2.0/8, 0.0/8, 1.0/8, 8.0/8, 6.0/8,
                              1.0/8, 2.0/8, 2.0/8, 6.0/8, 4.0/8,
                              2.0/8, 3.0/8, 6.0/8, 4.0/8,
                              0.0/8, 8.0/8, 6.0/8,
                              8.0/8, 6.0/8,
                              2.0/8]
        self.error_dna_asym = [1.0/8, 0.0/8, 1.0/8, 0.0/8.0, 0.0/8.0, 0.0/8.0, 0.0/8.0,
                               3.0/8, 2.0/8, 2.0/8, 3.0/8, 2.0/8, 2.0/8,
                               3.0/8, 2.0/8, 2.0/8, 2.0/8, 2.0/8,
                               4.0/8, 5.0/8, 4.0/8, 4.0/8,
                               8.0/8, 8.0/8, 8.0/8,
                               6.0/8, 6.0/8,
                               8.0/8]
        self.error_dna_mask = [1.0/6, 2.0/8, 3.0/6, 1.0, 0.0/2, 8.0/8, 6.0/8,
                               3.0/6, 2.0/6, 1.0, 1.0/2, 6.0/6, 4.0/6,
                               1.0/6, 1.0, 0.0/2, 6.0/8, 4.0/8,
                               1.0, 1.0/2, 4.0/6, 2.0/6,
                               1.0, 1.0, 1.0,
                               2.0/2, 2.0/2,
                               2.0/8]

        # Test amino acids sequences
        seq_aa = [Seq('PQRRRWQQ'),
                  Seq('PQXRRWQX'),
                  Seq('PQR--WQQ'),
                  Seq('PQX--WQX')]
        self.records_aa = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                           for i, s in enumerate(seq_aa, start=1)]
        self.weight_aa_mask = [8, 6, 8, 6]

        # Test DNA characters
        self.pairs_dna_chars = [('A', 'A'),
                                ('A', 'T'),
                                ('A', 'R'),
                                ('U', 'T'),
                                ('A', 'N'),
                                ('N', 'A'),
                                ('A', '-'),
                                ('-', 'A')]

        # Test amino acid characters
        self.pairs_aa_chars = [('P', 'P'),
                               ('P', 'Q'),
                               ('R', 'B'),
                               ('R', 'N'),
                               ('P', 'X'),
                               ('X', 'P'),
                               ('P', '-'),
                               ('-', 'P')]

        # Character pair scores for default case
        self.pairs_scores_def = [1, 0, 1, 1, 1, 1, 0, 0]
        # Character pair scores for symmetric case
        self.pairs_scores_sym = [1, 0, 1, 1, 1, 1, 1, 1]
        # Character pair scores for asymmetric case where N/gaps in character one are a mismatch
        self.pairs_scores_asym = [1, 0, 1, 1, 1, 0, 1, 0]

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    #@unittest.skip('-> weightSeq() skipped\n')
    def test_weightDNA(self):
        # DNA weights
        ignore_chars = set(['n', 'N'])
        weights = [weightSeq(x, ignore_chars=ignore_chars) for x in self.records_dna]
        print('DNA Weight>')
        for x, s in zip(self.records_dna, weights):
            print('  %s> %s' % (x.id, s))

        self.assertSequenceEqual(weights, self.weight_dna_mask)

        # Amino acid weights
        ignore_chars = set(['x', 'X'])
        weights = [weightSeq(x, ignore_chars=ignore_chars) for x in self.records_aa]
        print('AA Weight>')
        for x, s in zip(self.records_dna, weights):
            print('  %s> %s' % (x.id, s))

        self.assertSequenceEqual(weights, self.weight_aa_mask)

    #@unittest.skip('-> scoreSeqPair() skipped\n')
    def test_scoreSeqPair(self):
        # Default scoring
        scores = [scoreSeqPair(x, y) for x, y in self.seq_pairs]
        print('Default DNA Scores>')
        for (x, y), s in zip(self.seq_pairs, scores):
            print('    %s> %s' % (x.id, x.seq))
            print('    %s> %s' % (y.id, y.seq))
            print('   SCORE> %i' % s[0])
            print('  WEIGHT> %i' % s[1])
            print('   ERROR> %f\n' % s[2])

        self.assertSequenceEqual([round(s[2], 4) for s in scores],
                                 [round(s, 4) for s in self.error_dna_def])

        # Asymmetric scoring without position masking
        score_dict = getDNAScoreDict(mask_score=(0, 1), gap_score=(0, 1))
        scores = [scoreSeqPair(x, y, score_dict=score_dict) \
                  for x, y in self.seq_pairs]
        print('Asymmetric DNA Scores>')
        for (x, y), s in zip(self.seq_pairs, scores):
            print('    %s> %s' % (x.id, x.seq))
            print('    %s> %s' % (y.id, y.seq))
            print('   SCORE> %i' % s[0])
            print('  WEIGHT> %i' % s[1])
            print('   ERROR> %f\n' % s[2])

        self.assertSequenceEqual([round(s[2], 4) for s in scores],
                                 [round(s, 4) for s in self.error_dna_asym])

        # Symmetric scoring with N positions excluded
        ignore_chars = set(['n', 'N'])
        scores = [scoreSeqPair(x, y, ignore_chars=ignore_chars) \
                  for x, y in self.seq_pairs]
        print('Masked DNA Scores>')
        for (x, y), s in zip(self.seq_pairs, scores):
            print('    %s> %s' % (x.id, x.seq))
            print('    %s> %s' % (y.id, y.seq))
            print('   SCORE> %i' % s[0])
            print('  WEIGHT> %i' % s[1])
            print('   ERROR> %f\n' % s[2])

        self.assertSequenceEqual([round(s[2], 4) for s in scores],
                                 [round(s, 4) for s in self.error_dna_mask])

    #@unittest.skip('-> scoreDNA() skipped\n')
    def test_scoreDNA(self):
        scores = [scoreDNA(a, b) for a, b in self.pairs_dna_chars]
        print('Default DNA Scores>')
        for (a, b), s in zip(self.pairs_dna_chars, scores):
            print('  %s==%s> %s' % (a, b, s))

        self.assertSequenceEqual(self.pairs_scores_def, scores)

        scores = [scoreDNA(a, b, mask_score=(1, 1), gap_score=(1, 1)) \
                  for a, b in self.pairs_dna_chars]
        print('Symmetric DNA Scores>')
        for (a, b), s in zip(self.pairs_dna_chars, scores):
            print('  %s==%s> %s' % (a, b, s))

        self.assertSequenceEqual(self.pairs_scores_sym, scores)

        scores = [scoreDNA(a, b, mask_score=(0, 1), gap_score=(0, 1)) \
                  for a, b in self.pairs_dna_chars]
        print('Asymmetric DNA Scores>')
        for (a, b), s in zip(self.pairs_dna_chars, scores):
            print('  %s==%s> %s' % (a, b, s))

        self.assertSequenceEqual(self.pairs_scores_asym, scores)

    #@unittest.skip('-> calculateSetError() skipped\n')
    def test_calculateSetError(self):
        for r, e in zip(self.reference_dna, self.error_ref):
            result = calculateSetError(self.records_dna, r)
            print('REF>', r.seq, result)
            self.assertAlmostEqual(e, result, places=4)

    def test_meanQuality(self):
        for r, q in zip(self.phred_dna, self.qual_mean):
            result = meanQuality(r.letter_annotations['phred_quality'])
            print('RESULT>', result, r.letter_annotations['phred_quality'])
            self.assertEqual(q, result)

    def test_filterQuality(self):
        seq_data = [SeqData(i, x) for i, x in enumerate(self.phred_dna)]
        for r, q in zip(seq_data, self.qual_mean):
            result = filterQuality(r, min_qual=20, inner=False)
            print('RESULT>', result.valid, result.log['QUALITY'])
            self.assertEqual(q, result.log['QUALITY'])
            self.assertEqual(q >= 20, result.valid)

if __name__ == '__main__':
    unittest.main()