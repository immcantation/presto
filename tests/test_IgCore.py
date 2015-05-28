"""
Unit tests for IgCore
"""
__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.5'
__date__      = '2015.04.06'

# Imports
import itertools, time, unittest
import IgCore as mod
from collections import OrderedDict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class TestIgCore(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName
        # Test DNA
        seq_dna = [Seq('CGGCGTAA'),
                   Seq('CGNNGTAA'),
                   Seq('CGGC--AA'),
                   Seq('CGNN--AA'),
                   Seq('NNNNNNNN'),
                   Seq('NNNNNNAA'),
                   Seq('--------'),
                   Seq('CG------')]
        self.records_dna = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                            for i, s in enumerate(seq_dna, start=1)]
        self.weight_dna = [8, 6, 6, 4, 1, 2, 8, 8]
        # scoreSeqPair returns (score, minimum weight, error rate)
        # 1vs2, 1vs3, 1vs4, 1vs5, 1vs6, 1v7, 1v8
        # 2vs3, 2vs4, 2vs5, 2vs6, 2v7, 2v8,
        # 3vs4, 3vs5, 3vs6, 3vs7, 3v8,
        # 4vs5, 4vs6, 4vs7, 4vs8,
        # 5vs6, 5vs7, 5v8,
        # 6vs7, 6vs8,
        # 7vs8
        self.error_dna = [0.0/6, 2.0/8, 2.0/6, 1.0, 0.0/2, 8.0/8, 6.0/8,
                          2.0/6, 2.0/6, 1.0, 0.0/2, 6.0/6, 4.0/6,
                          0.0/6, 1.0, 0.0/2, 6.0/8, 4.0/8,
                          1.0, 0.0/2, 4.0/6, 2.0/6,
                          1.0, 1.0, 1.0,
                          1.0, 2.0/2,
                          2.0/8]

        #Test Amino Acids
        seq_aa = [Seq('PQRRRWQQ'),
                  Seq('PQXRRWQX'),
                  Seq('PQR--WQQ'),
                  Seq('PQX--WQX')]
        self.records_aa = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                           for i, s in enumerate(seq_aa, start=1)]
        self.weight_aa = [8, 6, 6, 4]

        # Annotation dictionaries
        self.ann_dict_1 = OrderedDict([('ID', 'SEQ1'), ('TEST1', 'A,B'), ('TEST2', [1, 2])])
        self.ann_dict_2 = OrderedDict([('ID', 'SEQ2'), ('TEST1', 'C,C'), ('TEST2', 3)])

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
        scores = [mod.weightDNA(x) for x in self.records_dna]
        print 'DNA Weight>'
        for x, s in zip(self.records_dna, scores):
            print '  %s> %s' % (x.id, s)

        self.assertEqual(scores, self.weight_dna)

    #@unittest.skip('-> weightAA() skipped\n')
    def test_weightAA(self):
        scores = [mod.weightAA(x) for x in self.records_aa]
        print 'AA Weight>'
        for x, s in zip(self.records_aa, scores):
            print '  %s> %s' % (x.id, s)

        self.assertEqual(scores, self.weight_aa)

    #@unittest.skip('-> scoreSeqPair() skipped\n')
    def test_scoreSeqPair(self):
        ignore_chars = set(['N', 'n'])
        pairs = list(itertools.combinations(self.records_dna, 2))
        scores = [mod.scoreSeqPair(x, y, ignore_chars=ignore_chars) for x, y in pairs]
        print 'DNA Scores>'
        for (x, y), s in zip(pairs, scores):
            print '  %s> %s' % (x.id, x.seq)
            print '  %s> %s' % (y.id, y.seq)
            print ' SCORE> %i' % s[0]
            print 'WEIGHT> %i' % s[1]
            print ' ERROR> %f\n' % s[2]

        self.assertEqual([round(s[2], 4) for s in scores],
                         [round(s, 4) for s in self.error_dna])

    #@unittest.skip('-> mergeAnnotation() skipped\n')
    def test_mergeAnnotation(self):
        result = mod.mergeAnnotation(self.ann_dict_1, self.ann_dict_2)
        print result
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('A,B,C,C', result['TEST1'])
        self.assertEqual('1,2,3', result['TEST2'])

        result = mod.mergeAnnotation(self.ann_dict_1, self.ann_dict_2, prepend=True)
        print result
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('C,C,A,B', result['TEST1'])
        self.assertEqual('3,1,2', result['TEST2'])

    #@unittest.skip('-> renameAnnotation() skipped\n')
    def test_renameAnnotation(self):
        result = mod.renameAnnotation(self.ann_dict_1, 'TEST1', 'RENAMED1')
        print result
        self.assertEqual('A,B', result['RENAMED1'])
        self.assertEqual(None, result.get('TEST1', None))

        result = mod.renameAnnotation(self.ann_dict_1, 'TEST1', 'RENAMED2')
        print result
        self.assertEqual('A,B', result['RENAMED2'])
        self.assertEqual(None, result.get('TEST1', None))

    #@unittest.skip('-> collapseAnnotation() skipped\n')
    def test_collapseAnnotation(self):
        result = mod.collapseAnnotation(self.ann_dict_1, 'first')
        print result
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('A', result['TEST1'])
        self.assertEqual(1, result['TEST2'])

        result = mod.collapseAnnotation(self.ann_dict_1, 'last')
        print result
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('B', result['TEST1'])
        self.assertEqual(2, result['TEST2'])

        result = mod.collapseAnnotation(self.ann_dict_1, 'set')
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual(['A' ,'B'], result['TEST1'])
        self.assertEqual([1, 2], result['TEST2'])
        print result

        result = mod.collapseAnnotation(self.ann_dict_1, 'cat')
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('AB', result['TEST1'])
        self.assertEqual('12', result['TEST2'])
        print result

        # TODO:  may be better for min/max/sum to return float/int
        result = mod.collapseAnnotation(self.ann_dict_1, 'min', fields='TEST2')
        self.assertEqual('1', result['TEST2'])
        print result

        result = mod.collapseAnnotation(self.ann_dict_1, 'max', fields='TEST2')
        self.assertEqual('2', result['TEST2'])
        print result

        result = mod.collapseAnnotation(self.ann_dict_1, 'sum', fields='TEST2')
        self.assertEqual('3', result['TEST2'])
        print result

     #@unittest.skip('-> scoreDNA() skipped\n')
    def test_scoreDNA(self):
        scores = [mod.scoreDNA(a, b) for a, b in self.pairs_dna_chars]
        print 'Default DNA Scores>'
        for (a, b), s in zip(self.pairs_dna_chars, scores):
            print '  %s==%s> %s' % (a, b, s)

        self.assertSequenceEqual(self.pairs_scores_def, scores)

        scores = [mod.scoreDNA(a, b, n_score=(1, 1), gap_score=(1, 1)) \
                  for a, b in self.pairs_dna_chars]
        print 'Symmetric DNA Scores>'
        for (a, b), s in zip(self.pairs_dna_chars, scores):
            print '  %s==%s> %s' % (a, b, s)

        self.assertSequenceEqual(self.pairs_scores_sym, scores)

        scores = [mod.scoreDNA(a, b, n_score=(0, 1), gap_score=(0, 1)) \
                  for a, b in self.pairs_dna_chars]
        print 'Asymmetric DNA Scores>'
        for (a, b), s in zip(self.pairs_dna_chars, scores):
            print '  %s==%s> %s' % (a, b, s)

        self.assertSequenceEqual(self.pairs_scores_asym, scores)