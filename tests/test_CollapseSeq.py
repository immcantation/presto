"""
Unit tests for CollapseSeq
"""
# Imports
import time
import unittest
from itertools import combinations
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from presto.Sequence import testSeqEqual

# Info
__author__ = 'Jason Anthony Vander Heiden'


class TestCollapseSeq(unittest.TestCase):

    def setUp(self):
        print '-> %s()' % self._testMethodName
        seq_dna = [Seq("CCACGTTTTAGTAATTAATA"),
                   Seq("CCACGTTTTAGTAATTAATA"),
                   Seq("CCACGTTTTACTAATTAATA"),
                   Seq("CCACGTTTTANTAATTAATA")]
        self.records_dna = [SeqRecord(s, id='SEQ%i|COUNT=%i' % (i, i), name='SEQ%i|COUNT=%i' % (i, i), description='')
                            for i, s in enumerate(seq_dna, start=1)]
        # Make sequence pairs
        self.pairs_dna = list(combinations(self.records_dna, 2))

        # Equality
        # 1vs2, 1vs3, 1vs4,
        # 2vs3, 2vs4,
        # 3vs4,
        self.equal_dna = [True, False, True,
                          False, True,
                          True]

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    #@unittest.skip("-> testSeqEqual() skipped\n")
    def test_testSeqEqual(self):
        results = [testSeqEqual(x, y) for x, y in self.pairs_dna]
        print 'DNA Equality>'
        for (x, y), s in zip(self.pairs_dna, results):
            print '  %s> %s' % (x.id, x.seq)
            print '  %s> %s' % (y.id, y.seq)
            print '         EQUAL> %s\n' % s

        self.assertSequenceEqual(self.equal_dna, results)
