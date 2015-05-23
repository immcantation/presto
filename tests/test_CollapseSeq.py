"""
Unit tests for CollapseSeq
"""
# Imports
import time, unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from presto.Sequence import testSeqEqual

# Info
__author__    = 'Jason Anthony Vander Heiden'


class TestCollapseSeq(unittest.TestCase):

    def setUp(self):
        print '-> %s()' % self._testMethodName
        self.seq_list = [Seq("CCACGTTTTAGTAATTAATA"),
                         Seq("CCACGTTTTAGTAATTAATA"),
                         Seq("CCACGTTTTACTAATTAATA"),
                         Seq("CCACGTTTTANTAATTAATA")]
        self.rec_list = [SeqRecord(s, id='SEQ%i|COUNT=%i' % (i, i), name='SEQ%i|COUNT=%i' % (i, i), description='')
                         for i, s in enumerate(self.seq_list, start=1)]
        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    #@unittest.skip("-> testSeqEqual() skipped\n")
    def test_testSeqEqual(self):
        result = testSeqEqual(self.rec_list[0], self.rec_list[1])
        print result
        self.fail()
