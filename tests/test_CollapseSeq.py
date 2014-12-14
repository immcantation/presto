"""
Unit tests for CollapseSeq
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.5'
__date__      = '2014.12.12'

# Imports
import unittest
import CollapseSeq as mod
import time
from timeit import timeit
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from IgCore import testSeqEqual


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
