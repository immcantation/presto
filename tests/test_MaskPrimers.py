"""
Unit tests for MaskPrimers
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.5'
__date__      = '2015.02.20'

# Imports
import time, unittest
import MaskPrimers as mod
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from IgCore import getScoreDict

class TestMaskPrimers(unittest.TestCase):

    def setUp(self):
        print '-> %s()' % self._testMethodName
        seq_list = [Seq('CCACGTTTTAGTAATTAATA'),
                    Seq('CCNCGTTTTAGTAATTAATA'),
                    Seq('GGGCGTTTTAGTAATTAATA'),
                    Seq('GGNNGTTTTACTAATTAATA')]
        self.records = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                         for i, s in enumerate(seq_list, start=1)]
        self.primers =  {1:'ACGTTT', 2:'GCCGTT'}
        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print '<- %s() %.3f' % (self._testMethodName, t)

    #@unittest.skip('-> scorePrimers() skipped\n')
    def test_scorePrimers(self):
        score_dict=getScoreDict(n_score=1, gap_score=0)
        results = [mod.scorePrimers(x, self.primers, start=1, score_dict=score_dict)
                   for x in self.records]
        for x in results:
            print '%s>' % x['seq'].id
            print '  SEQ> %s' % x['seq'].seq
            print '  ALN> %s' % x['align']
            print '  ERR> %f' % x['error']
        self.fail()

    #@unittest.skip('-> alignPrimers() skipped\n')
    def test_alignPrimers(self):
        score_dict=getScoreDict(n_score=1, gap_score=0)
        results = [mod.alignPrimers(x, self.primers, max_error=0.2, score_dict=score_dict)
                   for x in self.records]
        for x in results:
            print '%s>' % x['seq'].id
            print '  SEQ> %s' % x['seq'].seq
            print '  ALN> %s' % x['align']
            print '  ERR> %f' % x['error']
        self.fail()