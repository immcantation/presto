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
from IgCore import getDNAScoreDict

class TestMaskPrimers(unittest.TestCase):

    def setUp(self):
        print '-> %s()' % self._testMethodName
        # Test Ns
        seq_n = [Seq('CCACGTTTTAGTAATTAATA'),
                 Seq('CCNCGTTTTAGTAATTAATA'),
                 Seq('GGGCGTTTTAGTAATTAATA'),
                 Seq('GGNNGTTTTACTAATTAATA')]
        self.records_n = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                          for i, s in enumerate(seq_n, start=1)]
        self.primers_n =  {1:'ACGTTT', 2:'GCCGTT'}

        #Test indels
        seq_indel = [Seq('CGGATCTTCTACTCCATACGTCCGTCAGTCGTGGATCTGATCTAGCTGCGCCTTTTTCTCAG'),
                     Seq('CGGATCTTCTACTCAAAACCGTCCTCAGTCGTGGATGTGGTCTAGCTGGGGCTGTTTCCCTG'),
                     Seq('GGTTTAAGTTAAGATAATACGTCCGTCAGTCGTGATGTGTTTTACATGGGGGCATCACCCAG'),
                     Seq('CAACCACATCTGTCCTCTAGAGAATCCCCTGAGAGCTCCGTTCCTCACCATGGACTGGACCT'),
                     Seq('CAACCACATGGGTCCTCTAGAGAATCCCCTGAGAGCTCCGTTCCTCACCATGGACTGGACCT'),
                     Seq('AGGTGAAGAAGCCTGGGGCCTCCGTGAAGGTCTCCTGCTCGGCTTCTGGATACGCCTTCACC'),
                     Seq('A--TGAAGAAGCCTGGGGCCTCCGTGAAGGTCTCCTGCTCGGCTTCTGGATACGCCTTCACC'),
                     Seq('TCTGTCCTCTAGAGAATCCCCTGAGAGCTCCGTTCCTCACCATGGACTGGACCTCAACCACA'),
                     Seq('-CTGTCCTCTAGAGAATCCCCTGAGAGCTCCGTTCCTCACCATGGACTGGACCTCAACCACA')]
        self.records_indel = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                              for i, s in enumerate(seq_indel, start=1)]
        self.primers_indel =  {1:'AATACGTCCGTCAGTCGTGGATGT', 2:'CATCTGTCCTC', 3:'GTGAAGAGCCTGG'}

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print '<- %s() %.3f' % (self._testMethodName, t)

    @unittest.skip('-> scorePrimers() skipped\n')
    def test_scorePrimers(self):
        score_dict=getDNAScoreDict(n_score=1, gap_score=0)
        align = [mod.scorePrimers(x, self.primers_n, start=1, score_dict=score_dict)
                 for x in self.records_n]
        for x in align:
            print '%s>' % x.seq.id
            print '   IN> %s' % x.seq.seq
            print '  SEQ> %s' % x.align_seq
            print '   PR> %s' % x.align_primer
            print '  ERR> %f' % x.error
        self.fail()

    #@unittest.skip('-> alignPrimers() skipped\n')
    def test_alignPrimers(self):
        # print 'TEST Ns>'
        # score_dict=getScoreDict(n_score=1, gap_score=0)
        # align_n = [mod.alignPrimers(x, self.primers_n, max_error=0.2, score_dict=score_dict)
        #            for x in self.records_n]
        # for x in align_n:
        #     print '  %s>' % x.seq.id
        #     print '   IN> %s' % x.seq.seq
        #     print '  SEQ> %s' % x.align_seq
        #     print '   PR> %s' % x.align_primer
        #     print '  ERR> %f' % x.error

        print 'TEST INDELS>'
        align_indel = [mod.alignPrimers(x, self.primers_indel, max_error=0.2)
                       for x in self.records_indel]
        for x in align_indel:
            print '  %s>' % x.seq.id
            print '   IN> %s' % x.seq.seq
            print '  SEQ> %s' % x.align_seq
            print '   PR> %s' % x.align_primer
            print '  ERR> %f' % x.error

        self.fail()