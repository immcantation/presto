"""
Unit tests for MaskPrimers
"""
# Imports
import time, unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from presto.Sequence import getScoreDict
from bin import MaskPrimers as script

# Info
__author__ = 'Jason Anthony Vander Heiden'


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
        score_dict=getScoreDict(n_score=1, gap_score=0)
        align = [script.scorePrimers(x, self.primers_n, start=1, score_dict=score_dict)
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
        align_indel = [script.alignPrimers(x, self.primers_indel, max_error=0.2)
                       for x in self.records_indel]
        for x in align_indel:
            print '  %s>' % x.seq.id
            print '   IN> %s' % x.seq.seq
            print '  SEQ> %s' % x.align_seq
            print '   PR> %s' % x.align_primer
            print '  ERR> %f' % x.error

        self.fail()