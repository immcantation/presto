"""
Unit tests for MaskPrimers
"""
# Info
__author__    = 'Jason Anthony Vander Heiden'

# Imports
import os
import sys
import time
import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Sequence import getDNAScoreDict

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import MaskPrimers


class TestMaskPrimers(unittest.TestCase):

    def setUp(self):
        print '-> %s()' % self._testMethodName
        # Test Ns
        seq_n = [Seq('CCACGTTTTAGTAATTAATA'),
                 Seq('CCNCGTTTTAGTAATTAATA'),
                 Seq('GGGCGTTTTAGTAATTAATA'),
                 Seq('GGNNGTTTTACTAATTAATA'),
                 Seq('NNNNNNNNNACTAATTAATA'),
                 Seq('GGGATANNNACTAATTAATA'),
                 Seq('NNNNNNNNNNNNNNNNNNNN')]
        self.records_n = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                          for i, s in enumerate(seq_n, start=1)]
        self.primers_n =  {'PR1':'ACGTTT',
                           'PR2':'GCCGTT',
                           'PR3':'GGNATA'}

        # (primer name, error rate)
        self.align_n = [('PR1', 0.0/6),
                        ('PR1', 1.0/6),
                        ('PR1', 1.0/6),
                        ('PR1', 2.0/6),
                        ('PR3', 2.0/6),
                        ('PR3', 0.0/6),
                        (None, 1.0)]
        # Score primers with start=2
        self.score_n = [('PR1', 0.0/6),
                        ('PR1', 1.0/6),
                        ('PR1', 1.0/6),
                        ('PR1', 2.0/6),
                        ('PR1', 1.0),
                        ('PR3', 3.0/6),
                        ('PR1', 1.0)]

        #Test indels
        seq_indel = [Seq('CGGATCTTCTACTCCATACGTCCGTCAGTCGTGGATCTGATCTAGCTGCGCCTTTTTCTCAG'),
                     Seq('CGGATCTTCTACTCAAAACCGTCCTCAGTCGTGGATCTGGTCTAGCTGGGGCTGTTTCCCTG'),
                     Seq('GGTTTAAGTTAAGATAATACGTCCGTCAGTCGTGATCTGTTTTACATGGGGGCATCACCCAG'),
                     Seq('CAACCACATGGGTCCTCTAGAGAATCCCCTGAGAGCTCCGTTCCTCACCATGGACTGGACCT'),
                     Seq('CAACCACATGGTCCTCTAGAGAATCCCCTGAGAGCTCCGTTCCTCACCATGGACTGGACCTG'),
                     Seq('TCTGTCCTCTAGAGAATCCCCTGAGAGCTCCGTTCCTCACCATGGACTGGACCTCAACCACA'),
                     Seq('AGGTGAAGAAGCCTGGGGCCTCCGTGAAGGTCTCCTGCTCGGCTTCTGGATACGCCTTCACC'),
                     Seq('A--TGAAGAAGCCTGGGGCCTCCGTGAAGGTCTCCTGCTCGGCTTCTGGATACGCCTTCACC'),
                     Seq('ANNTGAAGAAGNNTGGGGCCTCCGTGAAGGTCTCCTGCTCGGCTTCTGGATACGCCTTCACC'),
                     Seq('--------------------------------------------------------------')]
        self.records_indel = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                              for i, s in enumerate(seq_indel, start=1)]
        self.primers_indel =  {'PR1':'AATACGTCCGTCAGTCGTGGATGT',
                               'PR2':'CATCTTCCTCTA',
                               'PR3':'GTGAAGAGCCTG'}

        # (primer name, error rate)
        self.align_indel = [('PR1', (2.0 + 0)/24),
                            ('PR1', (4.0 + 1)/24),
                            ('PR1', (2.0 + 1)/24),
                            ('PR2', (2.0 + 1)/12),
                            ('PR2', (2.0 + 0)/12),
                            ('PR2', (0.0 + 3)/12),
                            ('PR3', (0.0 + 1)/12),
                            ('PR3', (0.0 + 2)/12),
                            ('PR3', (3.0 + 1)/12),
                            (None, 1.0)]
        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print '<- %s() %.3f' % (self._testMethodName, t)

    #@unittest.skip('-> scorePrimers() skipped\n')
    def test_scorePrimers(self):
        score_dict=getDNAScoreDict(n_score=(0, 1), gap_score=(0, 0))
        align = [MaskPrimers.scorePrimers(x, self.primers_n, start=2, score_dict=score_dict)
                 for x in self.records_n]
        for x in align:
            print '  %s>' % x.seq.id
            print '      SEQ> %s' % x.seq.seq
            print '   PRIMER> %s' % x.primer
            print '  ALN-SEQ> %s' % x.align_seq
            print '   ALN-PR> %s' % x.align_primer
            print '    START> %s' % x.start
            print '      END> %s' % x.end
            print '     GAPS> %i' % x.gaps
            print '    ERROR> %f\n' % x.error

        self.assertSequenceEqual([(x, round(y, 4)) for x, y in self.score_n],
                                 [(x.primer, round(x.error, 4)) for x in align])

    #@unittest.skip('-> alignPrimers() skipped\n')
    def test_alignPrimers(self):
        score_dict=getDNAScoreDict(n_score=(0, 1), gap_score=(0, 0))

        # N character tests
        print 'TEST Ns>'
        align = [MaskPrimers.alignPrimers(x, self.primers_n, max_error=0.2, score_dict=score_dict)
                 for x in self.records_n]
        for x in align:
            print '  %s>' % x.seq.id
            print '      SEQ> %s' % x.seq.seq
            print '   PRIMER> %s' % x.primer
            print '  ALN-SEQ> %s' % x.align_seq
            print '   ALN-PR> %s' % x.align_primer
            print '    START> %s' % x.start
            print '      END> %s' % x.end
            print '     GAPS> %i' % x.gaps
            print '    ERROR> %f\n' % x.error

        self.assertSequenceEqual([(x, round(y, 4)) for x, y in self.align_n],
                                 [(x.primer, round(x.error, 4)) for x in align])

        # Indel tests
        print 'TEST INDELS>'
        align = [MaskPrimers.alignPrimers(x, self.primers_indel, max_error=0.2, gap_penalty=(1, 1))
                 for x in self.records_indel]
        for x in align:
            print '  %s>' % x.seq.id
            print '      SEQ> %s' % x.seq.seq
            print '   PRIMER> %s' % x.primer
            print '  ALN-SEQ> %s' % x.align_seq
            print '   ALN-PR> %s' % x.align_primer
            print '    START> %s' % x.start
            print '      END> %s' % x.end
            print '     GAPS> %i' % x.gaps
            print '    ERROR> %f\n' % x.error

        self.assertSequenceEqual([(x, round(y, 4)) for x, y in self.align_indel],
                                 [(x.primer, round(x.error, 4)) for x in align])


if __name__ == '__main__':
    unittest.main()