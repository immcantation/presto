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
from collections import OrderedDict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Sequence import getDNAScoreDict, localAlignment, scoreAlignment, extractAlignment, \
                            maskSeq, PrimerAlignment

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import MaskPrimers


class TestMaskPrimers(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        # Test Ns
        seq_n = [Seq('CCACGTTTTAGTAATTAATA'),
                 Seq('CCNCGTTTTAGTAATTAATA'),
                 Seq('GGGCGTTTTAGTAATTAATA'),
                 Seq('GGNNGTTTTACTAATTAATA'),
                 Seq('NNGCNNNNNACTAATTAATA'),
                 Seq('GGGATANNNACTAATTAATA'),
                 Seq('NNNNNNNNNNNNNNNNNNNN')]
        self.records_n = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                          for i, s in enumerate(seq_n, start=1)]
        self.primers_n =  OrderedDict([('PR1', 'CACGTTTT'),
                                       ('PR2', 'GGCGTTTT'),
                                       ('PR3', 'GGANATAA')])

        # (primer name, error rate)
        self.align_n = [('PR1', 0.0/8),
                        ('PR1', 1.0/8),
                        ('PR2', 0.0/8),
                        ('PR2', 2.0/8),
                        ('PR3', 3.0/8),
                        ('PR3', 3.0/8),
                        (None, 1.0)]
        # Score primers with start=1
        self.score_n = [('PR1', 0.0/8),
                        ('PR1', 1.0/8),
                        ('PR2', 0.0/8),
                        ('PR2', 2.0/8),
                        ('PR2', 6.0/8),
                        ('PR3', 3.0/8),
                        (None, 1.0)]
        # Extract primers with start=4 and length=8
        self.extract_n = ['GTTTTAGT',
                          'GTTTTAGT',
                          'GTTTTAGT',
                          'GTTTTACT',
                          'NNNNNACT',
                          'TANNNACT',
                          'NNNNNNNN']

        # Test indels
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
        self.primers_indel =  OrderedDict([('PR1', 'AATACGTCCGTCAGTCGTGGATGT'),
                                           ('PR2', 'CATCTTCCTCTA'),
                                           ('PR3', 'GTGAAGAGCCTG')])

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

        # Masking test case
        mask_seq = SeqRecord(Seq('CCACGTTTTAGTAATTAATA'), id='SEQ', description='SEQ')
        self.mask_primer = PrimerAlignment(mask_seq)
        self.mask_primer.primer = 'A'
        self.mask_primer.align_seq =    'CCACGTTTT'
        self.mask_primer.align_primer = '---CGTTTT'
        self.mask_primer.start = 3
        self.mask_primer.end = 9
        self.mask_primer.gaps = 0
        self.mask_primer.error = 0
        self.mask_primer.rev_primer = False
        self.mask_primer.valid = True

        self.mask_seq_cut = SeqRecord(Seq('AGTAATTAATA'), id='SEQ|PRIMER=A|BARCODE=CCA')
        self.mask_seq_mask = SeqRecord(Seq('NNNNNNAGTAATTAATA'), id='SEQ|PRIMER=A|BARCODE=CCA')
        self.mask_seq_trim = SeqRecord(Seq('CGTTTTAGTAATTAATA'), id='SEQ|PRIMER=A|BARCODE=CCA')
        self.mask_seq_tag = SeqRecord(Seq('CCACGTTTTAGTAATTAATA'), id='SEQ|PRIMER=A|BARCODE=CCA')
        self.mask_seq_cut_barcodelen = SeqRecord(Seq('AGTAATTAATA'), id='SEQ|PRIMER=A|BARCODE=CA')

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    #@unittest.skip('-> extractAlignment() skipped\n')
    def test_extractAlignment(self):
        align = [extractAlignment(x, start=4, length=8)
                 for x in self.records_n]
        for x in align:
            print('  %s>' % x.seq.id)
            print('      SEQ> %s' % x.seq.seq)
            print('  ALN-SEQ> %s' % x.align_seq)
            print('   ALN-PR> %s' % x.align_primer)
            print('   PRIMER> %s' % x.primer)
            print('    START> %s' % x.start)
            print('      END> %s' % x.end)
            print('     GAPS> %i' % x.gaps)
            print('    ERROR> %f\n' % x.error)

        self.assertListEqual([x for x in self.extract_n],
                             [str(x.primer) for x in align])

    #@unittest.skip('-> scoreAlignment() skipped\n')
    def test_scoreAlignment(self):
        score_dict = getDNAScoreDict(mask_score=(0, 1), gap_score=(0, 0))
        align = [scoreAlignment(x, self.primers_n, start=1, score_dict=score_dict)
                 for x in self.records_n]
        for x in align:
            print('  %s>' % x.seq.id)
            print('      SEQ> %s' % x.seq.seq)
            print('  ALN-SEQ> %s' % x.align_seq)
            print('   ALN-PR> %s' % x.align_primer)
            print('   PRIMER> %s' % x.primer)
            print('    START> %s' % x.start)
            print('      END> %s' % x.end)
            print('     GAPS> %i' % x.gaps)
            print('    ERROR> %f\n' % x.error)

        self.assertListEqual([(x, round(y, 4)) for x, y in self.score_n],
                             [(x.primer, round(x.error, 4)) for x in align])

    #@unittest.skip('-> localAlignment() skipped\n')
    def test_localAlignment(self):
        score_dict = getDNAScoreDict(mask_score=(0, 1), gap_score=(0, 0))

        # N character tests
        print('TEST Ns>')
        align = [localAlignment(x, self.primers_n, max_error=0.2, score_dict=score_dict)
                 for x in self.records_n]
        for x in align:
            print('  %s>' % x.seq.id)
            print('      SEQ> %s' % x.seq.seq)
            print('  ALN-SEQ> %s' % x.align_seq)
            print('   ALN-PR> %s' % x.align_primer)
            print('   PRIMER> %s' % x.primer)
            print('    START> %s' % x.start)
            print('      END> %s' % x.end)
            print('     GAPS> %i' % x.gaps)
            print('    ERROR> %f\n' % x.error)

        self.assertListEqual([(x, round(y, 4)) for x, y in self.align_n],
                             [(x.primer, round(x.error, 4)) for x in align])

        # Indel tests
        print('TEST INDELS>')
        align = [localAlignment(x, self.primers_indel, max_error=0.2, gap_penalty=(1, 1))
                 for x in self.records_indel]
        for x in align:
            print('  %s>' % x.seq.id)
            print('      SEQ> %s' % x.seq.seq)
            print('  ALN-SEQ> %s' % x.align_seq)
            print('   ALN-PR> %s' % x.align_primer)
            print('   PRIMER> %s' % x.primer)
            print('    START> %s' % x.start)
            print('      END> %s' % x.end)
            print('     GAPS> %i' % x.gaps)
            print('    ERROR> %f\n' % x.error)

        self.assertListEqual([(x, round(y, 4)) for x, y in self.align_indel],
                             [(x.primer, round(x.error, 4)) for x in align])

    #@unittest.skip('-> maskSeq() skipped\n')
    def test_maskSeq(self):
        print('TEST CUT>')
        result = maskSeq(self.mask_primer, mode='cut', barcode=True)
        print(' ID> %s' % result.id)
        print('SEQ> %s\n' % result.seq)
        self.assertEqual(str(self.mask_seq_cut.seq), str(result.seq))
        self.assertEqual(self.mask_seq_cut.id, result.id)

        print('TEST MASK>')
        result = maskSeq(self.mask_primer, mode='mask', barcode=True)
        print(' ID> %s' % result.id)
        print('SEQ> %s\n' % result.seq)
        self.assertEqual(str(self.mask_seq_mask.seq), str(result.seq))
        self.assertEqual(self.mask_seq_mask.id, result.id)

        print('TEST TRIM>')
        result = maskSeq(self.mask_primer, mode='trim', barcode=True)
        print(' ID> %s' % result.id)
        print('SEQ> %s\n' % result.seq)
        self.assertEqual(str(self.mask_seq_trim.seq), str(result.seq))
        self.assertEqual(self.mask_seq_trim.id, result.id)

        print('TEST TAG>')
        result = maskSeq(self.mask_primer, mode='tag', barcode=True)
        print(' ID> %s' % result.id)
        print('SEQ> %s\n' % result.seq)
        self.assertEqual(str(self.mask_seq_tag.seq), str(result.seq))
        self.assertEqual(self.mask_seq_tag.id, result.id)

        print('TEST CUT BARCODLEN>')
        result = maskSeq(self.mask_primer, mode='cut', barcode=True, barcode_length=2)
        print(' ID> %s' % result.id)
        print('SEQ> %s\n' % result.seq)
        self.assertEqual(str(self.mask_seq_cut_barcodelen.seq), str(result.seq))
        self.assertEqual(self.mask_seq_cut_barcodelen.id, result.id)

if __name__ == '__main__':
    unittest.main()