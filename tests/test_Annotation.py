"""
Unit tests for Annotation module
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import os
import sys
import time
import unittest
from collections import OrderedDict

# Presto imports
from presto.Annotation import collapseAnnotation, mergeAnnotation, renameAnnotation, getCoordKey

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))


class TestAnnotation(unittest.TestCase):
    """
    Notes:
      illumina  @MISEQ:132:000000000-A2F3U:1:1101:14340:1555 2:N:0:ATCACG
                @HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1
                @MS6_33112:1:1101:18371:1066/1
      sra       @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
                @SRR001666.1.2 1 length=250
      454       @000034_0199_0169 length=437 uaccno=GNDG01201ARRCR
      pesto     @AATCGGATTTGC|COUNT=2|PRIMER=IGHJ_RT|PRFREQ=1.0
    """
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        # Annotation dictionaries
        self.ann_dict_1 = OrderedDict([('ID', 'SEQ1'),
                                       ('TEST1', 'A,B'),
                                       ('TEST2', [1, 2]),
                                       ('TEST3', ['First','Second'])])
        self.ann_dict_2 = OrderedDict([('ID', 'SEQ2'),
                                       ('TEST1', 'C,C'),
                                       ('TEST2', 3)])

        # SRA, ENA
        self.desc_sra = ['SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36',
                         'SRR735691.1 GXGJ56Z01AE06X length=222',
                         'SRR1383326.1 1 length=250',
                         'SRR1383326.1.2 1 length=250',
                         'ERR346596.6 BS-DSFCONTROL04:4:000000000-A3F0Y:1:1101:13220:1649/1']
        self.id_sra = ['SRR001666.1',
                       'SRR735691.1',
                       'SRR1383326.1',
                       'SRR1383326.1',
                       'ERR346596.6']

        # 454
        self.desc_454 = ['000034_0199_0169 length=437 uaccno=GNDG01201ARRCR',
                         'GXGJ56Z01AE06X length=222']
        self.id_454 = ['000034_0199_0169',
                       'GXGJ56Z01AE06X']

        # Illumina, Solexa
        self.desc_illumina = ['MISEQ:132:000000000-A2F3U:1:1101:14340:1555 1:N:0:ATCACG',
                              'HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1',
                              'MS6_33112:1:1101:18371:1066/1']
        self.id_illumina = ['MISEQ:132:000000000-A2F3U:1:1101:14340:1555',
                            'HWI-EAS209_0006_FC706VJ:5:58:5894:21141',
                            'MS6_33112:1:1101:18371:1066']

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    #@unittest.skip('-> mergeAnnotation() skipped\n')
    def test_mergeAnnotation(self):
        result = mergeAnnotation(self.ann_dict_1, self.ann_dict_2)
        print(result)
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('A,B,C,C', result['TEST1'])
        self.assertEqual('1,2,3', result['TEST2'])

        result = mergeAnnotation(self.ann_dict_1, self.ann_dict_2, prepend=True)
        print(result)
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('C,C,A,B', result['TEST1'])
        self.assertEqual('3,1,2', result['TEST2'])

    #@unittest.skip('-> renameAnnotation() skipped\n')
    def test_renameAnnotation(self):
        result = renameAnnotation(self.ann_dict_1, 'TEST1', 'RENAMED1')
        print(result)
        self.assertEqual('A,B', result['RENAMED1'])
        self.assertEqual(None, result.get('TEST1', None))

        result = renameAnnotation(self.ann_dict_1, 'TEST1', 'RENAMED2')
        print(result)
        self.assertEqual('A,B', result['RENAMED2'])
        self.assertEqual(None, result.get('TEST1', None))

        result = renameAnnotation(self.ann_dict_1, 'TEST2', 'TEST1')
        print(result)
        self.assertEqual('A,B,1,2', result['TEST1'])
        self.assertEqual(None, result.get('TEST2', None))

        result = renameAnnotation(self.ann_dict_1, 'TEST2', 'TEST1')
        result = renameAnnotation(result, 'TEST3', 'TEST1')
        print(result)
        self.assertEqual('A,B,1,2,First,Second', result['TEST1'])
        self.assertEqual(None, result.get('TEST2', None))
        self.assertEqual(None, result.get('TEST3', None))

    #@unittest.skip('-> collapseAnnotation() skipped\n')
    def test_collapseAnnotation(self):
        result = collapseAnnotation(self.ann_dict_1, 'first')
        print(result)
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('A', result['TEST1'])
        self.assertEqual(1, result['TEST2'])

        result = collapseAnnotation(self.ann_dict_1, 'last')
        print(result)
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('B', result['TEST1'])
        self.assertEqual(2, result['TEST2'])

        result = collapseAnnotation(self.ann_dict_1, 'set')
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual(['A' ,'B'], result['TEST1'])
        self.assertEqual([1, 2], result['TEST2'])
        print(result)

        result = collapseAnnotation(self.ann_dict_1, 'cat')
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('AB', result['TEST1'])
        self.assertEqual('12', result['TEST2'])
        print(result)

        # TODO:  may be better for min/max/sum to return float/int
        result = collapseAnnotation(self.ann_dict_1, 'min', fields='TEST2')
        self.assertEqual('1', result['TEST2'])
        print(result)

        result = collapseAnnotation(self.ann_dict_1, 'max', fields='TEST2')
        self.assertEqual('2', result['TEST2'])
        print(result)

        result = collapseAnnotation(self.ann_dict_1, 'sum', fields='TEST2')
        self.assertEqual('3', result['TEST2'])
        print(result)

    #@unittest.skip('-> getCoordKey() skipped\n')
    def test_getCoordKey(self):
        result = [getCoordKey(x, coord_type='illumina') for x in self.desc_illumina]
        self.assertListEqual(result, self.id_illumina)
        print(result)

        result = [getCoordKey(x, coord_type='sra') for x in self.desc_sra]
        self.assertListEqual(result, self.id_sra)
        print(result)

        result = [getCoordKey(x, coord_type='454') for x in self.desc_454]
        self.assertListEqual(result, self.id_454)
        print(result)



if __name__ == '__main__':
    unittest.main()