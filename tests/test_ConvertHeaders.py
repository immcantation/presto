"""
Unit tests for ConvertHeaders
"""
# Future
from __future__ import absolute_import, division, print_function

# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import os
import sys
import time
import unittest

# Import script
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(test_dir, '..', 'bin'))
import ConvertHeaders


class TestConvertHeaders(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        # Test DNA
        self.desc_imgt = ['M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |',
						  'Z29978|IGHV1-69*07|Homo sapiens|F|V-REGION|1..233|233 nt|1| | | | |233+58=291|partial in 5\' and in 3\' | |',
						  'AB019439|IGHV(III)-22-2*01|Homo sapiens|P|L-PART1+V-EXON|168065..168119+168222..168262|96 nt|1| | | | |96+0=96| | |',
						  'AE000659|TRAV11*01|Homo sapiens|(P)|V-REGION|85121..85397|276 nt|1| | | | |276+42=318| | |']

        self.desc_genbank = ['gi|568336023|gb|CM000663.2| Homo sapiens chromosome 1, GRCh38 reference primary assembly',
                             'CM000663.2 Homo sapiens chromosome 1, GRCh38 reference primary assembly']

        self.desc_sra = ['SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36',
                         'SRR735691.1 GXGJ56Z01AE06X length=222',
                         'SRR1383326.1 1 length=250',
                         'SRR1383326.1.2 1 length=250']
        
        self.desc_454 = ['000034_0199_0169 length=437 uaccno=GNDG01201ARRCR',
                         'GXGJ56Z01AE06X length=222']

        self.desc_illumina = ['MISEQ:132:000000000-A2F3U:1:1101:14340:1555 1:N:0:ATCACG',
                              'HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1']

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    #@unittest.skip('-> convertGenericHeader() skipped\n')
    def test_convertGenericHeader(self):
        results = [ConvertHeaders.convertGenericHeader(x) for x in self.desc_genbank]
        for x in results:
            print('%s' % x)

        self.assertTrue(all([x is not None for x in results]))

        results = [ConvertHeaders.convertGenericHeader(x) for x in self.desc_imgt]
        for x in results:
            print('%s' % x)

        self.assertTrue(all([x is not None for x in results]))

    #@unittest.skip('-> convertIMGTHeader() skipped\n')
    def test_convertIMGTHeader(self):
        results = [ConvertHeaders.convertIMGTHeader(x) for x in self.desc_imgt]
        for x in results:
            print('%s' % x)

        self.assertTrue(all([x is not None for x in results]))

        results = [ConvertHeaders.convertIMGTHeader(x, simple=True) for x in self.desc_imgt]
        for x in results:
            print('%s' % x)

        self.assertTrue(all([x is not None for x in results]))

    #@unittest.skip('-> convertIlluminaHeader() skipped\n')
    def test_convertIlluminaHeader(self):
        results = [ConvertHeaders.convertIlluminaHeader(x) for x in self.desc_illumina]
        for x in results:
            print('%s' % x)

        self.assertTrue(all([x is not None for x in results]))

    #@unittest.skip('-> convertGenbankHeader() skipped\n')
    def test_convertGenbankHeader(self):
        results = [ConvertHeaders.convertGenbankHeader(x) for x in self.desc_genbank]
        for x in results:
            print('%s' % x)

        self.assertTrue(all([x is not None for x in results]))

    #@unittest.skip('-> convertSRAHeader() skipped\n')
    def test_convertSRAHeader(self):
        results = [ConvertHeaders.convertSRAHeader(x) for x in self.desc_sra]
        for x in results:
            print('%s' % x)

        self.assertTrue(all([x is not None for x in results]))

    #@unittest.skip('-> convert454Header() skipped\n')
    def test_convert454Header(self):
        results = [ConvertHeaders.convert454Header(x) for x in self.desc_454]
        for x in results:
            print('%s' % x)

        self.assertTrue(all([x is not None for x in results]))