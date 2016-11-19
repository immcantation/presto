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
from presto.IO import getFileType, readSeqFile

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


class TestIO(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        # Annotation dictionaries
        self.germ_file = os.path.join(data_path, 'human_igv_trv.fasta')
        self.missing_file = os.path.join(data_path, 'human_igv_trv.bob')

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    #@unittest.skip('-> getFileType() skipped\n')
    def test_getFileType(self):
        result = getFileType(self.germ_file)
        print(result)
        self.assertEqual('fasta', result)

        result = getFileType(self.missing_file)
        print(result)
        self.assertEqual('fasta', result)

    #@unittest.skip('-> readSeqFile() skipped\n')
    def test_getFileType(self):
        pass
        #result = readSeqFile(self.germ_file)
        #self.assertRaises(SystemExit, readSeqFile(self.missing_file))


if __name__ == '__main__':
    unittest.main()