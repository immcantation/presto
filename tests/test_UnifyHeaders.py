"""
Unit tests for UnifyHeaders
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import os
import sys
import time
import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Multiprocessing import SeqData
# Paths
test_path = os.path.dirname(os.path.realpath(__file__))

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import UnifyHeaders

class TestConvertHeaders(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        # Test records
        ann_list = [('AAAA', 'S1'), ('AAAA', 'S1'), ('AAAA', 'S2'), ('AAAA', 'S3'),
                    ('CCCC', 'S1'), ('CCCC', 'S2'),
                    ('GGGG', 'S2'), ('GGGG', 'S2'),
                    ('TTTT', 'S3')]
        id_list = ['SEQ%i|BARCODE=%s|SAMPLE=%s' % (i, x, y) for i, (x, y) in enumerate(ann_list, start=0)]
        seq_list = [SeqRecord(Seq('ACGTACGTACGT'), id=x, name=x, description=x) for x in id_list]

        # Data objects
        self.data_1 = SeqData('AAAA', seq_list[0:4])
        self.data_2 = SeqData('CCCC', seq_list[4:6])
        self.data_3 = SeqData('GGGG', seq_list[6:8])
        self.data_4 = SeqData('TTTT', seq_list[8])

        # Consensus results
        self.consensus_1 = ['SEQ0|BARCODE=AAAA|SAMPLE=S1',
                            'SEQ1|BARCODE=AAAA|SAMPLE=S1',
                            'SEQ2|BARCODE=AAAA|SAMPLE=S1',
                            'SEQ3|BARCODE=AAAA|SAMPLE=S1']
        self.consensus_2 = ['SEQ4|BARCODE=CCCC|SAMPLE=S1',
                            'SEQ5|BARCODE=CCCC|SAMPLE=S1']
        self.consensus_3 = ['SEQ6|BARCODE=GGGG|SAMPLE=S2',
                            'SEQ7|BARCODE=GGGG|SAMPLE=S2']

        # Delete results
        self.delete_1 = False
        self.delete_2 = False
        self.delete_3 = True

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    #@unittest.skip('-> consensusUnify() skipped\n')
    def test_consensusUnify(self):
        # data_1
        result = UnifyHeaders.consensusUnify(self.data_1, 'SAMPLE')
        for x, y in zip(self.data_1.data, result.results):
            print('ID> %s -> %s' % (x.id, y.id))
        self.assertListEqual(self.consensus_1, [x.id for x in result.results])

        # data_2
        result = UnifyHeaders.consensusUnify(self.data_2, 'SAMPLE')
        for x, y in zip(self.data_2.data, result.results):
            print('ID> %s -> %s' % (x.id, y.id))
        self.assertListEqual(self.consensus_2, [x.id for x in result.results])

        # data_3
        result = UnifyHeaders.consensusUnify(self.data_3, 'SAMPLE')
        for x, y in zip(self.data_3.data, result.results):
            print('ID> %s -> %s' % (x.id, y.id))
        self.assertListEqual(self.consensus_3, [x.id for x in result.results])


    #@unittest.skip('-> deletionUnify() skipped\n')
    def test_deletionUnify(self):
        # data_1
        result = UnifyHeaders.deletionUnify(self.data_1, 'SAMPLE')
        for x in self.data_1.data:
            print('ID> %s -> %s' % (x.id, result.valid))
        self.assertEqual(self.delete_1, result.valid)

        # data_2
        result = UnifyHeaders.deletionUnify(self.data_2, 'SAMPLE')
        for x in self.data_2.data:
            print('ID> %s -> %s' % (x.id, result.valid))
        self.assertEqual(self.delete_2, result.valid)

        # data_3
        result = UnifyHeaders.deletionUnify(self.data_3, 'SAMPLE')
        for x in self.data_3.data:
            print('ID> %s -> %s' % (x.id, result.valid))
        self.assertEqual(self.delete_3, result.valid)


if __name__ == '__main__':
    unittest.main()