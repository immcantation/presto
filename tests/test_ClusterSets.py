"""
Unit tests for ClusterSets
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

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import ClusterSets


class TestClusterSets(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        # Test DNA
        seq_clust = [Seq('CGGCGTAACGGCGTAATTTTTTTTTTCGGCGTAACGGCGTAACGGCGTAACGGCGTAA'),
                     Seq('CCGGGTAACGGCGTAATTTTTTTTTTCGGCGTAACGGCGTAACGGCGTAACGGCGTAA'),
                     Seq('CGGC--TACGGCGTAATTTTTTTTTTCGGCGTAACGGCGTAACGGCGTAACGGCGTAA'),
                     Seq('CGGG--TACGGCGTAACCCCCCCCCCCGGCGTAACGGCGTAACGGCGTAACGGCGTAA'),
                     Seq('CGGC--AACGGCGTAACCCCCCCCCCCGGCGTAACGGCGTAACGGCGTAACGGCGTAA')]
        self.records_clust = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='')
                              for i, s in enumerate(seq_clust, start=1)]
        self.results_clust = {1:['SEQ1', 'SEQ2', 'SEQ3'], 2:['SEQ4','SEQ5']}

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    @unittest.skip('-> runUClust() skipped\n')
    def test_runUClust(self):
        results = ClusterSets.runUClust(self.records_clust)
        print(results)

        self.assertEqual(sorted(self.results_clust), sorted(results))


if __name__ == '__main__':
    unittest.main()