"""
Unit tests for ClusterSets
"""
# Future
from __future__ import absolute_import, division, print_function

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.6'
__date__      = '2015.05.10'

# Imports
import time, unittest
from bin import ClusterSets as mod
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class TestBuildConsensus(unittest.TestCase):
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

    #@unittest.skip('-> runUClust() skipped\n')
    def test_runUClust(self):
        results = mod.runUClust(self.records_clust)
        print(results)

        self.assertEqual(sorted(self.results_clust), sorted(results))