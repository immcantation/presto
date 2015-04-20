"""
Unit tests for BuildConsensus
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.5'
__date__      = '2015.04.20'

# Imports
import itertools, time, unittest
import IgCore as core
import BuildConsensus as mod
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class TestBuildConsensus(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName
        # Test DNA
        seq_dna = [Seq('CGGCGTAA'),
                   Seq('CGNNGTAA'),
                   Seq('CGGC--AA'),
                   Seq('CGNN--AA'),
                   Seq('CGGC--AA')]
        seq_qual = [30, 40, 20, 20, 20, 20, 40, 40]
        self.records_dna = [SeqRecord(s, id='SEQ%i' % i, name='SEQ%i' % i, description='',
                                      letter_annotations={'phred_quality':seq_qual})
                            for i, s in enumerate(seq_dna, start=1)]

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print '<- %s() %.3f' % (self._testMethodName, t)

    #@unittest.skip('-> qualityConsensus() skipped\n')
    def test_qualityConsensus(self):
        result = core.qualityConsensus(self.records_dna)
        print result.seq

        result = core.qualityConsensus(self.records_dna, min_freq=0.2)
        print result.seq

        result = core.qualityConsensus(self.records_dna, min_qual=0)
        print result.seq

        result = core.qualityConsensus(self.records_dna, max_gap=0.5)
        print result.seq

        self.fail()

    #@unittest.skip('-> frequencyConsensus() skipped\n')
    def test_frequencyConsensus(self):
        result = core.frequencyConsensus(self.records_dna)
        print result.seq

        result = core.frequencyConsensus(self.records_dna, min_freq=0.2)
        print result.seq

        result = core.frequencyConsensus(self.records_dna, max_gap=0.5)
        print result.seq

        self.fail()
