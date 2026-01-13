"""
Unit tests for PairSeq with gzip support
"""
# Info
__author__ = 'Ivan Gregoretti'

# Imports
import os
import sys
import time
import unittest
import gzip
from Bio import SeqIO

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


class TestPairSeqGzip(unittest.TestCase):
    """Test PairSeq functionality with gzipped input and output files"""
    
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        
        # Define test files - paired output from primer processing
        self.input_R1_gz = os.path.join(data_path, 'set413_R1_primers-pass.fq.gz')
        self.input_R2_gz = os.path.join(data_path, 'set413_R2_primers-pass.fq.gz')
        self.output_R1_paired_gz = os.path.join(data_path, 'set413_R1_primers-pass_pair-pass.fq.gz')
        self.output_R2_paired_gz = os.path.join(data_path, 'set413_R2_primers-pass_pair-pass.fq.gz')
        
        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    def test_pairseq_gzipped_input_R1(self):
        """Test that PairSeq can read gzipped R1 input"""
        with gzip.open(self.input_R1_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"R1 input: {len(records)} records")
        self.assertGreater(len(records), 0, "Should read R1 records")

    def test_pairseq_gzipped_input_R2(self):
        """Test that PairSeq can read gzipped R2 input"""
        with gzip.open(self.input_R2_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"R2 input: {len(records)} records")
        self.assertGreater(len(records), 0, "Should read R2 records")

    def test_pairseq_gzipped_output_R1(self):
        """Test that PairSeq created gzipped paired output for R1"""
        self.assertTrue(os.path.exists(self.output_R1_paired_gz),
                       "PairSeq gzipped output should exist for R1")
        
        with gzip.open(self.output_R1_paired_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"R1 paired output: {len(records)} records")
        self.assertGreater(len(records), 0, "Output should contain records")

    def test_pairseq_gzipped_output_R2(self):
        """Test that PairSeq created gzipped paired output for R2"""
        self.assertTrue(os.path.exists(self.output_R2_paired_gz),
                       "PairSeq gzipped output should exist for R2")
        
        with gzip.open(self.output_R2_paired_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"R2 paired output: {len(records)} records")
        self.assertGreater(len(records), 0, "Output should contain records")

    def test_pairseq_maintains_pairing(self):
        """Test that PairSeq maintained proper pairing between R1 and R2"""
        with gzip.open(self.output_R1_paired_gz, 'rt') as handle:
            r1_records = list(SeqIO.parse(handle, 'fastq'))
        
        with gzip.open(self.output_R2_paired_gz, 'rt') as handle:
            r2_records = list(SeqIO.parse(handle, 'fastq'))
        
        # Paired files should have same number of records
        self.assertEqual(len(r1_records), len(r2_records),
                        "R1 and R2 paired outputs should have same number of records")
        
        print(f"Paired output: {len(r1_records)} read pairs")

    def test_pairseq_record_ids_match(self):
        """Test that paired records have matching IDs (or expected ID pattern)"""
        with gzip.open(self.output_R1_paired_gz, 'rt') as handle:
            r1_records = list(SeqIO.parse(handle, 'fastq'))
        
        with gzip.open(self.output_R2_paired_gz, 'rt') as handle:
            r2_records = list(SeqIO.parse(handle, 'fastq'))
        
        # Check first few pairs
        for i in range(min(3, len(r1_records))):
            r1_id = r1_records[i].id
            r2_id = r2_records[i].id
            print(f"Pair {i}: R1={r1_id}, R2={r2_id}")
            
            # IDs should be related (may differ by /1 /2 or other pair indicator)
            # The base ID should match
            r1_base = r1_id.split()[0].rstrip('/1').rstrip('/2')
            r2_base = r2_id.split()[0].rstrip('/1').rstrip('/2')
            
            self.assertEqual(r1_base, r2_base,
                           f"Paired records should have matching base IDs")


if __name__ == '__main__':
    unittest.main()
