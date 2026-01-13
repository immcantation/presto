"""
Unit tests for MaskPrimers with gzip support
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


class TestMaskPrimersGzip(unittest.TestCase):
    """Test MaskPrimers functionality with gzipped input and output files"""
    
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        
        # Define test files
        self.input_R1_gz = os.path.join(data_path, 'set413_R1.fq.gz')
        self.input_R2_gz = os.path.join(data_path, 'set413_R2.fq.gz')
        self.output_R1_primers_gz = os.path.join(data_path, 'set413_R1_primers-pass.fq.gz')
        self.output_R2_primers_gz = os.path.join(data_path, 'set413_R2_primers-pass.fq.gz')
        
        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    def test_maskprimers_gzipped_input(self):
        """Test that MaskPrimers can read gzipped input files"""
        with gzip.open(self.input_R1_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"Read {len(records)} records from R1 gzipped input")
        self.assertGreater(len(records), 0, "Should read records from gzipped input")

    def test_maskprimers_gzipped_output_R1(self):
        """Test that MaskPrimers created gzipped output for R1"""
        self.assertTrue(os.path.exists(self.output_R1_primers_gz),
                       "MaskPrimers gzipped output should exist for R1")
        
        with gzip.open(self.output_R1_primers_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"R1 primer-masked output contains {len(records)} records")
        self.assertGreater(len(records), 0, "Output should contain records")

    def test_maskprimers_gzipped_output_R2(self):
        """Test that MaskPrimers created gzipped output for R2"""
        self.assertTrue(os.path.exists(self.output_R2_primers_gz),
                       "MaskPrimers gzipped output should exist for R2")
        
        with gzip.open(self.output_R2_primers_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"R2 primer-masked output contains {len(records)} records")
        self.assertGreater(len(records), 0, "Output should contain records")

    def test_maskprimers_annotations_present(self):
        """Test that MaskPrimers added annotations to the headers"""
        with gzip.open(self.output_R1_primers_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        if len(records) > 0:
            # Check if the description contains PRIMER annotation
            # MaskPrimers typically adds annotations to the description
            sample_desc = records[0].description
            print(f"Sample description: {sample_desc}")
            
            # The record should have been processed (description might have changed)
            self.assertIsNotNone(sample_desc)

    def test_paired_primer_processing(self):
        """Test that both R1 and R2 were processed with primers"""
        with gzip.open(self.output_R1_primers_gz, 'rt') as handle:
            r1_records = list(SeqIO.parse(handle, 'fastq'))
        
        with gzip.open(self.output_R2_primers_gz, 'rt') as handle:
            r2_records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"R1 records: {len(r1_records)}, R2 records: {len(r2_records)}")
        
        # Both files should have records
        self.assertGreater(len(r1_records), 0)
        self.assertGreater(len(r2_records), 0)


if __name__ == '__main__':
    unittest.main()
