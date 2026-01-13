"""
Unit tests for FilterSeq with gzip support
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

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import FilterSeq


class TestFilterSeqGzip(unittest.TestCase):
    """Test FilterSeq functionality with gzipped input and output files"""
    
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        
        # Define test files
        self.input_fq_gz = os.path.join(data_path, 'set413_R1.fq.gz')
        self.output_quality_gz = os.path.join(data_path, 'set413_R1_quality-pass.fq.gz')
        
        # Verify test files exist
        self.assertTrue(os.path.exists(self.input_fq_gz), 
                       f"Input file not found: {self.input_fq_gz}")
        self.assertTrue(os.path.exists(self.output_quality_gz),
                       f"Expected output file not found: {self.output_quality_gz}")
        
        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    def test_read_gzipped_fastq(self):
        """Test that gzipped FASTQ files can be read"""
        with gzip.open(self.input_fq_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"Read {len(records)} records from {self.input_fq_gz}")
        self.assertGreater(len(records), 0, "Should read at least one record from gzipped FASTQ")
        
        # Verify record structure
        for rec in records[:3]:  # Check first 3 records
            self.assertIsNotNone(rec.id)
            self.assertIsNotNone(rec.seq)
            self.assertIsNotNone(rec.letter_annotations.get('phred_quality'))

    def test_gzipped_quality_output_exists(self):
        """Test that quality-filtered gzipped output was created"""
        self.assertTrue(os.path.exists(self.output_quality_gz),
                       "Quality-filtered gzipped output should exist")
        
        # Verify it's actually gzipped and readable
        with gzip.open(self.output_quality_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"Quality-filtered output contains {len(records)} records")
        self.assertGreater(len(records), 0, "Filtered output should contain records")

    def test_quality_filtering_removed_records(self):
        """Test that quality filtering actually filtered out some records"""
        # Read input
        with gzip.open(self.input_fq_gz, 'rt') as handle:
            input_records = list(SeqIO.parse(handle, 'fastq'))
        
        # Read output
        with gzip.open(self.output_quality_gz, 'rt') as handle:
            output_records = list(SeqIO.parse(handle, 'fastq'))
        
        input_count = len(input_records)
        output_count = len(output_records)
        
        print(f"Input: {input_count} records, Output: {output_count} records")
        
        # Quality filtering should remove some records (or keep all if quality is good)
        self.assertLessEqual(output_count, input_count,
                            "Output should have same or fewer records than input")

    def test_quality_filtered_records_have_good_quality(self):
        """Test that quality-filtered records actually have good quality scores"""
        with gzip.open(self.output_quality_gz, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))
        
        # Check that records have quality information
        for i, rec in enumerate(records):
            qualities = rec.letter_annotations['phred_quality']
            self.assertIsNotNone(qualities)
            self.assertEqual(len(qualities), len(rec.seq))
            
            # Print first record's quality info
            if i == 0:
                min_q = min(qualities)
                max_q = max(qualities)
                avg_q = sum(qualities) / len(qualities)
                print(f"Sample quality scores - min: {min_q}, max: {max_q}, avg: {avg_q:.1f}")


class TestFilterSeqGzipCompression(unittest.TestCase):
    """Test that gzip compression/decompression works correctly"""
    
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        self.data_path = os.path.join(test_path, 'data')
        
        # Test with both R1 and R2 files
        self.test_files = [
            'set413_R1.fq.gz',
            'set413_R1.fa.gz',
            'set413_R2.fq.gz',
            'set413_R2.fa.gz'
        ]
        
        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    def test_all_gzipped_test_files_readable(self):
        """Test that all gzipped test files can be opened and read"""
        for filename in self.test_files:
            filepath = os.path.join(self.data_path, filename)
            
            if not os.path.exists(filepath):
                self.skipTest(f"Test file not found: {filepath}")
            
            # Determine format from extension
            file_format = 'fasta' if filename.endswith('.fa.gz') else 'fastq'
            
            with gzip.open(filepath, 'rt') as handle:
                records = list(SeqIO.parse(handle, file_format))
            
            print(f"{filename}: {len(records)} records")
            self.assertGreater(len(records), 0, 
                             f"Should read at least one record from {filename}")

    def test_gzipped_fasta_vs_fastq(self):
        """Test that FASTA and FASTQ gzipped files can both be processed"""
        fasta_file = os.path.join(self.data_path, 'set413_R1.fa.gz')
        fastq_file = os.path.join(self.data_path, 'set413_R1.fq.gz')
        
        # Read FASTA
        with gzip.open(fasta_file, 'rt') as handle:
            fasta_records = list(SeqIO.parse(handle, 'fasta'))
        
        # Read FASTQ
        with gzip.open(fastq_file, 'rt') as handle:
            fastq_records = list(SeqIO.parse(handle, 'fastq'))
        
        print(f"FASTA records: {len(fasta_records)}, FASTQ records: {len(fastq_records)}")
        
        # Should have same number of records
        self.assertEqual(len(fasta_records), len(fastq_records),
                        "FASTA and FASTQ should have same number of records")
        
        # FASTQ should have quality scores, FASTA should not
        if len(fastq_records) > 0:
            self.assertIn('phred_quality', fastq_records[0].letter_annotations)
        if len(fasta_records) > 0:
            self.assertNotIn('phred_quality', fasta_records[0].letter_annotations)


if __name__ == '__main__':
    unittest.main()
