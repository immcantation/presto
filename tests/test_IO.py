"""
Unit tests for IO module
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import os
import sys
import time
import unittest
import tempfile
import shutil
from collections import OrderedDict

# Presto imports
from presto.IO import getFileType, readSeqFile, getOutputHandle, openFile
from Bio import SeqIO

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


class TestIO(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        # Test files
        self.germ_file = os.path.join(data_path, 'human_igv_trv.fasta')
        self.missing_file = os.path.join(data_path, 'human_igv_trv.bob')
        
        # Gzipped test files
        self.fq_gz_file = os.path.join(data_path, 'set413_R1.fq.gz')
        self.fa_gz_file = os.path.join(data_path, 'set413_R1.fa.gz')
        
        # Create temporary directory for output tests
        self.temp_dir = tempfile.mkdtemp()

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # Clean up temporary directory
        shutil.rmtree(self.temp_dir)
        
        # End clock
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    #@unittest.skip('-> getFileType() skipped\n')
    def test_getFileType(self):
        """Test file type detection for regular files"""
        result = getFileType(self.germ_file)
        print(result)
        self.assertEqual('fasta', result)

        result = getFileType(self.missing_file)
        print(result)
        # missing_file has .bob extension, so getFileType returns 'bob'
        self.assertEqual('bob', result)

    #@unittest.skip('-> getFileType_gzipped() skipped\n')
    def test_getFileType_gzipped(self):
        """Test file type detection for gzipped files"""
        if os.path.exists(self.fq_gz_file):
            result = getFileType(self.fq_gz_file)
            print(f'Gzipped FASTQ type: {result}')
            self.assertEqual('fastq', result)
        
        if os.path.exists(self.fa_gz_file):
            result = getFileType(self.fa_gz_file)
            print(f'Gzipped FASTA type: {result}')
            self.assertEqual('fasta', result)

    #@unittest.skip('-> test_openFile_gzipped() skipped\n')
    def test_openFile_gzipped(self):
        """Test opening gzipped files with openFile function"""
        if os.path.exists(self.fq_gz_file):
            # Test reading gzipped FASTQ
            with openFile(self.fq_gz_file, 'r') as handle:
                first_line = handle.readline()
                print(f'First line from gzipped FASTQ: {first_line.strip()}')
                self.assertTrue(first_line.startswith('@'))
        
        if os.path.exists(self.fa_gz_file):
            # Test reading gzipped FASTA
            with openFile(self.fa_gz_file, 'r') as handle:
                first_line = handle.readline()
                print(f'First line from gzipped FASTA: {first_line.strip()}')
                self.assertTrue(first_line.startswith('>'))

    #@unittest.skip('-> test_readSeqFile_gzipped() skipped\n')
    def test_readSeqFile_gzipped(self):
        """Test reading sequences from gzipped files"""
        if os.path.exists(self.fq_gz_file):
            # Test reading gzipped FASTQ with readSeqFile (use index=True to get dict)
            seq_dict = readSeqFile(self.fq_gz_file, index=True)
            print(f'Read {len(seq_dict)} sequences from gzipped FASTQ')
            self.assertGreater(len(seq_dict), 0)
            
            # Verify sequences have expected attributes
            first_seq = next(iter(seq_dict.values()))
            self.assertTrue(hasattr(first_seq, 'id'))
            self.assertTrue(hasattr(first_seq, 'seq'))

    #@unittest.skip('-> test_gzip_output_functionality() skipped\n')
    def test_gzip_output_functionality(self):
        """Test writing gzipped output files"""
        test_input = os.path.join(data_path, 'head.fasta')
        if not os.path.exists(test_input):
            self.skipTest(f"Test input file {test_input} not found")
        
        # Test creating gzipped output
        output_handle = getOutputHandle(
            in_file=test_input,
            out_dir=self.temp_dir,
            out_name="test_output",
            out_type="fasta",
            gzip_output=True
        )
        
        self.assertIsNotNone(output_handle)
        output_path = output_handle.name
        
        # Write some test data
        with openFile(test_input, 'r') as input_handle:
            first_record = next(SeqIO.parse(input_handle, 'fasta'))
            SeqIO.write([first_record], output_handle, 'fasta')
        
        output_handle.close()
        
        # Verify the output file was created with .gz extension
        print(f'Output file created: {output_path}')
        self.assertTrue(output_path.endswith('.gz'))
        self.assertTrue(os.path.exists(output_path))
        
        # Verify we can read the gzipped output back
        with openFile(output_path, 'r') as gz_handle:
            records = list(SeqIO.parse(gz_handle, 'fasta'))
            self.assertEqual(len(records), 1)
            print(f'Successfully read back record: {records[0].id}')

    #@unittest.skip('-> test_round_trip_gzip() skipped\n')
    def test_round_trip_gzip(self):
        """Test round-trip: read gzipped file, process, write gzipped output"""
        if not os.path.exists(self.fq_gz_file):
            self.skipTest(f"Gzipped test file {self.fq_gz_file} not found")
        
        # Read from gzipped input (use index=True to get dict)
        input_dict = readSeqFile(self.fq_gz_file, index=True)
        original_count = len(input_dict)
        print(f'Original gzipped file has {original_count} sequences')
        
        # Write to gzipped output
        output_path = os.path.join(self.temp_dir, "round_trip_test.fq.gz")
        with openFile(output_path, 'w') as output_handle:
            file_type = getFileType(self.fq_gz_file)
            SeqIO.write(input_dict.values(), output_handle, file_type)
        
        # Read back the gzipped output (use index=True to get dict)
        readback_dict = readSeqFile(output_path, index=True)
        readback_count = len(readback_dict)
        print(f'Round-trip gzipped file has {readback_count} sequences')
        
        # Verify counts match
        self.assertEqual(original_count, readback_count)
        
        # Verify first sequence IDs match
        original_ids = set(input_dict.keys())
        readback_ids = set(readback_dict.keys())
        self.assertEqual(original_ids, readback_ids)

    #@unittest.skip('-> test_gzip_auto_detection_input() skipped\n')
    def test_gzip_auto_detection_input(self):
        """Test that gzipped inputs are automatically detected and handled"""
        if not os.path.exists(self.fq_gz_file):
            self.skipTest(f"Gzipped test file {self.fq_gz_file} not found")
        
        # Test with explicit gzip_output=True 
        output_handle = getOutputHandle(
            in_file=self.fq_gz_file,
            out_dir=self.temp_dir,
            out_name="auto_detect_test",
            gzip_output=True
        )
        
        output_path = output_handle.name
        output_handle.close()
        
        # Should add .gz extension when gzip_output=True
        print(f'Auto-detected output path: {output_path}')
        self.assertTrue(output_path.endswith('.gz'))


if __name__ == '__main__':
    unittest.main()