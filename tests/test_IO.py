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

    #@unittest.skip('-> test_readSeqFile_gzipped_with_key_func() skipped\n')
    def test_readSeqFile_gzipped_with_key_func(self):
        """
        Test for issue #108: readSeqFile with .gz + index=True should respect key_func
        
        When using .gz FASTQ with index=True, readSeqFile should use key_func 
        to create keys instead of using rec.id directly. This is essential for 
        PairSeq --coord presto to work with gzipped consensus files.
        
        Issue: https://github.com/immcantation/presto/issues/108
        """
        # Test files with presto annotations (ID|CONSCOUNT=X|...)
        presto_r1 = os.path.join(data_path, 'test_presto_R1.fastq.gz')
        presto_r2 = os.path.join(data_path, 'test_presto_R2.fastq.gz')
        
        if not os.path.exists(presto_r1) or not os.path.exists(presto_r2):
            self.skipTest(f"Test files not found: {presto_r1} or {presto_r2}")
        
        # Import getCoordKey to test with presto coordinate extraction
        from presto.Annotation import getCoordKey
        
        # Define key function that extracts presto ID (strips annotations after |)
        def presto_key_func(header):
            return getCoordKey(header, coord_type='presto')
        
        # Read R1 file with index=True and key_func
        r1_dict = readSeqFile(presto_r1, index=True, key_func=presto_key_func)
        print(f'R1 sequences indexed: {len(r1_dict)}')
        print(f'R1 keys: {list(r1_dict.keys())}')
        
        # Read R2 file with index=True and key_func
        r2_dict = readSeqFile(presto_r2, index=True, key_func=presto_key_func)
        print(f'R2 sequences indexed: {len(r2_dict)}')
        print(f'R2 keys: {list(r2_dict.keys())}')
        
        # Verify that keys are clean IDs without annotations
        # (not full headers like "SEQUENCE001|CONSCOUNT=2|...")
        for key in r1_dict.keys():
            self.assertNotIn('|', key, 
                f"Key should not contain annotations: {key}")
            self.assertNotIn('CONSCOUNT', key,
                f"Key should not contain CONSCOUNT: {key}")
        
        # Verify that both R1 and R2 have the same keys (pairable sequences)
        self.assertEqual(set(r1_dict.keys()), set(r2_dict.keys()),
            "R1 and R2 should have matching keys when using presto key_func")
        
        # Expected keys based on our test data
        expected_keys = {'SEQUENCE001', 'SEQUENCE002', 'SEQUENCE003'}
        self.assertEqual(set(r1_dict.keys()), expected_keys,
            f"Keys should be clean sequence IDs: {expected_keys}")
        
        # Verify the original full headers are preserved in the records
        for key, record in r1_dict.items():
            self.assertIn('|', record.id,
                f"Original record.id should contain full header with annotations")
            self.assertTrue(record.id.startswith(key),
                f"Record ID should start with the key: {record.id} vs {key}")
        

if __name__ == '__main__':
    unittest.main()