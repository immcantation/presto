"""
Unit tests for AssemblePairs
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import os
import sys
import time
import unittest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.IO import readSeqFile

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')
usearch_exec = 'usearch'
blastn_exec = 'blastn'
blastdb_exec = 'makeblastdb'

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import AssemblePairs
import presto.Applications as Applications

# Test files
ref_file = os.path.join(data_path, 'human_igv_trv.fasta')
head_file = os.path.join(data_path, 'head.fasta')
tail_file = os.path.join(data_path, 'tail.fasta')

class TestAssemblePairs(unittest.TestCase):

    def setUp(self):
        print('-> %s()' % self._testMethodName)

        # Set pandas options
        pd.set_option('display.width', 200)

        # Load data
        self.head_list = [s.upper() for s in readSeqFile(head_file)]
        self.tail_list = [s.upper() for s in readSeqFile(tail_file)]
        self.ref_dict = {s.id: s.upper() for s in readSeqFile(ref_file)}

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    @unittest.skip("-> makeUSearchDb() skipped\n")
    def test_makeUSearchDb(self):
        database, db_handle = Applications.makeUBlastDb(ref_file, db_exec=usearch_exec)
        print('DATABASE>', database)
        print('  HANDLE>', db_handle.name)

        self.assertIsNotNone(db_handle)

    @unittest.skip("-> makeBlastnDb() skipped\n")
    def test_makeBlastnDb(self):
        database, db_handle = Applications.makeBlastnDb(ref_file, db_exec=blastdb_exec)
        print('DATABASE>', database)
        print('  HANDLE>', db_handle.name)
        print('  DIR>', os.listdir(db_handle.name))

        self.assertIsNotNone(db_handle)

    @unittest.skip("-> runUSearch() skipped\n")
    def test_runUSearch(self):
        # Build reference database
        database, db_handle = Applications.makeUBlastDb(ref_file, db_exec=usearch_exec)
        print('DB>', database)

        # Run alignment
        for head_seq, tail_seq in zip(self.head_list, self.tail_list):
            head_df = Applications.runUBlast(head_seq, database, evalue=1e-5,
                                             aligner_exec=usearch_exec)
            tail_df = Applications.runUBlast(tail_seq, database, evalue=1e-5,
                                             aligner_exec=usearch_exec)
            print('HEAD> %s' % head_seq.id)
            print(head_df)
            print('TAIL> %s' % tail_seq.id)
            print(tail_df)

        self.fail('TODO')

    @unittest.skip("-> runBlastn() skipped\n")
    def test_runBlastn(self):
        # Build reference database
        database, db_handle = Applications.makeBlastnDb(ref_file, db_exec=blastdb_exec)
        print('      DB>', database)

        # Run alignment
        for head_seq, tail_seq in zip(self.head_list, self.tail_list):
            head_df = Applications.runBlastn(head_seq, database, evalue=1e-5,
                                             aligner_exec=blastn_exec)
            tail_df = Applications.runBlastn(tail_seq, database, evalue=1e-5,
                                             aligner_exec=blastn_exec)
            print('HEAD> %s' % head_seq.id)
            print(head_df)
            print('TAIL> %s' % tail_seq.id)
            print(tail_df)

            self.fail('TODO')

    @unittest.skip("-> referenceAssembly() skipped\n")
    def test_referenceAssembly(self):
        stitch = AssemblePairs.referenceAssembly(self.head_rec, self.tail_rec, self.ref_dict, self.ref_file, fill=False)

        print('   REFID> %s' % stitch.ref_seq.id)
        print('  REFSEQ> %s' % (' ' * stitch.ref_pos[0] + stitch.ref_seq.seq))
        print('ASSEMBLY> %s' % stitch.seq.seq)
        print('     GAP> %s' % stitch.gap)
        print(' EVALUE1> %.4e' % stitch.evalue[0])
        print(' EVALUE2> %.4e' % stitch.evalue[1])
        print('IDENTITY> %.4f' % stitch.ident)

        stitch = AssemblePairs.referenceAssembly(self.head_rec, self.tail_rec, self.ref_dict, self.ref_file, fill=True)

        print('   REFID> %s' % stitch.ref_seq.id)
        print('  REFSEQ> %s' % (' ' * stitch.ref_pos[0] + stitch.ref_seq.seq))
        print('ASSEMBLY> %s' % stitch.seq.seq)
        print('     GAP> %s' % stitch.gap)
        print(' EVALUE1> %.4e' % stitch.evalue[0])
        print(' EVALUE2> %.4e' % stitch.evalue[1])
        print('IDENTITY> %.4f' % stitch.ident)

        #print tuple(stitch.evalue)
        self.fail()

    @unittest.skip("-> alignAssembly() skipped\n")
    def test_alignAssembly(self):
        head = SeqRecord(Seq("TTTCCGG"), id="HEAD",
                         letter_annotations={'phred_quality':[40,40,40,20,20,40,40]})
        tail = SeqRecord(Seq("CTGGAAA"), id="TAIL",
                         letter_annotations={'phred_quality':[40,20,40,40,40,40,40]})

        stitch = AssemblePairs.alignAssembly(head, tail, alpha=0.1)
        print('    HEAD> %s' % head.seq)
        print('    TAIL>    %s\n' % tail.seq)
        print('ASSEMBLY>', stitch.seq.seq)
        print('   ERROR>', stitch.error)
        print('  PVALUE>', stitch.pvalue, '\n')
        self.fail()


if __name__ == '__main__':
    unittest.main()