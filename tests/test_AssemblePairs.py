"""
Unit tests for AssemblePairs
"""
# Imports
import time
import unittest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from presto.IO import readSeqFile
from bin import AssemblePairs as script

# Info
__author__ = 'Jason Anthony Vander Heiden'


class TestAssemblePairs(unittest.TestCase):

    def setUp(self):
        print '-> %s()' % self._testMethodName

        # Set pandas options
        pd.set_option('display.width', 120)

        #self.ref_file = 'IMGT_Human_IGV.fasta'
        self.ref_file = 'data/IG_TR.V.human.F+ORF+infrP.ungapped.fasta'
        self.ref_dict = {s.id:s.upper() for s in readSeqFile(self.ref_file)}
        # MISEQ:121:000000000-A7VDM:1:1101:10041:1280
        #self.head_seq = Seq("CCACGTTTTAGTAATTAATACGGGAGCAAAAACCAGGGAAAGCCCCTAAGCTCCTGCTCTATGCTGCATCCACTTTGCAAAGTGTGGTCCCATCACGGTTCAGCGGCAGTGGATCTGGGACAGAATTCACTCTCACAATCAGCAGCCTGCAGCCTGAAGATTTTGCAACTTATTACTGTCAACAGCTTACTCCTTACCCTCCTACGTTCGCCCCAGGCCCCACGGTCGACCTCCACCCCCCTCTCGCTGCCCCCTCTCTCCCCTCCGACCCGCCTC")
        #self.tail_seq = Seq("GGCAAGAAGAAGAGGGGATACGAGAGCCAAGGGGGAGTGGAGTGGAGAGGTGTGCGCTTACGATCTACAAGTGTGAGTAATGAATACGGGAGCAAAAACGAGGGAAAGCCCCGAAGCTCCTGATCTATGCTGAATCCACTTTGCAAAGTGGGGTCCCATCAAGGTTCAGCGGCAGTGGATCTGGGACAGAATTCACTCTCACAATCAGCAGCCTGCAGCCTGAAGATTTTGCAACTTATTACTGTCAACAGCTTAATACTTACCCTCGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAACGAACTGTGGCTGCACCATCTGTC")

        # MISEQ:121:000000000-A7VDM:1:1101:13900:1538
        #self.head_seq = Seq("GGAGCTTGCATTAGCATCGATACGGGGAGCTGTGGGCTCAGAAGCAGAGTTCTGGGGTGTCTCCACCATGGCCTGGACCCCTCTCTGGCTCACTCTCCTCACTCTTTGCATAGGTTCTGTGGTTTCTTCTGAGCTGACTCAGGACCCTGCTGTGTCTGTGGCCTTGGGACAGACAGTCAGGATCACATGCCAAGGAGACCGCCTCAGAAGCTATTATGCAAGCTGGTCCCCGCAGAATCCAGGACAGGCCCCCCTCCTTCTCATCTATGGTATAAC")
        #self.tail_seq = Seq("CTTGGGACAGACAGTCAGGATCACAGGCCAAGGAGACAGGCTCAGAAGCTAGTATGCAAGCGGGTACCAGCAGAAGCCAGGACAGGCCCCTGTACTTGTCATCTATGGTAAAAACAACCGGCCCTCAGGGATCCCAGACCGATTCTCTGGCTTCAGCTCAGGAAACACAGCTTCCTTGACCATCACTGGGGCTCAGGCGGAAGATGAGGCTGACTATTACTGTAACTCCCGGGACAGCAGTGGTAACCATCCATTCGGCGGAGGGACCAAGCTGACCGTCCTAGGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCT")

        # Non-Overlaping MISEQ:121:000000000-A7VDM:1:1101:13900:1538
        #self.head_seq = Seq("GGAGCTTGCATTAGCATCGATACGGGGAGCTGTGGGCTCAGAAGCAGAGTTCTGGGGTGTCTCCACCATGGCCTGGACCCCTCTCTGGCTCACTCTCCTCACTCTTTGCATAGGTTCTGTGGTTTCTTCTGAGCTGACTCAGGACCCTGCTGTGTCTGTGGCCTTGGGACAGACAGTCAGGAC")
        #self.tail_seq = Seq("TTGTCATCTATGGTAAAAACAACCGGCCCTCAGGGATCCCAGACCGATTCTCTGGCTTCAGCTCAGGAAACACAGCTTCCTTGACCATCACTGGGGCTCAGGCGGAAGATGAGGCTGACTATTACTGTAACTCCCGGGACAGCAGTGGTAACCATCCATTCGGCGGAGGGACCAAGCTGACCGTCCTAGGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCT")

        # MISEQ:121:000000000-A7VDM:1:1101:22042:1446
        #self.head_seq = Seq("GTCTCAACAAAATACCATCTACGGGATGTGTATTGGTACCAGCAACTCCCAGGATCGGCCCCCAAGCTCCTCATTTATAGTAGTAGCCAGCCGCCCTCAGGGGTCCCTGACCGATTCTCTGGCTCCAAGTCTGGCACCTCAGCCTCCCTGGCCCTCCGTGTGCTCCGTTCCGAAGATGACGCCGCTTATTCCTGTTCAGCCTGGGCTTACCCCCTGATTCGTCTTTCCTTGTTCGCCCCACCGTCCCCTCCCCCCCTCCTACCTCCGCCCCACGCC")
        #self.tail_seq = Seq("GAGAGGTGGGAGGAGCCGATCTGGGTCAACAAAATAGGATAGACGGGAAGGGTATTGGTACAAGCAACTCCCAGGAGCGGAGCCCAAGCTCGTCATTTATAGGAGAAGCCAGGGGCCCTCAGGGGTACCTGACCGATTCTCTGGCTCCAAGTCTGGAAGCTAAGCCTCCCTGGCCATCAGTGGGCTCCGGTACGAAGATGAGGCTGATTATTACTGTGCAGCATGGGATGACAGCCTGAGTGGTCTTTGGGTGTTCGGCGGAGGGACCAGGCTGACCGTCCTAGGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCT")

        # TATATTATGTGCAGT_GCCTTC|PRCONS=TRB
        self.head_seq = Seq('gacacctctccccagagaaggtggtgtgagaccaccacggaagatgctgctgctactgctcctcctgtggataggctccgggcttggtgccgtcgtctctcaacatccgagccgggctatctgaaagtgtggaaccattgtcaaccacgagggccgttcaccggactttcggccccctacgaagttttggcaacctcagctcccgatacatggtcttttgctgctggcgacctccaaccggggctccac')
        self.tail_seq = Seq('gtcgagtgtcgttccctggactttcaagccacatcgaagtactggtatcggcagttcccgaaacagagtctaaggctggtggcaacatccaacgagtactcaaaggccacatacgagcaaggcgtcgagaaggacaagtatctcatcaaccatgcaagcctgaccttgtcctctcttacagtgaccagtgcccagcctgaagacagcagcttctacatctacagtgctgcaacgggtcaggggacgatcgagacccagtatttcgggcctggcacgcggctcctggtgctcgaggacctgaaaaac')
        #self.head_seq = self.head_seq.reverse_complement()

        self.head_rec = SeqRecord(self.head_seq, id='HEAD', name='HEAD', description='')
        self.tail_rec = SeqRecord(self.tail_seq, id='TAIL', name='TAIL', description='')
        #self.head_rec = SeqRecord(self.tail_seq, id='HEAD', name='HEAD', description='')
        #self.tail_rec = SeqRecord(self.head_seq, id='TAIL', name='TAIL', description='')

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    #@unittest.skip("-> getUBlastAlignment() skipped\n")
    def test_getUBlastAlignment(self):
        head_df = script.runUBlastAlignment(self.head_rec, self.ref_file)
        tail_df = script.runUBlastAlignment(self.tail_rec, self.ref_file)
        print 'HEAD SEQUENCE>'
        print head_df
        print 'TAIL SEQUENCE>'
        print tail_df
        self.fail()

    @unittest.skip("-> getBlastnAlignment() skipped\n")
    def test_getBlastnAlignment(self):
        head_df = script.runBlastnAlignment(self.head_rec, self.ref_file)
        tail_df = script.runBlastnAlignment(self.tail_rec, self.ref_file)
        print 'HEAD SEQUENCE>'
        print head_df
        print 'TAIL SEQUENCE>'
        print tail_df
        self.fail()

    #@unittest.skip("-> referenceAssembly() skipped\n")
    def test_referenceAssembly(self):
        stitch = script.referenceAssembly(self.head_rec, self.tail_rec, self.ref_dict, self.ref_file)

        print '   REFID> %s' % stitch.ref_seq.id
        print '  REFSEQ> %s' % (' ' * stitch.ref_pos[0] + stitch.ref_seq.seq)
        print 'ASSEMBLY> %s' % stitch.seq.seq
        print '     GAP> %s' % stitch.gap
        print ' EVALUE1> %.4e' % stitch.evalue[0]
        print ' EVALUE2> %.4e' % stitch.evalue[1]
        print 'IDENTITY> %.4f' % stitch.ident

        #print tuple(stitch.evalue)
        self.fail()

    @unittest.skip("-> alignAssembly() skipped\n")
    def test_alignAssembly(self):
        head = SeqRecord(Seq("TTTCCGG"), id="HEAD",
                         letter_annotations={'phred_quality':[40,40,40,20,20,40,40]})
        tail = SeqRecord(Seq("CTGGAAA"), id="TAIL",
                         letter_annotations={'phred_quality':[40,20,40,40,40,40,40]})

        stitch = script.alignAssembly(head, tail, alpha=0.1)
        print '    HEAD> %s' % head.seq
        print '    TAIL>    %s\n' % tail.seq
        print 'ASSEMBLY>', stitch.seq.seq
        print '   ERROR>', stitch.error
        print '  PVALUE>', stitch.pvalue, '\n'
        self.fail()