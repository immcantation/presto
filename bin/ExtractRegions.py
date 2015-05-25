#!/usr/bin/env python
"""
Extracts and annotated sequences regions based on reference alignments
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import csv
import os
import sys
import tempfile
from textwrap import dedent
from argparse import ArgumentParser
from collections import OrderedDict
from subprocess import CalledProcessError, check_output, STDOUT
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_delimiter, default_out_args, default_usearch_exec
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation
from presto.IO import printLog, readSeqFile, getFileType
from presto.Multiprocessing import SeqResult, manageProcesses, feedSeqQueue, \
                                   processSeqQueue, collectSeqQueue

# Defaults
default_action='cut'
default_min_ident = 0.5
default_evalue = 1e-5
default_region_field = 'REGION'

default_min_len = 1
default_max_len = 1000


def extractSeqRegion(seq, start, end, action=default_action):
    """
    Extracts a sub-sequence from a SeqRecord

    Arguments:
    seq = a SeqRecord object to modify
    start = sub-sequence start position
    end = sub-sequence end position
    action = action to take on input SeqRecord; one of 'cut', 'mask', 'trim', 'tag'

    Returns:
    a tuple of (modified input SeqRecord, extracted sub-sequence SeqRecord)
    """
    pass


def runUBlastAlignment(seq_file, ref_file, evalue=default_evalue,
                       usearch_exec=default_usearch_exec):
    """
    Aligns a sequence file against a reference database using the UBLAST algorithm of USEARCH

    Arguments:
    seq_file = the input sequence file to align
    ref_file = the reference sequence file file
    evalue = the E-value cut-off for ublast
    usearch_exec = the path to the usearch executable

    Returns:
    a DataFrame of alignment results
    """
    # Generate temporary FASTA input file if needed
    if getFileType(seq_file) != 'fasta':
        in_handle = tempfile.NamedTemporaryFile()
        seq_iter = readSeqFile(seq_file)
        SeqIO.write(seq_iter, in_handle, 'fasta')
        seq_file = in_handle.name
        in_handle.seek(0)
    else:
        in_handle = None

    # Open temporary output file
    out_handle = tempfile.NamedTemporaryFile()

    # Define usearch command
    cmd = [usearch_exec,
           '-ublast', seq_file,
           '-db', ref_file,
           '-strand', 'plus',
           '-evalue', str(evalue),
           '-maxhits', '1',
           '-output_no_hits',
           '-userout', out_handle.name,
           '-userfields', 'query+target+qlo+qhi+tlo+thi+alnlen+evalue+id',
           '-threads', '1']

    # Run usearch ublast algorithm
    stdout_str = check_output(cmd, stderr=STDOUT, shell=False)

    # Parse usearch output
    field_names = ['query', 'target', 'query_start', 'query_end', 'target_start', 'target_end',
                   'length', 'evalue', 'identity']
    align_df = pd.read_table(out_handle, header=None, names=field_names)
    # Convert to base-zero indices
    align_df[['query_start', 'query_end', 'target_start', 'target_end']] -= 1

    # Close temporary file handles
    if in_handle is not None:  in_handle.close()
    out_handle.close()

    return align_df


def extractRegions(seq_file, ref_file, action=default_action,
                   min_ident=default_min_ident, evalue=default_evalue,
                   region_field=default_region_field,
                   usearch_exec=default_usearch_exec,
                   out_args=default_out_args, nproc=None,
                   queue_size=None):
    """
    Extracts regions from sequences based on reference alignment

    Arguments:
    seq_file = the sample sequence file name
    ref_file = the path to the reference database file
    action = defines the action taken; one of 'cut', 'mask', 'tag', 'trim'
    region_field = the name of the field to use for region annotation
    min_ident = the minimum identity for a valid alignment
    evalue = the E-value cut-off for ublast
    usearch_exec = the path to the usearch executable
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc

    Returns:
    a list of (valid_file, invalid_file) names
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'ExtractRegions'
    log['FILE'] = os.path.basename(seq_file)
    log['REF_FILE'] = os.path.basename(ref_file)
    log['REGION_FIELD'] = region_field
    log['MIN_IDENT'] = min_ident
    log['EVALUE'] = evalue
    log['NPROC'] = nproc
    printLog(log)

    # Define alignment function and parameters
    align_func = runUBlastAlignment
    align_args = {'action':action,
                  'min_ident':min_ident,
                  'evalue':evalue,
                  'usearch_exec':usearch_exec}

    # Define feeder function and arguments
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file}
    # Define worker function and arguments
    work_func = processSeqQueue
    work_args = {'region_field': region_field,
                 'align_func': align_func,
                 'align_args': align_args,
                 'delimiter': out_args['delimiter']}
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'extract',
                    'out_args': out_args}

    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func,
                             feed_args, work_args, collect_args,
                             nproc, queue_size)

    # Print log
    log = OrderedDict()
    log['OUTPUT'] = result['log'].pop('OUTPUT')
    for k, v in result['log'].iteritems():  log[k] = v
    log['END'] = 'ExtractRegions'
    printLog(log)

    return result['out_files']


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments:
    None

    Returns:
    an ArgumentParser object
    """
    # Define output file names and header fields
    fields = dedent(
             '''
             output files:
               extract-pass          reads with regions extracted and annotations added.
               extract-fail          raw reads failing alignment against the reference.

             output annotation fields:
               REGION                the identifier of the reference alignment.

             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__),
                            parents=[getCommonArgParser(multiproc=True)],
                            formatter_class=CommonHelpFormatter)

    parser.add_argument('-r', action='store', dest='ref_file', required=True,
                        help='''A FASTA file containing the reference sequence database.''')
    parser.add_argument('--act', action='store', dest='action',
                        choices=('cut', 'mask', 'trim', 'tag'), default=default_action,
                        help='''Specifies the action to take with the extracted sequence.
                             The "cut" action will remove both the extracted region and
                             the preceding sequence. The "mask" action will replace the
                             extracted region with Ns and remove the preceding sequence.
                             The "trim" action will remove the sequence preceding the
                             extracted region, but leave the extracted region intact.
                             The "tag" action will leave the input sequence unmodified.''')
    parser.add_argument('--minident', action='store', dest='min_ident', type=float,
                        default=default_min_ident,
                        help='''Minimum identity of the aligned region required to
                             define an alignment as valid (between 0 and 1).''')
    parser.add_argument('--evalue', action='store', dest='evalue', type=float,
                        default=default_evalue,
                        help='''Minimum E-value of the ublast reference alignment
                             required to define a valid alignment.''')
    parser.add_argument('--rf', action='store', dest='region_field', type=str,
                        default=default_region_field,
                        help='''The name of the annotation field to add to the input
                             sequences that will contain the name of the extracted
                             region reference sequence.''')
    parser.add_argument('--exec', action='store', dest='usearch_exec',
                        default=default_usearch_exec,
                        help='''The path to the usearch executable file.''')

    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)

    # Convert fields to uppercase
    if 'region_field' in args_dict and args_dict['region_field'] is not None:
        args_dict['region_field'] = args_dict['region_field'].upper()

    # Check if a valid usearch executable was specified
    if not os.path.isfile(args.usearch_exec):
        parser.error('%s does not exist' % args.usearch_exec)


    # Call cluster for each input file
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        extractRegions(**args_dict)

