#!/usr/bin/env python3
"""
Removes primers and annotates sequences with primer and barcode identifiers
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import os
import sys
from argparse import ArgumentParser
from collections import OrderedDict

from textwrap import dedent

# Presto imports
from presto.Defaults import default_out_args, default_gap_penalty, default_max_error, \
                            default_max_len, default_start
from presto.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from presto.Sequence import alignPrimers, compilePrimers, getDNAScoreDict, maskSeq, \
                            reverseComplement, scorePrimers
from presto.IO import readPrimerFile, printLog
from presto.Multiprocessing import SeqResult, manageProcesses, feedSeqQueue, \
                                   collectSeqQueue


def processMPQueue(alive, data_queue, result_queue, align_func, align_args={}, 
                   mask_args={}, max_error=default_max_error):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
      alive : a multiprocessing.Value boolean controlling whether processing
              continues; when False function returns
      data_queue : a multiprocessing.Queue holding data to process
      result_queue : a multiprocessing.Queue to hold processed results
      align_func : the function to call for alignment
      align_args : a dictionary of arguments to pass to align_func
      mask_args : a dictionary of arguments to pass to maskSeq
      max_error : maximum acceptable error rate for a valid alignment
    
    Returns: 
      None
    """
    try:
        # Iterator over data queue until sentinel object reached
        while alive.value:
            # Get data from queue
            if data_queue.empty():  continue
            else:  data = data_queue.get()
            # Exit upon reaching sentinel
            if data is None:  break
            
            # Define result object for iteration
            in_seq = data.data
            result = SeqResult(in_seq.id, in_seq)
     
            # Align primers
            align = align_func(in_seq, **align_args)
            
            # Process alignment results
            if not align:
                # Update log if no alignment
                result.log['ALIGN'] = None
            else:
                # Create output sequence
                out_seq = maskSeq(align, **mask_args)
                result.results = out_seq
                result.valid = bool(align.error <= max_error) if len(out_seq) > 0 else False
                
                # Update log with successful alignment results
                result.log['SEQORIENT'] = out_seq.annotations['seqorient']
                result.log['PRIMER'] = out_seq.annotations['primer']
                result.log['PRORIENT'] = out_seq.annotations['prorient']
                result.log['PRSTART'] = out_seq.annotations['prstart']
                if 'barcode' in out_seq.annotations:  
                    result.log['BARCODE'] = out_seq.annotations['barcode']
                if not align.rev_primer:
                    align_cut = len(align.align_seq) - align.gaps
                    result.log['INSEQ'] = align.align_seq + \
                                          str(align.seq.seq[align_cut:])
                    result.log['ALIGN'] = align.align_primer
                    result.log['OUTSEQ'] = str(out_seq.seq).rjust(len(in_seq) + align.gaps)
                else:
                    align_cut = len(align.seq) - len(align.align_seq) + align.gaps
                    result.log['INSEQ'] = str(align.seq.seq[:align_cut]) + align.align_seq
                    result.log['ALIGN'] = align.align_primer.rjust(len(in_seq) + align.gaps)
                    result.log['OUTSEQ'] = str(out_seq.seq)
                result.log['ERROR'] = align.error
            
            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        sys.stderr.write('Error processing sequence with ID: %s.\n' % data.id)
        raise
    
    return None


def maskPrimers(seq_file, primer_file, mode, align_func, align_args={}, 
                max_error=default_max_error, barcode=False,
                out_args=default_out_args, nproc=None, queue_size=None):
    """
    Masks or cuts primers from sample sequences using local alignment

    Arguments: 
      seq_file : name of file containing sample sequences
      primer_file : name of the file containing primer sequences
      mode : defines the action taken; one of 'cut','mask','tag'
      align_func : the function to call for alignment
      align_arcs : a dictionary of arguments to pass to align_func
      max_error : maximum acceptable error rate for a valid alignment
      barcode : if True add sequence preceding primer to description
      out_args : common output argument dictionary from parseCommonArgs
      nproc : the number of processQueue processes;
              if None defaults to the number of CPUs
      queue_size : maximum size of the argument queue;
                   if None defaults to 2*nproc
                 
    Returns:
      list : a list of successful output file names
    """
    # Define subcommand label dictionary
    cmd_dict = {alignPrimers:'align', scorePrimers:'score'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MaskPrimers'
    log['COMMAND'] = cmd_dict.get(align_func, align_func.__name__)
    log['SEQ_FILE'] = os.path.basename(seq_file)
    log['PRIMER_FILE'] = os.path.basename(primer_file)
    log['MODE'] = mode
    log['BARCODE'] = barcode
    log['MAX_ERROR'] = max_error
    if 'start' in align_args: log['START_POS'] = align_args['start']
    if 'max_len' in align_args: log['MAX_LEN'] = align_args['max_len']
    if 'rev_primer' in align_args: log['REV_PRIMER'] = align_args['rev_primer']
    if 'skip_rc' in align_args: log['SKIP_RC'] = align_args['skip_rc']
    if 'gap_penalty' in align_args:
        log['GAP_PENALTY'] = ', '.join([str(x) for x in align_args['gap_penalty']])
    log['NPROC'] = nproc
    printLog(log)

    # Create dictionary of primer sequences to pass to maskPrimers
    primers = readPrimerFile(primer_file)
    if 'rev_primer' in align_args and align_args['rev_primer']:
        primers = {k: reverseComplement(v) for k, v in primers.items()}

    # Define alignment arguments and compile primers for align mode
    align_args['primers'] = primers 
    align_args['score_dict'] = getDNAScoreDict(mask_score=(0, 1), gap_score=(0, 0))
    if align_func is alignPrimers:
        align_args['max_error'] = max_error
        align_args['primers_regex'] = compilePrimers(primers)
    
    # Define sequence masking arguments
    mask_args = {'mode': mode, 
                 'barcode': barcode, 
                 'delimiter': out_args['delimiter']}

    # Define feeder function and arguments
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file}
    # Define worker function and arguments
    work_func = processMPQueue
    work_args = {'align_func': align_func, 
                 'align_args': align_args,
                 'mask_args': mask_args,
                 'max_error': max_error}
    
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'primers',
                    'out_args': out_args}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)

    # Print log
    result['log']['END'] = 'MaskPrimers'
    printLog(result['log'])
        
    return result['out_files']


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
      None
                      
    Returns: 
      argparse.ArgumentParser : argument parser.
    """
    # Define output file names and header fields
    fields = dedent(
             '''
             output files:
                 mask-pass
                     processed reads with successful primer matches.
                 mask-fail
                     raw reads failing primer identification.

             output annotation fields:
                 SEQORIENT
                     the orientation of the output sequence. Either F (input) or RC
                     (reverse complement of input).
                 PRIMER
                     name of the best primer match.
                 BARCODE
                     the sequence preceding the primer match. Only output when the
                     --barcode flag is specified.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter, add_help=False)
    group_help = parser.add_argument_group('help')
    group_help.add_argument('--version', action='version',
                            version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    group_help.add_argument('-h', '--help', action='help', help='show this help message and exit')
    subparsers = parser.add_subparsers(title='subcommands', metavar='',
                                       help='Alignment method')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True
    
    # Parent parser
    parent_parser = getCommonArgParser(multiproc=True)
    group_parent = parent_parser.add_argument_group('primer identification arguments')
    group_parent.add_argument('-p', action='store', dest='primer_file', required=True,
                              help='A FASTA or REGEX file containing primer sequences.')
    group_parent.add_argument('--mode', action='store', dest='mode',
                              choices=('cut', 'mask', 'trim', 'tag'), default='mask',
                              help='''Specifies the action to take with the primer sequence.
                                   The "cut" mode will remove both the primer region and
                                   the preceding sequence. The "mask" mode will replace the
                                   primer region with Ns and remove the preceding sequence.
                                   The "trim" mode will remove the region preceding the primer,
                                   but leave the primer region intact. The "tag" mode will
                                   leave the input sequence unmodified.''')
    group_parent.add_argument('--revpr', action='store_true', dest='rev_primer',
                              help='''Specify to match the tail-end of the sequence against the
                                   reverse complement of the primers. This also reverses the
                                   behavior of the --maxlen argument, such that the search
                                   window begins at the tail-end of the sequence.''')
    group_parent.add_argument('--barcode', action='store_true', dest='barcode',
                              help='''Specify to encode sequences with barcode sequences
                                   (unique molecular identifiers) found preceding the primer
                                   region.''')
    group_parent.add_argument('--maxerror', action='store', dest='max_error', type=float,
                              default=default_max_error, help='Maximum allowable error rate.')

    # Align mode argument parser
    parser_align = subparsers.add_parser('align', parents=[parent_parser],
                                         formatter_class=CommonHelpFormatter, add_help=False,
                                         help='Find primer matches using pairwise local alignment.',
                                         description='Find primer matches using pairwise local alignment.')
    group_align = parser_align.add_argument_group('alignment arguments')
    group_align.add_argument('--maxlen', action='store', dest='max_len', type=int,
                              default=default_max_len,
                              help='''Length of the sequence window to scan for primers.''')
    group_align.add_argument('--skiprc', action='store_true', dest='skip_rc',
                              help='Specify to prevent checking of sample reverse complement sequences.')
    group_align.add_argument('--gap', nargs=2, action='store', dest='gap_penalty',
                              type=float, default=default_gap_penalty,
                              help='''A list of two positive values defining the gap open
                                   and gap extension penalties for aligning the primers.
                                   Note: the error rate is calculated as the percentage
                                   of mismatches from the primer sequence with gap
                                   penalties reducing the match count accordingly; this may
                                   lead to error rates that differ from strict mismatch
                                   percentage when gaps are present in the alignment.''')
    parser_align.set_defaults(align_func=alignPrimers)
    

    # Score mode argument parser
    parser_score = subparsers.add_parser('score', parents=[parent_parser],
                                         formatter_class=CommonHelpFormatter, add_help=False,
                                         help='Find primer matches by scoring primers at a fixed position.',
                                         description='Find primer matches by scoring primers at a fixed position.')
    group_score = parser_score.add_argument_group('scoring arguments')
    group_score.add_argument('--start', action='store', dest='start', type=int, default=default_start,
                              help='The starting position of the primer')
    parser_score.set_defaults(align_func=scorePrimers)
    
    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    
    # Define align_args dictionary to pass to maskPrimers
    if args_dict['align_func'] is alignPrimers:
        args_dict['align_args'] = {'max_len':args_dict['max_len'],
                                   'rev_primer':args_dict['rev_primer'],
                                   'skip_rc':args_dict['skip_rc'],
                                   'gap_penalty':args_dict['gap_penalty']}
        del args_dict['max_len']
        del args_dict['rev_primer']
        del args_dict['skip_rc']
        del args_dict['gap_penalty']
    elif args_dict['align_func'] is scorePrimers:
        args_dict['align_args'] = {'start':args_dict['start'],
                                   'rev_primer':args_dict['rev_primer']}
        del args_dict['start']
        del args_dict['rev_primer']
    
    # Call maskPrimers for each sample file
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        maskPrimers(**args_dict)
    