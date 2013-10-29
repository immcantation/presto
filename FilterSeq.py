#!/usr/bin/env python
"""
Filters sequences or masks characters in FASTA/FASTQ files by missing nucleotides and quality scores
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.1'
__date__      = '2013.10.23'

# Imports
import os, sys
import multiprocessing as mp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from itertools import groupby
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_min_qual, default_out_args, default_missing_chars
from IgCore import getCommonParser, getFileType, parseCommonArgs, printLog
from IgCore import collectRecQueue, feedRecQueue

# Defaults
default_max_missing = 10
default_max_repeat = 15
default_min_length = 250
default_window = 10


def filterLength(seq, min_length=default_min_length, inner=True, 
                 missing_chars=''.join(default_missing_chars)):
    """
    Filters sequences by length
    
    Arguments: 
    seq = a SeqRecord object to process
    min_length = The minimum length allowed
    inner = if True exclude outer missing characters from calculation
    missing_chars = a string of missing character values
    
    Returns:
    a result dictionary containing {out_seq, pass, log}
    """
    # Remove outer missing characters if required
    if inner:  
        seq_str = str(seq.seq).strip(missing_chars)
        n = len(seq_str)
    else:
        n = len(seq)
    
    # Build result dictionary
    passed = (n >= min_length)
    result ={'out_seq':(seq if passed else None),
             'pass':passed,
             'log':OrderedDict([('SEQ', seq.seq),
                                ('LENGTH', n)])}
    
    return result


def filterMissing(seq, max_missing=default_max_missing, inner=True, 
                  missing_chars=''.join(default_missing_chars)):
    """
    Filters sequences by number of missing nucleotides
    
    Arguments: 
    seq = a SeqRecord object to process
    max_missing = The maximum number of allowed ambiguous characters
    inner = if True exclude outer missing characters from calculation
    missing_chars = a string of missing character values
    
    Returns:
    a result dictionary containing {out_seq, pass, log}
    """
    seq_str = str(seq.seq)
    # Remove outer missing character if required
    if inner:  seq_str = seq_str.strip(missing_chars)
    # Count missing characters
    n = len([c for c in seq_str if c in missing_chars])
    
    # Build result dictionary
    passed = (n <= max_missing)
    result ={'out_seq':(seq if passed else None),
             'pass':passed,
             'log':OrderedDict([('SEQ', seq.seq),
                                ('MISSING', n)])}
    
    return result


def filterRepeats(seq, max_repeat=default_max_repeat, include_missing=False, inner=True,
                  missing_chars=''.join(default_missing_chars)):
    """
    Filters sequences by fraction of ambiguous nucleotides
    
    Arguments: 
    seq = a SeqRecord object to process
    max_repeat = The maximum number of allowed repeating characters
    include_missing = if True count ambiguous character repeats;
                      if False do not consider ambiguous character repeats
    inner = if True exclude outer missing characters from calculation
    missing_chars = a string of missing character values
    
    Returns:
    a result dictionary containing {out_seq, pass, log}
    """
    seq_str = str(seq.seq)
    # Remove outer missing character if required
    if inner:  seq_str = seq_str.strip(missing_chars)
    # Remove missing characters if required
    if not include_missing:
        seq_str = ''.join([c for c in seq_str if c not in missing_chars])
    
    groups = ((c, len(list(g))) for c, g in groupby(seq_str))
    __, n = max(groups, key=lambda x: x[1])
    
    # Build result dictionary
    passed = (n <= max_repeat)
    result ={'out_seq':(seq if passed else None),
             'pass':passed,
             'log':OrderedDict([('SEQ', seq.seq),
                                ('REPEAT', n)])}
    
    return result


def filterQuality(seq, min_qual=default_min_qual, inner=True,
                  missing_chars=''.join(default_missing_chars)):
    """
    Filters sequences by quality score
    
    Arguments: 
    seq = a SeqRecord object to process
    min_qual = minimum mean quality score for retained sequences
    inner = if True exclude outer missing characters from calculation
    missing_chars = a string of missing character values
    
    Returns:
    a result dictionary containing {out_seq, pass, log}
    """
    if inner:  
        seq_str = str(seq.seq)
        seq_cut = seq_str.strip(missing_chars)
        s = seq_str.find(seq_cut)
        quals =  seq.letter_annotations['phred_quality'][s:s + len(seq_cut)]
    else:
        quals = seq.letter_annotations['phred_quality']
    
    q = sum(quals) / len(quals)

    # Build result dictionary
    passed = (q >= min_qual)
    result ={'out_seq':(seq if passed else None),
             'pass':passed,
             'log':OrderedDict([('SEQ', seq.seq),
                                ('QUALITY', q)])}
    
    return result


def trimQuality(seq, min_qual=default_min_qual, window=default_window, reverse=False):
    """
    Cuts sequences using a moving mean quality score
    
    Arguments: 
    seq = a SeqRecord object to process
    min_qual = minimum mean quality to define a cut point
    window = nucleotide window size
    reverse = if True cut the head of the sequence;
              if False cut the tail of the sequence
        
    Returns:
    a result dictionary containing {out_seq, pass, log}
    """
    quals = seq.letter_annotations['phred_quality']
    # Reverse quality scores if required
    if reverse:  quals = quals[::-1]
    
    # Scan across quality scores for first quality drop-off
    end = len(quals)
    for s in xrange(0, end, window):
        q_win = quals[s:s + window]
        q = sum(q_win) / len(q_win)
        if q < min_qual:
            end = s
            break

    # Define trimmed sequence
    if not reverse:
        trim_seq = seq[:end]
        out_str = str(trim_seq.seq)
    else:
        trim_seq = seq[len(seq) - end:]
        out_str =  ' ' * (len(seq) - end) + str(trim_seq.seq)
        
    # Build result dictionary
    passed = (len(trim_seq) > 0)
    result ={'out_seq':(trim_seq if passed else None),
             'pass':passed,
             'log':OrderedDict([('INSEQ', seq.seq),
                                ('OUTSEQ', out_str),
                                ('LENGTH', len(trim_seq))])}
    
    return result


def maskQuality(seq, min_qual=default_min_qual):
    """
    Masks characters by in sequence by quality score
    
    Arguments: 
    seq = a SeqRecord object to process
    min_qual = minimum quality for retained characters
        
    Returns:
    a result dictionary containing {out_seq, pass, log}
    """
    seq_str = str(seq.seq)
    quals = seq.letter_annotations['phred_quality']
    # Mask low quality nucleotides
    mask_chars = [seq_str[i] if q >= min_qual else 'N' for i, q in enumerate(quals)]
    # Define masked SeqRecord
    mask_seq = SeqRecord(Seq(''.join(mask_chars), IUPAC.ambiguous_dna), 
                         id=seq.id, 
                         name=seq.name, 
                         description=seq.description,
                         letter_annotations=seq.letter_annotations)

    # Build result dictionary
    result ={'out_seq':mask_seq,
             'pass':True,
             'log':OrderedDict([('INSEQ', seq.seq),
                                ('OUTSEQ', mask_seq.seq)])}
    
    return result


def processQueue(data_queue, result_queue, filter_func, filter_args={}):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    filter_func = the function to use for filtering sequences
    filter_args = a dictionary of arguments to pass to filter_func

    Returns: 
    None
    """
    # Iterator over data queue until sentinel object reached
    for args in iter(data_queue.get, None):
        in_seq = args['seq']        
        # Define result dictionary for iteration
        results = {'id':args['id'],
                   'in_seq':in_seq,
                   'out_seq':None,
                   'pass':False,
                   'log':OrderedDict([('ID', args['id'])])}
 
        # Perform filtering
        filter_result = filter_func(in_seq, **filter_args)
        results['out_seq'] = filter_result['out_seq']
        results['pass'] = filter_result['pass']
        for k, v in filter_result['log'].iteritems():  
            results['log'].update({k:v})
        
        # Feed results to result queue
        result_queue.put(results)

    return None


def filterSeq(seq_file, filter_func, filter_args={}, out_args=default_out_args, 
              nproc=None, queue_size=None):
    """
    Filters sequences by fraction of ambiguous nucleotides
    
    Arguments: 
    seq_file = the sequence file to filter
    filter_func = the function to use for filtering sequences
    filter_args = a dictionary of arguments to pass to filter_func
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                 
    Returns:
    a list of successful output file names
    """
    # Define number of processes and queue size
    if nproc is None:  nproc = mp.cpu_count()
    if queue_size is None:  queue_size = nproc * 2

    # Define output file label dictionary
    cmd_dict = {filterLength:'length', filterMissing:'missing', filterRepeats:'repeats', 
                  filterQuality:'quality', maskQuality:'maskqual', trimQuality:'trimqual'}

    # Print parameter info
    log = OrderedDict()
    log['START'] = 'FilterSeq'
    log['COMMAND'] = cmd_dict.get(filter_func, filter_func.__name__)
    log['FILE'] = os.path.basename(seq_file)
    for k in sorted(filter_args):  log[k.upper()] = filter_args[k]
    log['NPROC'] = nproc
    printLog(log)
    
    # Check input type
    in_type = getFileType(seq_file)
    if in_type != 'fastq' and filter_func in (filterQuality, maskQuality, trimQuality):
        sys.exit('ERROR:  Input file must be FASTQ for %s mode' % cmd_dict[filter_func])
    
    # Define shared data objects 
    manager = mp.Manager()
    collect_dict = manager.dict()
    data_queue = mp.Queue(queue_size)
    result_queue = mp.Queue(queue_size)
    
    # Initiate feeder process
    feeder = mp.Process(target=feedRecQueue, args=(data_queue, nproc, seq_file))
    feeder.start()

    # Initiate worker processes
    workers = []
    for __ in range(nproc):
        w = mp.Process(target=processQueue, args=(data_queue, result_queue, 
                                                  filter_func, filter_args))
        w.start()
        workers.append(w)

    # Initiate collector process
    collector = mp.Process(target=collectRecQueue, args=(result_queue, collect_dict, seq_file, 
                                                         cmd_dict.get(filter_func, 'filter'), 
                                                         out_args))
    collector.start()

    # Wait for feeder and worker processes to finish, add sentinel to result_queue
    feeder.join()
    for w in workers:  w.join()
    result_queue.put(None)
    
    # Wait for collector process to finish and shutdown manager
    collector.join()
    log = collect_dict['log']
    out_files = collect_dict['out_files']
    manager.shutdown()
    
    # Print log
    log['END'] = 'FilterSeq'
    printLog(log)
        
    return out_files


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
    None
                      
    Returns: 
    an ArgumentParser object
    """
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, version='%(prog)s:' + ' v%s-%s' %(__version__, __date__), 
                            formatter_class=ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()
    
    # Parent parser
    parser_parent = getCommonParser(annotation=False, log=True, multiproc=True)
    
    # Length filter mode argument parser
    parser_length = subparsers.add_parser('length', parents=[parser_parent],
                                          formatter_class=ArgumentDefaultsHelpFormatter, 
                                          help='Sequence length filtering mode')
    parser_length.add_argument('-n', action='store', dest='min_length', type=int, 
                               default=default_min_length, 
                               help='Minimum sequence length to retain')
    parser_length.add_argument('--inner', action='store_true', dest='inner',
                               help='If specified exclude consecutive missing characters at either end of the sequence')
    parser_length.set_defaults(filter_func=filterLength)
    
    # Missing character filter mode argument parser
    parser_missing = subparsers.add_parser('missing', parents=[parser_parent],
                                           formatter_class=ArgumentDefaultsHelpFormatter, 
                                           help='Missing nucleotide filtering mode')
    parser_missing.add_argument('-n', action='store', dest='max_missing', type=int, 
                                default=default_max_missing, 
                                help='Threshold for fraction of gap or N nucleotides')
    parser_missing.add_argument('--inner', action='store_true', dest='inner',
                                help='If specified exclude consecutive missing characters at either end of the sequence')
    parser_missing.set_defaults(filter_func=filterMissing)
    
    # Continuous repeating character filter mode argument parser
    parser_repeats = subparsers.add_parser('repeats', parents=[parser_parent],
                                           formatter_class=ArgumentDefaultsHelpFormatter, 
                                           help='Consencutive nucleotide repeating filtering mode')
    parser_repeats.add_argument('-n', action='store', dest='max_repeat', type=int, 
                                default=default_max_repeat, 
                                help='Threshold for fraction of repeating nucleotides')
    parser_repeats.add_argument('--missing', action='store_true', dest='include_missing',
                                help='If specified count consecutive gap and N characters in addition to {A,C,G,T}')
    parser_repeats.add_argument('--inner', action='store_true', dest='inner',
                                help='If specified exclude consecutive missing characters at either end of the sequence')
    parser_repeats.set_defaults(filter_func=filterRepeats)
    
    # Quality filter mode argument parser
    parser_quality = subparsers.add_parser('quality', parents=[parser_parent],
                                          formatter_class=ArgumentDefaultsHelpFormatter, 
                                          help='Quality filtering mode')
    parser_quality.add_argument('-q', action='store', dest='min_qual', type=float, 
                                default=default_min_qual, help='Quality score threshold')
    parser_quality.add_argument('--inner', action='store_true', dest='inner',
                                help='If specified exclude consecutive missing characters at either end of the sequence')
    parser_quality.set_defaults(filter_func=filterQuality)

    # Mask mode argument parser
    parser_maskqual = subparsers.add_parser('maskqual', parents=[parser_parent], 
                                        formatter_class=ArgumentDefaultsHelpFormatter,
                                        help='Character masking mode')
    parser_maskqual.add_argument('-q', action='store', dest='min_qual', type=float, 
                             default=default_min_qual, help='Quality score threshold')
    parser_maskqual.set_defaults(filter_func=maskQuality)

    # Trim mode argument parser
    parser_trimqual = subparsers.add_parser('trimqual', parents=[parser_parent], 
                                            formatter_class=ArgumentDefaultsHelpFormatter,
                                            help='Sequence trimming mode')
    parser_trimqual.add_argument('-q', action='store', dest='min_qual', type=float, 
                                 default=default_min_qual, help='Quality score threshold')
    parser_trimqual.add_argument('--win', action='store', dest='window', type=int, 
                                 default=default_window, 
                                 help='Nucleotide window size for moving average calculation')
    parser_trimqual.add_argument('--reverse', action='store_true', dest='reverse', 
                                 help='Specify to trim the head of the sequence rather than the tail')
    parser_trimqual.set_defaults(filter_func=trimQuality)
    
    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    
    # Create filter_args
    filter_keys = ['min_qual', 'max_repeat', 'max_missing', 'min_length', 'inner', 
                   'include_missing', 'window', 'reverse']
    args_dict['filter_args'] = dict((k, args_dict[k]) for k in args_dict if k in filter_keys)
    for k in args_dict['filter_args']:  del args_dict[k]
    
    # Calls quality processing function
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        filterSeq(**args_dict)
