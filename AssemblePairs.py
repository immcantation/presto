#!/usr/bin/env python
"""
Assembles paired-end reads into a single sequence
"""

__author__    = 'Jason Anthony Vander Heiden, Gur Yaari'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.1'
__date__      = '2013.10.12'

# Imports
import os, sys
import multiprocessing as mp
import numpy as np
import scipy.stats as stats
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from itertools import izip
from time import time
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_missing_chars, default_coord_choices, default_coord_type
from IgCore import default_delimiter, default_out_args
from IgCore import flattenAnnotation, mergeAnnotation, parseAnnotation
from IgCore import getCommonArgParser, parseCommonArgs
from IgCore import getFileType, getOutputHandle, printLog, printProgress
from IgCore import getScoreDict, reverseComplement, scoreSeqPair
from IgCore import getUnpairedIndex, indexSeqPairs, readSeqFile

# Defaults
default_alpha = 0.05
default_max_error = 0.2
default_min_len = 1
default_max_len = 1000
default_gap = 0


def getPMatrix(x):
    """
    Generates a matrix of p-values from a binomial distribution 

    Arguments: 
    x = maximum trials

    Returns:
    a numpy.array of successes by trials p-values 
    """
    p_matrix = np.ones([x, x], dtype=float) 
    a = np.array(range(x))
    for i in a:
        p_matrix[i, ] = 1 - stats.binom.cdf(i - 1, a, 0.25) - stats.binom.pmf(i, a, 0.25) / 2
    
    return p_matrix


def overlapConsensus(head_seq, tail_seq, ignore_chars=default_missing_chars):
    """
    Creates a consensus overlap sequences from two segments

    Arguments: 
    head_seq = the overlap head SeqRecord
    tail_seq = the overlap tail SeqRecord
    ignore_chars = list of characters which do not contribute to consensus
    
    Returns:
    A SeqRecord object with consensus characters and quality scores
    """
    # Initialize empty overlap character and quality score list
    seq_cons, score_cons = [], []
    # Define character and quality tuple iterators
    chars = zip(head_seq, tail_seq)
    quals = zip(head_seq.letter_annotations['phred_quality'], 
                 tail_seq.letter_annotations['phred_quality'])

    # Iterate over character and quality tuples and build consensus
    for c, q in izip(chars, quals):
        # Equivalent character case
        if c[0] == c[1]:
            c_cons = c[0]
            q_cons = max(q)
        # All ambiguous characters case
        elif all([x in ignore_chars for x in c]):
            c_cons = 'N'
            q_cons = max(q)
        # Some ambiguous characters case
        elif any([x in ignore_chars for x in c]):
            c_cons = [x for x in c if x not in ignore_chars][0]
            q_cons = q[c.index(c_cons)]
        # Conflicting character case        
        else:
            q_max = max(q)
            c_cons = c[q.index(q_max)]
            try:
                q_cons = q_max**2 / sum(q)
            except ZeroDivisionError:
                q_cons = 0
        # Append sequence and quality lists with consensus values
        seq_cons.append(c_cons)
        score_cons.append(q_cons)

    # Define overlap SeqRecord
    record = SeqRecord(Seq(''.join(seq_cons), IUPAC.ambiguous_dna), 
                       id='', 
                       name='', 
                       description='', 
                       letter_annotations={'phred_quality':score_cons})
    
    return record


def joinSeqPair(head_seq, tail_seq, gap=default_gap):
    """
    Concatenates two sequences 

    Arguments: 
    head_seq = the head SeqRecord
    tail_seq = the tail SeqRecord
    gap = number of gap characters to insert between head and tail
    
    Returns: 
    dictionary of {joined SeqRecord object, error rate, p-value, overlap length}
    """
    # Define joined ID
    join_id = head_seq.id if head_seq.id == tail_seq.id \
              else '+'.join([head_seq.id, tail_seq.id])
    join_seq = str(head_seq.seq) + '-' * gap + str(tail_seq.seq)
    
    # Define return record
    record = SeqRecord(Seq(join_seq, IUPAC.ambiguous_dna), 
                       id=join_id, 
                       name=join_id, 
                       description='')
    
    # Join quality score if present
    has_quality = hasattr(head_seq, 'letter_annotations') and \
                  hasattr(tail_seq, 'letter_annotations') and \
                  'phred_quality' in head_seq.letter_annotations and \
                  'phred_quality' in tail_seq.letter_annotations
    if has_quality:
        join_quality = head_seq.letter_annotations['phred_quality'] + \
                       [0] * gap + \
                       tail_seq.letter_annotations['phred_quality']
        record.letter_annotations = {'phred_quality':join_quality}
                               
    return {'seq':record, 'error':0, 'pval':0, 'len':-gap}


def alignSeqPair(head_seq, tail_seq, alpha=default_alpha, max_error=default_max_error, 
               min_len=default_min_len, max_len=default_max_len, 
               p_matrix=None, score_dict=getScoreDict(n_score=0, gap_score=0)):
    """
    Stitches two sequences together by aligning the ends

    Arguments: 
    head_seq = the head SeqRecord
    tail_seq = the tail SeqRecord
    alpha = maximum p-value for a valid stitch
    max_error = maximum error rate for a valid stitch
    min_len = minimum length of overlap to test
    max_len = maximum length of overlap to test
    p_matrix = optional successes by trials numpy.array of p-values 
    score_dict = optional dictionary of character scores in the 
                 form {(char1, char2): score}
                     
    Returns: 
    dictionary of {stitched SeqRecord object, error rate, p-value, overlap length}
    """
    # Define unpassed arguments
    if p_matrix is None:  p_matrix = getPMatrix(max_len + 1)
    
    # Define empty return dictionary
    best_dict = {'seq':None, 'error':1.0, 'pval':1.0, 'len':None}
    # Define general parameters
    head_len = len(head_seq)
    tail_len = len(tail_seq)
    max_len = min(head_len, tail_len, max_len)
    # Determine if quality scores are present
    has_quality = hasattr(head_seq, 'letter_annotations') and \
                  hasattr(tail_seq, 'letter_annotations') and \
                  'phred_quality' in head_seq.letter_annotations and \
                  'phred_quality' in tail_seq.letter_annotations


    # Iterate and score overlap segments
    pos = None
    head_str = str(head_seq.seq)
    tail_str = str(tail_seq.seq)
    for i in xrange(min_len, max_len + 1):
        score, weight, error = scoreSeqPair(head_str[-i:], tail_str[:i], max_error, i, 
                                            score_dict=score_dict)
        #print 'DEBUG> %s: score %i, weight %i' % (head_seq.id, score, weight)
        pval = p_matrix[score, weight]
        # Save stitch as optimal if p value improves and error criteria passed
        if pval <= alpha and pval < best_dict['pval'] and error <= max_error:
            pos = i
            best_dict['pval'] = pval
            best_dict['error'] = error
            
    # Build stitched sequences and assign best_dict values
    if pos is not None:
        best_dict['len'] = pos
        # Correct quality scores and resolve conflicts
        if has_quality:
            best_dict['seq'] = head_seq[:-pos] + \
                               overlapConsensus(head_seq[-pos:], tail_seq[:pos]) + \
                               tail_seq[pos:]
        # Assign head sequence to conflicts when no quality information is available
        else:
            best_dict['seq'] = head_seq + tail_seq[pos:]
        # Define best stitch ID
        best_dict['seq'].id = head_seq.id if head_seq.id == tail_seq.id \
                              else '+'.join([head_seq.id, tail_seq.id])
        best_dict['seq'].name = best_dict['seq'].id
        best_dict['seq'].description = ''

    return best_dict


def feedQueue(data_queue, nproc, seq_file_1, seq_file_2, index_dict):
    """
    Feeds the data queue with sequence pairs for processQueue processes

    Arguments: 
    data_queue = an multiprocessing.Queue to hold data for processing
    nproc = the number of processQueue processes
    seq_file_1 = the name of sequence file 1
    seq_file_2 = the name of sequence file 2
    index_dict = a dictionary returned by indexSeqPairs
         
    Returns: 
    None
    """
    # Open input files
    seq_dict_1 = readSeqFile(seq_file_1, index=True)
    seq_dict_2 = readSeqFile(seq_file_2, index=True)
    
    # Iterate over index_dict added sequences pairs to data_queue
    for pair_id, (key_1, key_2) in index_dict.iteritems():
        # Feed queue
        data_queue.put({'id':pair_id, 
                        'seq_1':seq_dict_1[key_1], 
                        'seq_2':seq_dict_2[key_2]})
        
    # Add sentinel object for each processQueue process
    for __ in range(nproc):
        data_queue.put(None)

    return None


def processQueue(data_queue, result_queue, assemble_func, assemble_args={}, rc=None, 
                 head_fields=None, tail_fields=None, delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    assemble_func = the function to use to assemble paired ends
    assemble_args = a dictionary of arguments to pass to the assembly function
    rc = Defines which sequences ('head','tail','both') to reverse complement before assembly;
         if None do not reverse complement sequences
    head_fields = list of annotations in head_file records to copy to assembled record;
                  if None do not copy an annotation
    tail_fields = list of annotations in tail_file records to copy to assembled record;
                  if None do not copy an annotation
    delimiter = a tuple of delimiters for (fields, values, value lists) 

    Returns: 
    None
    """
    # Iterator over data queue until sentinel object reached
    for args in iter(data_queue.get, None):
        # Reverse complement sequences if required
        if rc == 'both' or rc == 'head':  head_seq = reverseComplement(args['seq_1'])
        else:  head_seq = args['seq_1']
        if rc == 'both' or rc == 'tail':  tail_seq = reverseComplement(args['seq_2'])
        else:  tail_seq = args['seq_2']
        
        # Define result dictionary for iteration
        results = {'id':args['id'],
                   'in_seq_1':head_seq,
                   'in_seq_2':tail_seq,
                   'out_seq':None,
                   'error':None,
                   'pval':None,
                   'pass':False,
                   'log':OrderedDict([('ID', args['id'])])}
 
        # Assemble sequences
        stitch = assemble_func(head_seq, tail_seq, **assemble_args)
        
        # Define stitched sequence annotation
        stitch_ann = OrderedDict([('ID', args['id'])])                  
        if head_fields is not None:
            head_ann = parseAnnotation(head_seq.description, head_fields, delimiter=delimiter)
            stitch_ann = mergeAnnotation(stitch_ann, head_ann, delimiter=delimiter)
            results['log']['HEADFIELDS'] = '|'.join(['%s=%s' % (k, v) for k, v in head_ann.iteritems()])
        if tail_fields is not None:
            tail_ann = parseAnnotation(tail_seq.description, tail_fields, delimiter=delimiter)
            stitch_ann = mergeAnnotation(stitch_ann, tail_ann, delimiter=delimiter)
            results['log']['TAILFIELDS'] = '|'.join(['%s=%s' % (k, v) for k, v in tail_ann.iteritems()])
        
        # Define stitching log
        results['log']['HEADSEQ'] = head_seq.seq
        if stitch['seq'] is not None:
            out_seq = stitch['seq'] 
            # Update stitch annotation
            out_seq.id = flattenAnnotation(stitch_ann, delimiter=delimiter)
            out_seq.name = out_seq.id
            out_seq.description = '' 
            # Add stitch to results
            results['out_seq'] = out_seq
            results['pass'] = True
            # Update log
            results['log']['TAILSEQ'] = ' ' * (len(head_seq) - stitch['len']) + tail_seq.seq
            results['log']['ASSEMBLY'] = out_seq.seq
            if 'phred_quality' in out_seq.letter_annotations:
                results['log']['QUALITY'] = ''.join([chr(q+33) for q in out_seq.letter_annotations['phred_quality']])
            results['log']['LENGTH'] = len(out_seq)
            results['log']['OVERLAP'] = stitch['len']
        else:
            results['log']['TAILSEQ'] = ' ' * len(head_seq) + tail_seq.seq
            results['log']['ASSEMBLY'] = None
        results['log']['ERROR'] = stitch['error']
        results['log']['PVAL'] = stitch['pval']
            
        # Feed results to result queue
        result_queue.put(results)

    return None


def collectQueue(result_queue, collect_dict, seq_file_1, seq_file_2, out_args):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments: 
    result_queue = a multiprocessing.Queue holding processQueue results
    collect_dict = a multiprocessing.Manager.dict to store return values
    seq_file_1 = the first sequence file name
    seq_file_2 = the second sequence file name
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns: 
    None
    (adds 'log' and 'out_files' to collect_dict)
    """
    # Count records and define output format 
    out_type = getFileType(seq_file_1) if out_args['out_type'] is None \
               else out_args['out_type']
    result_count = collect_dict['result_count']
    
    # Defined valid assembly output handle
    pass_handle = getOutputHandle(seq_file_1, 
                                  'assemble-pass', 
                                  out_dir=out_args['out_dir'], 
                                  out_name=out_args['out_name'], 
                                  out_type=out_type)
    # Defined failed assembly output handles
    if out_args['clean']:
        fail_handle_1 = None
        fail_handle_2 = None
    else:
        # Define output name
        if out_args['out_name'] is None:
            out_name_1 = out_name_2 = None
        else: 
            out_name_1 = '%s-1' % out_args['out_name']
            out_name_2 = '%s-2' % out_args['out_name']
        fail_handle_1 = getOutputHandle(seq_file_1, 
                                        'assemble-fail', 
                                        out_dir=out_args['out_dir'], 
                                        out_name=out_name_1, 
                                        out_type=out_type)
        fail_handle_2 = getOutputHandle(seq_file_2, 
                                        'assemble-fail', 
                                        out_dir=out_args['out_dir'], 
                                        out_name=out_name_2, 
                                        out_type=out_type)

    # Define log handle
    if out_args['log_file'] is None:
        log_handle = None
    else:
        log_handle = open(out_args['log_file'], 'w')
        
    # Iterator over results queue until sentinel object reached
    start_time = time()
    pair_count = pass_count = fail_count = 0
    for result in iter(result_queue.get, None): 
        # Print progress for previous iteration
        printProgress(pair_count, result_count, 0.05, start_time)

        # Update counts for iteration
        pair_count += 1

        # Write log
        printLog(result['log'], handle=log_handle)

        # Write assembled sequences
        if result['out_seq'] is not None and result['pass']:
            pass_count += 1
            SeqIO.write(result['out_seq'], pass_handle, out_type)
        else:
            fail_count += 1
            if fail_handle_1 is not None and fail_handle_2 is not None:
                SeqIO.write(result['in_seq_1'], fail_handle_1, out_type)
                SeqIO.write(result['in_seq_2'], fail_handle_2, out_type)
 
    # Print total counts
    printProgress(pair_count, result_count, 0.05, start_time)

    # Update return values
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['PAIRS'] = pair_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    collect_dict['log'] = log
    collect_dict['out_files'] = [pass_handle.name]

    # Close file handles
    pass_handle.close()
    if fail_handle_1 is not None:  fail_handle_1.close()
    if fail_handle_2 is not None:  fail_handle_2.close()
    if log_handle is not None:  log_handle.close()

    return None


def assemblePairs(head_file, tail_file, assemble_func, assemble_args={}, 
                  coord_type=default_coord_type, rc=None, 
                  head_fields=None, tail_fields=None,  
                  out_args=default_out_args, nproc=None, queue_size=None):
    """
    Generates consensus sequences

    Arguments: 
    head_file = the head sequence file name
    tail_file = the tail sequence file name
    assemble_func = the function to use to assemble paired ends
    assemble_args = a dictionary of arguments to pass to the assembly function
    coord_type = the sequence header format
    rc = Defines which sequences ('head','tail','both') to reverse complement before assembly;
         if None do not reverse complement sequences
    head_fields = list of annotations in head_file records to copy to assembled record;
                  if None do not copy an annotation
    tail_fields = list of annotations in tail_file records to copy to assembled record;
                  if None do not copy an annotation
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
    
    # Define subcommand label dictionary
    cmd_dict = {alignSeqPair:'align', joinSeqPair:'join'}

    # Print parameter info
    log = OrderedDict()
    log['START'] = 'AssemblePairs'
    log['COMMAND'] = cmd_dict.get(assemble_func, assemble_func.__name__)
    log['HEAD_FILE'] = os.path.basename(head_file) 
    log['TAIL_FILE'] = os.path.basename(tail_file)
    log['COORD_TYPE'] = coord_type
    if 'alpha' in assemble_args:  log['ALPHA'] = assemble_args['alpha']
    if 'max_error' in assemble_args:  log['MAX_ERROR'] = assemble_args['max_error']
    if 'min_len' in assemble_args:  log['MIN_LEN'] = assemble_args['min_len']
    if 'max_len' in assemble_args:  log['MAX_LEN'] = assemble_args['max_len']
    if 'gap' in assemble_args:  log['GAP'] = assemble_args['gap']
    log['NPROC'] = nproc
    printLog(log)

    # Read input files
    head_type = getFileType(head_file)
    head_dict = readSeqFile(head_file, index=True)
    head_count = len(head_dict) 
    tail_type = getFileType(tail_file)
    tail_dict = readSeqFile(tail_file, index=True)
    tail_count = len(tail_dict)

    # Find paired sequences
    index_dict = indexSeqPairs(head_dict, tail_dict, coord_type, out_args['delimiter'])
    pair_count = len(index_dict)
    
    # Write unmatched entries to files
    if not out_args['clean']:
        # Define output type
        if out_args['out_type'] is not None:
            head_type = tail_type = out_args['out_type']        
        # Define output name
        if out_args['out_name'] is None:
            head_name = tail_name = None
        else: 
            head_name = '%s-1' % out_args['out_name']
            tail_name = '%s-2' % out_args['out_name']
        
        # Find unpaired sequences
        head_unpaired, tail_unpaired = getUnpairedIndex(head_dict, tail_dict, coord_type, 
                                                        out_args['delimiter'])
        # Write unpaired head records
        with getOutputHandle(head_file, 'assemble-unpaired', out_dir=out_args['out_dir'], 
                             out_name=head_name, out_type=head_type) as head_handle:
            for k in head_unpaired:
                SeqIO.write(head_dict[k], head_handle, head_type)

        # Write unpaired tail records
        with getOutputHandle(tail_file, 'assemble-unpaired', out_dir=out_args['out_dir'], 
                             out_name=tail_name, out_type=tail_type) as tail_handle:
            for k in tail_unpaired:
                SeqIO.write(tail_dict[k], tail_handle, tail_type)
    
    # Define shared data objects 
    manager = mp.Manager()
    collect_dict = manager.dict()
    data_queue = mp.Queue(queue_size)
    result_queue = mp.Queue(queue_size)
    
    # Initiate feeder process
    feeder = mp.Process(target=feedQueue, args=(data_queue, nproc, head_file, tail_file, 
                                                index_dict))
    feeder.start()

    # Initiate processQueue processes
    workers = []
    for __ in range(nproc):
        w = mp.Process(target=processQueue, args=(data_queue, result_queue, assemble_func,
                                                  assemble_args, rc, head_fields, tail_fields, 
                                                  out_args['delimiter']))
        w.start()
        workers.append(w)

    # Initiate collector process
    collect_dict['result_count'] = pair_count
    collector = mp.Process(target=collectQueue, args=(result_queue, collect_dict,
                                                      head_file, tail_file, out_args))
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
    log_end = OrderedDict()
    log_end['OUTPUT'] = log.pop('OUTPUT')
    log_end['SEQUENCES'] = '%i,%i' % (head_count, tail_count)
    log_end['UNPAIRED'] = '%i,%i' % (head_count - pair_count, 
                                     tail_count - pair_count)
    for k, v in log.iteritems():
        log_end[k] = v
    log_end['END'] = 'AssemblePairs'
    printLog(log_end)
        
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
    subparsers = parser.add_subparsers(dest='command', help='Assembly mode')
    
    # Parent parser    
    parser_parent = getCommonArgParser(paired=True, multiproc=True)
    parser_parent.add_argument('--coord', action='store', dest='coord_type', 
                               choices=default_coord_choices, default=default_coord_type,
                               help='The format of the sequence identifier which defines shared coordinate \
                                     information across paired ends')
    parser_parent.add_argument('--rc', action='store', dest='rc', choices=('head', 'tail', 'both'), 
                               default=None, help='Specify to reverse complement sequences before stitching')
    parser_parent.add_argument('--1f', nargs='+', action='store', dest='head_fields', type=str, default=None, 
                               help='Specify annotation fields to copy from head records into assembled record')
    parser_parent.add_argument('--2f', nargs='+', action='store', dest='tail_fields', type=str, default=None, 
                               help='Specify annotation fields to copy from tail records into assembled record')

    
    # Paired end overlap alignment mode argument parser
    parser_align = subparsers.add_parser('align', parents=[parser_parent],
                                         formatter_class=ArgumentDefaultsHelpFormatter,
                                         help='Assembled pairs by aligning ends')
    parser_align.add_argument('--alpha', action='store', dest='alpha', type=float,
                              default=default_alpha, help='Significance threshold for sequence assemble')
    parser_align.add_argument('--maxerror', action='store', dest='max_error', type=float,
                              default=default_max_error, help='Maximum allowable error rate')
    parser_align.add_argument('--minlen', action='store', dest='min_len', type=int,
                              default=default_min_len, help='Minimum sequence length to scan for overlap')
    parser_align.add_argument('--maxlen', action='store', dest='max_len', type=int,
                              default=default_max_len, help='Maximum sequence length to scan for overlap')
    parser_align.set_defaults(assemble_func=alignSeqPair)
    
    # Paired end concatenation mode argument parser
    parser_join = subparsers.add_parser('join', parents=[parser_parent],
                                         formatter_class=ArgumentDefaultsHelpFormatter,
                                         help='Assembled pairs by concatenating ends')
    parser_join.add_argument('--gap', action='store', dest='gap', type=int, default=default_gap, 
                             help='Number of gap characters to place between ends')
    parser_join.set_defaults(assemble_func=joinSeqPair)    
    
    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    
    # Convert case of fields
    if args_dict['head_fields']:  args_dict['head_fields'] = map(str.upper, args_dict['head_fields']) 
    if args_dict['tail_fields']:  args_dict['tail_fields'] = map(str.upper, args_dict['tail_fields']) 
    
    # Define assemble_args dictionary to pass to maskPrimers
    if args_dict['assemble_func'] is alignSeqPair:
        args_dict['assemble_args'] = {'alpha':args_dict['alpha'],
                                      'max_error':args_dict['max_error'],
                                      'min_len':args_dict['min_len'],
                                      'max_len':args_dict['max_len'],
                                      'p_matrix':getPMatrix(args_dict['max_len'] + 1)}
        del args_dict['alpha']
        del args_dict['max_error']
        del args_dict['min_len']
        del args_dict['max_len']
    elif args_dict['assemble_func'] is joinSeqPair:
        args_dict['assemble_args'] = {'gap':args_dict['gap']}
        del args_dict['gap']
    
    # Call assemblePairs for each sample file
    del args_dict['command']
    del args_dict['seq_files_1']
    del args_dict['seq_files_2']
    for head, tail in zip(args.__dict__['seq_files_1'], args.__dict__['seq_files_2']):
        args_dict['head_file'] = head
        args_dict['tail_file'] = tail
        assemblePairs(**args_dict)
        
        # Profiling
        #import cProfile, pstats
        #cProfile.run('assemblePairs(**args_dict)', 'profile.prof')
        #p = pstats.Stats('profile.prof')
        #p.strip_dirs().sort_stats('time').print_stats()       
