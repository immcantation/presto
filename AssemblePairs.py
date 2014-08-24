#!/usr/bin/env python
"""
Assembles paired-end reads into a single sequence
"""

__author__    = 'Jason Anthony Vander Heiden, Gur Yaari'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.4'
__date__      = '2014.6.10'

# Imports
import os, sys
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
from IgCore import manageProcesses, SeqData, SeqResult

# Defaults
default_alpha = 0.01
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
    # Define undefined arguments
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


def feedAPQueue(alive, data_queue, seq_file_1, seq_file_2, index_dict):
    """
    Feeds the data queue with sequence pairs for processQueue processes

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = an multiprocessing.Queue to hold data for processing
    seq_file_1 = the name of sequence file 1
    seq_file_2 = the name of sequence file 2
    index_dict = a dictionary returned by indexSeqPairs
         
    Returns: 
    None
    """
    try:
        # Open input files
        seq_dict_1 = readSeqFile(seq_file_1, index=True)
        seq_dict_2 = readSeqFile(seq_file_2, index=True)
        # Define data iterator
        data_iter = ((k, [seq_dict_1[i], seq_dict_2[j]]) \
                     for k, (i, j) in index_dict.iteritems())
    except:
        alive.value = False
        raise
    
    try:
        # Iterate over data_iter and feed data queue 
        while alive.value:
            # Get data from queue
            if data_queue.full():  continue
            else:  data = next(data_iter, None)
            # Exit upon reaching end of iterator
            if data is None:  break
            
            # Feed queue
            data_queue.put(SeqData(*data))
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        raise

    return None


def processAPQueue(alive, data_queue, result_queue, assemble_func, assemble_args={}, 
                   rc=None, fields_1=None, fields_2=None, 
                   delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    assemble_func = the function to use to assemble paired ends
    assemble_args = a dictionary of arguments to pass to the assembly function
    rc = Defines which sequences ('head','tail','both') to reverse complement 
         before assembly; if None do not reverse complement sequences
    fields_1 = list of annotations in head_file records to copy to assembled record;
               if None do not copy an annotation
    fields_2 = list of annotations in tail_file records to copy to assembled record;
               if None do not copy an annotation
    delimiter = a tuple of delimiters for (fields, values, value lists) 

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
            
            # Reverse complement sequences if required  
            head_seq = data.data[0] if rc not in ('head', 'both') \
                       else reverseComplement(data.data[0])
            tail_seq = data.data[1] if rc not in ('tail', 'both') \
                       else reverseComplement(data.data[1])

            # Define result object for iteration
            result = SeqResult(data.id, [head_seq, tail_seq])
            
            # Assemble sequences
            stitch = assemble_func(head_seq, tail_seq, **assemble_args)
            
            # Define stitched sequence annotation
            stitch_ann = OrderedDict([('ID', data.id)])                  
            if fields_1 is not None:
                head_ann = parseAnnotation(head_seq.description, fields_1, 
                                           delimiter=delimiter)
                stitch_ann = mergeAnnotation(stitch_ann, head_ann, delimiter=delimiter)
                result.log['HEADFIELDS'] = '|'.join(['%s=%s' % (k, v) 
                                                     for k, v in head_ann.iteritems()])
            if fields_2 is not None:
                tail_ann = parseAnnotation(tail_seq.description, fields_2, 
                                           delimiter=delimiter)
                stitch_ann = mergeAnnotation(stitch_ann, tail_ann, delimiter=delimiter)
                result.log['TAILFIELDS'] = '|'.join(['%s=%s' % (k, v) 
                                                     for k, v in tail_ann.iteritems()])
            
            # Define stitching log
            result.log['HEADSEQ'] = head_seq.seq
            if stitch['seq'] is not None:
                out_seq = stitch['seq'] 
                # Update stitch annotation
                out_seq.id = flattenAnnotation(stitch_ann, delimiter=delimiter)
                out_seq.name = out_seq.id
                out_seq.description = '' 
                # Add stitch to results
                result.results = out_seq
                result.valid = True
                # Update log
                result.log['TAILSEQ'] = ' ' * (len(head_seq) - stitch['len']) + tail_seq.seq
                result.log['ASSEMBLY'] = out_seq.seq
                if 'phred_quality' in out_seq.letter_annotations:
                    result.log['QUALITY'] = ''.join([chr(q+33) for q in out_seq.letter_annotations['phred_quality']])
                result.log['LENGTH'] = len(out_seq)
                result.log['OVERLAP'] = stitch['len']
            else:
                result.log['TAILSEQ'] = ' ' * len(head_seq) + tail_seq.seq
                result.log['ASSEMBLY'] = None
            result.log['ERROR'] = stitch['error']
            result.log['PVAL'] = stitch['pval']
                
            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        sys.stderr.write('Error processing sequence pair with ID: %s.\n' % data.id)
        raise
    
    return None


def collectAPQueue(alive, result_queue, collect_queue, result_count, seq_file_1, seq_file_2, 
                   out_args):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    result_queue = a multiprocessing.Queue holding worker results
    collect_queue = a multiprocessing.Queue holding collector return values
    result_count = the number of expected assembled sequences
    seq_file_1 = the first sequence file name
    seq_file_2 = the second sequence file name
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns: 
    None
    (adds a dictionary of {log: log object, out_files: output file names} to collect_queue)
    """
    try:
        # Count records and define output format 
        out_type = getFileType(seq_file_1) if out_args['out_type'] is None \
                   else out_args['out_type']
        
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
    except:
        alive.value = False
        raise
    
    try:
        # Iterator over results queue until sentinel object reached
        start_time = time()
        iter_count = pass_count = fail_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break

            # Print progress for previous iteration
            printProgress(iter_count, result_count, 0.05, start_time)
    
            # Update counts for iteration
            iter_count += 1
    
            # Write log
            printLog(result.log, handle=log_handle)
    
            # Write assembled sequences
            if result:
                pass_count += 1
                SeqIO.write(result.results, pass_handle, out_type)
            else:
                fail_count += 1
                if fail_handle_1 is not None and fail_handle_2 is not None:
                    SeqIO.write(result.data[0], fail_handle_1, out_type)
                    SeqIO.write(result.data[1], fail_handle_2, out_type)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
        
        # Print total counts
        printProgress(iter_count, result_count, 0.05, start_time)
    
        # Update return values
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name)
        log['PAIRS'] = iter_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count
        collect_dict = {'log':log, 'out_files': [pass_handle.name]}
        collect_queue.put(collect_dict)
        
        # Close file handles
        pass_handle.close()
        if fail_handle_1 is not None:  fail_handle_1.close()
        if fail_handle_2 is not None:  fail_handle_2.close()
        if log_handle is not None:  log_handle.close()
    except:
        alive.value = False
        raise
    
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
    # Define subcommand label dictionary
    cmd_dict = {alignSeqPair:'align', joinSeqPair:'join'}

    # Print parameter info
    log = OrderedDict()
    log['START'] = 'AssemblePairs'
    log['COMMAND'] = cmd_dict.get(assemble_func, assemble_func.__name__)
    log['FILE1'] = os.path.basename(head_file) 
    log['FILE2'] = os.path.basename(tail_file)
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
    
    # Define feeder function and arguments
    feed_func = feedAPQueue
    feed_args = {'seq_file_1': head_file,
                 'seq_file_2': tail_file,
                 'index_dict': index_dict}
    # Define worker function and arguments
    work_func = processAPQueue
    work_args = {'assemble_func': assemble_func, 
                 'assemble_args': assemble_args,
                 'rc': rc,
                 'fields_1': head_fields,
                 'fields_2': tail_fields,
                 'delimiter': out_args['delimiter']}
    # Define collector function and arguments
    collect_func = collectAPQueue
    collect_args = {'result_count': pair_count,
                    'seq_file_1': head_file,
                    'seq_file_2': tail_file,
                    'out_args': out_args}
                   
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    log = OrderedDict()
    log['OUTPUT'] = result['log'].pop('OUTPUT')
    log['SEQUENCES1'] = head_count
    log['SEQUENCES2'] = tail_count
    log['UNPAIRED1'] = head_count - pair_count 
    log['UNPAIRED2'] = tail_count - pair_count
    for k, v in result['log'].iteritems():  log[k] = v
    log['END'] = 'AssemblePairs'
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
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, version='%(prog)s:' + ' v%s-%s' %(__version__, __date__), 
                            formatter_class=ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', dest='command', help='Assembly method', metavar='')
    
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
    for head, tail in zip(args.__dict__['seq_files_1'], 
                          args.__dict__['seq_files_2']):
        args_dict['head_file'] = head
        args_dict['tail_file'] = tail
        assemblePairs(**args_dict)
            
            
            
