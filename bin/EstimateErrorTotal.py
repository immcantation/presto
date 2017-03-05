#!/usr/bin/env python3
"""
Calculates distance distribution from pairwise sequence comparisons within annotation sets 
"""
# Info
__author__ = 'Ruoyi Jiang, Jason Anthony Vander Heiden, Namita Gupta'
from presto import __version__, __date__

# Imports
import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import permutations
from textwrap import dedent
from time import time

# Presto imports
from presto.Defaults import default_barcode_field, default_missing_chars, default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.IO import getFileType, countSeqSets, getOutputHandle, printLog, printProgress
from presto.Sequence import getDNAScoreDict, indexSeqSets
from presto.Multiprocessing import SeqResult, manageProcesses, feedSeqQueue



#####NEED TO UPDATE API TO INCORPORATE CHANGES

from itertools import product
##Addition of invert flag
def scoreDNA(a, b, mask_score=None, gap_score=None, invert=False):
    """
    Returns the score for a pair of IUPAC Ambiguous Nucleotide characters

    Arguments:
      a : First characters
      b : Second character
      n_score : Tuple of length two defining scores for all matches against an N
                character for (a, b), with the score for character (a) taking precedence;
                if None score symmetrically according to IUPAC character identity
      gap_score : Tuple of length two defining score for all matches against a gap (-, .)
                  character for (a, b), with the score for character (a) taking precedence;
                  if None score symmetrically according to IUPAC character identity

    Returns:
      int : Score for the character pair
    """
    # Define ambiguous character translations
    IUPAC_trans = {'AGWSKMBDHV':'R', 'CTSWKMBDHV':'Y', 'CGKMBDHV':'S', 'ATKMBDHV':'W', 'GTBDHV':'K',
                   'ACBDHV':'M', 'CGTDHV':'B', 'AGTHV':'D', 'ACTV':'H', 'ACG':'V', 'ABCDGHKMRSTVWY':'N',
                   '-.':'.'}
    # Create list of tuples of synonymous character pairs
    IUPAC_matches = [p for k, v in IUPAC_trans.items() for p in list(product(k, v))]
    # Check gap and N-value conditions, prioritizing score for first character
    if gap_score is not None and a in '-.':
        return gap_score[0]
    elif mask_score is not None and a in 'nN':
        return mask_score[0]
    elif gap_score is not None and b in '-.':
        return gap_score[1]
    elif mask_score is not None and b in 'nN':
        return mask_score[1]
    # Return symmetric and reflexive score for IUPAC match conditions
    if not invert:
        if a == b or (a, b) in IUPAC_matches or (b, a) in IUPAC_matches:
            return 1
        else:
            return 0
    else:
        if a == b or (a, b) in IUPAC_matches or (b, a) in IUPAC_matches:
            return 0
        else:
            return 1

from itertools import product
def getDNAScoreDict(mask_score=None, gap_score=None, invert =False):
    """
    Generates a score dictionary

    Arguments:
      mask_score : Tuple of length two defining scores for all matches against an N
                   character for (a, b), with the score for character (a) taking precedence;
                   if None score symmetrically according to IUPAC character identity
      gap_score : Tuple of length two defining score for all matches against a [-, .]
                  character for (a, b), with the score for character (a) taking precedence;
                  if None score symmetrically according to IUPAC character identity

    Returns:
      dict : Score dictionary with keys (char1, char2) mapping to scores
    """
    chars = '-.ACGTRYSWKMBDHVN'
    score_dict = {k:scoreDNA(*k, mask_score=mask_score, gap_score=gap_score, invert=invert)
                  for k in product(chars, repeat=2)}
    return score_dict


default_score_dict = getDNAScoreDict(invert=True)


def scoreHamPair(seq1, seq2, score_dict=default_score_dict):
    """
    Simple hamming calculator derived from scoreSeqPair

    Arguments:
        seq1 : string representing an nt sequence with valid chars
        seq2 : string representing an nt sequence with valid chars
        score_dict : optional dictionary of alignment scores

    Returns:
          dict : Score dictionary with keys (char1, char2) mapping to scores
    """
    nts = zip(seq1, seq2)
    score = sum([score_dict[(c1, c2)] for c1, c2 in nts])
    
    return score


from itertools import combinations

def calcDistancesPairwise(sequences, score_dict=default_score_dict):
    """
    Calculate pairwise distances between input sequences (currently hamming only)
    
    Arguments:
        sequences: List of sequences for which to calculate pairwise distances
        score_dict : optional dictionary of alignment scores
    
    Returns:
        ndarray: numpy matrix of pairwise distances between input sequneces
    """
    #Initialize output distance matrix 
    dists = np.zeros((len(sequences), len(sequences)))
    
    #Iterate over combinations of input sequences
    for j,k in combinations(list(range(len(sequences))), 2):
    
        seq1 = str(sequences[j].seq)
        seq2 = str(sequences[k].seq)
        
        min_len = min(len(seq1), len(seq2))
        
        #Calculate distances
        try:
        	#TODO: hamming calculator ignores differences in length
            dists[j, k] = dists[k, j] = scoreHamPair(seq1, seq2, score_dict) 
        except (KeyError):
            raise KeyError('Unrecognized character in sequence.')
        
    return dists


def histogramDistMatrix(array, max_dist, triangle=True):
    """
    Bins the output distance matrix from pairwise distance calculations
    
    Arguments:
        array: a numpy matrix of pairwise distances between input sequences
        max_dist: the maximum plausible hamming distance within the pairwise comparisons
        triangle: True/False for whether to flatten the whole matrix or just the upper triangular
    
    Returns:
        output_hist: a histogram/count distribution of hamming distances 
    """
	#we assume the input array is a square matrix/array
    n = len(array)
	#check that the array is not empty
    if not n:
        raise Exception('barcode/annotation group with no sequences encountered')
	#flattens whole or upper triangular part of matrix (based on triangle flag)
    output = []
    if triangle:
    	for i in range(n):
    		for j in range(i+1, n):
    			output.append(array[i,j])
    else:
    	for i in range(n):
    		for j in range(n):
    			output.append(array[i,j])
	#check that the hamming distance resulted in plausible distances
    if max(output) > max_dist:
        raise Exception('distance observed between sequences is greater than maximum plausible distance')
	#generate the histogram from the flattened matrix
    output_hist = list(np.histogram(output, bins=list(range(max_dist)))[0])
    return output_hist


def logHistogram(hist):
    """
    Prints a concise string summarizing histogram/count distribution (hist, a list)
    
    Arguments:
        hist: a list of ints corresponding to a histogram/count distribution
    
    Returns:
        output_string: a string representing the simplified form of hist
    """
    index = 0
    output_string = ''
    for count in hist:
        if count:
            output_string += str(index) + ':' + str(count) + ' '
        index += 1
    return output_string








def processEETQueue(alive, data_queue, result_queue):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    cons_func = the function to use for consensus generation 
    cons_args = a dictionary of optional arguments for the consensus function
    min_count = threshold number of sequences to retain a set
    max_diversity = the minimum diversity score to retain a set;
                    if None do not calculate diversity
                        
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
            
            # Define sequences set
            seq_list = data.data
            seq_id = data.id
            seq_count = len(seq_list)
            # Gets the length of the shortest seq in the list
            seq_min_len = min([len(seq) for seq in seq_list])

            result = SeqResult(seq_id, seq_list)

            # Update log
            result.log['SET'] = seq_id
            result.log['SEQCOUNT'] = seq_count
            result.log['MAXLEN'] = seq_min_len
            for i, s in enumerate(seq_list):
                result.log['SEQ%i' % (i + 1)] = str(s.seq)

            #Calculate within UMI group distance matrix and bin into a density distribution
            if not seq_count:
            	raise Exception('barcode/annotation group with no associated SeqRecords encountered.') 
            elif seq_count is 1:
                result.log['HIST'] = 0
            else:
                try: 
                    hist = histogramDistMatrix(calcDistancesPairwise(seq_list), seq_min_len)
                except:
                    raise
                #Update log with distance distribution
                str_hist = logHistogram(hist)
                result.log['HIST'] = str_hist
                
                #Add the result of the histogram
                result.results = {'id': seq_id, 'hist': hist}
                
                
            # Update results and feed result queue (including results without hists)
            result.valid = (seq_count > 0)
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        sys.stderr.write('Error processing sequence set with ID: %s.\n' % data.id)
        raise

    return None



def collectEETQueue(alive, result_queue, collect_queue, seq_file, out_args, set_field):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    result_queue = a multiprocessing.Queue holding worker results
    collect_queue = a multiprocessing.Queue to store collector return values
    seq_file = the sample sequence file name
    out_args = common output argument dictionary from parseCommonArgs
    set_field = the field defining set membership

    Returns:
    None
    (adds a dictionary of {log: log object, out_files: output file names} to collect_queue)
    """
    try:
        # Count sets in file
        #This takes a very long time
        result_count = countSeqSets(seq_file, set_field, out_args['delimiter'])
        
        ##### Create an empty one row dataframe (with variable col)

        empty_df = pd.DataFrame(None, dtype = int)
        hist_df = empty_df

        # # Define empty DataFrames to store assembled results
        # pos_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'],
        #                       dtype=float)
        # qual_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'],
        #                        dtype=float)
        # set_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'],
        #                       dtype=float)
        # #nuc_pairs = list(permutations(['A', 'C', 'G', 'T'], 2))
        # #nuc_index = pd.MultiIndex.from_tuples(nuc_pairs, names=['obs', 'ref'])
        # nuc_index = pd.MultiIndex(levels=[[], []], labels=[[], []], names=['obs', 'ref'])
        # nuc_df = pd.DataFrame(None, index=nuc_index,
        #                       columns=['mismatch', 'q_sum', 'total'],
        #                       dtype=float)

        # Open log file
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
        #########@@
        set_count = seq_count = pass_count = fail_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  
            	break

            # Print progress for previous iteration
            printProgress(set_count, result_count, 0.05, start_time)
            
            # Update counts for progress logs
            set_count += 1
            seq_count += result.data_count
            
            # only add results that contain a hist to the eventual output hist_df
            if result.results is not None:
                result_id = result.results['id']
                result_hist = result.results['hist']
                row_df = pd.DataFrame(np.array([result_hist]), index = [result_id], dtype=int)

                try:
                	hist_df = pd.concat([hist_df, row_df])
                	pass_count += 1
                except:
                	raise
            else:
                # fail if a result is a single sequence UMI group (no histogram)
                fail_count += 1

            # Write log
            printLog(result.log, handle=log_handle)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
        
        ###Exit MP loop

        # Print final progress
        printProgress(set_count, result_count, 0.05, start_time)
        
        #Eliminate NA from hist_df
        hist_df = hist_df.fillna(0)

        #Collapse a version of hist_df -> sum_df. Transpose the df.
       	sum_df = pd.DataFrame(hist_df.sum(), columns=['Sum']).transpose()

        #Check if threshold can be found from sum_df. Generate input for the function (list)
        threshold_calc_input = [i for i in sum_df.as_matrix()[0]]

        #TODO: Find threshold statistics on the log scale version of hist_df

        # Generate log
        log = OrderedDict()
        for i in range(2): 
            log['OUTPUT%i' % (i + 1)] = None
        log['SETS'] = set_count
        log['SEQUENCES'] = seq_count

        log['PASS'] = pass_count
        log['FAIL'] = fail_count

        #log['THRESHOLD'] = None 
    
        # Build results dictionary
        assembled = {'hist':hist_df, 'sum':sum_df}

        # Write assembled error counts to output files
        out_files = writeResults(assembled, seq_file, out_args)
        
        # Update log
        for i, f in enumerate(out_files):
            log['OUTPUT%i' % (i + 1)] = os.path.basename(f)

        #log['THRESHOLD'] = None 

        # Update collector results
        collect_dict = {'log':log, 'out_files': out_files}
        collect_queue.put(collect_dict)
    except:
        alive.value = False
        raise
    
    return None


def writeResults(results, seq_file, out_args):
    """
    Formats results and writes to output files

    Arguments: 
    results = assembled results dictionary
    seq_file = the sample sequence file name
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    a tuple of (position error, nucleotide pairwise error, quality error, sequence set) file names
    """
    hist_df = results['hist']
    sum_df = results['sum']

    # Type conversion to int of all values
    hist_df = hist_df.astype(int) 
    sum_df = sum_df.astype(int) 

    
    # Write to tab delimited files
    file_args = {'out_dir':out_args['out_dir'], 'out_name':out_args['out_name'], 'out_type':'tab'}
    
    with getOutputHandle(seq_file, 'dist-all', **file_args) as hist_handle, \
            getOutputHandle(seq_file, 'dist-sum', **file_args) as sum_handle: \

        hist_df.to_csv(hist_handle, sep='\t', index=True, header = True)
        sum_df.to_csv(sum_handle, sep='\t', index = False, header = False) 

    return (hist_handle.name, sum_handle.name)





def estimateErrorTotal(seq_file,  
                  set_field=default_barcode_field,  
                  out_args=default_out_args, nproc=None, queue_size=None):
    """
    Calculates error rates of sequence sets

    Arguments: 
    seq_file = the sample sequence file name
    cons_func = the function to use for consensus generation 
    cons_args = a dictionary of arguments for the consensus function
    set_field = the annotation field containing set IDs
    min_count = threshold number of sequences to consider a set
    max_diversity = a threshold defining the average pairwise error rate required to retain a read group;
                    if None do not calculate diversity
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                    
    Returns: 
    a list of tuples of (position error, quality error, nucleotide pairwise error) output file names
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'EstimateErrorTotal'
    log['FILE'] = os.path.basename(seq_file)
    log['SET_FIELD'] = set_field
    log['NPROC'] = nproc
    printLog(log)
    
    # Check input file type
    in_type = getFileType(seq_file)
    #if in_type != 'fastq':  sys.exit('ERROR:  Input file must be FASTQ')
    
    # Define feeder function and arguments
    index_args = {'field': set_field, 'delimiter': out_args['delimiter']}
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets, 
                 'index_args': index_args}
    # Define worker function and arguments
    work_func = processEETQueue
    work_args = {}
    # Define collector function and arguments
    collect_func = collectEETQueue
    collect_args = {'seq_file': seq_file,
                    'out_args': out_args,
                    'set_field': set_field}

    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    result['log']['END'] = 'EstimateErrorTotal'
    printLog(result['log'])
        
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
                 error-position
                     estimated error by read position.
                 error-quality
                     estimated error by the quality score assigned within the input file.
                 error-nucleotide
                     estimated error by nucleotide.
                 error-set
                     estimated error by barcode read group size.

             output fields:
                 POSITION
                     read position with base zero indexing.
                 Q
                     Phred quality score.
                 OBSERVED
                     observed nucleotide value.
                 REFERENCE
                     consensus nucleotide for the barcode read group.
                 SET_COUNT
                     barcode read group size.
                 REPORTED_Q
                     mean Phred quality score reported within the input file for the given
                     position, quality score, nucleotide or read group.
                 MISMATCHES
                     count of observed mismatches from consensus for the given position,
                     quality score, nucleotide or read group.
                 OBSERVATIONS
                     total count of observed values for each position, quality score,
                     nucleotide or read group size.
                 ERROR
                     estimated error rate.
                 EMPIRICAL_Q
                     estimated error rate converted to a Phred quality score.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(seq_out=False,
                                                        failed=False,
                                                        multiproc=True)],
                            formatter_class=CommonHelpFormatter, add_help=False)

    # Error profiling arguments
    group_error = parser.add_argument_group('error profiling arguments')
    group_error.add_argument('-f', action='store', dest='set_field', type=str, default=default_barcode_field,
                             help='The name of the annotation field to group sequences by')

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
    if args_dict['set_field']:  args_dict['set_field'] = args_dict['set_field'].upper()
    
    # # Define cons_func and cons_args
    # if args_dict['mode'] == 'freq':
    #     args_dict['cons_func'] = frequencyConsensus
    #     args_dict['cons_args'] = {'min_freq':args_dict['min_freq']}
    # elif args_dict['mode'] == 'qual':
    #     args_dict['cons_func'] = qualityConsensus
    #     args_dict['cons_args'] = {'min_qual':args_dict['min_qual'],
    #                               'min_freq':args_dict['min_freq'],
    #                               'dependent':False}
    # del args_dict['mode']
    # if 'min_freq' in args_dict:  del args_dict['min_freq']
    # if 'min_qual' in args_dict:  del args_dict['min_qual']
    
    # Call estimateError for each sample file    
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        estimateErrorTotal(**args_dict)