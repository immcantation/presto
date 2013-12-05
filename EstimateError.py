#!/usr/bin/env python
"""
Calculates annotation set error rates
"""

__author__    = 'Jason Anthony Vander Heiden, Namita Gupta'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.1'
__date__      = '2013.10.12'

# Imports
import os, sys
import multiprocessing as mp
import numpy as np
import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from itertools import permutations
from time import time

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_missing_chars, default_barcode_field, default_out_args
from IgCore import default_min_freq, default_min_qual
from IgCore import getCommonArgParser, parseCommonArgs
from IgCore import getOutputHandle, printLog, printProgress, getFileType
from IgCore import feedSetQueue, getScoreDict, calculateDiversity, countSeqSets
from IgCore import frequencyConsensus, qualityConsensus

# Defaults
default_min_count = 10
        

def countMismatches(seq_list, ref_seq, ignore_chars=default_missing_chars, 
                    score_dict=getScoreDict(n_score=1, gap_score=1)):
    """
    Counts the occurrence of nucleotide mismatches in a set of sequences

    Arguments: 
    seq_list = a list of SeqRecord objects with aligned sequences
    ref_seq = a SeqRecord object containing the reference sequence to match against
    ignore_chars = list of characters to exclude from mismatch counts
    score_dict = optional dictionary of alignment scores as {(char1, char2): score}

    Returns: 
    a dictionary of pandas.DataFrame objects containing [mismatch, qsum, total] counts  
    for {pos:sequence position, nuc:nucleotide pairs, qual:quality score, set:sequence set} 
    """
    # Define position mismatch DataFrame
    pos_max = max([len(s) for s in seq_list])
    pos_df = pd.DataFrame(0, index=range(pos_max), 
                          columns=['mismatch', 'q_sum', 'total'], dtype=float)
    # Define nucleotide mismatch DataFrame
    nuc_pairs = list(permutations(['A', 'C', 'G', 'T'], 2))
    nuc_df = pd.DataFrame(0, index=pd.MultiIndex.from_tuples(nuc_pairs, names=['obs', 'ref']), 
                          columns=['mismatch', 'q_sum', 'total'], dtype=float)
    # Define quality mismatch DataFrame
    qual_df = pd.DataFrame(0, index=range(94), 
                           columns=['mismatch', 'q_sum', 'total'], dtype=float)    

    # Iterate over seq_list and count mismatches
    for seq in seq_list:
        qual = seq.letter_annotations['phred_quality']
        for i, b in enumerate(seq):
            a = ref_seq[i]
            q = qual[i]
            # Update total counts and qualities
            if {a, b}.isdisjoint(ignore_chars):  
                pos_df.ix[i, 'total'] += 1
                pos_df.ix[i, 'q_sum'] += q
                nuc_df.ix[b]['total'] += 1
                nuc_df.ix[b]['q_sum'] += q
                qual_df.ix[q, 'total'] += 1
                qual_df.ix[q, 'q_sum'] += q
            # Update mismatch counts
            if score_dict[(a, b)] == 0:
                pos_df.ix[i, 'mismatch'] += 1
                nuc_df.ix[(a, b), 'mismatch'] += 1
                qual_df.ix[q, 'mismatch'] += 1

    # Define set total mismatch DataFrame
    set_df = pd.DataFrame([pos_df.sum(axis=0)], index=[len(seq_list)], 
                          columns=['mismatch', 'q_sum', 'total'], dtype=float)

    return {'pos':pos_df, 'nuc':nuc_df, 'qual':qual_df, 'set':set_df}

    
def processQueue(data_queue, result_queue, cons_func, cons_args={}, 
                 min_count=default_min_count, max_diversity=None):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
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
    # Iterator over data queue until sentinel object reached
    for args in iter(data_queue.get, None):        
        seq_list = args['seq_list']
        seq_count = len(seq_list)
        # Define result dictionary for iteration
        results = {'id':args['id'],
                   'seq_count':seq_count,
                   'diversity':None,
                   'pos':None,
                   'nuc':None,
                   'qual':None,
                   'set':None,
                   'pass':False,
                   'log':OrderedDict()}
        # Update log
        results['log']['SET'] = args['id']
        results['log']['SEQCOUNT'] = seq_count
        for i, s in enumerate(seq_list):
            results['log']['SEQ%i' % (i + 1)] = str(s.seq)
        
        # Check count threshold and continue if failed
        if seq_count < min_count:
            result_queue.put(results)
            continue
            
        # Calculate average pairwise error rate
        if max_diversity is not None:
            diversity = calculateDiversity(args['seq_list'])
            results['diversity'] = diversity
            results['log']['DIVERSITY'] = diversity
            # Check diversity threshold and continue if failed
            if diversity > max_diversity:
                result_queue.put(results)
                continue
            
        # Define reference sequence by consensus
        ref_seq = cons_func(seq_list, **cons_args)
        
        # Count mismatches against consensus
        mismatch = countMismatches(seq_list, ref_seq)
        
        # Calculate average reported and observed error
        reported_q = mismatch['set']['q_sum'].sum() / mismatch['set']['total'].sum()
        error_rate = mismatch['set']['mismatch'].sum() / mismatch['set']['total'].sum()

        # Update log
        results['log']['REFERENCE'] = str(ref_seq.seq)
        results['log']['MISMATCH'] = ''.join(['*' if x > 0 else ' ' for x in mismatch['pos']['mismatch']])
        results['log']['ERROR'] = '%.6f' % error_rate
        results['log']['REPORTED_Q'] = '%.2f' % reported_q
        results['log']['EMPIRICAL_Q'] = '%.2f' % (-10 * np.log10(error_rate))
            
        # Update results and feed result queue
        results['pass'] = True
        results.update(mismatch)
        result_queue.put(results)


    return None


def collectQueue(result_queue, collect_dict, seq_file, field, out_args):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments: 
    result_queue = a multiprocessing.Queue holding processQueue results
    collect_dict = a multiprocessing.Manager.dict to store return values
    seq_file = the sample sequence file name
    field = the field defining set membership
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns:
    None
    (adds 'log' and 'out_files' to collect_dict)
    """
    # Count sets
    result_count = countSeqSets(seq_file, field, out_args['delimiter'])
    
    # Define empty DataFrames to store assembled results
    pos_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'], dtype=float)
    qual_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'], dtype=float)
    nuc_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'], dtype=float)
    set_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'], dtype=float)
    
    # Open log file
    if out_args['log_file'] is None:
        log_handle = None
    else:
        log_handle = open(out_args['log_file'], 'w')

    # Iterator over results queue until sentinel object reached
    start_time = time()
    set_count = seq_count = pass_count = fail_count = 0
    for result in iter(result_queue.get, None):
        # Print progress for previous iteration
        printProgress(set_count, result_count, 0.05, start_time)
        
        # Update counts for iteration
        set_count += 1
        seq_count += result['seq_count']
        
        # Sum results
        if result['pass']:
            pass_count += 1
            pos_df = pos_df.add(result['pos'], fill_value=0)
            qual_df = qual_df.add(result['qual'], fill_value=0)
            nuc_df = nuc_df.add(result['nuc'], fill_value=0)
            set_df = set_df.add(result['set'], fill_value=0)
        else:
            fail_count += 1
            
        # Write log
        printLog(result['log'], handle=log_handle)
    
    # Print final progress
    printProgress(set_count, result_count, 0.05, start_time)
    
    # Generate log
    log = OrderedDict()
    for i in xrange(4): 
        log['OUTPUT%i' % (i + 1)] = None
    log['SETS'] = set_count
    log['SEQUENCES'] = seq_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['POSITION_ERROR'] = None 
    log['NUCLEOTIDE_ERROR'] = None
    log['QUALITY_ERROR'] = None
    log['SET_ERROR'] = None
 
    # Return if no mismatch data
    if pass_count == 0:
        collect_dict['out_files'] = None
        collect_dict['log'] = log
        return None

    # Calculate error rates
    pos_df['error'] = pos_df['mismatch'] / pos_df['total'] 
    nuc_df['error'] = nuc_df['mismatch'] / nuc_df['total']
    qual_df['error'] = qual_df['mismatch'] / qual_df['total']
    set_df['error'] = set_df['mismatch'] / set_df['total']
    
    # Convert error to empirical quality score
    pos_df['emp_q'] = -10 * np.log10(pos_df['error'])
    nuc_df['emp_q'] = -10 * np.log10(nuc_df['error'])
    qual_df['emp_q'] = -10 * np.log10(qual_df['error'])
    set_df['emp_q'] = -10 * np.log10(set_df['error'])

    # Calculate reported quality means
    pos_df['rep_q'] = pos_df['q_sum'] / pos_df['total'] 
    nuc_df['rep_q'] = nuc_df['q_sum'] / nuc_df['total']
    qual_df['rep_q'] = qual_df['q_sum'] / qual_df['total']
    set_df['rep_q'] = set_df['q_sum'] / set_df['total']
        
    # Calculate overall error rate
    pos_error = pos_df['mismatch'].sum() / pos_df['total'].sum() 
    qual_error = qual_df['mismatch'].sum() / qual_df['total'].sum() 
    nuc_error = nuc_df['mismatch'].sum() / nuc_df['total'].groupby(level='obs').mean().sum()
    set_error = set_df['mismatch'].sum() / set_df['total'].sum() 

    # Build results dictionary
    assembled = {'pos':pos_df, 'qual':qual_df, 'nuc':nuc_df, 'set':set_df}
    
    # Write assembled error counts to output files
    out_files = writeResults(assembled, seq_file, out_args)
    
    # Update log
    for i, f in enumerate(out_files): 
        log['OUTPUT%i' % (i + 1)] = os.path.basename(f)
    log['POSITION_ERROR'] = pos_error 
    log['NUCLEOTIDE_ERROR'] = nuc_error 
    log['QUALITY_ERROR'] = qual_error
    log['SET_ERROR'] = set_error
    
    # Update collector results
    collect_dict['out_files'] = out_files
    collect_dict['log'] = log

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
    pos_df = results['pos']
    nuc_df = results['nuc']
    qual_df = results['qual']
    set_df = results['set']

    # Type conversion to int of mismatch and total columns
    pos_df[['mismatch', 'total']] = pos_df[['mismatch', 'total']].astype(int) 
    nuc_df[['mismatch', 'total']] = nuc_df[['mismatch', 'total']].astype(int) 
    qual_df[['mismatch', 'total']] = qual_df[['mismatch', 'total']].astype(int) 
    set_df[['mismatch', 'total']] = set_df[['mismatch', 'total']].astype(int) 
    
    # Write to tab delimited files
    file_args = {'out_dir':out_args['out_dir'], 'out_name':out_args['out_name'], 'out_type':'tab'}
    with getOutputHandle(seq_file, 'error-position', **file_args) as pos_handle, \
            getOutputHandle(seq_file, 'error-quality', **file_args) as qual_handle, \
            getOutputHandle(seq_file, 'error-nucleotide', **file_args) as nuc_handle, \
            getOutputHandle(seq_file, 'error-set', **file_args) as set_handle:
        pos_df.to_csv(pos_handle, sep='\t', na_rep='NA', index_label='POSITION', 
                      cols=['rep_q', 'mismatch', 'total', 'error', 'emp_q'],
                      header=['REPORTED_Q', 'MISMATCHES', 'OBSERVATIONS', 'ERROR', 'EMPIRICAL_Q'], float_format='%.6f')
        nuc_df.to_csv(nuc_handle, sep='\t', na_rep='NA', index_label=['OBSERVED', 'REFERENCE'], 
                      cols=['rep_q', 'mismatch', 'total', 'error', 'emp_q'],
                      header=['REPORTED_Q', 'MISMATCHES', 'OBSERVATIONS', 'ERROR', 'EMPIRICAL_Q'], float_format='%.6f')
        qual_df.to_csv(qual_handle, sep='\t', na_rep='NA', index_label='Q',
                       cols=['rep_q', 'mismatch', 'total', 'error', 'emp_q'], 
                       header=['REPORTED_Q', 'MISMATCHES', 'OBSERVATIONS', 'ERROR', 'EMPIRICAL_Q'], float_format='%.6f')
        set_df.to_csv(set_handle, sep='\t', na_rep='NA', index_label='SET_COUNT', 
                      cols=['rep_q', 'mismatch', 'total', 'error', 'emp_q'],
                      header=['REPORTED_Q', 'MISMATCHES', 'OBSERVATIONS', 'ERROR', 'EMPIRICAL_Q'], float_format='%.6f')

    return (pos_handle.name, qual_handle.name, nuc_handle.name, set_handle.name)


def estimateError(seq_file, cons_func=frequencyConsensus, cons_args={}, 
                  set_field=default_barcode_field, min_count=default_min_count, max_diversity=None, 
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
    # Define number of processes and queue size
    if nproc is None:  nproc = mp.cpu_count()
    if queue_size is None:  queue_size = nproc * 2
    
    # Define subcommand label dictionary
    cmd_dict = {frequencyConsensus:'freq', qualityConsensus:'qual'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'EstimateError'
    log['FILE'] = os.path.basename(seq_file)
    log['MODE'] = cmd_dict.get(cons_func, cons_func.__name__)
    log['SET_FIELD'] = set_field
    log['MIN_COUNT'] = min_count
    log['MAX_DIVERSITY'] = max_diversity
    log['NPROC'] = nproc
    printLog(log)
    
    # Check input file type
    in_type = getFileType(seq_file)
    if in_type != 'fastq':  sys.exit('ERROR:  Input file must be FASTQ')
    
    # Define shared data objects 
    manager = mp.Manager()
    collect_dict = manager.dict()
    data_queue = mp.Queue(queue_size)
    result_queue = mp.Queue(queue_size)
    
    # Initiate feeder process
    feeder = mp.Process(target=feedSetQueue, args=(data_queue, nproc, seq_file, set_field, 
                                                   out_args['delimiter']))
    feeder.start()
    
    # Initiate processQueue processes
    workers = []
    for __ in range(nproc):
        w = mp.Process(target=processQueue, args=(data_queue, result_queue, cons_func, cons_args, 
                                                  min_count, max_diversity))
        w.start()
        workers.append(w)

    # Initiate collector process
    collector = mp.Process(target=collectQueue, args=(result_queue, collect_dict, seq_file, 
                                                      set_field, out_args))
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
    log['END'] = 'EstimateError'
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
                            parents=[getCommonArgParser(seq_out=False, multiproc=True)], 
                            formatter_class=ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-f', action='store', dest='set_field', type=str, default=default_barcode_field, 
                        help='The name of the annotation field to group sequences by')
    parser.add_argument('-n', action='store', dest='min_count', type=int, default=default_min_count,
                        help='The minimum number of sequences needed to consider a set')
    parser.add_argument('--mode', action='store', dest='mode', choices=('freq', 'qual'), default='freq', 
                        help='Specifies which method to use to determin the consensus sequence; \
                              either by pure frequency (freq) or with consideration of quality score (qual)')
    parser.add_argument('-q', action='store', dest='min_qual', type=float, default=default_min_qual,
                        help='Consensus quality score cut-off under which an ambiguous character is assigned')
    parser.add_argument('--freq', action='store', dest='min_freq', type=float, default=default_min_freq,
                        help='Fraction of character occurrences under which an ambiguous character is assigned')
    parser.add_argument('--maxdiv', action='store', dest='max_diversity', type=float, default=None,
                        help='Specify to calculate the nucleotide diversity of each read group \
                              (average pairwise error rate) and remove groups exceeding the given diversity threshold')
    
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
    
    # Define cons_func and cons_args
    if args_dict['mode'] == 'freq':
        args_dict['cons_func'] = frequencyConsensus
        args_dict['cons_args'] = {'min_freq':args_dict['min_freq']}
    elif args_dict['mode'] == 'qual':
        args_dict['cons_func'] = qualityConsensus
        args_dict['cons_args'] = {'min_qual':args_dict['min_qual'], 'dependent':False}
    del args_dict['mode']
    if 'min_freq' in args_dict:  del args_dict['min_freq']
    if 'min_qual' in args_dict:  del args_dict['min_qual']
    
    # Call estimateError for each sample file    
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        estimateError(**args_dict)

        # Profiling
        #import cProfile, pstats
        #cProfile.run('estimateError(**args_dict)', 'profile.prof')
        #p = pstats.Stats('profile.prof')
        #p.strip_dirs().sort_stats('time').print_stats()