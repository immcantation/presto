#!/usr/bin/env python
"""
Builds a consensus sequence for each set of input sequences
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.2'
__date__      = '2014.3.5'

# Imports
import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict

# IgPipeline imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_barcode_field, default_delimiter, default_out_args
from IgCore import flattenAnnotation, mergeAnnotation
from IgCore import getCommonArgParser, parseCommonArgs, printLog, getFileType
from IgCore import annotationConsensus, frequencyConsensus, qualityConsensus
from IgCore import calculateDiversity, indexSeqSets, subsetSeqSet
from IgCore import collectSeqQueue, feedSeqQueue
from IgCore import manageProcesses, SeqResult

# Defaults
default_min_count = 1
default_min_qual = 0
default_min_freq = 0.3


def processBCQueue(alive, data_queue, result_queue, cons_func, cons_args={}, 
                   min_count=default_min_count, primer_field=None, primer_freq=None, 
                   max_diversity=None, delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    cons_func = the function to use for consensus generation 
    cons_args = a dictionary of optional arguments for the consensus function
    min_count = threshold number of sequences to define a consensus 
    primer_field = the annotation field containing primer names;
                   if None do not annotate with primer names
    primer_freq = the maximum primer frequency that must be meet to build a consensus;
                  if None do not filter by primer frequency
    max_diversity = the minimum diversity score to retain a set;
                    if None do not calculate diversity
    delimiter = a tuple of delimiters for (annotations, field/values, value lists) 

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
            
            # Define result dictionary for iteration
            result = SeqResult(data.id, data.data)
            result.log['BARCODE'] = data.id
            result.log['SEQCOUNT'] = len(data)
    
            # Define primer annotations and consensus primer if applicable
            if primer_field is None:
                primer_ann = None
                seq_list = data.data
            else:
                # Calculate consensus primer            
                primer_ann = OrderedDict()
                prcons = annotationConsensus(data.data, primer_field, delimiter=delimiter)
                result.log['PRIMER'] = ','.join(prcons['set'])
                result.log['PRCOUNT'] = ','.join([str(c) for c in prcons['count']])
                result.log['PRCONS'] = prcons['cons']
                result.log['PRFREQ'] = prcons['freq']
                if primer_freq is None:
                    # Retain full sequence set if not in primer consensus mode
                    seq_list = data.data
                    primer_ann = mergeAnnotation(primer_ann, {'PRIMER':prcons['set']}, 
                                                 delimiter=delimiter)
                    primer_ann = mergeAnnotation(primer_ann, {'PRCOUNT':prcons['count']}, 
                                                 delimiter=delimiter)
                elif prcons['freq'] >= primer_freq:
                    # Define consensus subset
                    seq_list = subsetSeqSet(data.data, primer_field, prcons['cons'], 
                                            delimiter=delimiter)
                    primer_ann = mergeAnnotation(primer_ann, {'PRCONS':prcons['cons']}, 
                                                 delimiter=delimiter)
                    primer_ann = mergeAnnotation(primer_ann, {'PRFREQ':prcons['freq']}, 
                                                 delimiter=delimiter)
                else:
                    # If set fails primer consensus, feed result queue and continue
                    result_queue.put(result)
                    continue
    
            # Update log
            cons_count = len(seq_list)
            result.log['CONSCOUNT'] = cons_count
            if cons_count < min_count:
                # If set fails count threshold, feed result queue and continue
                result_queue.put(result)
                continue
                
            # Calculate average pairwise error rate
            if max_diversity is not None:
                diversity = calculateDiversity(seq_list)
                result.log['DIVERSITY'] = diversity
                if diversity > max_diversity:
                    # If diversity exceeds threshold, feed result queue and continue
                    for i, s in enumerate(seq_list):
                        result.log['INSEQ%i' % (i + 1)] = str(s.seq)
                    result_queue.put(result)
                    continue
    
            # If primer and diversity filters pass, generate consensus sequence
            consensus = cons_func(seq_list, **cons_args)
    
            # Update log
            for i, s in enumerate(seq_list):
                result.log['INSEQ%i' % (i + 1)] = str(s.seq)
            result.log['CONSENSUS'] = str(consensus.seq)
            if 'phred_quality' in consensus.letter_annotations:
                result.log['QUALITY'] = ''.join([chr(c+33) for c in consensus.letter_annotations['phred_quality']])
            
            # Define annotation for consensus sequence
            cons_ann = OrderedDict([('ID', data.id),
                                    ('CONSCOUNT', cons_count)])
            if primer_ann is not None:
                cons_ann = mergeAnnotation(cons_ann, primer_ann, delimiter=delimiter)
            consensus.id = consensus.name = flattenAnnotation(cons_ann, delimiter=delimiter)
            consensus.description = ''
            result.results = consensus
            result.valid = True
            
            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        raise
    
    return None


def buildConsensus(seq_file, barcode_field=default_barcode_field, 
                   min_count=default_min_count, min_freq=default_min_freq,
                   min_qual=default_min_qual, max_miss=None, primer_field=None, 
                   primer_freq=None, max_diversity=None, dependent=False, 
                   out_args=default_out_args, nproc=None, queue_size=None):
    """
    Generates consensus sequences

    Arguments: 
    seq_file = the sample sequence file name
    barcode_field = the annotation field containing set IDs
    min_count = threshold number of sequences to define a consensus 
    min_freq = the frequency cutoff to assign a base 
    min_qual = the quality cutoff to assign a base
    max_miss = the maximum frequency of (., -, N) characters allowed before 
               deleting a position; if None do not delete positions 
    primer_field = the annotation field containing primer tags;
                   if None do not annotate with primer tags
    primer_freq = the maximum primer tag frequency that must be meet to build a consensus;
                  if None do not filter by primer frequency
    max_diversity = a threshold defining the average pairwise error rate required to retain a read group;
                    if None do not calculate diversity
    dependent = if False treat barcode group sequences as independent data
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                    
    Returns: 
    a list of successful output file names
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'BuildConsensus'
    log['FILE'] = os.path.basename(seq_file)
    log['BARCODE_FIELD'] = barcode_field
    log['MIN_COUNT'] = min_count
    log['MIN_FREQUENCY'] = min_freq
    log['MIN_QUALITY'] = min_qual
    log['MAX_MISSING'] = max_miss
    log['PRIMER_FIELD'] = primer_field
    log['PRIMER_FREQUENCY'] = primer_freq
    log['MAX_DIVERSITY'] = max_diversity
    log['DEPENDENT'] = dependent
    log['NPROC'] = nproc
    printLog(log)
    
    # Set consensus building function
    in_type = getFileType(seq_file)
    if in_type == 'fastq':
        cons_func = qualityConsensus
        cons_args = {'min_qual': min_qual, 
                     'min_freq': min_freq,
                     'max_miss': max_miss,
                     'dependent': dependent}
    elif in_type == 'fasta':  
        cons_func = frequencyConsensus
        cons_args = {'min_freq': min_freq,
                     'max_miss': max_miss}
    else:
        sys.exit('ERROR:  Input file must be FASTA or FASTQ')
    
    # Define feeder function and arguments
    index_args = {'field': barcode_field, 'delimiter': out_args['delimiter']}
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets, 
                 'index_args': index_args}
    # Define worker function and arguments
    work_func = processBCQueue
    work_args = {'cons_func': cons_func, 
                 'cons_args': cons_args,
                 'min_count': min_count,
                 'primer_field': primer_field,
                 'primer_freq': primer_freq,
                 'max_diversity': max_diversity,
                 'delimiter': out_args['delimiter']}
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'consensus',
                    'out_args': out_args,
                    'index_field': barcode_field}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    result['log']['END'] = 'BuildConsensus'
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
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, version='%(prog)s:' + ' v%s-%s' %(__version__, __date__), 
                            parents=[getCommonArgParser(multiproc=True)], 
                            formatter_class=ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-n', action='store', dest='min_count', type=int, default=default_min_count,
                        help='The minimum number of sequences needed to define a valid consensus')
    parser.add_argument('--bf', action='store', dest='barcode_field', type=str,
                        default=default_barcode_field, 
                        help='Position of description barcode field to group sequences by')
    parser.add_argument('-q', action='store', dest='min_qual', type=float, default=default_min_qual,
                        help='Consensus quality score cut-off under which an ambiguous character is assigned; \
                              does not apply when quality scores are unavailable')
    parser.add_argument('--freq', action='store', dest='min_freq', type=float, default=default_min_freq,
                        help='Fraction of character occurrences under which an ambiguous character is assigned.')
    parser.add_argument('--maxmiss', action='store', dest='max_miss', type=float, default=None,
                        help='If specified, this defines a cut-off for the frequency of allowed \
                              missing values for each position. Positions exceeding the threshold \
                              are deleted from the consensus. If not defined, positions are always \
                              retained.')
    parser.add_argument('--pf', action='store', dest='primer_field', type=str, default=None, 
                        help='Specifies the field name of the primer annotations')
    parser.add_argument('--prcons', action='store', dest='primer_freq', type=float, default=None, 
                        help='Specify to define a minimum primer frequency required to assign a consensus primer, \
                              and filter out sequences with minority primers from the consensus building step')
    parser.add_argument('--maxdiv', action='store', dest='max_diversity', type=float, default=None,
                        help='Specify to calculate the nucleotide diversity of each read group \
                              (average pairwise error rate) and remove groups exceeding the given diversity threshold')
    parser.add_argument('--dep', action='store_true', dest='dependent',
                        help='Specify to calculate consensus quality with a non-independence assumption')
    
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
    if args_dict['barcode_field']:  args_dict['barcode_field'] = args_dict['barcode_field'].upper()
    if args_dict['primer_field']:  args_dict['primer_field'] = args_dict['primer_field'].upper()
    
    # Check prcons argument dependencies
    if args.primer_freq and not args.primer_field:
        parser.error('You must define a primer field with --prf to use the --prcons option')
        
    # Call buildConsensus for each sample file    
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        buildConsensus(**args_dict)
