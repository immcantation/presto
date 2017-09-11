#!/usr/bin/env python3
"""
Collects annotation fields based on grouping scheme
"""
# Info
__author__ = 'Ruoyi Jiang'
from presto import __version__, __date__

# Imports
import os
import sys
from argparse import ArgumentParser
from collections import OrderedDict
from textwrap import dedent

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field, \
                            default_primer_field, default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation, annotationConsensus
from presto.IO import printLog, readSeqFile
from presto.Sequence import indexSeqSets
from presto.Multiprocessing import SeqData, SeqResult, manageProcesses, feedSeqQueue, \
                                   collectSeqQueue

from Bio.SeqRecord import SeqRecord

default_collect_field='SAMPLE'
default_bool_del_func=False

def processCollectMajorityQueue(alive, data_queue, result_queue, collect_field,
                  collect_args={}, delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    cluster_args = a dictionary of optional arguments for the clustering function
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

            # Define result object
            result = SeqResult(data.id, data.data)
            result.log['CLUSTER'] = clust = data.id
            result.log['SEQCOUNT'] = len(data)
            
            # Get number of seqs in the cluster
            result.log['SEQS_IN_CLUST'] = len(data.data)

            cons_dict = annotationConsensus(data.data, collect_field)
            result.log['COLLECT_FIELD'] = cons_dict['cons']

            # Update sequence annotations with cluster assignments
            results = list()
            for i,seq in enumerate(data.data):
                header = parseAnnotation(seq.description, delimiter=delimiter)
                header[collect_field] = cons_dict['cons']
                seq.id = seq.name = flattenAnnotation(header, delimiter=delimiter)
                seq.description = ''

                result.log['SEQ: %i' % (i)] = str(seq.seq)
                results.append(seq)
                
            # Check results
            result.results = results
            result.valid = True
            # Feed results to result queue
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

def processCollectDeleteQueue(alive, data_queue, result_queue, collect_field,
                  collect_args={}, delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    cluster_args = a dictionary of optional arguments for the clustering function
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

            # Define result object
            result = SeqResult(data.id, data.data)
            result.log['CLUSTER'] = clust = data.id
            result.log['SEQCOUNT'] = len(data)
            
            # Get number of seqs in the cluster
            result.log['SEQS_IN_CLUST'] = len(data.data)

            #[parseAnnotation(seq.description, delimiter=delimiter)[collect_field] for seq in data.data]

            #if the number of unique identities in the annotation field is not 1, then the group is invalid and should be removed
            if len(set(parseAnnotation(seq.description, delimiter=delimiter)[collect_field] for seq in data.data)) == 1:
                result.valid = True
            else:
                result.valid = False

            results = data.data
            result.results = results

            # Feed results to result queue
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

def collectFields(seq_file, barcode_field=default_barcode_field,
                collect_field=default_collect_field, bool_del_func=default_bool_del_func,
                out_args=default_out_args, nproc=None,
                queue_size=None):
    """
    Performs clustering on sets of sequences

    Arguments:
    seq_file = the sample sequence file name
    barcode_field = the annotation containing set IDs
    ident = the identity threshold for clustering sequences
    seq_start = the start position to trim sequences at before clustering
    seq_end = the end position to trim sequences at before clustering
    cluster_exec = the path to the executable for usearch
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc

    Returns:
    the clustered output file name
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'CollectFields'
    log['FILE'] = os.path.basename(seq_file)
    log['BARCODE_FIELD'] = barcode_field
    log['COLLECT_FIELD'] = collect_field
    log['NPROC'] = nproc
    printLog(log)

    # Define feeder function and arguments
    index_args = {'field': barcode_field, 'delimiter': out_args['delimiter']}
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets,
                 'index_args': index_args}
    # Define worker function and arguments
    work_func = processCollectDeleteQueue if bool_del_func else processCollectMajorityQueue
    work_args = {'collect_field': collect_field,
                'delimiter': out_args['delimiter']}
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'collect',
                    'out_args': out_args,
                    'index_field': barcode_field}

    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func,
                             feed_args, work_args, collect_args,
                             nproc, queue_size)

    # Print log
    log = OrderedDict()
    log['OUTPUT'] = result['log'].pop('OUTPUT')
    for k, v in result['log'].items():  log[k] = v
    log['END'] = 'CollectFields'
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
                 collect-pass
                    clustered reads with consensus annotation.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(multiproc=True)],
                            formatter_class=CommonHelpFormatter, add_help=False)

    # Clustering arguments
    group_collect = parser.add_argument_group('collect field arguments')
    group_collect.add_argument('-f', action='store', dest='barcode_field', type=str,
                             default=default_barcode_field,
                             help='''The annotation field containing annotations, such as UID
                                  barcode, for sequence grouping.''')
    group_collect.add_argument('-k', action='store', dest='collect_field', type=str,
                             default=default_collect_field,
                             help='''The name of the annotation field to find a consensus for
                                per each barcode group.''')
    group_collect.add_argument('--del', action='store_true', dest='bool_del_func', 
                         default=default_bool_del_func,
                         help='''if set, delete collect function will be used.''')

    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """

#     #Read in each uid group
# #Find the majority of a particular group and assign everything to that group
# CollapseFields.py -s SEQFILE -k GROUP -f FIELD --maj 




    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)

    # Convert fields to uppercase
    if 'barcode_field' in args_dict and args_dict['barcode_field'] is not None:
        args_dict['barcode_field'] = args_dict['barcode_field'].upper()
    if 'collect_field' in args_dict and args_dict['collect_field'] is not None:
        args_dict['collect_field'] = args_dict['collect_field'].upper()

    # Call cluster for each input file
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        collectFields(**args_dict)

