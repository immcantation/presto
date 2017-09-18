#!/usr/bin/env python3
"""
Collects annotation fields based on grouping scheme
"""
# Info
__author__ = 'Ruoyi Jiang, Jason Vander Heiden'
from presto import __version__, __date__

# Imports
import os
import sys
from argparse import ArgumentParser
from collections import OrderedDict
from textwrap import dedent

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field, default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, annotationConsensus
from presto.IO import printLog
from presto.Sequence import indexSeqSets
from presto.Multiprocessing import SeqResult, manageProcesses, collectSeqQueue, feedSeqQueue, \
                                   processSeqQueue

# Defaults
default_filter_field = 'SAMPLE'


def consensusField(data, field, delimiter=default_delimiter):
    """
    Reassigns all annotations to the consensus annotation in group

    Returns:

    """
    # Define result object
    result = SeqResult(data.id, data.data)
    result.log['CLUSTER'] = data.id
    result.log['SEQCOUNT'] = len(data)

    # Get number of seqs in the cluster
    result.log['SEQS_IN_CLUST'] = len(data.data)

    cons_dict = annotationConsensus(data.data, field)
    result.log['FILTER_FIELD'] = cons_dict['cons']

    # Update sequence annotations with cluster assignments
    results = list()
    for i, seq in enumerate(data.data):
        header = parseAnnotation(seq.description, delimiter=delimiter)
        header[field] = cons_dict['cons']
        seq.id = seq.name = flattenAnnotation(header, delimiter=delimiter)
        seq.description = ''

        result.log['SEQ: %i' % (i)] = str(seq.seq)
        results.append(seq)

    # Check results
    result.results = results
    result.valid = True

    return result


def deleteField(data, field, delimiter=default_delimiter):
    """
    Removes all sequences with differing fields in a group

    Returns:

    """
    # Define result object
    result = SeqResult(data.id, data.data)
    result.log['CLUSTER'] = data.id
    result.log['SEQCOUNT'] = len(data)

    # Get number of seqs in the cluster
    result.log['SEQS_IN_CLUST'] = len(data.data)

    # [parseAnnotation(seq.description, delimiter=delimiter)[filter_field] for seq in data.data]

    # if the number of unique identities in the annotation field is not 1, then the group is invalid and should be removed
    if len(set(parseAnnotation(seq.description, delimiter=delimiter)[field] for seq in data.data)) == 1:
        result.valid = True
    else:
        result.valid = False

    results = data.data
    result.results = results

    return None


def minorityField(data, field, delimiter=default_delimiter):
    """
    Removes all sequences in the minority annotation

    Returns:

    """
    pass


def processCollectMajorityQueue(alive, data_queue, result_queue, filter_field,
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

            cons_dict = annotationConsensus(data.data, filter_field)
            result.log['FILTER_FIELD'] = cons_dict['cons']

            # Update sequence annotations with cluster assignments
            results = list()
            for i,seq in enumerate(data.data):
                header = parseAnnotation(seq.description, delimiter=delimiter)
                header[filter_field] = cons_dict['cons']
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


def processCollectDeleteQueue(alive, data_queue, result_queue, filter_field,
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

            #[parseAnnotation(seq.description, delimiter=delimiter)[filter_field] for seq in data.data]

            #if the number of unique identities in the annotation field is not 1, then the group is invalid and should be removed
            if len(set(parseAnnotation(seq.description, delimiter=delimiter)[filter_field] for seq in data.data)) == 1:
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


def filterFields(seq_file, filter_func, set_field=default_barcode_field,
                  filter_field=default_filter_field,
                  out_args=default_out_args, nproc=None,
                  queue_size=None):
    """
    Filters and modifies fields within annotation groups

    Arguments:
      seq_file : the sample sequence file name
      filter_func : the function to use for filtering fields
      set_field : the annotation containing set IDs
      filter_field : the field for collection criteria
      out_args : common output argument dictionary from parseCommonArgs
      nproc : the number of processQueue processes;
              if None defaults to the number of CPUs
      queue_size : maximum size of the argument queue;
                   if None defaults to 2*nproc

    Returns:
      str : output file name
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'CollectFields'
    log['FILE'] = os.path.basename(seq_file)
    log['SET_FIELD'] = set_field
    log['FILTER_FIELD'] = filter_field
    log['NPROC'] = nproc
    printLog(log)

    # Define feeder function and arguments
    index_args = {'field': set_field,
                  'delimiter': out_args['delimiter']}
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets,
                 'index_args': index_args}
    # Define worker function and arguments
    filter_args = {'field': filter_field,
                   'delimiter': out_args['delimiter']}
    work_func = processSeqQueue
    work_args = {'process_func': filter_func,
                 'process_args': filter_args}
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'collect',
                    'out_args': out_args,
                    'index_field': set_field}

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
                            formatter_class=CommonHelpFormatter, add_help=False)
    group_help = parser.add_argument_group('help')
    group_help.add_argument('--version', action='version',
                            version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    group_help.add_argument('-h', '--help', action='help', help='show this help message and exit')
    subparsers = parser.add_subparsers(title='subcommands', metavar='',
                                       help='Filtering operation')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Parent parser
    parser_parent = getCommonArgParser(annotation=False, log=True, multiproc=True)
    group_parent = parser_parent.add_argument_group('annotation arguments')
    group_parent.add_argument('-f', action='store', dest='set_field', type=str,
                             default=default_barcode_field,
                             help='''The annotation field containing annotations, such as the UMI
                                  barcode, for sequence grouping.''')
    group_parent.add_argument('-k', action='store', dest='filter_field', type=str,
                             default=default_filter_field,
                             help='''The name of the annotation field to find a consensus for
                                  per each sequence group.''')
    # Consensus arguments
    parser_cons = subparsers.add_parser('consensus', parents=[parser_parent],
                                          formatter_class=CommonHelpFormatter, add_help=False,
                                          help='Reassign fields to consensus values.',
                                          description='Reassign fields to consensus values.')
    parser_cons.set_defaults(filter_func=consensusField)

    # Deletion arguments
    parser_del = subparsers.add_parser('delete', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter, add_help=False,
                                        help='Delete fields with differen values.',
                                        description='Delete fields with differen values.')
    parser_del.set_defaults(filter_func=deleteField)

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
    if 'set_field' in args_dict and args_dict['set_field'] is not None:
        args_dict['set_field'] = args_dict['set_field'].upper()
    if 'filter_field' in args_dict and args_dict['filter_field'] is not None:
        args_dict['filter_field'] = args_dict['filter_field'].upper()

    # Call cluster for each input file
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        filterFields(**args_dict)

