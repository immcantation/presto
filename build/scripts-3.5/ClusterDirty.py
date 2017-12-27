#!/bin/python3
"""
Cluster sequences by group
"""
# Info
__author__ = 'Christopher Bolen, Jason Anthony Vander Heiden, Ruoyi Jiang'
from presto import __version__, __date__

# Imports
import os
import shutil
import sys
from argparse import ArgumentParser
from collections import OrderedDict
from textwrap import dedent
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field, \
                            default_cluster_field, default_out_args, \
                            default_usearch_exec
from presto.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation
from presto.Applications import runUClust
from presto.IO import countSeqFile, getFileType, getOutputHandle, printLog, printMessage, \
                      printProgress, readSeqFile
from presto.Sequence import indexSeqSets
from presto.Multiprocessing import SeqResult, manageProcesses, feedSeqQueue, \
                                   collectSeqQueue

# Defaults
default_cluster_prefix=''
default_cluster_exec = default_usearch_exec
default_ident = 0.9


def processCSQueue(alive, data_queue, result_queue,
                   cluster_field=default_cluster_field,
                   cluster_prefix=default_cluster_prefix,
                   cluster_args={}, delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
      alive : a multiprocessing.Value boolean controlling whether processing
              continues; when False function returns
      data_queue : a multiprocessing.Queue holding data to process
      result_queue : a multiprocessing.Queue to hold processed results
      cluster_field : string defining the output cluster field name.
      cluster_prefix : string defining a prefix for the cluster identifier.
      cluster_args : a dictionary of optional arguments for the clustering function
      delimiter : a tuple of delimiters for (annotations, field/values, value lists)

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
            result.log['BARCODE'] = data.id
            result.log['SEQCOUNT'] = len(data)

            # Perform clustering
            cluster_dict = runUClust(data.data, **cluster_args)

            # Process failed result
            if cluster_dict is None:
                # Update log
                result.log['CLUSTERS'] = 0
                for i, seq in enumerate(data.data, start=1):
                    result.log['CLUST0-%i' % i] = str(seq.seq)

                # Feed results queue and continue
                result_queue.put(result)
                continue

            # Get number of clusters
            result.log['CLUSTERS'] = len(cluster_dict)

            # Update sequence annotations with cluster assignments
            results = list()
            seq_dict = {s.id: s for s in data.data}
            for cluster, id_list in cluster_dict.items():
                for i, seq_id in enumerate(id_list, start=1):
                    # Add cluster annotation
                    seq = seq_dict[seq_id]
                    label = '%s%i' % (cluster_prefix, cluster)
                    header = parseAnnotation(seq.description, delimiter=delimiter)
                    header = mergeAnnotation(header, {cluster_field: label}, delimiter=delimiter)
                    seq.id = seq.name = flattenAnnotation(header, delimiter=delimiter)
                    seq.description = ''

                    # Update log and results
                    result.log['CLUST%i-%i' % (cluster, i)] = str(seq.seq)
                    results.append(seq)

            # Check results
            result.results = results
            result.valid = (len(results) == len(seq_dict))
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


def clusterSets(seq_file, set_field=default_barcode_field,
                cluster_field=default_cluster_field,
                cluster_prefix=default_cluster_prefix,
                ident=default_ident, seq_start=0, seq_end=None,
                cluster_exec=default_cluster_exec,
                out_args=default_out_args, nproc=None,
                queue_size=None):
    """
    Performs clustering on sets of sequences

    Arguments:
      seq_file : the sample sequence file name
      set_field : the annotation containing set IDs
      cluster_field : the name of the output cluster field
      cluster_prefix : string defining a prefix for the cluster identifier.
      ident : the identity threshold for clustering sequences
      seq_start : the start position to trim sequences at before clustering
      seq_end : the end position to trim sequences at before clustering
      cluster_exec : the path to the executable for usearch
      nproc : the number of processQueue processes;
              if None defaults to the number of CPUs
      queue_size : maximum size of the argument queue;
                   if None defaults to 2*nproc

    Returns:
      str : the clustered output file name
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'ClusterSets'
    log['FILE'] = os.path.basename(seq_file)
    log['SET_FIELD'] = set_field
    log['CLUSTER_FIELD'] = cluster_field
    log['CLUSTER_PREFIX'] = cluster_prefix
    log['IDENTITY'] = ident
    log['SEQUENCE_START'] = seq_start
    log['SEQUENCE_END'] = seq_end
    log['NPROC'] = nproc
    printLog(log)

    # Define cluster function parameters
    cluster_args = {'cluster_exec':cluster_exec,
                    'ident':ident,
                    'seq_start':seq_start,
                    'seq_end':seq_end}

    # Define feeder function and arguments
    index_args = {'field': set_field, 'delimiter': out_args['delimiter']}
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets,
                 'index_args': index_args}
    # Define worker function and arguments
    work_func = processCSQueue
    work_args = {'cluster_field': cluster_field,
                 'cluster_prefix': cluster_prefix,
                 'cluster_args': cluster_args,
                 'delimiter': out_args['delimiter']}
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'cluster',
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
    log['END'] = 'ClusterSets'
    printLog(log)

    return result['out_files']


def clusterAll(seq_file, cluster_field=default_cluster_field,
               cluster_prefix=default_cluster_prefix,
               ident=default_ident, seq_start=0, seq_end=None,
               cluster_exec=default_cluster_exec,
               out_args=default_out_args, nproc=None):
    """
    Performs clustering on sets of sequences

    Arguments:
      seq_file : the sample sequence file name
      cluster_field : the name of the output cluster field
      cluster_prefix : string defining a prefix for the cluster identifier.
      ident : the identity threshold for clustering sequences
      seq_start : the start position to trim sequences at before clustering
      seq_end : the end position to trim sequences at before clustering
      cluster_exec : the path to the executable for usearch
      nproc : the number of processQueue processes;
              if None defaults to the number of CPUs

    Returns:
      str : the clustered output file name
    """
    # Function to modify SeqRecord header with cluster identifier
    def _header(seq, cluster, field=cluster_field, prefix=cluster_prefix,
                delimiter=out_args['delimiter']):
        label = '%s%i' % (prefix, cluster)
        header = parseAnnotation(seq.description, delimiter=delimiter)
        header = mergeAnnotation(header, {field: label}, delimiter=delimiter)
        seq.id = seq.name = flattenAnnotation(header, delimiter=delimiter)
        seq.description = ''
        return seq

    # Print parameter info
    log = OrderedDict()
    log['START'] = 'ClusterSets'
    log['COMMAND'] = 'all'
    log['FILE'] = os.path.basename(seq_file)
    log['CLUSTER_FIELD'] = cluster_field
    log['CLUSTER_PREFIX'] = cluster_prefix
    log['IDENTITY'] = ident
    log['SEQUENCE_START'] = seq_start
    log['SEQUENCE_END'] = seq_end
    log['NPROC'] = nproc
    printLog(log)

    # Count sequence file and parse into a list of SeqRecords
    result_count = countSeqFile(seq_file)
    seq_iter = readSeqFile(seq_file)

    # Perform clustering
    start_time = time()
    printMessage('Running uclust', start_time=start_time, width=25)
    cluster_dict = runUClust(seq_iter, ident=ident,
                             seq_start=seq_start, seq_end=seq_end,
                             threads=nproc, cluster_exec=cluster_exec)
    printMessage('Done', start_time=start_time, end=True, width=25)

    # Determine file type
    if out_args['out_type'] is None:
        out_args['out_type'] = getFileType(seq_file)

    # Open output file handles
    pass_handle = getOutputHandle(seq_file,
                                  'cluster-pass',
                                  out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'],
                                  out_type=out_args['out_type'])

    # Open indexed sequence file
    seq_dict = readSeqFile(seq_file, index=True)

    # Iterate over sequence records and update header with cluster annotation
    start_time = time()
    rec_count = pass_count = 0
    for cluster, id_list in cluster_dict.items():
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += len(id_list)

        # Define output sequences
        seq_output = [_header(seq_dict[x], cluster) for x in id_list]

        # Write output
        pass_count += len(seq_output)
        SeqIO.write(seq_output, pass_handle, out_args['out_type'])

    # Update progress
    printProgress(rec_count, result_count, 0.05, start_time)

    # Print log
    log = OrderedDict()
    log['SEQUENCES'] = result_count
    log['PASS'] = pass_count
    log['CLUSTERS'] = len(cluster_dict)
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['END'] = 'ClusterSets'
    printLog(log)

    # Close handles
    pass_handle.close()

    return pass_handle.name


def clusterBarcodes(seq_file, barcode_field=default_barcode_field,
                    cluster_field=default_cluster_field,
                    cluster_prefix=default_cluster_prefix,
                    ident=default_ident, cluster_exec=default_cluster_exec,
                    out_args=default_out_args, nproc=None):
    """
    Performs clustering on sets of sequences

    Arguments:
      seq_file : the sample sequence file name
      barcode_field : the annotation field containing barcode sequences
      cluster_field : the name of the output cluster field
      cluster_prefix : string defining a prefix for the cluster identifier.
      ident : the identity threshold for clustering sequences
      seq_start : the start position to trim sequences at before clustering
      seq_end : the end position to trim sequences at before clustering
      cluster_exec : the path to the executable for usearch
      nproc : the number of processQueue processes;
              if None defaults to the number of CPUs

    Returns:
      str : the clustered output file name
    """
    g_vec = 16 *'G'
    
    # Function to modify SeqRecord header with cluster identifier
    def _header(seq, cluster, field=cluster_field, prefix=cluster_prefix,
                delimiter=out_args['delimiter']):
        label = '%s%i' % (prefix, cluster)
        header = parseAnnotation(seq.description, delimiter=delimiter)
        header = mergeAnnotation(header, {field: label}, delimiter=delimiter)
        seq.id = seq.name = flattenAnnotation(header, delimiter=delimiter)
        seq.description = ''
        return seq

    # Function to extract to make SeqRecord object from a barcode annotation
    def _barcode(seq, field=barcode_field, delimiter=out_args['delimiter']):
        header = parseAnnotation(seq.description, delimiter=delimiter)
        return SeqRecord(Seq(header[field] + g_vec), id=seq.id)

    # Print parameter info
    log = OrderedDict()
    log['START'] = 'ClusterSets'
    log['COMMAND'] = 'barcode'
    log['FILE'] = os.path.basename(seq_file)
    log['BARCODE_FIELD'] = barcode_field
    log['CLUSTER_FIELD'] = cluster_field
    log['CLUSTER_PREFIX'] = cluster_prefix
    log['IDENTITY'] = ident
    log['NPROC'] = nproc
    printLog(log)

    # Count sequence file and parse into a list of SeqRecords
    result_count = countSeqFile(seq_file)
    barcode_iter = (_barcode(x) for x in readSeqFile(seq_file))

    # Perform clustering
    start_time = time()
    printMessage('Running uclust', start_time=start_time, width=25)
    cluster_dict = runUClust(barcode_iter, ident=ident, seq_start=0, seq_end=None,
                             threads=nproc, cluster_exec=cluster_exec)
    printMessage('Done', start_time=start_time, end=True, width=25)

    # Determine file type
    if out_args['out_type'] is None:
        out_args['out_type'] = getFileType(seq_file)

    # Open output file handles
    pass_handle = getOutputHandle(seq_file,
                                  'cluster-pass',
                                  out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'],
                                  out_type=out_args['out_type'])

    # Open indexed sequence file
    seq_dict = readSeqFile(seq_file, index=True)

    # Iterate over sequence records and update header with cluster annotation
    start_time = time()
    rec_count = pass_count = 0
    for cluster, id_list in cluster_dict.items():
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += len(id_list)

        # Define output sequences
        seq_output = [_header(seq_dict[x], cluster) for x in id_list]

        # Write output
        pass_count += len(seq_output)
        SeqIO.write(seq_output, pass_handle, out_args['out_type'])

    # Update progress
    printProgress(rec_count, result_count, 0.05, start_time)

    # Print log
    log = OrderedDict()
    log['SEQUENCES'] = result_count
    log['PASS'] = pass_count
    log['CLUSTERS'] = len(cluster_dict)
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['END'] = 'ClusterSets'
    printLog(log)

    # Close handles
    pass_handle.close()

    return pass_handle.name


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
                 cluster-pass
                    clustered reads.
                 cluster-fail
                    raw reads failing clustering.

             output annotation fields:
                 CLUSTER
                    a numeric cluster identifier defining the within-group cluster.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter, add_help=False)
    group_help = parser.add_argument_group('help')
    group_help.add_argument('--version', action='version',
                            version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    group_help.add_argument('-h', '--help', action='help', help='show this help message and exit')
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Clustering method')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Parent parser
    parent_parser = getCommonArgParser(multiproc=True)

    # Sequence set clustering arguments
    parser_set = subparsers.add_parser('set', parents=[parent_parser],
                                       formatter_class=CommonHelpFormatter, add_help=False,
                                       help='Cluster sequences within annotation sets.',
                                       description='Cluster sequences within annotation sets.')
    group_set = parser_set.add_argument_group('clustering arguments')
    group_set.add_argument('-f', action='store', dest='set_field', type=str,
                           default=default_barcode_field,
                           help='''The annotation field containing annotations, such as UMI
                                barcode, for sequence grouping.''')
    group_set.add_argument('-k', action='store', dest='cluster_field', type=str,
                           default=default_cluster_field,
                           help='''The name of the output annotation field to add with the
                                cluster information for each sequence.''')
    group_set.add_argument('--ident', action='store', dest='ident', type=float,
                           default=default_ident,
                           help='The sequence identity threshold for the uclust algorithm.')
    group_set.add_argument('--start', action='store', dest='seq_start', type=int, default=0,
                           help='''The start of the region to be used for clustering.
                                Together with --end, this parameter can be used to specify a
                                subsequence of each read to use in the clustering algorithm.''')
    group_set.add_argument('--end', action='store', dest='seq_end', type=int,
                           help='The end of the region to be used for clustering.')
    group_set.add_argument('--prefix', action='store', dest='cluster_prefix', type=str,
                           default=default_cluster_prefix,
                           help='''A string to use as the prefix for each cluster identifier.
                                By default, cluster identifiers will be numeric values only.''')
    group_set.add_argument('--exec', action='store', dest='cluster_exec',
                           default=default_cluster_exec,
                           help='The name or location of the usearch or vsearch executable.')
    parser_set.set_defaults(func=clusterSets)

    # Total sequence clustering arguments
    parser_all = subparsers.add_parser('all', parents=[parent_parser],
                                       formatter_class=CommonHelpFormatter, add_help=False,
                                       help='Cluster all sequences regardless of annotation.',
                                       description='Cluster all sequences regardless of annotation.')
    group_all = parser_all.add_argument_group('clustering arguments')
    group_all.add_argument('-k', action='store', dest='cluster_field', type=str,
                           default=default_cluster_field,
                           help='''The name of the output annotation field to add with the
                                cluster information for each sequence.''')
    group_all.add_argument('--ident', action='store', dest='ident', type=float,
                           default=default_ident,
                           help='The sequence identity threshold for the uclust algorithm.')
    group_all.add_argument('--start', action='store', dest='seq_start', type=int,
                           help='''The start of the region to be used for clustering.
                                Together with --end, this parameter can be used to specify a
                                subsequence of each read to use in the clustering algorithm.''')
    group_all.add_argument('--end', action='store', dest='seq_end', type=int,
                           help='The end of the region to be used for clustering.')
    group_all.add_argument('--prefix', action='store', dest='cluster_prefix', type=str,
                           default=default_cluster_prefix,
                           help='''A string to use as the prefix for each cluster identifier.
                                By default, cluster identifiers will be numeric values only.''')
    group_all.add_argument('--exec', action='store', dest='cluster_exec',
                           default=default_cluster_exec,
                           help='The name or location of the usearch or vsearch executable.')
    parser_all.set_defaults(func=clusterAll)

    # Sequence set clustering arguments
    parser_barcode = subparsers.add_parser('barcode', parents=[parent_parser],
                                           formatter_class=CommonHelpFormatter, add_help=False,
                                           help='Cluster reads by clustering barcode sequences.',
                                           description='Cluster reads by clustering barcode sequences.')
    group_barcode = parser_barcode.add_argument_group('clustering arguments')
    group_barcode.add_argument('-f', action='store', dest='barcode_field', type=str,
                               default=default_barcode_field,
                               help='''The annotation field containing barcode sequences.''')
    group_barcode.add_argument('-k', action='store', dest='cluster_field', type=str,
                               default=default_cluster_field,
                               help='''The name of the output annotation field to add with the
                                    cluster information for each sequence.''')
    group_barcode.add_argument('--ident', action='store', dest='ident', type=float,
                               default=default_ident,
                               help='The sequence identity threshold for the uclust algorithm.')
    group_barcode.add_argument('--prefix', action='store', dest='cluster_prefix', type=str,
                               default=default_cluster_prefix,
                               help='''A string to use as the prefix for each cluster identifier.
                                    By default, cluster identifiers will be numeric values only.''')
    group_barcode.add_argument('--exec', action='store', dest='cluster_exec',
                               default=default_cluster_exec,
                               help='The name or location of the usearch or vsearch executable.')
    parser_barcode.set_defaults(func=clusterBarcodes)

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

    # Convert fields to uppercase
    if 'set_field' in args_dict and args_dict['set_field'] is not None:
        args_dict['set_field'] = args_dict['set_field'].upper()
    if 'cluster_field' in args_dict and args_dict['cluster_field'] is not None:
        args_dict['cluster_field'] = args_dict['cluster_field'].upper()
    
    # Check if a valid usearch executable was specified
    if not shutil.which(args.cluster_exec):
        parser.error('%s does not exist' % args.cluster_exec)

    # Check for valid start and end input
    if ('seq_start' in args_dict and 'seq_end' in args_dict) and \
            args_dict['seq_start'] is not None and args_dict['seq_end'] is not None and \
            args_dict['seq_start'] >= args_dict['seq_end']:
        parser.error('--start must be less than --end')

    # Call cluster main function for each input file
    del args_dict['seq_files']
    del args_dict['func']
    del args_dict['command']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        args.func(**args_dict)