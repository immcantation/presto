#!/usr/bin/env python
"""
Cluster sequences by group
"""

__author__    = 'Christopher Bolen, Jason Vander Heiden'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.6'
__date__      = '2015.05.13'

# Imports
import csv, os, sys, tempfile, textwrap, time
from argparse import ArgumentParser
from collections import OrderedDict
from subprocess import CalledProcessError, check_output, STDOUT
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_delimiter, default_out_args
from IgCore import default_barcode_field
from IgCore import parseAnnotation, mergeAnnotation, flattenAnnotation
from IgCore import CommonHelpFormatter, getCommonArgParser, parseCommonArgs, printLog
from IgCore import indexSeqSets, collectSeqQueue, feedSeqQueue
from IgCore import manageProcesses, SeqResult

# Defaults
default_cluster_field = r'CLUSTER'
default_usearch_exec = r'/usr/local/bin/usearch'
default_ident = 0.9


def runUClust(seq_list, ident=default_ident, seq_start=0, seq_end=None,
              usearch_exec=default_usearch_exec):
    """
    Cluster a set of sequences using the UCLUST algorithm from USEARCH

    Arguments: 
    seq_list = a list of SeqRecord objects to align
    ident = the sequence identity cutoff to be passed to usearch
    seq_start = the start position to trim sequences at before clustering
    seq_end = the end position to trim sequences at before clustering
    usearch_exec = the path to the usearch executable

    Returns: 
    a dictionary object containing {sequence id: cluster id}
    """
    # Return sequence if only one sequence in seq_list
    if len(seq_list) < 2:
        #return {seq_list[0].id:0}
        return {1:[seq_list[0].id]}

    # Format sequences and make a copy so we don't mess up original sequences
    short_list = list()
    for rec in seq_list:
        seq = rec.seq[seq_start:seq_end]
        seq = seq.ungap('-')
        seq = seq.ungap('.')
        short_list.append(SeqRecord(seq, id=rec.id, name=rec.name,
                                    description=rec.description))

    # Open temporary files
    in_handle = tempfile.NamedTemporaryFile()
    out_handle = tempfile.NamedTemporaryFile()

    # Define usearch command
    cmd = [usearch_exec,
           '-cluster_fast', in_handle.name,
           '-uc', out_handle.name,
           '-id', str(ident),
           '-threads', '1']

    # Write usearch input fasta file
    SeqIO.write(short_list, in_handle, 'fasta')
    in_handle.seek(0)

    # Run usearch uclust algorithm
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False)
        #check_call(cmd, stderr=STDOUT, shell=False)
    except CalledProcessError:
        group_dict = None
    else:
        # TODO:  unsure about this return object.
        # Parse the results of usearch
        # Output columns for the usearch 'uc' output format
        #   0 = entry type -- S: centroid seq, H: hit, C: cluster record (redundant with S)
        #   1 = group the sequence is assigned to
        #   8 = the id of the sequence
        group_dict = {}
        for row in csv.reader(out_handle, delimiter='\t'):
            if row[0] in ('S', 'H'):
                key = int(row[1]) + 1
                group = group_dict.setdefault(key, [])
                group.append(row[8])
        #out_list = [r for r in csv.reader(out_handle, delimiter='\t')]
        #group_dict = {r[8]: int(r[1]) + 1 for r in out_list if r[0] in ('S', 'H')}

    return group_dict if group_dict else None


def processCSQueue(alive, data_queue, result_queue, cluster_field,
                  cluster_args={}, delimiter=default_delimiter):
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
            result.log['BARCODE'] = data.id
            result.log['SEQCOUNT'] = len(data)
    
            # Perform clustering
            cluster_dict = runUClust(data.data, **cluster_args)

            # Process failed result
            if cluster_dict is None:
                # Update log
                result.log['CLUSTERS'] = 0
                for seq in enumerate(data.data, start=1):
                    result.log['CLUST0-%i' % i] = str(seq.seq)

                # Feed results queue and continue
                result_queue.put(result)
                continue

            # Get number of clusters
            result.log['CLUSTERS'] = len(cluster_dict)

            # Update sequence annotations with cluster assignments
            results = list()
            seq_dict = {s.id: s for s in data.data}
            for clust, id_list in cluster_dict.iteritems():
                for i, seq_id in enumerate(id_list, start=1):
                    # Add cluster annotation
                    seq = seq_dict[seq_id]
                    header = parseAnnotation(seq.description, delimiter=delimiter)
                    header = mergeAnnotation(header, {cluster_field:clust}, delimiter=delimiter)
                    seq.id = seq.name = flattenAnnotation(header, delimiter=delimiter)
                    seq.description = ''

                    # Update log and results
                    result.log['CLUST%i-%i' % (clust, i)] = str(seq.seq)
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


def clusterSets(seq_file, barcode_field=default_barcode_field,
                 cluster_field=default_cluster_field, 
                 ident=default_ident, seq_start=None, seq_end=None,
                 usearch_exec=default_usearch_exec,
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
    usearch_exec = the path to the executable for usearch
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                      
    Returns: 
    the clustered output file name
    """    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'ClusterSets'
    log['FILE'] = os.path.basename(seq_file)
    log['BARCODE_FIELD'] = barcode_field
    log['CLUSTER_FIELD'] = cluster_field
    log['IDENTITY'] = ident
    log['SEQUENCE_START'] = seq_start
    log['SEQUENCE_END'] = seq_end
    log['NPROC'] = nproc
    printLog(log)

    # Define cluster function parameters
    cluster_args = {'usearch_exec':usearch_exec,
                    'ident':ident,
                    'seq_start':seq_start,
                    'seq_end':seq_end}
 
    # Define feeder function and arguments
    index_args = {'field': barcode_field, 'delimiter': out_args['delimiter']}
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets,
                 'index_args': index_args}
    # Define worker function and arguments
    work_func = processCSQueue
    work_args = {'cluster_field': cluster_field, 
                 'cluster_args': cluster_args,
                 'delimiter': out_args['delimiter']}
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'cluster',
                    'out_args': out_args,
                    'index_field': barcode_field}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    log = OrderedDict()
    log['OUTPUT'] = result['log'].pop('OUTPUT')
    for k, v in result['log'].iteritems():  log[k] = v
    log['END'] = 'ClusterSets'
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
    fields = textwrap.dedent(
         '''
         output files:
           cluster-pass          clustered reads.
           cluster-fail          raw reads failing clustering.

         output annotation fields:
           CLUSTER               a numeric cluster identifier defining the within-group
                                 cluster.

         ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__),
                            parents=[getCommonArgParser(multiproc=True)], 
                            formatter_class=CommonHelpFormatter)
      
    parser.add_argument('--bf', action='store', dest='barcode_field', type=str,
                        default=default_barcode_field,
                        help='''The annotation field containing barcode labels for
                             sequence grouping.''')
        
    # usearch arguments
    parser.add_argument('--cf', action='store', dest='cluster_field', type=str, 
                        default=default_cluster_field,
                        help='''The name of the annotation field to add with the cluster
                             information for each sequence.''')
    parser.add_argument('--id', action='store', dest='ident', type=float,
                        default=default_ident,
                        help='The sequence identity threshold for the usearch algorithm.')
    parser.add_argument('--start', action='store', dest='seq_start', type=int,
                        help='''The start of the region to be used for clustering.
                             Together with --end, this parameter can be used to specify a
                             subsequence of each read to use in the clustering algorithm.''')
    parser.add_argument('--end', action='store', dest='seq_end', type=int,
                        help='The end of the region to be used for clustering.')
    parser.add_argument('--exec', action='store', dest='usearch_exec',
                        default=default_usearch_exec,
                        help='The location of the USEARCH executable.')
    
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
    if 'barcode_field' in args_dict and args_dict['barcode_field'] is not None:
        args_dict['barcode_field'] = args_dict['barcode_field'].upper()
    if 'cluster_field' in args_dict and args_dict['cluster_field'] is not None:
        args_dict['cluster_field'] = args_dict['cluster_field'].upper()
    
    # Check if a valid usearch executable was specified
    if not os.path.isfile(args.usearch_exec):
        parser.error('%s does not exist' % args.usearch_exec)
    
        
    # Call cluster for each input file
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        clusterSets(**args_dict)

