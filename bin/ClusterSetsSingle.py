#!/usr/bin/env python3
"""
Cluster Sequences (Wrapper for UCLUST)
"""
# Info
__author__ = 'Ruoyi Jiang, Christopher Bolen, Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import os
import shutil
import sys
import tempfile
from argparse import ArgumentParser
from collections import OrderedDict
from textwrap import dedent

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field, \
                            default_cluster_field, default_out_args, \
                            default_usearch_exec
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation
from presto.Applications import runUClust
from presto.IO import printLog, readSeqFile
from presto.Sequence import indexSeqSets
from presto.Multiprocessing import SeqData, SeqResult, manageProcesses, feedSeqQueue, \
                                   collectSeqQueue

from Bio.SeqRecord import SeqRecord


# Defaults
default_cluster_exec = default_usearch_exec
default_ident = 0.9
default_cluster_ident = 0.9



#Temporary. runUClust uses a different set of run flags as of Feb 2017. Should be updated in Applications API. 
#In particular, not the absence of a thread specification.

import tempfile
import csv
from subprocess import CalledProcessError, check_output, PIPE, Popen, STDOUT
from Bio.Seq import Seq
from Bio import SeqIO

def runUClust(seq_list, ident=default_cluster_ident, seq_start=0, seq_end=None,
              cluster_exec=default_usearch_exec):
    """
    Cluster a set of sequences using the UCLUST algorithm from USEARCH

    Arguments:
      seq_list : a list of SeqRecord objects to align.
      ident : the sequence identity cutoff to be passed to usearch.
      seq_start : the start position to trim sequences at before clustering.
      seq_end : the end position to trim sequences at before clustering.
      cluster_exec : the path to the usearch executable.

    Returns:
      dict : {sequence id: cluster id}.
    """
    # Function to trim and mask sequences
    gap_trans = str.maketrans({'-': 'N', '.': 'N'})
    def _clean(rec, i, j):
        seq = str(rec.seq[i:j])
        seq = seq.translate(gap_trans)
        return SeqRecord(Seq(seq), id=rec.id, name=rec.name, description=rec.description)
    
    # Return sequence if only one sequence in seq_list
    if len(seq_list) < 2:
        return {1:[seq_list[0].id]}
    
    # Make a trimmed and masked copy of each sequence so we don't mess up originals
    seq_trimmed = [_clean(x, seq_start, seq_end) for x in seq_list]
    
    # If there are any empty sequences after trimming return None
    if any([len(x.seq) == 0 for x in seq_trimmed]):
        return None
    
    # Open temporary files
    in_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    out_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    
    # Define usearch command
    #/usr/local/bin/usearch --input /tmp/tmpg24_0tae --uc /tmp/tmp2na2_on8 --id 0.9 --minlen 1 --usersort
    cmd = [cluster_exec,
           '--input', in_handle.name,
           '--uc', out_handle.name,
           '--id', str(ident),
           '--minlen', '1',
           '--usersort']
    
    # Write usearch input fasta file
    SeqIO.write(seq_trimmed, in_handle, 'fasta')
    in_handle.seek(0)
    
    # Run usearch uclust algorithm
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                              universal_newlines=True)
        #check_call(cmd, stderr=STDOUT, shell=False)
    except CalledProcessError:
        group_dict = None
    
    # TODO:  unsure about this return object.
    # Parse the results of usearch
    # Output columns for the usearch 'uc' output format
    #   0 = entry type -- S: centroid seq, H: hit, C: cluster record (redundant with S)
    #   1 = group the sequence is assigned to
    #   8 = the id of the sequence
    #   9 = id of the centroid for cluster
    cluster_dict = {}
    for row in csv.reader(out_handle, delimiter='\t'):
        if row[0] in ('H', 'S'):
            # Trim sequence label to portion before space for usearch v9 compatibility
            key = int(row[1]) + 1
            # Trim sequence label to portion before space for usearch v9 compatibility
            hit = row[8].split()[0]
            # Update cluster dictionary
            cluster = cluster_dict.setdefault(key, [])
            cluster.append(hit)
    
    return cluster_dict if cluster_dict else None
    


def feedSingleQueue(alive, data_queue, seq_file, index_func=None, index_args={}, cluster_args={}):
    """
    Runs UClust and feeds the data queue with SeqRecord objects

    Arguments:
      alive : multiprocessing.Value boolean controlling whether processing
              continues; when False function returns
      data_queue : multiprocessing.Queue to hold data for processing
      seq_file : Sequence file to read input from
      index_func : Function to use to define sequence sets
                   if None do not index sets and feed individual records
      index_args : Dictionary of arguments to pass to index_func

    Returns:
      None
    """
    def _returnSeqRecordsfromDict(dict):
        output = []
        for UID in dict:
            for seq_id in dict[UID]: 
                seq_record = SeqRecord(UID)
                seq_record.name = seq_record.id = seq_id
                output.append(seq_record)
        return output
    
    try:
        seq_dict = readSeqFile(seq_file, index=True) #seq_id + seq_record
        index_dict = index_func(seq_dict, **index_args) #UID + seq_id
        seq_list_umi = _returnSeqRecordsfromDict(index_dict) # list of SeqRecord(seq=UID, name/id = list(seq_id))
        
        #Runs Uclust (No threads specified.)
        try:
            clust_dict = runUClust(seq_list_umi, **cluster_args)
        except:
            raise
        
        #Separate raise because runUClust is used in other bin tools
        if clust_dict is None:
            sys.stderr.write('Uclust failure. No clusters generated or IO error.')
            raise

        #Generates a data_iter from UClust
        data_iter = ((k, [seq_dict[i] for i in v]) \
                         for k, v in clust_dict.items())

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



def processSingleQueue(alive, data_queue, result_queue, cluster_field,
                  cluster_args={}, delimiter=default_delimiter):
    """
    Pulls from data queue, and feeds results queue. 

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
            
            # Get number of clusters
            result.log['SEQS_IN_CLUST'] = len(data.data)

            # Update sequence annotations with cluster assignments
            
            results = list()
            for i,seq in enumerate(data.data):
                header = parseAnnotation(seq.description, delimiter=delimiter)
                header = mergeAnnotation(header, {cluster_field: clust}, delimiter=delimiter)
                seq.id = seq.name = flattenAnnotation(header, delimiter=delimiter)
                seq.description = ''
                result.log['SEQ:%i' % (i)] = str(seq.seq)
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



def clusterSetsSingle(seq_file, barcode_field=default_barcode_field,
                cluster_field=default_cluster_field,
                ident=default_ident, seq_start=None, seq_end=None,
                cluster_exec=default_cluster_exec,
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
    log['START'] = 'ClusterSetsSingle'
    log['FILE'] = os.path.basename(seq_file)
    log['BARCODE_FIELD'] = barcode_field
    log['CLUSTER_FIELD'] = cluster_field
    log['IDENTITY'] = ident
    log['SEQUENCE_START'] = seq_start
    log['SEQUENCE_END'] = seq_end
    log['NPROC'] = nproc
    printLog(log)

    # Define cluster function parameters
    cluster_args = {'ident':ident,
                    'seq_start':seq_start,
                    'seq_end':seq_end,
                    'cluster_exec':cluster_exec,}

    # Define feeder function and arguments
    index_args = {'field': barcode_field, 'delimiter': out_args['delimiter']}
    feed_func = feedSingleQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets,
                 'index_args': index_args}
    # Define worker function and arguments
    work_func = processSingleQueue
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
    for k, v in result['log'].items():  log[k] = v
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
                            parents=[getCommonArgParser(multiproc=True)],
                            formatter_class=CommonHelpFormatter, add_help=False)

    # Clustering arguments
    group_clust = parser.add_argument_group('clustering arguments')
    group_clust.add_argument('-f', action='store', dest='barcode_field', type=str,
                             default=default_barcode_field,
                             help='''The annotation field containing annotations, such as UID
                                  barcode, for sequence grouping.''')
    group_clust.add_argument('-k', action='store', dest='cluster_field', type=str,
                             default=default_cluster_field,
                             help='''The name of the output annotation field to add with the
                                  cluster information for each sequence.''')
    group_clust.add_argument('--id', action='store', dest='ident', type=float,
                             default=default_ident,
                             help='The sequence identity threshold for the uclust algorithm.')
    group_clust.add_argument('--start', action='store', dest='seq_start', type=int,
                             help='''The start of the region to be used for clustering.
                                  Together with --end, this parameter can be used to specify a
                                  subsequence of each read to use in the clustering algorithm.''')
    group_clust.add_argument('--end', action='store', dest='seq_end', type=int,
                             help='The end of the region to be used for clustering.')
    group_clust.add_argument('--exec', action='store', dest='cluster_exec',
                             default=default_cluster_exec,
                             help='The name or location of the usearch or vsearch executable.')

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
    if not shutil.which(args.cluster_exec):
        parser.error('%s does not exist' % args.cluster_exec)

    # Check for valid start and end input
    if ('seq_start' in args_dict and 'seq_end' in args_dict) and \
            args_dict['seq_start'] is not None and args_dict['seq_end'] is not None and \
            args_dict['seq_start'] >= args_dict['seq_end']:
        parser.error('--start must be less than --end')
    
        
    # Call cluster for each input file
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        clusterSetsSingle(**args_dict)

