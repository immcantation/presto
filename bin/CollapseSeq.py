#!/usr/bin/env python3
"""
Removes duplicate sequences from FASTA/FASTQ files
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import os
import re
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import chain
from textwrap import dedent
from time import time
from Bio import SeqIO

# Presto imports
from presto.Defaults import default_delimiter, default_out_args
from presto.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation, \
                              collapseAnnotation
from presto.Sequence import checkSeqEqual
from presto.IO import getFileType, readSeqFile, getOutputHandle, printLog, printProgress

# Default parameters
default_max_missing = 0

class DuplicateSet:
    """
    A class defining unique sequence sets

    Attributes:
      seq : SeqRecord.
      missing : missing character count of SeqRecord.
      annotation: copy field annotation dictionary.
      keys : list of sequence identifier.
      count : duplicate count.
    """
    # Instantiation
    def __init__(self, seq, key, missing, annotations):
        self.seq = seq
        self.missing = missing
        self.annotations = annotations
        self.keys = [key]
        self.count = 1

    # Set length evaluation to number of duplicates
    def __len__(self):
        return len(self.keys)


def findUID(uid, search_dict, score=False):
    """
    Checks if a unique identifier is already present in a unique dictionary

    Arguments: 
    uid = the unique identifier key to check
    search_dict = a dictionary to search for key matches in
    score = if True score sequence element of the uid against each sequence in search_list

    Returns: 
    uid of match if found; None otherwise
    """
    match = None
    # Check for exact matches
    if not score:
        match = uid if uid in search_dict else None
    # Check for ambiguous matches
    else:
        for key in search_dict:
            if uid[1:] == key[1:] and checkSeqEqual(uid[0], key[0]):
                match = key 
                break
    
    # Return search boolean
    return match


def findUniqueSeq(uniq_dict, search_keys, seq_dict, max_missing=default_max_missing, 
                  uniq_fields=None, copy_fields=None, max_field=None, min_field=None, 
                  inner=False, delimiter=default_delimiter):
    """
    Finds unique sequences 

    Arguments: 
    uniq_dict = a dictionary of unique sequences generated by findUniqueSeq()
    search_keys = a list containing the subset of dictionary keys to be checked
    seq_dict = a SeqRecords dictionary generated by SeqIO.index()
    max_missing = the number of missing characters to allow in a unique sequences
    uniq_fields = a list of annotations that define a sequence as unique if they differ
    copy_fields = a list of annotations to copy into unique sequence annotations
    max_field = a numeric field whose maximum value determines the retained sequence
    min_field = a numeric field whose minimum value determines the retained sequence
    inner = if True exclude consecutive outer ambiguous characters from iterations and matching
    delimiter = description field delimiter
    
    Returns: 
    a tuple of (uniq_dict, search_keys, dup_keys) modified from passed values
    """
    # Define local variables
    ambig_re = re.compile(r'[\.\-N]')
    score = (max_missing > 0)
    dup_keys = []
    to_remove = []
    
    start_time = time()
    result_count = len(search_keys)
    # Iterate over search keys and update uniq_dict and dup_keys
    for idx, key in enumerate(search_keys):
        # Print progress of previous iteration
        printProgress(idx, result_count, 0.05, start_time, task='%i missing' % max_missing)
        
        # Define sequence to process
        seq = seq_dict[key]
        seq_str = str(seq.seq)
        if inner:  seq_str = seq_str.strip('.-N')
        
        # Skip processing of ambiguous sequences over max_missing threshold 
        ambig_count = len(ambig_re.findall(seq_str))
        if ambig_count > max_missing:  continue
        
        # Parse annotation and define unique identifiers (uid)
        if uniq_fields is not None:
            ann = parseAnnotation(seq_dict[key].description, uniq_fields, delimiter=delimiter)
            uid = tuple(chain([seq_str], list(ann.values())))             
        else:
            uid = (seq_str, None)

        # Parse annotation and define copied identifiers (cid)        
        if copy_fields is not None:
            ann = parseAnnotation(seq.description, copy_fields, delimiter=delimiter)
            cid = {k:[ann.get(k)] for k in copy_fields}
        else:
            cid = {}

        # Store new unique sequences and process duplicates
        match = findUID(uid, uniq_dict, score)
        if match is None:
            uniq_dict[uid] = DuplicateSet(seq, key=key, missing=ambig_count, annotations=cid)
        else:
            # Updated sequence, count, ambiguous character count, and count sets
            dup_key = key
            uniq_dict[match].count += 1
            uniq_dict[match].keys.append(key)
            for k, v in cid.items():
                uniq_dict[match].annotations[k].extend(v)
            # Check whether to replace previous unique sequence with current sequence
            if ambig_count <= uniq_dict[match].missing:
                swap = False
                seq_last = uniq_dict[match].seq
                if max_field is not None:
                    swap = float(parseAnnotation(seq.description, delimiter=delimiter)[max_field]) > \
                           float(parseAnnotation(seq_last.description, delimiter=delimiter)[max_field])
                elif min_field is not None:
                    swap = float(parseAnnotation(seq.description, delimiter=delimiter)[min_field]) > \
                           float(parseAnnotation(seq_last.description, delimiter=delimiter)[min_field])
                # TODO:  quality evaluation is a bottleneck
                else:
                    if hasattr(seq, 'letter_annotations') and 'phred_quality' in seq.letter_annotations:
                        q_this = float(sum(seq.letter_annotations['phred_quality'])) / len(seq)
                        q_last = float(sum(seq_last.letter_annotations['phred_quality'])) / len(seq_last)
                        swap = q_this > q_last
                # Replace old sequence if criteria passed
                if swap:
                    dup_key = seq_last.id
                    uniq_dict[match].seq = seq
                    uniq_dict[match].missing = ambig_count

            # Update duplicate list
            dup_keys.append(dup_key)

        # Mark seq for removal from later steps
        to_remove.append(idx)
        
    # Remove matched sequences from search_keys
    for j in reversed(to_remove):  del search_keys[j]

    # Update progress
    printProgress(result_count, result_count, 0.05, start_time, task='%i missing' % max_missing)
        
    return (uniq_dict, search_keys, dup_keys)


def collapseSeq(seq_file, max_missing=default_max_missing, uniq_fields=None,
                copy_fields=None, copy_actions=None, max_field=None, min_field=None, 
                inner=False, keep_missing=False, out_args=default_out_args):
    """
    Removes duplicate sequences from a file

    Arguments: 
    seq_file = filename of the sequence file to sample from
    max_missing = number of ambiguous characters to allow in a unique sequence
    uniq_fields = a list of annotations that define a sequence as unique if they differ
    copy_fields = a list of annotations to copy into unique sequence annotations
    copy_actions = the list of collapseAnnotation actions to take on copy_fields 
    max_field = a numeric field whose maximum value determines the retained sequence
    min_field = a numeric field whose minimum value determines the retained sequence
    inner = if True exclude consecutive outer ambiguous characters from iterations and matching
    keep_missing = if True retain sequences with more ambiguous characters than max_missing as unique
    out_args = common output argument dictionary from parseCommonArgs
              
    Returns: 
    the collapsed output file name
    """
    log = OrderedDict()
    log['START'] = 'CollapseSeq'
    log['FILE'] = os.path.basename(seq_file)
    log['MAX_MISSING'] = max_missing
    log['UNIQ_FIELDS'] = ','.join([str(x) for x in uniq_fields]) \
                         if uniq_fields is not None else None
    log['COPY_FIELDS'] = ','.join([str(x) for x in copy_fields]) \
                         if copy_fields is not None else None
    log['COPY_ACTIONS'] = ','.join([str(x) for x in copy_actions]) \
                          if copy_actions is not None else None
    log['MAX_FIELD'] = max_field
    log['MIN_FIELD'] = min_field
    log['INNER'] = inner
    log['KEEP_MISSING'] = keep_missing
    printLog(log)
    
    # TODO:  storing all sequences in memory is faster
    # Read input file
    in_type = getFileType(seq_file)
    #seq_dict = readSeqFile(seq_file, index=True)
    seq_dict = SeqIO.to_dict(readSeqFile(seq_file, index=False))
    if out_args['out_type'] is None:  out_args['out_type'] = in_type

    # Count total sequences
    rec_count = len(seq_dict)

    # Define log handle
    if out_args['log_file'] is None:  
        log_handle = None
    else:  
        log_handle = open(out_args['log_file'], 'w')

    # Find sequences with duplicates
    uniq_dict = {}
    # Added list typing for compatibility issue with Python 2.7.5 on OS X
    # TypeError: object of type 'dictionary-keyiterator' has no len()
    search_keys = list(seq_dict.keys())
    dup_keys = []
    for n in range(0, max_missing + 1):
        # Find unique sequences
        uniq_dict, search_keys, dup_list = findUniqueSeq(uniq_dict, search_keys, seq_dict, n, 
                                                         uniq_fields, copy_fields,
                                                         max_field, min_field, inner, 
                                                         out_args['delimiter'])

        # Update list of duplicates
        dup_keys.extend(dup_list)
                
        # Break if no keys to search remain
        if len(search_keys) == 0:  break
    
    # Write unique sequences
    with getOutputHandle(seq_file, 'collapse-unique', out_dir=out_args['out_dir'], 
                         out_name=out_args['out_name'], out_type=out_args['out_type']) \
            as uniq_handle:
        for val in uniq_dict.values():
            # Define output sequence
            out_seq = val.seq
            out_ann = parseAnnotation(out_seq.description, delimiter=out_args['delimiter'])
            out_app = OrderedDict()
            if copy_fields  is not None and copy_actions is not None:
                for f, a in zip(copy_fields, copy_actions):
                    x = collapseAnnotation(val.annotations, a, f, delimiter=out_args['delimiter'])
                    out_app[f] = x[f]
                    out_ann.pop(f, None)
            out_app['DUPCOUNT'] = val.count
            out_ann = mergeAnnotation(out_ann, out_app, delimiter=out_args['delimiter'])
            out_seq.id = out_seq.name = flattenAnnotation(out_ann, delimiter=out_args['delimiter'])
            out_seq.description = ''
            # Write unique sequence
            SeqIO.write(out_seq, uniq_handle, out_args['out_type'])

            # Update log
            log = OrderedDict()
            log['HEADER'] = out_seq.id
            log['DUPCOUNT'] = val.count
            for i, k in enumerate(val.keys, start=1):
                log['ID%i' % i] = k
            for i, k in enumerate(val.keys, start=1):
                log['SEQ%i' % i] = str(seq_dict[k].seq)
            printLog(log, handle=log_handle)

        # Write sequence with high missing character counts
        if keep_missing:
            for k in search_keys:
                out_seq = seq_dict[k]
                out_ann = parseAnnotation(out_seq.description, delimiter=out_args['delimiter'])
                out_ann = mergeAnnotation(out_ann, {'DUPCOUNT':1}, delimiter=out_args['delimiter'])
                out_seq.id = out_seq.name = flattenAnnotation(out_ann, delimiter=out_args['delimiter'])
                out_seq.description = ''
                SeqIO.write(out_seq, uniq_handle, out_args['out_type'])

    # Write sequence with high missing character counts
    if out_args['failed'] and not keep_missing:
        with getOutputHandle(seq_file, 'collapse-undetermined', out_dir=out_args['out_dir'],
                             out_name=out_args['out_name'], out_type=out_args['out_type']) \
                as missing_handle:
            for k in search_keys:
                SeqIO.write(seq_dict[k], missing_handle, out_args['out_type'])

    if out_args['failed']:
        # Write duplicate sequences 
        with getOutputHandle(seq_file, 'collapse-duplicate', out_dir=out_args['out_dir'], 
                             out_name=out_args['out_name'], out_type=out_args['out_type']) \
                as dup_handle:
            for k in dup_keys:
                SeqIO.write(seq_dict[k], dup_handle, out_args['out_type'])

    # Print log
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(uniq_handle.name)
    log['SEQUENCES'] = rec_count
    log['UNIQUE'] = len(uniq_dict)
    log['DUPLICATE'] = len(dup_keys)
    log['UNDETERMINED'] = len(search_keys)
    log['END'] = 'CollapseSeq'
    printLog(log)
        
    # Close file handles
    if log_handle is not None:  log_handle.close()
    
    return uniq_handle.name
    

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
                 collapse-unique
                     unique sequences. Contains one representative from each set of
                     duplicate sequences. The retained representative is determined by
                     user defined criteria.
                 collapse-duplicate
                     raw reads which are duplicates of the sequences retained in the
                     collapse-unique file.
                 collapse-undetermined
                     raw reads which were excluded from consideration due to having too
                     many N characters in the sequence.

             output annotation fields:
                 DUPCOUNT
                     total number of sequences within the set of duplicates for each
                     retained unique sequence. Meaning, the copy number of each unique
                     sequence within the data file.
                 <user defined>
                     annotation fields specified by the --cf parameter.
             ''')

    # TODO: add exact mode which does hash table search only
    # TODO: write better algorithm for ambiguous character mode
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser()], 
                            formatter_class=CommonHelpFormatter, add_help=False)

    # Collapse arguments
    group_dedup = parser.add_argument_group('collapse arguments')
    group_dedup.add_argument('-n', action='store', dest='max_missing', type=int, default=default_max_missing,
                             help='''Maximum number of missing nucleotides to consider for collapsing
                                   sequences. A sequence will be considered undetermined if it contains too
                                   many missing nucleotides.''')
    group_dedup.add_argument('--uf', nargs='+', action='store', dest='uniq_fields', type=str, default=None,
                             help='''Specifies a set of annotation fields that must match for sequences
                                  to be considered duplicates.''')
    group_dedup.add_argument('--cf', nargs='+', action='store', dest='copy_fields', type=str, default=None,
                             help='''Specifies a set of annotation fields to copy into the unique
                                  sequence output.''')
    group_dedup.add_argument('--act', nargs='+', action='store', dest='copy_actions', default=None,
                             choices=['min', 'max', 'sum', 'set'],
                             help='''List of actions to take for each copy field which defines how
                                  each annotation will be combined into a single value. The actions
                                  "min", "max", "sum" perform the corresponding mathematical
                                  operation on numeric annotations. The action "set" collapses
                                  annotations into a comma delimited list of unique values.''')
    group_dedup.add_argument('--inner', action='store_true', dest='inner',
                             help='''If specified, exclude consecutive missing characters at either end of
                                  the sequence.''')
    group_dedup.add_argument('--keepmiss', action='store_true', dest='keep_missing',
                             help='''If specified, sequences with more missing characters than the
                                  threshold set by the -n parameter will be written to the unique
                                  sequence output file with a DUPCOUNT=1 annotation. If not specified,
                                  such sequences will be written to a separate file.''')

    # Mutually exclusive argument group
    group_field = group_dedup.add_mutually_exclusive_group()
    group_field.add_argument('--maxf', action='store', dest='max_field', type=str, default=None,
                             help='''Specify the field whose maximum value determines the retained sequence;
                                  mutually exclusive with --minf.''')
    group_field.add_argument('--minf', action='store', dest='min_field', type=str, default=None,
                             help='''Specify the field whose minimum value determines the retained sequence;
                                   mutually exclusive with --minf.''')

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
    
    # Convert case of fields
    if 'uniq_fields' in args_dict and args_dict['uniq_fields']:  
        args_dict['uniq_fields'] = list(map(str.upper, args_dict['uniq_fields'])) 
    if 'copy_fields' in args_dict and args_dict['copy_fields']:
        args_dict['copy_fields'] = list(map(str.upper, args_dict['copy_fields']))
    if 'copy_actions' in args_dict and args_dict['copy_actions']:
        args_dict['copy_actions'] = list(map(str.lower, args_dict['copy_actions']))
    if 'max_field' in args_dict and args_dict['max_field']:  
        args_dict['max_field'] = args_dict['max_field'].upper() 
    if 'min_field' in args_dict and args_dict['min_field']:  
        args_dict['min_field'] = args_dict['min_field'].upper()
    
    # Check copy field and action arguments
    if bool(args_dict['copy_fields']) ^ bool(args_dict['copy_actions']) or \
       len((args_dict['copy_fields'] or '')) != len((args_dict['copy_actions'] or '')):
            parser.error('You must specify exactly one copy action (--act) per copy field (--cf)')
    
    # Call appropriate function for each sample file
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        collapseSeq(**args_dict)
