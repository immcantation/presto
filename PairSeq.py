#!/usr/bin/env python
"""
Sorts and matches sequence records with matching coordinates across files
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.5'
__date__      = '2014.10.2'

# Imports
import os, sys
from argparse import ArgumentParser
from collections import OrderedDict
from time import time
from Bio import SeqIO

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_coord_choices, default_coord_type, default_out_args
from IgCore import flattenAnnotation, mergeAnnotation, parseAnnotation
from IgCore import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from IgCore import getOutputHandle, printLog, printProgress, readSeqFile, getFileType
from IgCore import indexSeqPairs, getUnpairedIndex


def pairSeq(seq_file_1, seq_file_2, fields=None, coord_type=default_coord_type, 
            out_args=default_out_args):
    """
    Generates consensus sequences

    Arguments: 
    seq_file_1 = the file containing the grouped sequences and annotations
    seq_file_2 = the file to assign annotations to from seq_file_1
    fields = list of annotations in seq_file_1 records to copy to seq_file_2 records;
             if None do not copy an annotation
    coord_type = the sequence header format
    out_args = common output argument dictionary from parseCommonArgs
                    
    Returns: 
    a list of tuples holding successfully paired filenames for (seq_file_1, seq_file_2)
    """
    log = OrderedDict()
    log['START'] = 'PairSeq'
    log['FILE1'] = os.path.basename(seq_file_1) 
    log['FILE2'] = os.path.basename(seq_file_2)
    log['COORD_TYPE'] = coord_type
    printLog(log)
    
    # Read input files and open output files
    in_type_1 = getFileType(seq_file_1)
    seq_dict_1 = readSeqFile(seq_file_1, index=True)
    in_type_2 = getFileType(seq_file_2)
    seq_dict_2 = readSeqFile(seq_file_2, index=True)

    # Define output type
    if out_args['out_type'] is None:
        out_type_1 = in_type_1
        out_type_2 = in_type_2
    else: 
        out_type_1 = out_type_2 = out_args['out_type']

    # Define output name
    if out_args['out_name'] is None:
        out_name_1 = out_name_2 = None
    else: 
        out_name_1 = '%s-1' % out_args['out_name']
        out_name_2 = '%s-2' % out_args['out_name']
    
    # Find paired sequences
    index_dict = indexSeqPairs(seq_dict_1, seq_dict_2, coord_type, delimiter=out_args['delimiter'])
    
    # Write unmatched entries to files
    if not out_args['clean']:
        unpaired_1, unpaired_2 = getUnpairedIndex(seq_dict_1, seq_dict_2, coord_type, 
                                                  delimiter=out_args['delimiter'])
        # Write unpaired seq_file_1 records
        with getOutputHandle(seq_file_1, 'pair-fail', out_dir=out_args['out_dir'], 
                             out_name=out_name_1, out_type=out_type_1) as out_handle_1:
            for k in unpaired_1:
                SeqIO.write(seq_dict_1[k], out_handle_1, out_type_1)

        # Write unpaired seq_file_2 records
        with getOutputHandle(seq_file_2, 'pair-fail', out_dir=out_args['out_dir'], 
                             out_name=out_name_2, out_type=out_type_2) as out_handle_2:
            for k in unpaired_2:
                SeqIO.write(seq_dict_2[k], out_handle_2, out_type_2)
    
    # Open paired output file handles
    pass_handle_1 = getOutputHandle(seq_file_1, 'pair-pass', out_args['out_dir'], 
                                    out_name=out_name_1, out_type=out_type_1)
    pass_handle_2 = getOutputHandle(seq_file_2, 'pair-pass', out_args['out_dir'], 
                                    out_name=out_name_2, out_type=out_type_2)
    
    # Iterate over pairs and write to output files
    start_time = time()
    result_count = len(index_dict)
    pair_count = 0
    for key_1, key_2 in index_dict.itervalues():
        # Print progress for previous iteration
        printProgress(pair_count, result_count, 0.05, start_time)
        
        # Get records
        rec_1 = seq_dict_1[key_1]
        rec_2 = seq_dict_2[key_2]
        
        # Update counts for iteration
        pair_count += 1
        
        # Copy annotations from rec_1 to rec_2
        if fields is not None:
            copy_ann = parseAnnotation(rec_1.description, fields, delimiter=out_args['delimiter'])
            ann_2 = parseAnnotation(rec_2.description, delimiter=out_args['delimiter'])
            ann_2 = mergeAnnotation(ann_2, copy_ann, delimiter=out_args['delimiter'])
            rec_2.id = flattenAnnotation(ann_2, delimiter=out_args['delimiter']) 
            rec_2.description = ''
        
        # Write paired records
        SeqIO.write(rec_1, pass_handle_1, out_type_1)
        SeqIO.write(rec_2, pass_handle_2, out_type_2)

    # Print log
    printProgress(pair_count, result_count, 0.05, start_time)
    count_1, count_2 = len(seq_dict_1), len(seq_dict_2)
    log = OrderedDict()
    log['OUTPUT1'] = os.path.basename(pass_handle_1.name) 
    log['OUTPUT2'] = os.path.basename(pass_handle_2.name) 
    log['SEQUENCES1'] = count_1
    log['SEQUENCES2'] = count_2
    log['FAIL1'] = count_1 - pair_count
    log['FAIL2'] = count_2 - pair_count
    log['PASS'] = pair_count
    log['END'] = 'PairSeq'
    printLog(log)
   
    # Close file handles
    pass_handle_1.close()
    pass_handle_2.close()

    return [(pass_handle_1.name, pass_handle_2.name)]


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
                            parents=[getCommonArgParser(paired=True, log=False)], 
                            formatter_class=CommonHelpFormatter)
    
    parser.add_argument('-f', nargs='+', action='store', dest='fields', type=str, default=None, 
                        help='Specify the annotation fields to copy from file 1 records into file 2 records')
    parser.add_argument('--coord', action='store', dest='coord_type', 
                        choices=default_coord_choices, default=default_coord_type,
                        help='The format of the sequence identifier which defines shared coordinate \
                              information across paired ends')
    
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
    if args_dict['fields']:  args_dict['fields'] = map(str.upper, args_dict['fields'])

    # Call pairSeq
    del args_dict['seq_files_1']
    del args_dict['seq_files_2']
    for file_1, file_2 in zip(args.__dict__['seq_files_1'], 
                              args.__dict__['seq_files_2']):
        args_dict['seq_file_1'] = file_1
        args_dict['seq_file_2'] = file_2
        pairSeq(**args_dict)
