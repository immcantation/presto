#!/usr/bin/env python
"""
Parses records in the console log of pRESTO modules
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2013.9.27'

# Imports
import csv, os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from time import time

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_out_args
from IgCore import getCommonParser, parseCommonArgs
from IgCore import getOutputHandle, printLog, printProgress


def parseLogRecord(record):
    """
    Parses an pRESTO log record

    Arguments: 
    record = a string of lines representing a log record
                    
    Returns: 
    an OrderedDict of {field names: values}
    """
    record_dict = OrderedDict()
    for line in record.split('\n'):
        line = line.strip().split('> ')
        if len(line) == 2:
            record_dict[line[0]] = line[1]
    
    return record_dict


def tableLog(record_file, fields, out_args=default_out_args):
    """
    Converts a pRESTO log to a table of annotations

    Arguments: 
    record_file = the log file name
    fields = the list of fields to output
    out_args = common output argument dictionary from parseCommonArgs
                    
    Returns: 
    the output table file name
    """
    log = OrderedDict()
    log['START'] = 'ParseLog'
    log['FILE'] = os.path.basename(record_file)
    printLog(log)
    
    # Open file handles
    log_handle = open(record_file, 'rU')
    out_handle = getOutputHandle(record_file, 'table', out_dir=out_args['out_dir'], 
                                 out_name=out_args['out_name'], out_type='tab')
    
    # Open csv writer and write header
    out_writer = csv.DictWriter(out_handle, extrasaction='ignore', restval='', 
                                delimiter='\t', fieldnames=fields)
    out_writer.writeheader()
    
    # Iterate over log records
    start_time = time()
    record = ''
    rec_count = 0
    for line in log_handle:
        if line.strip() == '' and record:
            # At end of record write to csv and reset record
            record_dict = parseLogRecord(record)
            record = ''
            if record_dict and any([f in fields for f in record_dict]):
                # Print progress for previous iteration
                printProgress(rec_count, None, 1e5, start_time)
                # Update count and write record
                rec_count += 1
                out_writer.writerow(record_dict)
        else:
            # Append to record
            record += line
    else:
        # Write final record
        if record: 
            record_dict = parseLogRecord(record)
            if record_dict:
                out_writer.writerow(record_dict)                
                rec_count += 1
    
    # Print counts
    printProgress(rec_count, None, 1e5, start_time, end=True)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(out_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ParseLog'
    printLog(log)

    # Close file handles
    log_handle.close()
    out_handle.close()
 
    return log_handle.name


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
                            parents=[getCommonParser(seq_in=False, seq_out=False, log=False)], 
                            formatter_class=ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-l', nargs='+', action='store', dest='record_file', required=True,
                        help='List of log files to parse')
    parser.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                        help='List of fields to collect')
    
    return parser

    
if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args, 'record_file')
    # Convert case of fields
    if args_dict['fields']:  args_dict['fields'] = map(str.upper, args_dict['fields']) 
    
    # Call parseLog for each log file
    for k in args_dict:
        print k, args_dict[k]
    del args_dict['record_file']
    for f in args.__dict__['record_file']:
        args_dict['record_file'] = f
        tableLog(**args_dict)

        # Profiling
        #import cProfile, pstats
        #cProfile.run('parseLog(**args_dict)', 'profile.prof')
        #p = pstats.Stats('profile.prof')
        #p.strip_dirs().sort_stats('time').print_stats() 