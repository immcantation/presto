"""
File I/O and logging functions
"""
# Future
from __future__ import division, absolute_import, print_function

# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import math
import os
import re
import sys
from collections import OrderedDict
from time import time, strftime
from Bio import SeqIO
from Bio.Alphabet import IUPAC

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field
from presto.Annotation import parseAnnotation
from presto.Sequence import translateIUPAC


def readPrimerFile(primer_file):
    """
    Processes primer sequences from file

    Arguments:
    primer_file = name of file containing primer sequences

    Returns:
    dictionary of {primer id, primer sequence}
    """

    # Parse primers from .fasta and .regex files
    ext_name = os.path.splitext(primer_file)[1]
    if ext_name == '.fasta':
        with open(primer_file, 'rU') as primer_handle:
            primer_iter = SeqIO.parse(primer_handle, 'fasta', IUPAC.ambiguous_dna)
            primers = {p.description: str(p.seq).upper()
                        for p in primer_iter}
    elif ext_name == '.regex':
        with open(primer_file, 'rU') as primer_handle:
            primer_list = [a.split(':') for a in primer_handle]
            primers = {a[0].strip(): re.sub(r'\[([^\[^\]]+)\]', translateIUPAC, a[1].strip().upper())
                        for a in primer_list}
    else:
        sys.exit('ERROR:  The primer file %s is not a supported type' % ext_name.upper())
        primers = None

    return primers


def getFileType(filename):
    """
    Determines the type of a file by file extension

    Arguments:
    filename = a filename

    Returns:
    a string defining the sequence type for SeqIO operations
    """
    # Read and check file
    try:
        file_type = os.path.splitext(filename)[1].lower().lstrip('.')
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % filename)
    except Exception as e:
        sys.exit('ERROR:  File %s is invalid with exception %s' % (filename, e))
    else:
        if file_type not in ['fasta', 'fastq', 'embl', 'gb', 'tab']:
            sys.exit('ERROR:  File %s has an unrecognized type' % filename)

    return file_type


def readSeqFile(seq_file, index=False, key_func=None):
    """
    Reads FASTA/FASTQ files

    Arguments:
    seq_file = a FASTA or FASTQ file containing sample sequences
    index = if True return a dictionary from SeqIO.index();
            if False return an iterator from SeqIO.parse()
    key_func = the key_function argument to pass to SeqIO.index if
               index=True

    Returns:
    a tuple of (input file type, sequence record object)
    """
    # Read and check file
    try:
        seq_type = getFileType(seq_file)
        if index:
            seq_records = SeqIO.index(seq_file, seq_type,
                                      alphabet=IUPAC.ambiguous_dna,
                                      key_function=key_func)
        else:
            seq_records = SeqIO.parse(seq_file, seq_type,
                                      alphabet=IUPAC.ambiguous_dna)
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % seq_file)
    except Exception as e:
        sys.exit('ERROR:  File %s is invalid with exception %s' % (seq_file, e))

    return seq_records


def countSeqFile(seq_file):
    """
    Counts the records in FASTA/FASTQ files

    Arguments:
    seq_file = a FASTA or FASTQ file containing sample sequences

    Returns:
    the count of records in the sequence file
    """
    # Count records and check file
    try:
        result_count = len(readSeqFile(seq_file, index=True))
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % seq_file)
    except Exception as e:
        sys.exit('ERROR:  File %s is invalid with exception %s' % (seq_file, e))
    else:
        if result_count == 0:  sys.exit('ERROR:  File %s is empty' % seq_file)

    return result_count


def countSeqSets(seq_file, field=default_barcode_field, delimiter=default_delimiter):
    """
    Identifies sets of sequences with the same ID field

    Arguments:
    seq_file = a FASTA or FASTQ file containing sample sequences
    field = the annotation field containing set IDs
    delimiter = a tuple of delimiters for (fields, values, value lists)

    Returns:
    the count of unit set IDs in the sequence file
    """
    # Count records and check file
    try:
        id_set = set()
        for seq in readSeqFile(seq_file):
            id_set.add(parseAnnotation(seq.description, delimiter=delimiter)[field])
        result_count = len(id_set)
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % seq_file)
    except Exception as e:
        sys.exit('ERROR:  File %s is invalid with exception %s' % (seq_file, e))
    else:
        if result_count == 0:  sys.exit('ERROR:  File %s is empty' % seq_file)

    return result_count


def getOutputHandle(in_file, out_label=None, out_dir=None, out_name=None, out_type=None):
    """
    Opens an output file handle

    Arguments:
    in_file = the input filename
    out_label = text to be inserted before the file extension;
                if None do not add a label
    out_type = the file extension of the output file;
               if None use input file extension
    out_dir = the output directory;
              if None use directory of input file
    out_name = the short filename to use for the output file;
               if None use input file short name

    Returns:
    a file handle
    """
    # Get in_file components
    dir_name, file_name = os.path.split(in_file)
    short_name, ext_name = os.path.splitext(file_name)

    # Define output directory
    if out_dir is None:
        out_dir = dir_name
    else:
        out_dir = os.path.abspath(out_dir)
        if not os.path.exists(out_dir):  os.mkdir(out_dir)
    # Define output file prefix
    if out_name is None:  out_name = short_name
    # Define output file extension
    if out_type is None:  out_type = ext_name.lstrip('.')

    # Define output file name
    if out_label is None:
        out_file = os.path.join(out_dir, '%s.%s' % (out_name, out_type))
    else:
        out_file = os.path.join(out_dir, '%s_%s.%s' % (out_name, out_label, out_type))

    # Open and return handle
    try:
        return open(out_file, 'wb')
    except:
        #raise
        sys.exit('ERROR:  File %s cannot be opened' % out_file)


def printLog(record, handle=sys.stdout, inset=None):
    """
    Formats a dictionary into an IgPipeline log string

    Arguments:
    record = a dict or OrderedDict of {field names: values}
    handle = the file handle to write the log to;
             if None do not write to file
    inset = minimum field name inset;
            if None automatically space field names

    Returns:
    a formatted multi-line string in IgPipeline log format
    """
    # Return empty string if empty dictionary passed
    if not record:
        return ''

    # Determine inset
    max_len = max(map(len, record))
    inset = max(max_len, inset)

    # Assemble log string
    record_str = ''
    if isinstance(record, OrderedDict):
        key_list = record.keys()
    else:
        key_list = sorted(record)
    for key in key_list:
        record_str += '%s%s> %s\n' % (' ' * (inset - len(key)), key, record[key])

    # Write log record
    if handle is not None:
        try:
            handle.write('%s\n' % record_str)
        except IOError as e:
            sys.stderr.write('I/O error writing to log file: %s\n' % e.strerror)

    return record_str


def printMessage(message, start_time=None, end=False, width=20):
    """
    Prints a progress message to standard out

    Arguments:
    message = the current task message
    start_time = task start time returned by time.time();
                 if None do not add run time to progress
    end = if True print final message (add newline)
    width = maximum number of characters for messages

    Returns:
    None
    """
    # Define progress bar
    bar = 'PROGRESS> %s [%s]' % (strftime('%H:%M:%S'), message.ljust(width))

    # Add run time to bar if start_time is specified
    if start_time is not None:
        bar = '%s %.1f min' % (bar, (time() - start_time)/60)

    # Print progress bar
    if end:
        print('\r%s\n' % bar)
        sys.stdout.flush()
    else:
        print('\r%s' % bar, end='')
        sys.stdout.flush()

    return None


# TODO:  probably better to split this into printProgress and printCount (like printMessage)
# TODO:  might be worth adding task argument for the name of the current operation.
def printProgress(current, total=None, step=None, start_time=None, end=False):
    """
    Prints a progress bar to standard out

    Arguments:
    current = the count of completed tasks
    total = the total task count;
            if None do not print percentage
    step = a float defining the fractional progress increment to print if total is defined;
           an int defining the progress increment to print at if total is not defined;
           if None always output the progress
    start_time = task start time returned by time.time();
                 if None do not add run time to progress
    end = if True print final log (add newline)

    Returns:
    None
    """
    try:
        # Check update condition
        if total is None:
            update = (current % step == 0)
        else:
            update = (current % math.ceil(step*total) == 0)
    except:
        # Return on modulo by zero error
        return None
    else:
        # Check end condition and return if no update needed
        if current == total:
            end = True
        if not end and not update:
            return None

    # Define progress bar
    if total is not None and total != 0:
        p = float(current) / total
        c = format(current, "%i,d" % len(format(total, ",d")))
        bar = 'PROGRESS> %s [%-20s] %3.0f%% (%s)' \
              % (strftime('%H:%M:%S'), '#' * int(p*20), p*100, c)
    else:
        bar = 'PROGRESS> %s (%s)' % (strftime('%H:%M:%S'), current)

    # Add run time to bar if start_time is specified
    if start_time is not None:
        bar = '%s %.1f min' % (bar, (time() - start_time)/60)

    # Print progress bar
    if current == 0:
        print('%s' % bar, end='')
        sys.stdout.flush()
    elif end:
        print('\r%s\n' % bar)
        sys.stdout.flush()
    else:
        print('\r%s' % bar, end='')
        sys.stdout.flush()

    return None
