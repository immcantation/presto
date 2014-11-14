#!/usr/bin/env python
"""
Core functions shared by pRESTO modules
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.5'
__date__      = '2014.10.2'

# Imports
import ctypes, math, os, re, signal, sys
import multiprocessing as mp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter
from itertools import izip, izip_longest, product
from collections import OrderedDict
from time import time, strftime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


# Defaults
default_delimiter = ('|', '=', ',')
default_separator = default_delimiter[2]
default_coord_choices = ['illumina', 'solexa', 'sra', '454', 'presto']
default_action_choices = ['min', 'max', 'sum', 'first', 'last', 'set']
default_coord_type = 'presto'
default_barcode_field = 'BARCODE'
default_primer_field = 'PRIMER'
default_missing_chars = ['-', '.', 'N']
default_min_freq = 0.6
default_min_qual = 20
default_out_args = {'log_file':None, 
                    'delimiter':default_delimiter,
                    'separator':default_separator,
                    'out_dir':None,
                    'out_name':None,
                    'out_type':None,
                    'clean':False}

# Constants
TERMINATION_SENTINEL = None
EXCEPTION_SENTINEL = None


class CommonHelpFormatter(RawDescriptionHelpFormatter, ArgumentDefaultsHelpFormatter):
    """
    Custom argparse.HelpFormatter
    """
    pass
    

class SeqData:
    """
    A class defining sequence data objects for worker processes
    """
    # Instantiation
    def __init__(self, key, records):
        self.id = key
        self.data = records
        self.valid = (key is not None and records is not None)
        
    # Set boolean evaluation to valid value
    def __nonzero__(self): 
        return self.valid
    
    # Set length evaluation to number of data records
    def __len__(self):
        if isinstance(self.data, SeqRecord) or isinstance(self.data, Seq):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


class SeqResult:
    """
    A class defining sequence result objects for collector processes
    """
    # Instantiation
    def __init__(self, key, records):
        self.id = key
        self.data = records
        self.results = None
        self.valid = False
        self.log = OrderedDict([('ID', key)])
    
    # Set boolean evaluation to valid value
    def __nonzero__(self): 
        return self.valid

    # Set length evaluation to number of results
    def __len__(self):
        if isinstance(self.results, SeqRecord) or isinstance(self.results, Seq):
            return 1
        elif self.results is None:
            return 0
        else:
            return len(self.results)
    
    # Set data_count to number of data records
    @property
    def data_count(self):
        if isinstance(self.data, SeqRecord) or isinstance(self.data, Seq):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


def parseAnnotation(record, fields=None, delimiter=default_delimiter):
    """
    Extracts annotations from a FASTA/FASTQ sequence description

    Arguments:
    record = description string to extract annotations from
    fields = a list of fields to subset the return dictionary to
             if None return all fields
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a OrderedDict of {field:value}
    """
    annotation = record.split(delimiter[0])
    field_dict = OrderedDict([('ID', annotation.pop(0))])
    for a in annotation:
        vals = a.split(delimiter[1])
        field_dict[vals[0].upper()] = vals[1] 

    # Subset return dictionary to requested fields
    if fields is not None:
        if not isinstance(fields, list):  fields = [fields]
        for f in set(field_dict).difference(fields):  del field_dict[f]
    
    return field_dict


def flattenAnnotation(ann_dict, delimiter=default_delimiter):
    """
    Converts annotations from a dictionary to a FASTA/FASTQ sequence description

    Arguments: 
    ann_dict = a dictionary of field/value pairs
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a formatted sequence description string
    """
    annotation = ann_dict.get('ID', 'NONE')
    for k, v in ann_dict.iteritems():
        if k.upper() == 'ID':  continue
        if isinstance(v, list):
            v = delimiter[2].join([str(x) for x in v])
        annotation += '%s%s%s%s' % (delimiter[0], k.upper(), delimiter[1], v) 

    return annotation


def mergeAnnotation(ann_dict_1, ann_dict_2, delimiter=default_delimiter):
    """
    Merges non-ID field annotations from one field dictionary into another

    Arguments:
    ann_dict_1 = a dictionary of field/value pairs to add to
    ann_dict_2 = a dictionary of field/value pairs to merge with fields_2
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a modified fields_1 dictonary of {field:value}
    """
    for k, v in ann_dict_2.iteritems():
        if k.upper() == 'ID':  continue
        if k in ann_dict_1:  
            if isinstance(ann_dict_1[k], list):
                ann_dict_1[k] = delimiter[2].join([str(x) for x in ann_dict_1[k]])
            if isinstance(v, list):
                v = delimiter[2].join([str(x) for x in v])
            ann_dict_1[k] += '%s%s' % (delimiter[2], v)
        else:  
            ann_dict_1[k.upper()] = v

    return ann_dict_1


def renameAnnotation(ann_dict, old_field, new_field, delimiter=default_delimiter):
    """
    Renames an annotation and merges annotations if the new name already exists

    Arguments:
    ann_dict = a dictionary of field/value pairs
    old_field = the old field name
    new_field = the new field name
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a modified fields dictonary of {field:value}
    """
    if new_field in ann_dict:
        new_ann = ann_dict.copy()
        del new_ann[old_field]
        mergeAnnotation(new_ann, {new_field:ann_dict[old_field]}, delimiter=delimiter)
    else:
        new_ann = OrderedDict([(new_field, v) if k == old_field else (k, v) \
                               for k, v in ann_dict.iteritems()])

    return new_ann


def collapseAnnotation(ann_dict, action, fields=None, delimiter=default_delimiter):
    """
    Convert annotations from a dictionary to a FASTA/FASTQ sequence description

    Arguments: 
    ann_dict = a dictionary of field/value pairs
    action = the _collapse action to take;
             one of {min, max, sum, first, last, set}
    fields = a subset of ann_dict to _collapse;
             if None _collapse all but the ID field
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a modified field dictionary of {field:value}
    """
    # Define _collapse action            
    if action == 'set':
        def _collapse(value):  return list(set(value))
    elif action == 'first':
        def _collapse(value):  return value[0]
    elif action == 'last':
        def _collapse(value):  return value[-1]
    elif action == 'min':
        def _collapse(value):  return '%.12g' % min([float(x or 0) for x in value])
    elif action == 'max':
        def _collapse(value):  return '%.12g' % max([float(x or 0) for x in value])
    elif action == 'sum':
        def _collapse(value):  return '%.12g' % sum([float(x or 0) for x in value])
    else:
        def _collapse(value):  return value

    # Collapse fields        
    for k, v in ann_dict.iteritems():
        if k.upper() == 'ID':
            continue
        if fields is None or k in fields:
            # Convert field to list
            if not isinstance(v, list) and isinstance(v, str):
                v = v.split(delimiter[2])
            elif not isinstance(v, list):
                v = [v]
            # Perform _collapse and reassign field
            ann_dict[k] = _collapse(v)       

    return ann_dict


def getAnnotationValues(seq_iter, field, unique=False, delimiter=default_delimiter):
    """
    Gets the set of unique annotation values in a sequence set

    Arguments: 
    seq_iter = an iterator or list of SeqRecord objects
    field = the annotation field to retrieve values for
    unique = if True return a list of only the unique values;
             if Falser return a list of all values
    delimiter = a tuple of delimiters for (fields, values, value lists) 

    Returns: 
    the list of values for the field
    """
    # Parse annotations from seq_list records
    ann_iter = (parseAnnotation(s.description, delimiter=delimiter) for s in seq_iter)
    values = [a[field] for a in ann_iter]

    return list(set(values)) if unique else values


def annotationConsensus(seq_iter, field, delimiter=default_delimiter):
    """
    Calculate a consensus annotation for a set of sequences
    
    Arguments:
    seq_iter = an iterator or list of SeqRecord objects
    field = the annotation field to take a consensus of
    delimiter = a tuple of delimiters for (annotations, field/values, value lists) 
    
    Returns:
    a dictionary {set:list of annotation values, count:annotation counts,
                  cons:consensus annotation,freq:majority annotation frequency)
    """
    # Define return dictionary
    cons_dict = {'set':None, 'count':None, 'cons':None, 'freq':None}
    
    # Parse annotations from seq_list records
    val_list = getAnnotationValues(seq_iter, field, delimiter=delimiter)
    
    # Define annotation set and counts
    cons_dict['set'] = sorted(set(val_list))
    cons_dict['count'] = [val_list.count(v) for v in cons_dict['set']]

    # Define consensus annotation
    i = cons_dict['count'].index(max(cons_dict['count']))
    cons_dict['cons'] = cons_dict['set'][i]
    cons_dict['freq'] = float(cons_dict['count'][i]) / len(val_list)
                
    return cons_dict


def subsetSeqSet(seq_iter, field, values, delimiter=default_delimiter):
    """
    Subsets a sequence set by annotation value
    
    Arguments:
    seq_iter = an iterator or list of SeqRecord objects
    field = the annotation field to select by
    values = a list of annotation values that define the retained sequences
    delimiter = a tuple of delimiters for (annotations, field/values, value lists) 
    
    Returns:
    a modified list of SeqRecord objects
    """
    # Parse annotations from seq_list records
    ann_list = [parseAnnotation(s.description, delimiter=delimiter) for s in seq_iter]
    
    # Subset seq_list by annotation
    if not isinstance(values, list):  values = [values]
    seq_subset = [seq_iter[i] for i, a in enumerate(ann_list) if a[field] in values]

    return seq_subset


def subsetSeqIndex(seq_dict, field, values, delimiter=default_delimiter):
    """
    Subsets a sequence set by annotation value
    
    Arguments:
    seq_dict = a dictionary index of sequences returned from SeqIO.index()
    field = the annotation field to select keys by
    values = a list of annotation values that define the retained keys
    delimiter = a tuple of delimiters for (annotations, field/values, value lists) 
    
    Returns:
    a list of keys
    """
    # Parse annotations from seq_dict and subset keys
    key_subset = [k for k in seq_dict \
                  if parseAnnotation(seq_dict[k].description, delimiter=delimiter)[field] \
                  in values]

    return key_subset


def compilePrimers(primers):
    """
    Translates IUPAC Ambiguous Nucleotide characters to regular expressions and compiles them

    Arguments: 
    key = a dictionary of sequences to translate

    Returns:
    dictionary of compiled regular expressions
    """
    
    primers_regex = {k: re.compile(re.sub(r'([RYSWKMBDHVN])', translateIUPAC, v)) 
                     for k, v in primers.iteritems()}  
    
    return primers_regex


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


def readSeqFile(seq_file, index=False):
    """
    Reads FASTA/FASTQ files

    Arguments: 
    seq_file = a FASTA or FASTQ file containing sample sequences
    index = if True return a dictionary from SeqIO.index();
            if False return an iterator from SeqIO.parse()
    Returns: 
    a tuple of (input file type, sequence record object)
    """
    # Read and check file
    try:
        seq_type = getFileType(seq_file)
        if index:  
            seq_records = SeqIO.index(seq_file, seq_type, IUPAC.ambiguous_dna)
        else:  
            seq_records = SeqIO.parse(seq_file, seq_type, IUPAC.ambiguous_dna)
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


def translateIUPAC(key):
    """
    Translates IUPAC Ambiguous Nucleotide characters to or from character sets

    Arguments: 
    key = a string or re.search object containing the character set to translate

    Returns:
    character translation
    """
    
    # Define valid characters and character translations
    IUPAC_uniq = '-.ACGT'
    IUPAC_ambig = 'BDHKMNRSVWY'
    IUPAC_trans = {'AG':'R', 'CT':'Y', 'CG':'S', 'AT':'W', 'GT':'K', 'AC':'M',
                  'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACG':'V', 'ABCDGHKMRSTVWY':'N',
                  '-.':'.'}

    # Convert passed regular expression match to a string
    if hasattr(key, 'group'): key = key.group(1)
    
    # Sort character in alphabetic order
    key = ''.join(sorted(key))
    
    # Return input character if no translation needed
    if len(key) == 1 and key in IUPAC_uniq:
        return key
    # Return regular expression string for ambiguous single character
    elif len(key) == 1 and key in IUPAC_ambig:
        return ['[' + k + ']' for k, v in IUPAC_trans.iteritems() if v == key][0]
    # Return single ambiguous character for character set 
    elif key in IUPAC_trans:
        return IUPAC_trans[key]
    else:
        return 'N'


def scoreDNA(a, b, n_score=None, gap_score=None):
    """
    Returns the score for a pair of IUPAC Ambiguous Nucleotide characters

    Arguments: 
    a = first characters
    b = second character
    n_score = score for all matches against an N character;
              if None score according to IUPAC character identity
    gap_score = score for all matches against a [-, .] character;
                if None score according to IUPAC character identity

    Returns:
    score for the character pair
    """
    # Define ambiguous character translations    
    IUPAC_trans = {'AG':'R', 'CT':'Y', 'CG':'S', 'AT':'W', 'GT':'K', 'AC':'M',
                   'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACG':'V', 'ABCDGHKMRSTVWY':'N',
                   '-.':'.'}
    # Create list of tuples of synonymous character pairs
    IUPAC_matches = [p for k, v in IUPAC_trans.iteritems() for p in list(product(k, v))]

    # Check gap condition
    if gap_score is not None and (a in '-.' or b in '-.'):
        return gap_score

    # Check N-value condition
    if n_score is not None and (a == 'N' or b == 'N'):
        return n_score

    # Determine and return score for IUPAC match conditions
    # Symmetric and reflexive
    if a == b:
        return 1
    elif (a, b) in IUPAC_matches:
        return 1
    elif (b, a) in IUPAC_matches:
        return 1
    else:
        return 0


def scoreAA(a, b, x_score=None, gap_score=None):
    """
    Returns the score for a pair of IUPAC Extended Protein characters

    Arguments: 
    a = first character
    b = second character
    x_score = score for all matches against an X character;
              if None score according to IUPAC character identity
    gap_score = score for all matches against a [-, .] character;
                if None score according to IUPAC character identity

    Returns:
    score for the character pair
    """
    # Define ambiguous character translations    
    IUPAC_trans = {'RN':'B', 'EQ':'Z', 'LI':'J', 'ABCDEFGHIJKLMNOPQRSTUVWYZ':'X',
                   '-.':'.'}
    # Create list of tuples of synonymous character pairs
    IUPAC_matches = [p for k, v in IUPAC_trans.iteritems() for p in list(product(k, v))]

    # Check gap condition
    if gap_score is not None and (a in '-.' or b in '-.'):
        return gap_score

    # Check X-value condition
    if x_score is not None and (a == 'X' or b == 'X'):
        return x_score

    # Determine and return score for IUPAC match conditions
    # Symmetric and reflexive
    if a == b:
        return 1
    elif (a, b) in IUPAC_matches:
        return 1
    elif (b, a) in IUPAC_matches:
        return 1
    else:
        return 0


def getScoreDict(n_score=None, gap_score=None, alphabet='dna'):
    """
    Generates a score dictionary

    Arguments: 
    n_score = score for all matches against an N character;
              if None score according to IUPAC character identity
    gap_score = score for all matches against a [-, .] character;
                if None score according to IUPAC character identity
    alphabet = the type of score dictionary to generate;
               one of [dna, aa] for DNA and amino acid characters
                
    Returns:
    a score dictionary of the form {(char1, char2) : score}
    """
    if alphabet=='dna':
        IUPAC_chars = '-.ACGTRYSWKMBDHVN'
        IUPAC_dict = {k:scoreDNA(*k, n_score=n_score, gap_score=gap_score) 
                      for k in product(IUPAC_chars, repeat=2)}
    elif alphabet=='aa':
        IUPAC_chars = '-.*ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        IUPAC_dict = {k:scoreAA(*k, x_score=n_score, gap_score=gap_score) 
                      for k in product(IUPAC_chars, repeat=2)}
    else:
        sys.stderr.write('ERROR:  The alphabet %s is not a recognized type\n' % alphabet)

    return IUPAC_dict


def reverseComplement(seq):
    """
    Takes the reverse complement of a sequence

    Arguments: 
    seq = a SeqRecord object, Seq object or string to reverse complement
    
    Returns:
    a object of the same type as the input with the reverse complement sequence 
    """
    
    if isinstance(seq, SeqRecord):
        new_record = seq.reverse_complement(id=True, name=True, description=True, 
                                            features=True, annotations=True, 
                                            letter_annotations=True)
        #new_record.annotations['orientation'] = 'RC'
    elif isinstance(seq, Seq):
        new_record = seq.reverse_complement()
    elif isinstance(seq, str):
        new_record = str(Seq(seq, IUPAC.ambiguous_dna).reverse_complement())
    else:
        sys.exit('ERROR:  Invalid record type passed to reverseComplement')
        new_record = None
    
    return new_record


def testSeqEqual(seq1, seq2, ignore_chars=default_missing_chars):
    """
    Determine if two sequences are equal, excluding missing positions
    
    Arguments: 
    seq1 = a SeqRecord object
    seq2 = a SeqRecord object
    ignore_chars = a list of characters to ignore
    
    Returns:
    True if the sequences are equal
    """
    equal = True
    for a, b in izip(seq1.upper(), seq2.upper()):
        if a != b and a not in ignore_chars and b not in ignore_chars:
            equal = False
            break
    
    return equal
 

def weightSeq(seq):
    """
    Returns a score for a single sequence
    
    Arguments: 
    seq = a SeqRecord object
    
    Returns:
    The sum of the character scores for the sequence
    """
    nuc_score = sum([c in 'ACGTRYSWKMBDHV' for c in seq.upper()])
    gap_score = 0
      
    return max(nuc_score + gap_score, 1)


def scoreSeqPair(seq1, seq2, max_error=None, max_weight=None, 
                 score_dict=getScoreDict(n_score=0, gap_score=0)):
    """
    Determine the error rate for a pair of sequences
    
    Arguments: 
    seq1 = a SeqRecord object
    seq2 = a SeqRecord object
    max_error = the maximum error rate; once reached return (0, 0, 1.0)
                if None always return accurate score, weight and error
    max_weight = the maximum weight to use when checking the max_error break condition;
                 if None use the minimum length of seq1,seq2
    score_dict = optional dictionary of alignment scores as {(char1, char2): score}
    
    Returns:
    A tuple of the (score, minimum weight, error rate) for the pair of sequences
    """
    # Determine score
    if max_error is None:
        # Return accurate values when max_error is undefined
        chars = izip(seq1.upper(), seq2.upper())
        score = sum([score_dict[(a, b)] for a, b in chars])
        weight = min(weightSeq(seq1), weightSeq(seq2))
        error = 1.0 - float(score) / weight
    else:
        # If max_error defined return when worst case reach
        score = 0
        if not max_weight:  max_weight = min(len(seq1), len(seq2))
        for i, (a, b) in enumerate(izip(seq1, seq2)):
            score += score_dict[(a, b)]
            if (i - float(score)) / max_weight > max_error:
                score, weight, error = 0, 0, 1.0
                break
        else:
            weight = min(weightSeq(seq1), weightSeq(seq2))
            error = 1.0 - float(score) / weight

    return (score, weight, error)


def calculateDiversity(seq_list, score_dict=getScoreDict(n_score=0, gap_score=0)):
    """
    Determine the average pairwise error rate for a list of sequences

    Arguments: 
    seq_list = a list of SeqRecord objects to score
    score_dict = optional dictionary of alignment scores as {(char1, char2): score}
        
    Returns:
    The average pairwise error rate for the list of sequences
    """
    # Return 0 if less than 2 sequences
    if len(seq_list) <= 1:
        return 0
    
    scores = []
    for i, seq1 in enumerate(seq_list):
        for seq2 in seq_list[(i + 1):]:
            scores.append(scoreSeqPair(seq1, seq2, score_dict=score_dict)[2])
    
    return sum(scores) / len(scores)

# TODO: apply min_qual to single sequence case
def qualityConsensus(seq_list, min_qual=default_min_qual, 
                     min_freq=default_min_freq, max_miss=None,
                     dependent=False):
    """
    Builds a consensus sequence from a set of sequences

    Arguments: 
    seq_list = a list of SeqRecord objects
    min_qual = the quality cutoff to assign a base
    min_freq = the frequency cutoff to assign a base
    max_miss = the maximum frequency of (., -, N) characters allowed before 
               deleting a position; if None do not delete positions 
    dependent = if False assume sequences are independent for quality calculation
    
    Returns: 
    a consensus SeqRecord object
    """
    # Return a copy of the input SeqRecord upon singleton
    if len(seq_list) == 1:
        return seq_list[0].upper()
    
    # Define characters to be ignored in consensus building
    ignore_set = set(['.', '-', 'N'])

    # Create sequence and annotation iterators
    # Pad unequal length sequences with character '-' and quality 0 
    seq_str = [str(s.seq) for s in seq_list]
    seq_iter = izip_longest(*seq_str, fillvalue='-')
    ann_list = [s.letter_annotations['phred_quality'] for s in seq_list]
    ann_iter = izip_longest(*ann_list, fillvalue=0)
    
    # Build consensus
    consensus_seq = []
    consensus_qual = []
    for chars, quals in izip(seq_iter, ann_iter):
        # Remove position if max_gap is defined and positions contains too many gaps
        if max_miss is not None:
            gap_count = sum([chars.count(c) for c in ignore_set])
            gap_freq = float(gap_count) / len(chars)
            if gap_freq > max_miss:  continue
            
        # Define set of non-missing characters
        char_set = set(chars).difference(ignore_set)
        
        # Assign N if no non-N/gap characters and proceed to next position
        if not char_set:
            consensus_seq.append('N')
            consensus_qual.append(0)
            continue
        
        # Define non-missing character frequencies
        char_count = float(len([c for c in chars if c in char_set]))
        char_freq = {c: chars.count(c) / char_count for c in char_set}

        # Create per character quality sets and quality sums 
        qual_total = float(sum(quals))
        qual_set, qual_sum = {}, {}
        for c in char_set:
            qual_set[c] = [q for i, q in enumerate(quals) if chars[i] == c]
            qual_sum[c] = sum(qual_set[c])
        
        # Calculate per character consensus quality scores
        if dependent:
            qual_cons = {c:int(max(qual_set[c]) * qual_sum[c] / qual_total) for c in qual_set}
        else:
            qual_cons = {c:int(qual_sum[c] * qual_sum[c] / qual_total) for c in qual_set}
            
        # Select character with highest consensus quality
        cons = [(c, min(q, 93)) for c, q in qual_cons.iteritems() \
                if q == max(qual_cons.values())][0]
        # Assign N if consensus quality or frequency threshold is failed
        if cons[1] < min_qual or char_freq[cons[0]] < min_freq:  
            cons = ('N', 0)
        
        # Assign consensus base and quality
        consensus_seq.append(cons[0])
        consensus_qual.append(cons[1])
    
    # Define return SeqRecord
    record = SeqRecord(Seq(''.join(consensus_seq), IUPAC.ambiguous_dna),
                       id='consensus',
                       name='consensus',
                       description='',
                       letter_annotations={'phred_quality':consensus_qual})

    return record


def frequencyConsensus(seq_list, min_freq=default_min_freq, max_miss=None):
    """
    Builds a consensus sequence from a set of sequences

    Arguments: 
    set_seq = a list of SeqRecord objects
    min_freq = the frequency cutoff to assign a base
    max_miss = the maximum frequency of (., -, N) characters allowed before 
               deleting a position; if None do not delete positions 
    
    Returns: 
    a consensus SeqRecord object
    """
    # Return a copy of the input SeqRecord upon singleton
    if len(seq_list) == 1:
        return seq_list[0].upper()
    
    # Define gap characters
    ignore_set = set(['.', '-', 'N'])
    
    # Build consensus
    seq_str = [str(s.seq) for s in seq_list]
    consensus_seq = []
    for chars in izip_longest(*seq_str, fillvalue='-'):
        # Delete position if number of gap characters exceeds max_gap threshold
        if max_miss is not None:
            gap_count = sum([chars.count(c) for c in ignore_set])
            gap_freq = float(gap_count) / len(chars)
            if gap_freq > max_miss:  continue
            
        # Define set of non-missing characters
        char_set = set(chars).difference(ignore_set)
        
        # Assign N if no non-N/gap characters and proceed to next position
        if not char_set:
            consensus_seq.append('N')
            continue
        
        # Define non-missing character frequencies
        char_count = float(len([c for c in chars if c in char_set]))
        char_freq = {c: chars.count(c) / char_count for c in char_set}
        freq_max = max(char_freq.values())

        # Assign consensus as most frequent character
        cons = [c if char_freq[c] >= min_freq else 'N' \
                for c in char_set if char_freq[c] == freq_max][0]
        consensus_seq.append(cons)

    # Define return SeqRecord
    record = SeqRecord(Seq(''.join(consensus_seq), IUPAC.ambiguous_dna), 
                       id='consensus',
                       name='consensus',
                       description='')
        
    return record


def indexSeqSets(seq_dict, field=default_barcode_field, delimiter=default_delimiter):
    """
    Identifies sets of sequences with the same ID field

    Arguments: 
    seq_dict = a dictionary index of sequences returned from SeqIO.index()
    field = the annotation field containing set IDs
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a dictionary of {set name:[record names]}
    """
    set_dict = {}
    for key, rec in seq_dict.iteritems():
        tag = parseAnnotation(rec.description, delimiter=delimiter)[field]
        set_dict.setdefault(tag, []).append(key)

    return set_dict


def indexSeqPairs(seq_dict_1, seq_dict_2, coord_type=default_coord_type, 
                  delimiter=default_delimiter):
    """
    Identifies sequence pairs by coordination information

    Arguments: 
    seq_dict_1 = a dictionary index of sequences returned from SeqIO.index()
    seq_dict_2 = a dictionary index of sequences returned from SeqIO.index()
    coord_type = the sequence header format
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a dictionary of {coordinate:(seq_dict_1 key, seq_dict_2 key)}
    """
    # Get coordinate dictionaries
    coord_1 = getSeqCoord(seq_dict_1, coord_type, delimiter=delimiter)
    coord_2 = getSeqCoord(seq_dict_2, coord_type, delimiter=delimiter)
    
    # Form sets of keys
    set_1 = set(coord_1)
    set_2 = set(coord_2)
    
    # Find matching entries in key sets
    match_set = set_1.intersection(set_2)
    index_dict = {n:(coord_1[n], coord_2[n]) for n in match_set}
    
    return index_dict


def getUnpairedIndex(seq_dict_1, seq_dict_2, coord_type=default_coord_type, 
                     delimiter=default_delimiter):
    """
    Identifies sequence without a mate pair by coordination information

    Arguments: 
    seq_dict_1 = a dictionary index of sequences returned from SeqIO.index()
    seq_dict_2 = a dictionary index of sequences returned from SeqIO.index()
    coord_type = the sequence header format
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a tuple of lists containing unique ([seq_dict_1 keys], [seq_dict_2 keys])
    """
    # Get coordinate dictionaries
    coord_1 = getSeqCoord(seq_dict_1, coord_type, delimiter=delimiter)
    coord_2 = getSeqCoord(seq_dict_2, coord_type, delimiter=delimiter)
    
    # Form sets of coordinate keys
    keys_1 = set(coord_1)
    keys_2 = set(coord_2)
    
    # Find unmatched entries in coordinate key sets
    uniq_keys_1 = keys_1.difference(keys_2)
    uniq_keys_2 = keys_2.difference(keys_1)
    
    # Translate unmatched coordinate keys to sequences dictionary keys
    unpaired_1 = [coord_1[n] for n in uniq_keys_1]
    unpaired_2 = [coord_2[n] for n in uniq_keys_2]
    
    return unpaired_1, unpaired_2


def getSeqCoord(seq_dict, coord_type=default_coord_type, delimiter=default_delimiter):
    """
    Create a sequence ID translation using delimiter truncation
    
    Arguments: 
    seq_dict = a sequence dictionary generate by SeqIO.index()
    coord_type = the sequence header format; 
                 one of ['illumina', 'solexa', 'sra', '454', 'presto'];
                 if unrecognized type or None create ID dictionary without translation
    delimiter = a tuple of delimiters for (fields, values, value lists) 
                    
    Returns: 
    a dictionary of {coordinate: full ID}
    """
    # Sequencing identifier formats supported
    # illumina:  @MISEQ:132:000000000-A2F3U:1:1101:14340:1555 2:N:0:ATCACG
    #            @HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1
    # sra:       @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
    # 454:       @000034_0199_0169 length=437 uaccno=GNDG01201ARRCR
    # pesto      @AATCGGATTTGC|COUNT=2|PRIMER=IGHJ_RT|PRFREQ=1.0
     
    # Define _getCoord function by coordinate format
    if coord_type in ('illumina', 'solexa'):
        def _getCoord(s):  return s.split()[0].split('#')[0]
    elif coord_type in ('sra', '454'):
        def _getCoord(s):  return s.split()[0]
    elif coord_type == 'presto':
        def _getCoord(s):  return parseAnnotation(s, delimiter=delimiter)['ID']
    else:
        def _getCoord(s):  return s

    # Create a seq_dict ID translation using IDs truncated by _getCoord()
    ids = {}
    for seq_id in seq_dict:
        try:
            id_key = _getCoord(seq_id)
            ids.update({id_key:seq_id})
        except:
            sys.exit('ERROR:  No ID for sequence %s using coordinate format %s' \
                     % (seq_id, coord_type))
    
    return ids


def manageProcesses(feed_func, work_func, collect_func, 
                    feed_args={}, work_args={}, collect_args={}, 
                    nproc=None, queue_size=None):
    """
    Manages feeder, worker and collector processes
    
    Arguments:
    feed_func = the data Queue feeder function
    work_func = the worker function
    collect_func = the result Queue collector function
    feed_args = a dictionary of arguments to pass to feed_func
    work_args = a dictionary of arguments to pass to work_func
    collect_args = a dictionary of arguments to pass to collect_func
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc    
    
    Returns:
    a dictionary of collector results
    """
    # Define signal handler that raises KeyboardInterrupt
    def _signalHandler(s, f):
        raise SystemExit

    # Define function to terminate child processes
    def _terminate():
        sys.stderr.write('Terminating child processes...')
        # Terminate feeders
        feeder.terminate()
        feeder.join()
        # Terminate workers
        for w in workers:
            w.terminate()
            w.join()
        # Terminate collector
        collector.terminate()
        collector.join
        sys.stderr.write('  Done.\n')
    
    # Raise SystemExit upon termination signal
    signal.signal(signal.SIGTERM, _signalHandler)
        
    # Define number of processes and queue size
    if nproc is None:  nproc = mp.cpu_count()
    if queue_size is None:  queue_size = nproc * 2
    
    # Define shared child process keep alive flag
    alive = mp.Value(ctypes.c_bool, True)
    
    # Define shared data queues
    data_queue = mp.Queue(queue_size)
    result_queue = mp.Queue(queue_size)
    collect_queue = mp.queues.SimpleQueue()
    # Initiate manager and define shared data objects

    try:  
        # Initiate feeder process
        feeder = mp.Process(target=feed_func, 
                            args=(alive, data_queue), 
                            kwargs=feed_args)
        feeder.start()

        # Initiate worker processes
        workers = []
        for __ in range(nproc):
            w = mp.Process(target=work_func, 
                           args=(alive, data_queue, result_queue), 
                           kwargs=work_args) 
            w.start()
            workers.append(w)
    
        # Initiate collector process
        collector = mp.Process(target=collect_func, 
                               args=(alive, result_queue, collect_queue), 
                               kwargs=collect_args)
        collector.start()
    
        # Wait for feeder to finish and add sentinel objects to data_queue
        feeder.join()
        for __ in range(nproc):  data_queue.put(None)
        
        # Wait for worker processes to finish and add sentinel to result_queue
        for w in workers:  w.join()
        result_queue.put(None)
        
        # Wait for collector process to finish and add sentinel to collect_queue
        collector.join()
        collect_queue.put(None)

        # Get collector return values
        collected = collect_queue.get()
    except (KeyboardInterrupt, SystemExit):
        sys.stderr.write('Exit signal received\n')
        _terminate()
        sys.exit()
    except Exception as e:
        sys.stderr.write('ERROR:  %s\n' % e)
        _terminate()
        sys.exit()
    except:
        sys.stderr.write('ERROR:  Exiting with unknown exception\n')
        _terminate()
        sys.exit()
    else:
        if not alive.value:
            sys.stderr.write('ERROR:  Exiting due to child process error\n')
            _terminate()
            sys.exit()
    
    return collected


def feedSeqQueue(alive, data_queue, seq_file, index_func=None, index_args={}):
    """
    Feeds the data queue with SeqRecord objects

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = a multiprocessing.Queue to hold data for processing
    seq_file = the sequence file to read input from
    index_func = the function to use to define sequence sets
                 if None do not index sets and feed individual records
    index_args = a dictionary of arguments to pass to index_func
    
    Returns: 
    None
    """
    try:
        # Read input file and index sequence sets if required
        if index_func is None:
            seq_iter = readSeqFile(seq_file)
            data_iter = ((s.id, s) for s in seq_iter)
        else:
            seq_dict = readSeqFile(seq_file, index=True)
            index_dict = index_func(seq_dict, **index_args)
            data_iter = ((k, [seq_dict[i] for i in v]) \
                         for k, v in index_dict.iteritems()) 
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


def processSeqQueue(alive, data_queue, result_queue, process_func, process_args={}):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    process_func = the function to use for filtering sequences
    process_args = a dictionary of arguments to pass to process_func

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

            # Perform work
            result = process_func(data, **process_args)

            #import cProfile
            #prof = cProfile.Profile()
            #result = prof.runcall(process_func, data, **process_args)
            #prof.dump_stats('worker-%d.prof' % os.getpid())

            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        sys.stderr.write('Error processing sequence with ID: %s.\n' % data.id)
        raise
    
    return None


def collectSeqQueue(alive, result_queue, collect_queue, seq_file, 
                    task_label, out_args, index_field=None):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    result_queue = a multiprocessing.Queue holding worker results
    collect_queue = a multiprocessing.Queue to store collector return values
    seq_file = the sample sequence file name
    task_label = the task label used to tag the output files
    out_args = common output argument dictionary from parseCommonArgs
    index_field = the field defining set membership for sequence sets
                  if None data queue contained individual records
    
    Returns: 
    None
    (adds a dictionary of {log: log object, out_files: output file names} to collect_queue)
    """
    try:
        # Count records 
        if index_field is None:
            result_count = countSeqFile(seq_file)
        else:
            result_count = countSeqSets(seq_file, index_field, out_args['delimiter'])

        # Define output format
        out_type = getFileType(seq_file) if out_args['out_type'] is None \
                   else out_args['out_type']
            
        # Defined valid alignment output handle
        pass_handle = getOutputHandle(seq_file, 
                                      '%s-pass' % task_label, 
                                      out_dir=out_args['out_dir'], 
                                      out_name=out_args['out_name'], 
                                      out_type=out_type)
        # Defined failed alignment output handle
        if out_args['clean']:  
            fail_handle = None
        else:  
            fail_handle = getOutputHandle(seq_file, 
                                          '%s-fail'  % task_label, 
                                          out_dir=out_args['out_dir'], 
                                          out_name=out_args['out_name'], 
                                          out_type=out_type)
        # Define log handle
        if out_args['log_file'] is None:  
            log_handle = None
        else:  
            log_handle = open(out_args['log_file'], 'w')
    except:
        alive.value = False
        raise
    
    try:
        # Iterator over results queue until sentinel object reached
        start_time = time()
        set_count = seq_count = pass_count = fail_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break
            
            # Print progress for previous iteration
            printProgress(set_count, result_count, 0.05, start_time) 
            
            # Update counts for current iteration
            set_count += 1
            seq_count += result.data_count
            
            # Write log
            printLog(result.log, handle=log_handle)
    
            # Write alignments
            if result:
                pass_count += 1
                SeqIO.write(result.results, pass_handle, out_type)
            else:
                fail_count += 1
                if fail_handle is not None:
                    SeqIO.write(result.data, fail_handle, out_type)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
        
        # Print total counts
        printProgress(set_count, result_count, 0.05, start_time)
    
        # Update return values    
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name)
        log['SEQUENCES'] = seq_count
        if index_field is not None:
            log['SETS'] = set_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count
        collect_dict = {'log':log, 'out_files': [pass_handle.name]}
        collect_queue.put(collect_dict)
        
        # Close file handles
        pass_handle.close()
        if fail_handle is not None:  fail_handle.close()
        if log_handle is not None:  log_handle.close()
    except:
        alive.value = False
        raise
    
    return None


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
        print '%s' % bar,
        sys.stdout.flush()
    elif end:
        print '\r%s\n' % bar
        sys.stdout.flush()
    else:
        print '\r%s' % bar,
        sys.stdout.flush()
    
    return None
        
        
def getCommonArgParser(seq_in=True, seq_out=True, paired=False, db_in=False, db_out=False,
                       annotation=True, log=True, multiproc=False):
    """
    Defines an ArgumentParser object with common pRESTO arguments

    Arguments: 
    seq_in = if True include sequence input arguments
    seq_out = if True include sequence output arguments
    paired = if True defined paired-end sequence input and output arguments
    db_in = if True include tab delimited database input arguments
    db_out = if True include tab delimited database output arguments
    annotation = if True include annotation arguments
    log = if True include log arguments
    multiproc = if True include multiprocessing arguments
    
    Returns:
    an ArgumentParser object
    """
    parser = ArgumentParser(add_help=False, formatter_class=CommonHelpFormatter)

    # Database arguments
    if db_in:
        parser.add_argument('-d', nargs='+', action='store', dest='db_files', required=True,
                        help='A list of tab delimited database files.')
    
    if db_out:
        parser.add_argument('--clean', action='store_true', dest='clean', 
                            help='If specified do not create files containing records that \
                                  fail processing.')
            
    # Sequence arguments
    if seq_in and not paired:
        parser.add_argument('-s', nargs='+', action='store', dest='seq_files', required=True,
                            help='A list of FASTA/FASTQ files containing sequences to process.')
    elif seq_in and paired:
        parser.add_argument('-1', nargs='+', action='store', dest='seq_files_1', required=True,
                            help='An ordered list of FASTA/FASTQ files containing head/primary \
                                  sequences.')
        parser.add_argument('-2', nargs='+', action='store', dest='seq_files_2', required=True,
                            help='An ordered list of FASTA/FASTQ files containing tail/secondary \
                                  sequences.')
    if seq_out:
        parser.add_argument('--fasta', action='store_const', dest='out_type', const='fasta',
                            help='Specify to force output as FASTA rather than FASTQ.')
        parser.add_argument('--clean', action='store_true', dest='clean', 
                            help='If specified do not create files containing sequences that \
                                  fail processing.')
            
    # Annotation arguments
    if annotation:
        parser.add_argument('--delim', nargs=3, action='store', dest='delimiter', 
                            type=str, default=default_delimiter, 
                            help='A list of the three delimiters that separate annotation blocks, \
                                  field names and values, and values within a field, respectively.')

    # Log arguments
    if log:
        parser.add_argument('--log', action='store', dest='log_file', default=None,
                            help='Specify to write verbose logging to a file. May not be \
                                  specified with multiple input files.')
    
    # Multiprocessing arguments
    if multiproc:
        parser.add_argument('--nproc', action='store', dest='nproc', type=int, default=mp.cpu_count(),
                            help='The number of simultaneous computational processes to execute \
                                  (CPU cores to utilized).')
    
    # Universal arguments
    parser.add_argument('--outdir', action='store', dest='out_dir', default=None,
                        help='Specify to changes the output directory to the location specified. \
                              The input file directory is used if this is not specified.')
    parser.add_argument('--outname', action='store', dest='out_name', default=None,
                        help='Changes the prefix of the successfully processed output file \
                              to the string specified. May not be specified with multiple \
                              input files.')
    
    return parser


def parseCommonArgs(args, in_arg=None, in_types=None):
    """
    Checks common arguments from getCommonArgParser and transforms output options to a dictionary

    Arguments: 
    args = argument Namespace defined by ArgumentParser.parse_args
    in_arg = a string defining a non-standard input file argument to verify; by default
             ['db_files', 'seq_files', 'seq_files_1', 'seq_files_2', 'primer_file'] 
             are supported in that order
    in_types = a list of types (file extensions as strings) to allow for files in file_arg
               if None do not check type
                    
    Returns:
    a dictionary copy of args with output arguments embedded in the dictionary out_args
    """ 
    db_types = ['.tab']
    seq_types = ['.fasta', '.fastq']
    primer_types = ['.fasta', '.regex']
    if in_types is not None:  in_types = [f.lower for f in in_types]
    args_dict = args.__dict__.copy()
    
    # Count input files
    if 'db_files' in args_dict:      
        input_count = len(args_dict['db_files'])
    elif 'seq_files' in args_dict:       
        input_count = len(args_dict['seq_files'])
    elif 'seq_files_1' in args_dict and 'seq_files_2' in args_dict:   
        input_count = len(args_dict['seq_files_1'])
    elif 'primer_file' in args_dict:   
        input_count = 1
    elif in_arg is not None and in_arg in args_dict: 
        input_count = len(args_dict[in_arg])
    else:
        sys.exit('ERROR:  Cannot determine input file argument')
    
    # Exit if output names or log files are specified with multiple input files    
    if args_dict.get('out_name', None) is not None and input_count > 1:
        sys.exit('ERROR:  The --outname argument may not be specified with multiple input files')
    if args_dict.get('log_file', None) is not None and input_count > 1:
        sys.exit('ERROR:  The --log argument may not be specified with multiple input files')
        
    # Verify database files
    if 'db_files' in args_dict:
        for f in args_dict['db_files']:
            if not os.path.isfile(f):
                sys.exit('ERROR:  Database file %s does not exist' % f)
            if os.path.splitext(f)[-1].lower() not in db_types:
                sys.exit('ERROR:  Database file %s is not a supported type. Must be one: %s' \
                         % (','.join(db_types), f))       
    
    # Verify single-end sequence files
    if 'seq_files' in args_dict:
        for f in args_dict['seq_files']:
            if not os.path.isfile(f):
                sys.exit('ERROR:  Sequence file %s does not exist' % f)
            if os.path.splitext(f)[-1].lower() not in seq_types:
                sys.exit('ERROR:  Sequence file %s is not a supported type. Must be one: %s' \
                         % (','.join(seq_types), f))
    
    # Verify paired-end sequence files
    if 'seq_files_1' and 'seq_files_2' in args_dict:
        if len(args_dict['seq_files_1']) != len(args_dict['seq_files_2']):
            sys.exit('ERROR:  The -1 and -2 arguments must contain the same number of files')
        for f1, f2 in zip(args_dict['seq_files_1'], args_dict['seq_files_2']):
            if os.path.splitext(f1)[-1].lower() != os.path.splitext(f2)[-1].lower():
                sys.exit('ERROR:  Each pair of files in the -1 and -2 arguments must be the same file type')
        for f in (args_dict['seq_files_1'] + args_dict['seq_files_2']):
            if not os.path.isfile(f):
                sys.exit('ERROR:  Sequence file %s does not exist' % f)
            if os.path.splitext(f)[-1].lower() not in seq_types:
                sys.exit('ERROR:  Sequence file %s is not a supported type. Must be one: %s' \
                         % (','.join(seq_types), f))
                    
    # Verify primer file
    if 'primer_file' in args_dict:
        primer_file = args_dict['primer_file']
        if not os.path.isfile(primer_file):
            sys.exit('ERROR:  Primer file %s does not exist' % primer_file)
        if os.path.splitext(primer_file)[-1].lower() not in primer_types:
            sys.exit('ERROR:  Primer file %s is not a supported type. Must be one: %s' \
                     % (','.join(primer_types), primer_file))
    
    # Verify non-standard input files
    if in_arg is not None and in_arg in args_dict:
        files = args_dict[in_arg] if isinstance(args_dict[in_arg], list) \
                else [args_dict[in_arg]]
        for f in files:
            if not os.path.exists(f):
                sys.exit('ERROR:  Input %s does not exist' % f)
            if in_types is not None and os.path.splitext(f)[-1].lower() not in in_types:
                sys.exit('ERROR:  Input %s is not a supported type. Must be one: %s' \
                         % (','.join(in_types), f))
    
    # Verify output directory
    if 'out_dir' in args_dict and args_dict['out_dir'] is not None:
        if os.path.exists(args_dict['out_dir']) and not os.path.isdir(args_dict['out_dir']):
            sys.exit('ERROR:  Directory %s exists but it is not a directory' % args_dict['out_dir'])

    # Redefine common output options as out_args dictionary
    out_args = ['log_file', 'delimiter', 'separator', 
                'out_dir', 'out_name', 'out_type', 'clean']
    args_dict['out_args'] = {k:args_dict.setdefault(k, None) for k in out_args}
    for k in out_args: del args_dict[k]
    
    return args_dict


if __name__ == '__main__':
    """
    Print module information
    """
    print 'Version: %s %s %s' % (os.path.basename(__file__), __version__, __date__)
    print 'Location: %s' % os.path.dirname(os.path.realpath(__file__))
    #print 'Parent Dir: %s' % path.join(path.dirname(path.realpath(__file__)), path.pardir)
    