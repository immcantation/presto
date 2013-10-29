#!/usr/bin/env python
"""
Removes primers and annotates sequences with primer and barcode identifiers
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.1'
__date__      = '2013.10.29'

# Imports
import os, sys
import multiprocessing as mp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from itertools import izip
from Bio import pairwise2

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_delimiter, default_out_args
from IgCore import flattenAnnotation, parseAnnotation, mergeAnnotation
from IgCore import getCommonParser, parseCommonArgs, printLog
from IgCore import getScoreDict, reverseComplement, weightSeq
from IgCore import compilePrimers, readPrimerFile
from IgCore import collectRecQueue, feedRecQueue

# Defaults
default_max_error = 0.2
default_max_len = 50
default_start = 0


def alignPrimers(seq_record, primers, primers_regex=None, max_error=default_max_error, 
                 max_len=default_max_len, reverse=False, skip_rc=False, 
                 score_dict=getScoreDict(n_score=0, gap_score=0)):
    """
    Performs pairwise local alignment of a list of short sequences against a long sequence

    Arguments: 
    seq_record = a SeqRecord object to align primers against
    primers = dictionary of {names: short IUPAC ambiguous sequence strings}
    primers_regex = optional dictionary of {names: compiled primer regular expressions}
    max_error = maximum acceptable error rate before aligning reverse complement
    max_len = maximum length of sample sequence to align
    reverse = if True align with the tail end of the sequence
    skip_rc = if True do not check reverse complement sequences
    score_dict = optional dictionary of alignment scores as {(char1, char2): score}

    Returns:
    dictionary of {seq:input SeqRecord object, align:primer alignment, 
                   primer:primer id, start:alignment start, 
                   end:alignment end, offset:input sequence alignment offset,
                   error:error rate, reverse:True if alignment is tail-ended}
    """
    # Defined undefined parameters
    if primers_regex is None:  primers_regex = compilePrimers(primers)
    
    # Create empty return dictionary
    align_dict = {'seq':seq_record, 'align':None, 'primer':None,
                  'start':None, 'end':None, 'error':1, 
                  'reverse':reverse}

    # Define undefined arguments
    rec_len = len(seq_record)
    max_len = min(rec_len, max_len)
    
    # Define sequences to align and assign orientation tags
    if not skip_rc:
        seq_list = [seq_record, reverseComplement(seq_record)]
        seq_list[0].annotations['seqorient'] = 'F'
        seq_list[1].annotations['seqorient'] = 'RC'
    else:
        seq_list = [seq_record]
        seq_list[0].annotations['seqorient'] = 'F'
    
    # Assign primer orientation tags
    for rec in seq_list:
        rec.annotations['prorient'] = 'F' if not reverse else 'RC' 
    
    # Attempt regular expression match first
    for rec in seq_list:
        scan_rec = rec[:max_len] if not reverse else rec[-max_len:]
        for adpt_id, adpt_regex in primers_regex.iteritems():
            adpt_match = adpt_regex.search(str(scan_rec.seq))
            # Parse matches
            if adpt_match:
                align_dict['seq'] = rec
                align_dict['align'] = '-' * adpt_match.start(0) + primers[adpt_id] + \
                                      '-' * (max_len - adpt_match.end(0))
                align_dict['primer'] = adpt_id
                align_dict['seq'].annotations['primer'] = adpt_id
                align_dict['error'] = 0
                align_dict['offset'] = 0
                if not reverse:
                    align_dict['start'] = adpt_match.start(0)
                    align_dict['end'] = adpt_match.end(0)
                else:
                    align_dict['start'] = adpt_match.start(0) + rec_len - max_len
                    align_dict['end'] = adpt_match.end(0) + rec_len - max_len

                return align_dict
    
    # Perform local alignment
    best_align, best_rec, best_adpt, best_err = None, None, None, 1
    for rec in seq_list:
        scan_rec = rec[:max_len] if not reverse else rec[-max_len:]
        this_align = {}
        for adpt_id, adpt_seq in primers.iteritems():
            pw2_align = pairwise2.align.localds(scan_rec.seq, adpt_seq, score_dict, 
                                                -1, -1, one_alignment_only=True)
            if pw2_align:  this_align.update({adpt_id: pw2_align[0]})
        if not this_align:  continue
        
        # Determine alignment with lowest error rate
        for adpt, algn in this_align.iteritems():
            adpt_err = 1 - algn[2] / weightSeq(primers[adpt])
            if adpt_err < best_err: 
                best_align = this_align
                best_rec = rec
                best_adpt = adpt
                best_err = adpt_err
        
        # Skip reverse complement if forward sequence error within defined threshold
        if best_err <= max_error:  break

    # Set return dictionary to lowest error rate alignment
    if best_align:
        align_offset = len(best_align[best_adpt][0]) - max_len
        align_dict['seq'] = best_rec
        align_dict['align'] = best_align[best_adpt][1]
        align_dict['primer'] = best_adpt
        align_dict['error'] = best_err
        align_dict['offset'] = align_offset
        if not reverse:
            align_dict['start'] = max([best_align[best_adpt][3] - align_offset, 0])
            align_dict['end'] = best_align[best_adpt][4] - align_offset
        else:
            align_dict['start'] = best_align[best_adpt][3] + rec_len - max_len
            align_dict['end'] = best_align[best_adpt][4] + rec_len - max_len
    
    return align_dict


def scorePrimers(seq_record, primers, start=default_start, reverse=False, 
                 score_dict=getScoreDict(n_score=0, gap_score=0)):
    """
    Performs simple alignment of primers with a fixed starting position, 
    no reverse complement alignment, and no tail alignment option

    Arguments: 
    seq_record = a SeqRecord object to align primers against
    primers = dictionary of {names: short IUPAC ambiguous sequence strings}
    start = position where primer alignment starts
    reverse = if True align with the tail end of the sequence
    score_dict = optional dictionary of {(char1, char2): score} alignment scores
    
    Returns:
    dictionary of {seq:input SeqRecord object, align:primer alignment, 
                   primer:primer id, start:alignment start, 
                   end:alignment end, offset:input sequence alignment offset,
                   error:error rate, reverse:True if alignment is tail-ended}
    """
    # Create empty return dictionary
    align_dict = {'seq':seq_record, 'align':None, 'primer':None,
                  'start':None, 'end':None, 'offset':0, 'error':1,
                  'reverse':reverse}

    # Define orientation variables
    seq_record.annotations['seqorient'] = 'F'
    seq_record.annotations['prorient'] = 'F' if not reverse else 'RC'

    this_align = {}
    rec_len = len(seq_record)
    if reverse:  end = rec_len - start
    for adpt_id, adpt_seq in primers.iteritems():
        if reverse:  start = end - len(adpt_seq)
        else:  end = start + len(adpt_seq)
        chars = izip(adpt_seq, seq_record[start:end])
        score = sum([score_dict[(c1, c2)] for c1, c2 in chars])
        this_align.update({adpt_id: (score, start, end)})
        
    # Determine alignment with lowest error rate
    best_align, best_adpt, best_err = None, None, 1
    for adpt, algn in this_align.iteritems():
        adpt_err = 1.0 - float(algn[0]) / weightSeq(primers[adpt])
        if adpt_err < best_err:
            best_align = algn
            best_adpt = adpt
            best_err = adpt_err

    # Set return dictionary to lowest error rate alignment
    if best_align:
        align_dict['align'] = '-' * best_align[1] + primers[best_adpt] if not reverse \
                              else primers[best_adpt] + '-' * (rec_len - best_align[2]) 
        align_dict['primer'] = best_adpt
        align_dict['start'] = best_align[1]
        align_dict['end'] = best_align[2]
        align_dict['error'] = best_err
    
    return align_dict


def getMaskedSeq(align, mode='mask', barcode=False, delimiter=default_delimiter):
    """
    Create an output sequence with primers masked or cut

    Arguments: 
    align = an alignment dictionary from alignPrimers or scorePrimers
    mode = defines the action taken; one of ['cut','mask','tag','trim']
    barcode = if True add sequence preceding primer to description
    delimiter = a tuple of delimiters for (annotations, field/values, value lists) 

    Returns:
    output SeqRecord object
    """
    seq = align['seq']
    reverse = align['reverse']
    
    # Build output sequence
    if mode == 'tag' or not align['align']:
        out_seq = seq
    elif mode == 'trim':
        if not reverse:  
            out_seq = seq[align['start']:]
        else:  
            out_seq = seq[:align['end']]
    elif mode == 'cut':
        if not reverse:  
            out_seq = seq[align['end']:]
        else: 
            out_seq = seq[:align['start']]
    elif mode == 'mask':
        if not reverse:
            out_seq = 'N' * (align['end'] - align['start']) + seq[align['end']:]
            if hasattr(seq, 'letter_annotations') and 'phred_quality' in seq.letter_annotations:
                out_seq.letter_annotations['phred_quality'] = seq.letter_annotations['phred_quality'][align['start']:]
        else:
            out_seq = seq[:align['start']] + 'N' * (min(align['end'], len(seq)) - align['start'])
            if hasattr(seq, 'letter_annotations') and 'phred_quality' in seq.letter_annotations:
                out_seq.letter_annotations['phred_quality'] = seq.letter_annotations['phred_quality'][:len(out_seq)]
            
    # Add alignment annotations to output SeqRecord
    out_seq.annotations = seq.annotations    
    out_seq.annotations['primer'] = align['primer']
    out_seq.annotations['prstart'] = align['start']
    out_seq.annotations['error'] = align['error']

    # Parse seq annotation and create output annotation
    seq_ann = parseAnnotation(seq.description, delimiter=delimiter)
    out_ann = OrderedDict([('SEQORIENT', align['seq'].annotations['seqorient']),
                           ('PRIMER', align['primer'])])
    
    # Add ID sequence to description
    if barcode:
        seq_code = seq[:align['start']].seq if not reverse \
                   else seq[align['end']:].seq
        out_seq.annotations['barcode'] = seq_code
        out_ann['BARCODE'] = seq_code
    
    out_ann = mergeAnnotation(seq_ann, out_ann, delimiter=delimiter)
    out_seq.id = flattenAnnotation(out_ann, delimiter=delimiter)
    out_seq.description = ''

    #return out_seq if len(out_seq) > 0 else None
    return out_seq


def processQueue(data_queue, result_queue, primers, mode, align_func, align_args={}, 
                 max_error=default_max_error, barcode=False, delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    primers = dictionary of primer sequences
    mode = defines the action taken; one of ['cut','mask','tag','trim']
    align_func = the function to call for alignment
    align_args = a dictionary of arguments to pass to align_func
    max_error = maximum acceptable error rate for a valid alignment
    barcode = if True add sequence preceding primer to description
    delimiter = a tuple of delimiters for (fields, values, value lists) 

    Returns: 
    None
    """
    # Iterator over data queue until sentinel object reached
    for args in iter(data_queue.get, None):
        in_seq = args['seq']
        
        # Define result dictionary for iteration
        results = {'id':args['id'],
                   'in_seq':in_seq,
                   'out_seq':None,
                   'error': None,
                   'pass':False,
                   'log':OrderedDict([('ID', args['id'])])}
 
        # Align primers
        align = align_func(in_seq, primers, **align_args)
        
        # Process alignment results
        if align['align'] is None:
            # Update log if no alignment
            results['log']['INSEQ'] = in_seq.seq
            results['log']['ALIGN'] = None
        else:
            # Create output sequence
            out_seq = getMaskedSeq(align, mode, barcode, delimiter)        
            results['out_seq'] = out_seq
            results['error'] = align['error']
            results['pass'] = (align['error'] <= max_error) if len(out_seq) > 0 else False
            
            # Update log with successful alignment results
            results['log']['SEQORIENT'] = out_seq.annotations['seqorient']
            results['log']['PRIMER'] = out_seq.annotations['primer']
            results['log']['PRORIENT'] = out_seq.annotations['prorient']
            results['log']['PRSTART'] = out_seq.annotations['prstart']
            if 'barcode' in out_seq.annotations:  
                results['log']['BARCODE'] = out_seq.annotations['barcode']
            if not align['reverse']:
                results['log']['INSEQ'] = '-' * align['offset'] + str(align['seq'].seq)
                results['log']['ALIGN'] = align['align']
                results['log']['OUTSEQ'] = str(out_seq.seq).rjust(len(in_seq) + align['offset'])
            else:
                results['log']['INSEQ'] = str(align['seq'].seq) + '-' * align['offset']
                results['log']['ALIGN'] = align['align'].rjust(len(in_seq) + align['offset'])
                results['log']['OUTSEQ'] = str(out_seq.seq)
            results['log']['ERROR'] = align['error']
        
        # Feed results to result queue
        result_queue.put(results)

    return None


def maskPrimers(seq_file, primer_file, mode, align_func, align_args={}, 
                max_error=default_max_error, barcode=False,
                out_args=default_out_args, nproc=None, queue_size=None):
    """
    Masks or cuts primers from sample sequences using local alignment

    Arguments: 
    seq_file = name of file containing sample sequences
    primer_file = name of the file containing primer sequences
    mode = defines the action taken; one of 'cut','mask','tag'
    align_func = the function to call for alignment
    align_arcs = a dictionary of arguments to pass to align_func
    max_error = maximum acceptable error rate for a valid alignment
    barcode = if True add sequence preceding primer to description
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                 
    Returns:
    a list of successful output file names
    """
    # Define number of processes and queue size
    if nproc is None:  nproc = mp.cpu_count()
    if queue_size is None:  queue_size = nproc * 2
    
    # Define subcommand label dictionary
    cmd_dict = {alignPrimers:'align', scorePrimers:'score'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MaskPrimers'
    log['COMMAND'] = cmd_dict.get(align_func, align_func.__name__)
    log['SEQ_FILE'] = os.path.basename(seq_file)
    log['PRIMER_FILE'] = os.path.basename(primer_file)
    log['MODE'] = mode
    log['BARCODE'] = barcode
    log['MAX_ERROR'] = max_error
    if 'start' in align_args: log['START_POS'] = align_args['start']
    if 'max_len' in align_args: log['MAX_LEN'] = align_args['max_len']
    if 'reverse' in align_args: log['REVERSE'] = align_args['reverse']
    if 'skip_rc' in align_args: log['SKIP_RC'] = align_args['skip_rc']
    log['NPROC'] = nproc
    printLog(log)

    # Create dictionary of primer sequences to pass to maskPrimers
    primers = readPrimerFile(primer_file)
    if 'reverse' in align_args and align_args['reverse']:
        primers = {k: reverseComplement(v) for k, v in primers.iteritems()}
    
    # Define alignment arguments and compile primers for align mode
    align_args['score_dict'] = getScoreDict(n_score=0, gap_score=0)
    if align_func is alignPrimers:
        align_args['max_error'] = max_error
        align_args['primers_regex'] = compilePrimers(primers)

    # Define shared data objects
    manager = mp.Manager()
    collect_dict = manager.dict()
    data_queue = mp.Queue(queue_size)
    result_queue = mp.Queue(queue_size)
    
    # Initiate feeder process
    feeder = mp.Process(target=feedRecQueue, args=(data_queue, nproc, seq_file))
    feeder.start()

    # Initiate processQueue processes
    workers = []
    for __ in range(nproc):
        w = mp.Process(target=processQueue, args=(data_queue, result_queue, primers, mode, 
                                                  align_func, align_args, max_error, barcode, 
                                                  out_args['delimiter']))
        w.start()
        workers.append(w)

    # Initiate collector process
    collector = mp.Process(target=collectRecQueue, args=(result_queue, collect_dict, 
                                                         seq_file, 'primers', out_args))
    collector.start()

    # Wait for feeder and worker processes to finish, add sentinel to result_queue
    feeder.join()
    for w in workers:  w.join()
    result_queue.put(None)
    
    # Wait for collector process to finish and shutdown manager
    collector.join()
    log = collect_dict['log']
    out_files = collect_dict['out_files']
    manager.shutdown()
    
    # Print log
    log['END'] = 'MaskPrimers'
    printLog(log)
        
    return out_files


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
                            formatter_class=ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()
    
    # Parent parser
    parser_parent = getCommonParser(multiproc=True)
    parser_parent.add_argument('-p', action='store', dest='primer_file', required=True, 
                               help='A FASTA or REGEX file containing primer sequences')
    parser_parent.add_argument('--mode', action='store', dest='mode', 
                               choices=('cut', 'mask', 'tag', 'trim'), default='mask', 
                               help='Specifies whether or mask primers, cut primers, \
                                     trim region preceding primer, or only tag sample sequences')
    parser_parent.add_argument('--maxerror', action='store', dest='max_error', type=float,
                               default=default_max_error, help='Maximum allowable error rate')
    parser_parent.add_argument('--reverse', action='store_true', dest='reverse', 
                              help='Specify to reverse complement primer sequences and scan/cut from the end')
    parser_parent.add_argument('--barcode', action='store_true', dest='barcode', 
                               help='Specify to encode sequences with barcode sequences preceding primer matches')
    
    # Align mode argument parser
    parser_align = subparsers.add_parser('align', parents=[parser_parent],
                                         formatter_class=ArgumentDefaultsHelpFormatter,
                                         help='Find primer matches using pairwise local alignment')
    parser_align.add_argument('--maxlen', action='store', dest='max_len', type=int,
                              default=default_max_len, help='Maximum sequence length to scan for primers')
    parser_align.add_argument('--skiprc', action='store_true', dest='skip_rc', 
                              help='Specify to prevent checking of sample reverse complement sequences')
    #parser_align.set_defaults(start=None)
    parser_align.set_defaults(align_func=alignPrimers)
    

    # Score mode argument parser
    parser_score = subparsers.add_parser('score', parents=[parser_parent], 
                                         formatter_class=ArgumentDefaultsHelpFormatter,
                                         help='Find primer matches by scoring primers at a fixed position')
    parser_score.add_argument('--start', action='store', dest='start', type=int, default=default_start, 
                              help='The starting position of the primer')

    parser_score.set_defaults(align_func=scorePrimers)
    
    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    
    # Define align_args dictionary to pass to maskPrimers
    if args_dict['align_func'] is alignPrimers:
        args_dict['align_args'] = {'max_len':args_dict['max_len'],
                                   'reverse':args_dict['reverse'],
                                   'skip_rc':args_dict['skip_rc']}
        del args_dict['max_len']
        del args_dict['reverse']
        del args_dict['skip_rc']
    elif args_dict['align_func'] is scorePrimers:
        args_dict['align_args'] = {'start':args_dict['start'],
                                   'reverse':args_dict['reverse']}
        del args_dict['start']
        del args_dict['reverse']
    
    # Call maskPrimers for each sample file
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        maskPrimers(**args_dict)

        # Profiling
        #import cProfile, pstats
        #cProfile.run('maskPrimers(**args_dict)', 'profile.prof')
        #p = pstats.Stats('profile.prof')
        #p.strip_dirs().sort_stats('time').print_stats()      