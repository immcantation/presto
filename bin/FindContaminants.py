#!/usr/bin/env python3
"""
Identifies reads shared across samples by UMI barcode and sequence
"""
# Info
__author__ = 'Gisela Gabernet'
from presto import __version__, __date__

# Imports
import os
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import chain, combinations
from textwrap import dedent
from time import time
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field, default_out_args, \
                            default_missing_chars, choices_coord, default_coord
from presto.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation, getCoordKey
from presto.IO import openFile, getFileType, getOutputHandle, printLog, printMessage, \
                      printError, printWarning

# Defaults
default_contam_max_error = 0.0
default_shared_field = 'SHARED'


def getSampleLabel(seq_file):
    """
    Defines a default sample label from a file name

    Arguments:
      seq_file : sequence file name.

    Returns:
      str : the file name without directory, .gz extension and file extension.
    """
    name = os.path.basename(seq_file)
    if name.endswith('.gz'):  name = name[:-3]
    return os.path.splitext(name)[0]


def iterSeqFile(seq_file):
    """
    Iterates over the raw records of a FASTA or FASTQ file without building SeqRecords

    Arguments:
      seq_file : FASTA or FASTQ file name, optionally gzip compressed.

    Returns:
      iter : an iterator of (title, sequence, quality) string tuples;
             quality is None for FASTA files.
    """
    seq_type = getFileType(seq_file)
    with openFile(seq_file, 'r') as handle:
        if seq_type == 'fastq':
            for title, seq, qual in FastqGeneralIterator(handle):
                yield title, seq, qual
        elif seq_type == 'fasta':
            for title, seq in SimpleFastaParser(handle):
                yield title, seq, None
        else:
            printError('File %s must be FASTA or FASTQ.' % seq_file)


def matchSeq(seq_1, seq_2, max_error=default_contam_max_error,
             ignore_chars=default_missing_chars):
    """
    Checks whether two read sequences match

    Sequences of different lengths never match. Positions holding a missing character
    (N, -, .) in either sequence are ignored.

    Arguments:
      seq_1 : first sequence as an uppercase string.
      seq_2 : second sequence as an uppercase string.
      max_error : the maximum mismatch rate allowed over the compared positions.
      ignore_chars : set of characters to ignore.

    Returns:
      bool : True if the sequences have the same length, the mismatch rate is at most
             max_error and at least one position was compared.
    """
    if len(seq_1) != len(seq_2):
        return False
    if seq_1 == seq_2:
        return any(c not in ignore_chars for c in seq_1)

    # The compared length never exceeds the sequence length, so exceeding this bound means failure
    max_mismatch = max_error * len(seq_1)
    mismatch = weight = 0
    for a, b in zip(seq_1, seq_2):
        if a in ignore_chars or b in ignore_chars:
            continue
        weight += 1
        if a != b:
            mismatch += 1
            if mismatch > max_mismatch:
                return False

    return weight > 0 and mismatch <= max_error * weight


def collectSampleBarcodes(seq_file_1, seq_file_2, field=default_barcode_field,
                          delimiter=default_delimiter):
    """
    Collects the set of UMI barcodes annotated in a paired-end sample

    Arguments:
      seq_file_1 : read 1 sequence file.
      seq_file_2 : read 2 sequence file.
      field : the annotation field containing the UMI barcode.
      delimiter : a tuple of delimiters for (annotations, field/values, value lists).

    Returns:
      set : barcodes annotated in either file.
    """
    barcodes = set()
    for seq_file in (seq_file_1, seq_file_2):
        for title, _, _ in iterSeqFile(seq_file):
            tag = parseAnnotation(title, delimiter=delimiter).get(field, None)
            if tag is not None:  barcodes.add(tag)

    return barcodes


def collectSampleReads(seq_file_1, seq_file_2, candidates=None, field=default_barcode_field,
                       coord_type=default_coord, delimiter=default_delimiter):
    """
    Groups the read pairs of a paired-end sample by UMI barcode

    The barcode of a read pair is taken from the read 1 annotation when present,
    otherwise from the read 2 annotation. Mates are matched by their coordinate keys,
    so the barcode only needs to be annotated in one of the two files. Reads without
    a mate are not returned.

    Arguments:
      seq_file_1 : read 1 sequence file.
      seq_file_2 : read 2 sequence file.
      candidates : set of barcodes to keep; if None keep all barcodes.
      field : the annotation field containing the UMI barcode.
      coord_type : the sequence header format used to match mates.
      delimiter : a tuple of delimiters for (annotations, field/values, value lists).

    Returns:
      tuple : (groups, conflicts) where groups maps each barcode to a list of
              (coordinate, read 1 sequence, read 2 sequence) tuples, and conflicts is the
              number of read pairs found with different barcodes in read 1 and read 2.
    """
    def _parse(title):
        key = getCoordKey(title.split(None, 1)[0], coord_type=coord_type, delimiter=delimiter)
        tag = parseAnnotation(title, delimiter=delimiter).get(field, None)
        return key, tag

    def _keep(tag):
        return tag is not None and (candidates is None or tag in candidates)

    # Read 1 reads carrying a barcode
    reads = {}
    r1_untagged = 0
    for title, seq, _ in iterSeqFile(seq_file_1):
        key, tag = _parse(title)
        if tag is None:
            r1_untagged += 1
        elif _keep(tag):
            reads[key] = [tag, seq.upper(), None]

    # Add mates from read 2. Read 2 barcodes are only needed for untagged read 1 mates.
    pending = {}
    conflicts = 0
    for title, seq, _ in iterSeqFile(seq_file_2):
        key, tag = _parse(title)
        rec = reads.get(key, None)
        if rec is not None:
            rec[2] = seq.upper()
            if tag is not None and tag != rec[0]:  conflicts += 1
        elif r1_untagged > 0 and _keep(tag):
            pending[key] = [tag, None, seq.upper()]

    # Fetch untagged read 1 mates of read 2 reads carrying the barcode
    if pending:
        for title, seq, _ in iterSeqFile(seq_file_1):
            key, tag = _parse(title)
            rec = pending.get(key, None)
            if rec is None:
                continue
            if tag is None:
                rec[1] = seq.upper()
            else:
                # Read 1 barcode takes precedence and was not kept
                conflicts += 1
                del pending[key]
        reads.update(pending)

    groups = {}
    for key, (tag, seq_1, seq_2) in reads.items():
        if seq_1 is not None and seq_2 is not None:
            groups.setdefault(tag, []).append((key, seq_1, seq_2))

    return groups, conflicts


def findSharedReads(sample_reads, max_error=default_contam_max_error,
                    ignore_chars=default_missing_chars):
    """
    Identifies read pairs matching a read pair from another sample within one barcode

    Uses the duplicate search strategy of CollapseSeq. Identical read pairs are collapsed
    within each sample. When max_error is 0, read pairs without missing characters are then
    matched across samples by hashing. Read pairs containing missing characters are compared
    position by position against read pairs of the same read lengths only. When max_error is
    above 0, all read pairs of the same read lengths are compared position by position.
    Read pairs whose read 1 or read 2 lengths differ never match.

    Arguments:
      sample_reads : OrderedDict mapping sample label to a list of
                     (coordinate, read 1 sequence, read 2 sequence) tuples.
      max_error : the maximum mismatch rate allowed for both read 1 and read 2.
      ignore_chars : set of missing characters ignored when comparing sequences.

    Returns:
      dict : maps each sample label to a dictionary of the matched read coordinates and
             the list of other samples each read matched.
    """
    labels = list(sample_reads)

    # Collapse identical read pairs within each sample
    unique = OrderedDict()
    for label, reads in sample_reads.items():
        seqs = OrderedDict()
        for key, seq_1, seq_2 in reads:
            seqs.setdefault((seq_1, seq_2), []).append(key)
        unique[label] = seqs

    matched = {label: {} for label in labels}
    def _record(a, pair_a, b, pair_b):
        matched[a].setdefault(pair_a, set()).add(b)
        matched[b].setdefault(pair_b, set()).add(a)

    def _compare(pair_a, pair_b):
        return matchSeq(pair_a[0], pair_b[0], max_error, ignore_chars) and \
               matchSeq(pair_a[1], pair_b[1], max_error, ignore_chars)

    # Group read pairs by read lengths, since only equal lengths can match
    len_dict = OrderedDict()
    for label in labels:
        for pair in unique[label]:
            len_dict.setdefault((len(pair[0]), len(pair[1])), []).append((label, pair))

    for entries in len_dict.values():
        if max_error > 0:
            # Mismatches are allowed, so every read pair from different samples must be scored
            for i, (label, pair) in enumerate(entries):
                for label_b, pair_b in entries[i + 1:]:
                    if label_b != label and _compare(pair, pair_b):
                        _record(label, pair, label_b, pair_b)
            continue

        # Split read pairs into those with and without missing characters
        clean, ambig = [], []
        for label, pair in entries:
            if not pair[0] or not pair[1] or \
                    any(c in ignore_chars for c in chain(pair[0], pair[1])):
                ambig.append((label, pair))
            else:
                clean.append((label, pair))

        # Exact matching of clean read pairs by hashing
        index = {}
        for label, pair in clean:
            index.setdefault(pair, []).append(label)
        for pair, pair_labels in index.items():
            for a, b in combinations(pair_labels, 2):
                _record(a, pair, b, pair)

        # Score read pairs with missing characters against all other read pairs
        for i, (label, pair) in enumerate(ambig):
            for label_b, pair_b in chain(clean, ambig[i + 1:]):
                if label_b != label and _compare(pair, pair_b):
                    _record(label, pair, label_b, pair_b)

    shared = {label: {} for label in labels}
    for label in labels:
        for pair, others in matched[label].items():
            others = sorted(others, key=labels.index)
            for key in unique[label][pair]:
                shared[label][key] = others

    return shared


def writeSampleReads(seq_file, flagged, shared_field=default_shared_field,
                     coord_type=default_coord, out_args=default_out_args):
    """
    Writes the reads of one file to pass and fail files

    Arguments:
      seq_file : input sequence file.
      flagged : dictionary mapping flagged read coordinates to the list of samples shared with.
      shared_field : the annotation field listing the samples a read is shared with.
      coord_type : the sequence header format used to match mates.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      tuple : (pass file name, pass count, fail count).
    """
    delimiter = out_args['delimiter']
    out_type = out_args['out_type'] if out_args['out_type'] else getFileType(seq_file)
    gzip_output = out_args.get('gzip_output', False)

    def _write(handle, title, seq, qual):
        if out_type == 'fastq':
            handle.write('@%s\n%s\n+\n%s\n' % (title, seq, qual))
        else:
            handle.write('>%s\n%s\n' % (title, seq))

    pass_handle = getOutputHandle(seq_file, 'contam-pass', out_dir=out_args['out_dir'],
                                  out_type=out_type, gzip_output=gzip_output)
    fail_handle = None
    if out_args['failed']:
        fail_handle = getOutputHandle(seq_file, 'contam-fail', out_dir=out_args['out_dir'],
                                      out_type=out_type, gzip_output=gzip_output)

    pass_count = fail_count = 0
    for title, seq, qual in iterSeqFile(seq_file):
        key = getCoordKey(title.split(None, 1)[0], coord_type=coord_type, delimiter=delimiter)
        others = flagged.get(key, None)
        if others is None:
            pass_count += 1
            _write(pass_handle, title, seq, qual)
        else:
            fail_count += 1
            if fail_handle is not None:
                ann = parseAnnotation(title, delimiter=delimiter)
                ann = mergeAnnotation(ann, {shared_field: others}, delimiter=delimiter)
                _write(fail_handle, flattenAnnotation(ann, delimiter=delimiter), seq, qual)

    pass_handle.close()
    if fail_handle is not None:  fail_handle.close()

    return pass_handle.name, pass_count, fail_count


def findContaminants(seq_files_1, seq_files_2, labels=None,
                     barcode_field=default_barcode_field, coord_type=default_coord,
                     max_error=default_contam_max_error,
                     shared_field=default_shared_field, out_args=default_out_args):
    """
    Identifies reads shared across samples by UMI barcode and sequence

    A read pair is shared when another sample contains a read pair with the same UMI
    barcode whose read 1 and read 2 sequences both match it. Shared read pairs are
    removed from every sample containing them.

    To limit memory usage the input files are read in three passes. The first pass
    collects only the barcodes of each sample. The second keeps only the read pairs whose
    barcode occurs in multiple samples. The third writes the output files.

    Arguments:
      seq_files_1 : list of read 1 sequence files, one per sample.
      seq_files_2 : list of read 2 sequence files, in the same sample order.
      labels : list of sample labels; if None use the read 1 file names.
      barcode_field : the annotation field containing the UMI barcode.
      coord_type : the sequence header format used to match mates.
      max_error : the maximum read 1 and read 2 mismatch rate allowed between two
                  read pairs to consider them shared.
      shared_field : the annotation field listing the samples a flagged read is shared with.
      out_args : common output argument dictionary from parseCommonArgs. The log_file entry
                 defines the contaminant log; if None it is written to contam_log.txt in the
                 output directory, or in the directory of the first read 1 file.

    Returns:
      list : a list of tuples holding the (read 1, read 2) pass file names for each sample.
    """
    delimiter = out_args['delimiter']
    if labels is None:
        labels = [getSampleLabel(f) for f in seq_files_1]
    sample_files = OrderedDict(zip(labels, zip(seq_files_1, seq_files_2)))

    # Print parameter info
    log = OrderedDict()
    log['START'] = 'FindContaminants'
    log['FILES1'] = ','.join(os.path.basename(f) for f in seq_files_1)
    log['FILES2'] = ','.join(os.path.basename(f) for f in seq_files_2)
    log['SAMPLES'] = ','.join(labels)
    log['BARCODE_FIELD'] = barcode_field
    log['COORD_TYPE'] = coord_type
    log['MAX_ERROR'] = max_error
    printLog(log)

    # Pass 1: find barcodes annotated in multiple samples
    start_time = time()
    printMessage('Collecting barcodes', start_time=start_time)
    barcode_count = {}
    for f1, f2 in sample_files.values():
        for tag in collectSampleBarcodes(f1, f2, field=barcode_field, delimiter=delimiter):
            barcode_count[tag] = barcode_count.get(tag, 0) + 1
    candidates = {t for t, n in barcode_count.items() if n > 1}
    del barcode_count
    printMessage('Done', start_time=start_time, end=True)

    # Pass 2: group paired reads of candidate barcodes
    start_time = time()
    printMessage('Collecting reads', start_time=start_time)
    sample_groups = OrderedDict()
    for label, (f1, f2) in sample_files.items():
        groups, conflicts = collectSampleReads(f1, f2, candidates=candidates, field=barcode_field,
                                               coord_type=coord_type, delimiter=delimiter)
        if conflicts > 0:
            printWarning('%i read pairs in sample %s have different barcodes in read 1 and read 2. '
                         'The read 1 barcode was used.' % (conflicts, label))
        sample_groups[label] = groups
    printMessage('Done', start_time=start_time, end=True)

    # Find barcodes with paired reads in multiple samples
    barcode_samples = OrderedDict()
    for label, groups in sample_groups.items():
        for tag in groups:
            barcode_samples.setdefault(tag, []).append(label)
    multi_barcodes = [t for t, s in barcode_samples.items() if len(s) > 1]

    # Compare reads sharing a barcode
    start_time = time()
    printMessage('Comparing reads', start_time=start_time)
    log_file = out_args['log_file']
    if log_file is None:
        log_dir = out_args['out_dir'] if out_args['out_dir'] is not None \
                  else os.path.dirname(seq_files_1[0])
        log_file = os.path.join(log_dir, 'contam_log.txt')
    log_dir = os.path.dirname(log_file)
    if log_dir:  os.makedirs(log_dir, exist_ok=True)

    flagged = {label: {} for label in labels}
    contam_count = 0
    with open(log_file, 'w') as log_handle:
        for tag in multi_barcodes:
            tag_reads = OrderedDict((s, sample_groups[s][tag]) for s in barcode_samples[tag])
            shared = findSharedReads(tag_reads, max_error=max_error)
            sharing = [s for s in tag_reads if shared[s]]
            if not sharing:
                continue
            contam_count += 1
            for s in sharing:
                flagged[s].update(shared[s])

            # Write a log record with the read counts of each contaminated sample
            rec_log = OrderedDict()
            rec_log['BARCODE'] = tag
            rec_log['SAMPLES'] = ','.join(sharing)
            rec_log['SEQCOUNT'] = ','.join(str(len(shared[s])) for s in sharing)
            rec_log['GROUPCOUNT'] = ','.join(str(len(tag_reads[s])) for s in sharing)
            printLog(rec_log, handle=log_handle)
    del sample_groups
    printMessage('Done', start_time=start_time, end=True)

    # Pass 3: write output files
    start_time = time()
    printMessage('Writing output', start_time=start_time)
    out_files = []
    out_log = OrderedDict()
    for label, files in sample_files.items():
        pass_names = []
        pass_total = fail_total = 0
        for seq_file in files:
            name, pass_count, fail_count = writeSampleReads(seq_file, flagged[label],
                                                            shared_field=shared_field,
                                                            coord_type=coord_type,
                                                            out_args=out_args)
            pass_names.append(name)
            pass_total += pass_count
            fail_total += fail_count
        out_files.append(tuple(pass_names))
        out_log['%s_PASS' % label] = pass_total
        out_log['%s_FAIL' % label] = fail_total
    printMessage('Done', start_time=start_time, end=True)

    # Print log
    log = OrderedDict()
    log['OUTPUT'] = ','.join(os.path.basename(f) for p in out_files for f in p)
    log['MULTISAMPLE'] = len(multi_barcodes)
    log['CONTAMINATED'] = contam_count
    log['LOG'] = os.path.basename(log_file)
    log.update(out_log)
    log['END'] = 'FindContaminants'
    printLog(log)

    return out_files


def getArgParser():
    """
    Defines the ArgumentParser

    Returns:
      argparse.ArgumentParser: argument parser object.
    """
    # Define output file names and header fields
    fields = dedent(
             '''
             A read pair is considered contamination when another sample contains a read
             pair with the same UMI barcode and matching read 1 and read 2 sequences.
             Such read pairs are removed from every sample containing them. Reads are
             compared individually, so other reads with the same barcode are retained.
             To limit memory usage, the input files are read in three passes: the first
             collects only barcodes, the second keeps only reads whose barcode occurs in
             multiple samples, and the third writes the output.

             output files:
                 contam-pass
                     reads not shared with another sample, one file per input file.
                 contam-fail
                     reads shared with another sample, one file per input file.
                     Only written with --failed.

             output annotation fields:
                 SHARED
                     added to contam-fail reads. A comma delimited list of the other
                     samples containing a matching read pair with the same UMI barcode.

             log fields:
                 The contaminant log contains one record per UMI barcode with shared
                 read pairs.
                 BARCODE
                     the UMI barcode.
                 SAMPLES
                     the samples containing shared read pairs with this barcode.
                 SEQCOUNT
                     the number of shared read pairs removed from each of those samples,
                     in the same order.
                 GROUPCOUNT
                     the total number of read pairs with this barcode in each of those
                     samples, in the same order.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(seq_paired=True, out_file=False, log=False)],
                            formatter_class=CommonHelpFormatter, add_help=False)

    group_contam = parser.add_argument_group('contaminant identification arguments')
    group_contam.add_argument('--log', action='store', dest='contam_log', default=None,
                              help='''The contaminant log file, with one record per UMI barcode
                                   with shared read pairs. If unspecified, the log is written to
                                   contam_log.txt in the output directory, or in the directory
                                   of the first read 1 file.''')
    group_contam.add_argument('--labels', nargs='+', action='store', dest='labels', default=None,
                              help='''Sample labels, one per pair of input files and in the same
                                   order. If unspecified, the read 1 file names without extension
                                   are used.''')
    group_contam.add_argument('--bf', action='store', dest='barcode_field', type=str,
                              default=default_barcode_field,
                              help='''The annotation field containing the UMI barcode. It may be
                                   present in read 1, read 2 or both.''')
    group_contam.add_argument('--coord', action='store', dest='coord_type',
                              choices=choices_coord, default=default_coord,
                              help='''The format of the sequence identifier which defines shared
                                   coordinate information across mate pairs.''')
    group_contam.add_argument('--maxerror', action='store', dest='max_error', type=float,
                              default=default_contam_max_error,
                              help='''The maximum mismatch rate between two read pairs, for both
                                   read 1 and read 2, to consider them shared. Ambiguous and gap
                                   positions are ignored. Reads of different lengths never
                                   match. The default of 0 requires identical sequences and
                                   matches reads by hashing, as CollapseSeq does. Values above 0
                                   compare every pair of reads sharing a barcode and are slower.''')
    group_contam.add_argument('--sf', action='store', dest='shared_field', type=str,
                              default=default_shared_field,
                              help='''The annotation field added to flagged reads listing the
                                   samples sharing the read pair.''')

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

    # Check inputs
    n_samples = len(args_dict['seq_files_1'])
    if n_samples < 2:
        parser.error('At least two samples (two pairs of -1 and -2 files) are required.')
    if args_dict['out_args']['out_name'] is not None:
        parser.error('The --outname argument is not supported because multiple samples are processed.')
    if len(set(getFileType(f) for f in args_dict['seq_files_1'] + args_dict['seq_files_2'])) > 1:
        parser.error('All input files must have the same file type.')
    labels = args_dict['labels']
    if labels is None:
        labels = [getSampleLabel(f) for f in args_dict['seq_files_1']]
    elif len(labels) != n_samples:
        parser.error('The --labels argument requires one label per sample.')
    if len(set(labels)) != len(labels):
        parser.error('Sample labels must be unique: %s. Use --labels to set them.' % ','.join(labels))
    delimiter = args_dict['out_args']['delimiter']
    for l in labels:
        if any(d in l for d in delimiter) or ':' in l:
            parser.error('Sample label %s must not contain annotation delimiters or ":".' % l)
    args_dict['labels'] = labels

    # Convert case of fields
    args_dict['barcode_field'] = args_dict['barcode_field'].upper()
    args_dict['shared_field'] = args_dict['shared_field'].upper()

    # Move log file into output arguments
    args_dict['out_args']['log_file'] = args_dict.pop('contam_log')

    findContaminants(**args_dict)
