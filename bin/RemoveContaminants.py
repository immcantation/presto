#!/usr/bin/env python3
"""
Identifies reads shared across samples by UMI barcode and sequence
"""
# Info
__author__ = 'Gisela Gabernet'
from presto import __version__, __date__

# Imports
import heapq
import io
import os
import pickle
import shutil
import tempfile
import zlib
import multiprocessing as mp
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
default_contam_shards = 16


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


def hammingDist(seq_1, seq_2, ignore_chars=default_missing_chars):
    """
    Calculates the Hamming distance between two sequences of equal length

    Arguments:
      seq_1 : first sequence as an uppercase string.
      seq_2 : second sequence as an uppercase string, of the same length.
      ignore_chars : set of missing characters; positions holding one in either
                     sequence are not counted.

    Returns:
      int : the number of mismatched positions.
    """
    return sum(1 for a, b in zip(seq_1, seq_2)
               if a != b and a not in ignore_chars and b not in ignore_chars)


def writeGroupDistances(handle, tag, sample_reads, hist=None, ignore_chars=default_missing_chars):
    """
    Writes the pairwise Hamming distances of all read pairs within each UMI and length group

    Read pairs with the same barcode are grouped by read 1 and read 2 length. Within each
    group, one row is written for every pair of read pairs, from the same or different
    samples. Distances are calculated once per pair of distinct sequences.

    Arguments:
      handle : open output handle of the tab-delimited distance table.
      tag : the UMI barcode.
      sample_reads : OrderedDict mapping sample label to a list of
                     (coordinate, read 1 sequence, read 2 sequence) tuples.
      hist : if a dictionary, distance counts are added to it, keyed by
             ('R1' or 'R2', True if both read pairs are from the same sample)
             and then by distance.
      ignore_chars : set of missing characters not counted as mismatches.

    Returns:
      int : the number of rows written.
    """
    len_dict = OrderedDict()
    for label, reads in sample_reads.items():
        for key, seq_1, seq_2 in reads:
            len_dict.setdefault((len(seq_1), len(seq_2)), []).append((label, key, seq_1, seq_2))

    rows = 0
    for (len_1, len_2), entries in len_dict.items():
        # Index distinct read pairs so each distance is only calculated once
        seq_index = {}
        ids = [seq_index.setdefault((e[2], e[3]), len(seq_index)) for e in entries]
        cache = {}
        for i, (label_a, key_a, a_1, a_2) in enumerate(entries):
            for j in range(i + 1, len(entries)):
                label_b, key_b, b_1, b_2 = entries[j]
                pair = (min(ids[i], ids[j]), max(ids[i], ids[j]))
                dist = cache.get(pair, None)
                if dist is None:
                    dist = (hammingDist(a_1, b_1, ignore_chars), hammingDist(a_2, b_2, ignore_chars))
                    cache[pair] = dist
                handle.write('%s\t%i\t%i\t%s\t%s\t%s\t%s\t%i\t%i\n' %
                             (tag, len_1, len_2, label_a, key_a, label_b, key_b, dist[0], dist[1]))
                rows += 1
                if hist is not None:
                    same = label_a == label_b
                    for read, d in (('R1', dist[0]), ('R2', dist[1])):
                        counts = hist.setdefault((read, same), {})
                        counts[d] = counts.get(d, 0) + 1

    return rows


def plotDistanceHistogram(hist, plot_file):
    """
    Plots histograms of the pairwise Hamming distances to a PNG file

    Draws a grid of panels with one column per read and one row per comparison type:
    pairs of read pairs from different samples and from the same sample. All panels
    share a log-scaled count axis and a distance axis.

    Arguments:
      hist : dictionary of distance counts from writeGroupDistances.
      plot_file : output PNG file name.

    Returns:
      str : the plot file name.
    """
    # Imported here so matplotlib is only loaded when a plot is requested
    from matplotlib.figure import Figure

    # Colors from a colorblind-validated categorical palette
    surface, text_1, text_2, muted = '#fcfcfb', '#0b0b0b', '#52514e', '#898781'
    grid, baseline = '#e1e0d9', '#c3c2b7'
    rows = [(False, 'different samples', '#2a78d6'),
            (True, 'same sample', '#eb6834')]
    cols = [('R1', 'Read 1'), ('R2', 'Read 2')]

    fig = Figure(figsize=(10, 6), dpi=150, facecolor=surface)
    axes = fig.subplots(2, 2, sharex=True, sharey=True)
    max_dist = max([d for counts in hist.values() for d in counts] or [0])
    max_count = max([n for counts in hist.values() for n in counts.values()] or [1])

    for r, (same, row_label, color) in enumerate(rows):
        for c, (read, col_label) in enumerate(cols):
            ax = axes[r][c]
            ax.set_facecolor(surface)
            counts = hist.get((read, same), {})
            if counts:
                dists = sorted(counts)
                ax.bar(dists, [counts[d] for d in dists], width=0.9, color=color,
                       linewidth=0, log=True)
                # Label zero-distance pairs across samples, the signature of contamination
                if not same and 0 in counts:
                    ax.annotate('{:,} pairs at distance 0'.format(counts[0]), xy=(0, counts[0]),
                                xytext=(6, 0), textcoords='offset points',
                                va='center', ha='left', color=text_1, fontsize=8)
            else:
                ax.text(0.5, 0.5, 'No pairs', transform=ax.transAxes,
                        ha='center', va='center', color=text_2, fontsize=9)
            ax.set_title('%s, %s' % (col_label, row_label), loc='left',
                         color=text_1, fontsize=10)
            ax.set_yscale('log')
            ax.grid(axis='y', color=grid, linewidth=0.6)
            ax.set_axisbelow(True)
            ax.tick_params(colors=muted, labelcolor=text_2, labelsize=8)
            ax.tick_params(axis='y', which='minor', left=False)
            for side in ('top', 'right', 'left'):
                ax.spines[side].set_visible(False)
            ax.spines['bottom'].set_color(baseline)
            if r == len(rows) - 1:
                ax.set_xlabel('Hamming distance (mismatched positions)', color=text_2, fontsize=9)
            if c == 0:
                ax.set_ylabel('Pairs of read pairs', color=text_2, fontsize=9)

    axes[0][0].set_xlim(-1, max_dist + 1)
    axes[0][0].set_ylim(0.8, max_count * 3)
    fig.suptitle('Pairwise Hamming distances within UMI and length groups (log scale)',
                 x=0.01, ha='left', color=text_1, fontsize=12)
    fig.tight_layout()
    fig.savefig(plot_file, facecolor=surface)

    return plot_file


def countSharedPairs(flagged, labels):
    """
    Counts shared read pairs between every pair of samples

    Arguments:
      flagged : dictionary mapping each sample label to a dictionary of flagged read
                coordinates and the list of other samples each read matched.
      labels : ordered list of sample labels.

    Returns:
      list : a square matrix as a list of lists, where row i, column j holds the number
             of read pairs of sample i that matched a read pair in sample j. Matching read
             pairs are removed from both samples, so cell (j, i) holds the matching read
             pairs of sample j. A read pair matching several samples is counted once for
             each of them.
    """
    index = {label: i for i, label in enumerate(labels)}
    matrix = [[0] * len(labels) for _ in labels]
    for label in labels:
        row = matrix[index[label]]
        for others in flagged.get(label, {}).values():
            for other in others:
                row[index[other]] += 1

    return matrix


def plotSharedHeatmap(matrix, labels, plot_file, totals=None, max_annotate=20):
    """
    Plots a heatmap of the read pairs shared between every pair of samples to a PNG file

    Arguments:
      matrix : square count matrix from countSharedPairs.
      labels : ordered list of sample labels.
      plot_file : output PNG file name.
      totals : optional list of the total number of reads of each sample, written in the
               diagonal cells. Diagonal cells are not part of the color scale.
      max_annotate : the maximum number of samples for which counts are written in the cells.

    Returns:
      str : the plot file name.
    """
    # Imported here so matplotlib is only loaded when a plot is requested
    from matplotlib.figure import Figure
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
    from matplotlib.patches import Rectangle

    # Sequential single-hue ramp and neutral colors from a validated palette
    ramp = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
    surface, text_1, text_2, muted, zero_color = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#f0efec'
    baseline = '#c3c2b7'
    cmap = LinearSegmentedColormap.from_list('shared', ramp)
    cmap.set_bad(zero_color)

    # Keep the matrix at most 10 inches wide; cells shrink as samples are added
    n = len(labels)
    side = min(max(4.5, 0.75 * n), 10.0)
    cell_pt = side * 72.0 / n
    # Sample labels never go below 6 points; label every k-th sample if cells are smaller
    label_size = max(6.0, min(10.0, cell_pt * 0.75))
    label_step = max(1, int(-(-label_size * 1.15 // cell_pt)))
    gap = min(2.0, cell_pt * 0.15)
    fig = Figure(figsize=(side + 3.2, side + 2.2), dpi=150, facecolor=surface)
    ax = fig.subplots()
    ax.set_facecolor(surface)

    # Mask zero counts and the diagonal so they are not drawn on the color scale
    values = [[float('nan') if i == j or matrix[i][j] == 0 else matrix[i][j]
               for j in range(n)] for i in range(n)]
    positive = [matrix[i][j] for i in range(n) for j in range(n) if i != j and matrix[i][j] > 0]
    norm = LogNorm(vmin=1, vmax=max(positive)) if positive and max(positive) > 1 else None
    image = ax.imshow(values, cmap=cmap, norm=norm, aspect='equal',
                      vmin=None if norm else 0, vmax=None if norm else 1)

    # Diagonal cells are drawn as blank surface with a hairline border, off the color scale
    for i in range(n):
        ax.add_patch(Rectangle((i - 0.5, i - 0.5), 1, 1, facecolor=surface, edgecolor='none'))
        ax.add_patch(Rectangle((i - 0.45, i - 0.45), 0.9, 0.9, facecolor='none',
                               edgecolor=baseline, linewidth=min(0.8, cell_pt * 0.08)))

    # Cell borders as surface-colored gaps
    ax.set_xticks([x - 0.5 for x in range(1, n)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n)], minor=True)
    ax.grid(which='minor', color=surface, linewidth=gap)
    ax.tick_params(which='minor', length=0)

    def _cell_text(i, j, text, color):
        # Shrink long numbers so they fit within the cell width
        size = min(10.0, 0.8 * cell_pt / (0.6 * len(text)))
        if size >= 5:
            ax.text(j, i, text, ha='center', va='center', fontsize=size, color=color)

    if n <= max_annotate:
        for i in range(n):
            for j in range(n):
                if i == j:
                    if totals is not None:
                        _cell_text(i, j, '{:,}'.format(totals[i]), text_2)
                    continue
                count = matrix[i][j]
                if count > 0 and norm is not None:
                    dark = norm(count) > 0.5
                else:
                    dark = count > 0
                _cell_text(i, j, '{:,}'.format(count), '#ffffff' if dark else text_1)

    ticks = list(range(0, n, label_step))
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    # Vertical labels for many samples, so neighbouring labels never overlap
    rotation = 90 if n > 12 else (45 if n > 3 else 0)
    ax.set_xticklabels([labels[t] for t in ticks], rotation=rotation,
                       ha='right' if rotation else 'center',
                       va='center' if rotation == 90 else 'top', rotation_mode='anchor')
    ax.set_yticklabels([labels[t] for t in ticks])
    ax.tick_params(colors=muted, labelcolor=text_2, labelsize=label_size, length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    title = 'Pairwise sample contaminants' if positive else \
            'Pairwise sample contaminants: none found'
    ax.set_title(title, loc='left', color=text_1, fontsize=14, pad=10)

    if norm is not None:
        cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03, shrink=0.8)
        cbar.set_label('Number of contaminants removed', color=text_2, fontsize=11)
    elif positive:
        cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03, shrink=0.8, ticks=[1])
        cbar.set_label('Number of contaminants removed', color=text_2, fontsize=11)
    else:
        cbar = None
    if cbar is not None:
        if norm is not None:
            # Label powers of ten and the maximum as plain numbers
            vmax = int(norm.vmax)
            ticks = [10 ** k for k in range(len(str(vmax))) if 2 * 10 ** k <= vmax] + [vmax]
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(['{:,}'.format(t) for t in ticks])
        cbar.ax.minorticks_off()
        cbar.outline.set_visible(False)
        cbar.ax.tick_params(colors=muted, labelcolor=text_2, labelsize=10)

    fig.tight_layout()
    fig.savefig(plot_file, facecolor=surface)

    return plot_file


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


def runTasks(func, tasks, nproc=1, initializer=None, initargs=(), chunksize=1):
    """
    Applies a function to tasks in order, optionally in a pool of processes

    Arguments:
      func : a module level function taking a single task argument.
      tasks : an iterable of tasks.
      nproc : the number of processes; if 1 or less, tasks run in the current process.
      initializer : optional function called once in each process before its tasks.
      initargs : arguments for initializer.
      chunksize : the number of tasks sent to a process at a time.

    Returns:
      iter : the results, in the order of the tasks.
    """
    if nproc <= 1:
        if initializer is not None:  initializer(*initargs)
        for task in tasks:
            yield func(task)
    else:
        with mp.Pool(nproc, initializer=initializer, initargs=initargs) as pool:
            for result in pool.imap(func, tasks, chunksize):
                yield result


def _barcodeTask(task):
    """
    Pass 1 task: collects the barcodes of one sample
    """
    seq_file_1, seq_file_2, field, delimiter = task
    return collectSampleBarcodes(seq_file_1, seq_file_2, field=field, delimiter=delimiter)


# Candidate barcodes shared with pass 2 tasks, set once per process
_candidates = None

def _setCandidates(candidates):
    """
    Sets the candidate barcodes used by pass 2 tasks in the current process
    """
    global _candidates
    _candidates = candidates


def getShard(tag, shards):
    """
    Assigns a UMI barcode to a shard

    Uses a checksum instead of the built-in hash, which differs between processes.

    Arguments:
      tag : the UMI barcode.
      shards : the number of shards.

    Returns:
      int : the shard index, from 0 to shards - 1.
    """
    return zlib.crc32(tag.encode('utf-8')) % shards


def _shardFile(tmp_dir, shard, sample):
    """
    Returns the file name holding the reads of one sample in one shard
    """
    return os.path.join(tmp_dir, 'shard%i_sample%i.pkl' % (shard, sample))


def _distFile(tmp_dir, shard):
    """
    Returns the file name holding the distance table rows of one shard
    """
    return os.path.join(tmp_dir, 'shard%i_dist.tab' % shard)


def _readTask(task):
    """
    Pass 2 task: collects the reads of candidate barcodes in one sample and writes them to shards

    Each shard file holds a list of (barcode, rank, reads) tuples, where rank is the position
    of the barcode within the sample, used to restore the barcode order after comparison.
    Only the conflict count is returned, so the reads are never sent between processes.
    """
    sample, seq_file_1, seq_file_2, field, coord_type, delimiter, tmp_dir, shards = task
    groups, conflicts = collectSampleReads(seq_file_1, seq_file_2, candidates=_candidates,
                                           field=field, coord_type=coord_type,
                                           delimiter=delimiter)
    shard_groups = {}
    for rank, (tag, reads) in enumerate(groups.items()):
        shard_groups.setdefault(getShard(tag, shards), []).append((tag, rank, reads))
    del groups
    for shard, entries in shard_groups.items():
        with open(_shardFile(tmp_dir, shard, sample), 'wb') as handle:
            pickle.dump(entries, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return conflicts


def compareShard(task):
    """
    Compares the read pairs of each barcode of one shard found in multiple samples

    Each barcode is compared across all samples containing it. Barcodes are ordered by the
    first sample containing them and then by their position in that sample, which is the
    order of the barcodes when all samples are read sequentially.

    Arguments:
      task : a tuple of (shard, tmp_dir, labels, max_error, distances), where distances is
             True to also write the pairwise Hamming distances to the shard distance file,
             each row prefixed by the barcode order.

    Returns:
      tuple : (multi_count, records, dist_rows, hist), where multi_count is the number of
              barcodes found in multiple samples, records is a list of
              (order, barcode, sharing samples, shared reads by sample, group sizes) tuples for
              the barcodes with shared read pairs, dist_rows is the number of distance rows,
              and hist is the distance histogram counts, or None when distances is False.
    """
    shard, tmp_dir, labels, max_error, distances = task

    # Combine the reads of each barcode across samples
    barcode_reads, order = {}, {}
    for i, label in enumerate(labels):
        path = _shardFile(tmp_dir, shard, i)
        if not os.path.exists(path):
            continue
        with open(path, 'rb') as handle:
            entries = pickle.load(handle)
        for tag, rank, reads in entries:
            sample_reads = barcode_reads.get(tag, None)
            if sample_reads is None:
                sample_reads = barcode_reads[tag] = OrderedDict()
                order[tag] = (i, rank)
            sample_reads[label] = reads
        del entries
    multi = sorted((t for t, r in barcode_reads.items() if len(r) > 1), key=order.get)
    barcode_reads = {t: barcode_reads[t] for t in multi}

    records = []
    dist_rows, hist, dist_handle = 0, None, None
    if distances:
        hist = {}
        dist_handle = open(_distFile(tmp_dir, shard), 'w')
    try:
        for tag in multi:
            sample_reads = barcode_reads.pop(tag)
            if dist_handle is not None:
                buffer = io.StringIO()
                dist_rows += writeGroupDistances(buffer, tag, sample_reads, hist=hist)
                prefix = '%i\t%i\t' % order[tag]
                for line in buffer.getvalue().splitlines(True):
                    dist_handle.write(prefix + line)
            shared = findSharedReads(sample_reads, max_error=max_error)
            sharing = [s for s in sample_reads if shared[s]]
            if sharing:
                records.append((order[tag], tag, sharing, {s: shared[s] for s in sharing},
                                [len(sample_reads[s]) for s in sharing]))
    finally:
        if dist_handle is not None:  dist_handle.close()

    return len(multi), records, dist_rows, hist


def _writeTask(task):
    """
    Pass 3 task: writes the pass and fail files of one input file
    """
    seq_file, flagged, shared_field, coord_type, out_args = task
    return writeSampleReads(seq_file, flagged, shared_field=shared_field,
                            coord_type=coord_type, out_args=out_args)


def removeContaminants(seq_files_1, seq_files_2, labels=None,
                     barcode_field=default_barcode_field, coord_type=default_coord,
                     max_error=default_contam_max_error, dist_file=None,
                     shared_field=default_shared_field, out_args=default_out_args, nproc=1,
                     shards=None):
    """
    Identifies reads shared across samples by UMI barcode and sequence

    A read pair is shared when another sample contains a read pair with the same UMI
    barcode whose read 1 and read 2 sequences both match it. Shared read pairs are
    removed from every sample containing them.

    To limit memory usage the input files are read in three passes. The first pass
    collects only the barcodes of each sample. The second keeps only the read pairs whose
    barcode occurs in multiple samples. The third writes the output files.

    The second pass writes the kept read pairs to temporary shard files, split by barcode,
    in a directory created in the output directory, or in the directory of the first read 1
    file. Each shard is then compared separately, so all read pairs of a barcode are compared
    together while the main process never holds the reads. The directory is removed at the
    end of the run.

    With multiple processes, the samples are read in parallel in the first two passes, the
    shards are compared in parallel, and the files are written in parallel in the third pass.
    Results are restored to the barcode order of a sequential run, so they do not depend on
    the number of processes or shards.

    Arguments:
      seq_files_1 : list of read 1 sequence files, one per sample.
      seq_files_2 : list of read 2 sequence files, in the same sample order.
      labels : list of sample labels; if None use the read 1 file names.
      barcode_field : the annotation field containing the UMI barcode.
      coord_type : the sequence header format used to match mates.
      max_error : the maximum read 1 and read 2 mismatch rate allowed between two
                  read pairs to consider them shared.
      dist_file : if specified, write the pairwise Hamming distances of all read pairs
                  within each UMI and length group to this tab-delimited file. Only barcodes
                  with read pairs in multiple samples are included. A histogram of the
                  distances is plotted to a PNG file with the same name and a .png extension.
      shared_field : the annotation field listing the samples a flagged read is shared with.
      out_args : common output argument dictionary from parseCommonArgs. The log_file entry
                 defines the contaminant log; if None it is written to contam_log.txt in the
                 output directory, or in the directory of the first read 1 file. A heatmap
                 of the read pairs shared between samples is plotted next to it, with the
                 same name followed by _heatmap.png.
      nproc : the number of processes to use.
      shards : the number of shards the candidate read pairs are split into; if None use
               the larger of 16 and four times nproc.

    Returns:
      list : a list of tuples holding the (read 1, read 2) pass file names for each sample.
    """
    delimiter = out_args['delimiter']
    if labels is None:
        labels = [getSampleLabel(f) for f in seq_files_1]
    sample_files = OrderedDict(zip(labels, zip(seq_files_1, seq_files_2)))

    # Print parameter info
    log = OrderedDict()
    log['START'] = 'RemoveContaminants'
    log['FILES1'] = ','.join(os.path.basename(f) for f in seq_files_1)
    log['FILES2'] = ','.join(os.path.basename(f) for f in seq_files_2)
    log['SAMPLES'] = ','.join(labels)
    log['BARCODE_FIELD'] = barcode_field
    log['COORD_TYPE'] = coord_type
    log['MAX_ERROR'] = max_error
    log['DIST_FILE'] = os.path.basename(dist_file) if dist_file is not None else None
    log['NPROC'] = nproc
    printLog(log)

    # Pass 1: find barcodes annotated in multiple samples
    start_time = time()
    printMessage('Collecting barcodes', start_time=start_time)
    barcode_count = {}
    tasks = [(f1, f2, barcode_field, delimiter) for f1, f2 in sample_files.values()]
    for barcodes in runTasks(_barcodeTask, tasks, nproc=min(nproc, len(tasks))):
        for tag in barcodes:
            barcode_count[tag] = barcode_count.get(tag, 0) + 1
    candidates = {t for t, n in barcode_count.items() if n > 1}
    del barcode_count
    printMessage('Done', start_time=start_time, end=True)

    # Candidate reads are passed between passes in shard files, split by barcode, so the
    # main process never holds them
    if shards is None:
        shards = max(default_contam_shards, 4 * nproc)
    tmp_parent = out_args['out_dir'] if out_args['out_dir'] is not None \
                 else (os.path.dirname(seq_files_1[0]) or '.')
    os.makedirs(tmp_parent, exist_ok=True)
    tmp_dir = tempfile.mkdtemp(prefix='contam_tmp_', dir=tmp_parent)
    try:
        # Pass 2: group paired reads of candidate barcodes into shards
        start_time = time()
        printMessage('Collecting reads', start_time=start_time)
        tasks = [(i, f1, f2, barcode_field, coord_type, delimiter, tmp_dir, shards)
                 for i, (f1, f2) in enumerate(sample_files.values())]
        results = runTasks(_readTask, tasks, nproc=min(nproc, len(tasks)),
                           initializer=_setCandidates, initargs=(candidates,))
        for label, conflicts in zip(sample_files, results):
            if conflicts > 0:
                printWarning('%i read pairs in sample %s have different barcodes in read 1 and '
                             'read 2. The read 1 barcode was used.' % (conflicts, label))
        _setCandidates(None)
        del candidates
        printMessage('Done', start_time=start_time, end=True)

        # Compare reads sharing a barcode, one shard per task
        start_time = time()
        printMessage('Comparing reads', start_time=start_time)
        log_file = out_args['log_file']
        if log_file is None:
            log_dir = out_args['out_dir'] if out_args['out_dir'] is not None \
                      else os.path.dirname(seq_files_1[0])
            log_file = os.path.join(log_dir, 'contam_log.txt')
        log_dir = os.path.dirname(log_file)
        if log_dir:  os.makedirs(log_dir, exist_ok=True)

        distances = dist_file is not None
        dist_hist = {} if distances else None
        multi_count = dist_count = 0
        records = []
        tasks = [(shard, tmp_dir, labels, max_error, distances) for shard in range(shards)]
        for shard_multi, shard_records, dist_rows, hist in \
                runTasks(compareShard, tasks, nproc=min(nproc, shards)):
            multi_count += shard_multi
            records.extend(shard_records)
            if distances:
                dist_count += dist_rows
                for key, counts in hist.items():
                    merged = dist_hist.setdefault(key, {})
                    for d, c in counts.items():
                        merged[d] = merged.get(d, 0) + c

        # Restore the barcode order of a sequential run
        records.sort(key=lambda r: r[0])
        flagged = {label: {} for label in labels}
        with open(log_file, 'w') as log_handle:
            for _, tag, sharing, shared, group_counts in records:
                for s in sharing:
                    flagged[s].update(shared[s])

                # Write a log record with the read counts of each contaminated sample
                rec_log = OrderedDict()
                rec_log['BARCODE'] = tag
                rec_log['SAMPLES'] = ','.join(sharing)
                rec_log['SEQCOUNT'] = ','.join(str(len(shared[s])) for s in sharing)
                rec_log['GROUPCOUNT'] = ','.join(str(n) for n in group_counts)
                printLog(rec_log, handle=log_handle)
        contam_count = len(records)
        del records

        # Merge the sorted distance rows of all shards
        plot_file = None
        if distances:
            plot_file = os.path.splitext(dist_file)[0] + '.png'
            dist_dir = os.path.dirname(dist_file)
            if dist_dir:  os.makedirs(dist_dir, exist_ok=True)
            shard_handles = [open(_distFile(tmp_dir, shard)) for shard in range(shards)]
            try:
                with open(dist_file, 'w') as dist_handle:
                    dist_handle.write('BARCODE\tR1_LENGTH\tR2_LENGTH\tSAMPLE1\tID1\tSAMPLE2'
                                      '\tID2\tR1_DIST\tR2_DIST\n')
                    def _order(line):
                        i, rank, _ = line.split('\t', 2)
                        return int(i), int(rank)
                    for line in heapq.merge(*shard_handles, key=_order):
                        dist_handle.write(line.split('\t', 2)[2])
            finally:
                for handle in shard_handles:  handle.close()
            plotDistanceHistogram(dist_hist, plot_file)

        printMessage('Done', start_time=start_time, end=True)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    # Pass 3: write output files
    start_time = time()
    printMessage('Writing output', start_time=start_time)
    out_files = []
    out_log = OrderedDict()
    read_totals = []
    tasks = [(seq_file, flagged[label], shared_field, coord_type, out_args)
             for label, files in sample_files.items() for seq_file in files]
    results = list(runTasks(_writeTask, tasks, nproc=min(nproc, len(tasks))))
    for k, (label, files) in enumerate(sample_files.items()):
        pass_names = []
        pass_total = fail_total = 0
        for i, seq_file in enumerate(files):
            name, pass_count, fail_count = results[2 * k + i]
            pass_names.append(name)
            pass_total += pass_count
            fail_total += fail_count
            # Read 1 count, one read per read pair, for the heatmap diagonal
            if i == 0:  read_totals.append(pass_count + fail_count)
        out_files.append(tuple(pass_names))
        out_log['%s_PASS' % label] = pass_total
        out_log['%s_FAIL' % label] = fail_total
    printMessage('Done', start_time=start_time, end=True)

    # Plot read pairs shared between each pair of samples
    heatmap_file = os.path.splitext(log_file)[0] + '_heatmap.png'
    plotSharedHeatmap(countSharedPairs(flagged, labels), labels, heatmap_file,
                      totals=read_totals)

    # Print log
    log = OrderedDict()
    log['OUTPUT'] = ','.join(os.path.basename(f) for p in out_files for f in p)
    log['MULTISAMPLE'] = multi_count
    log['CONTAMINATED'] = contam_count
    log['LOG'] = os.path.basename(log_file)
    log['HEATMAP'] = os.path.basename(heatmap_file)
    if dist_file is not None:
        log['DISTANCES'] = dist_count
        log['DIST_PLOT'] = os.path.basename(plot_file)
    log.update(out_log)
    log['END'] = 'RemoveContaminants'
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
             multiple samples, and the third writes the output. The kept reads are stored
             in a temporary directory within the output directory, split into shards by
             barcode, and each shard is compared separately. With --nproc, samples are
             read, shards are compared and files are written in parallel.

             output files:
                 contam-pass
                     reads not shared with another sample, one file per input file.
                 contam-fail
                     reads shared with another sample, one file per input file.
                     Only written with --failed.
                 <log>_heatmap.png
                     heatmap of the read pairs shared between every pair of samples,
                     written next to the contaminant log. The cell in row i and column j
                     is the number of read pairs of sample i that matched a read pair in
                     sample j. Matching read pairs are removed from both samples, and
                     cell (j, i) holds those of sample j. The diagonal shows the total
                     number of read 1 reads of each sample and is not part of the color
                     scale.

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

             distance table fields:
                 The --distlog table contains one row for every pair of read pairs within
                 a UMI and length group, that is, read pairs with the same barcode, read 1
                 length and read 2 length, from the same or different samples. Only
                 barcodes with read pairs in multiple samples are included.
                 BARCODE, R1_LENGTH, R2_LENGTH
                     the UMI barcode and the read lengths defining the group.
                 SAMPLE1, ID1, SAMPLE2, ID2
                     the sample label and read identifier of the two read pairs.
                 R1_DIST, R2_DIST
                     the read 1 and read 2 Hamming distances. Positions with N or gap
                     characters in either read are not counted.
                 A histogram of the distances, split into pairs from the same sample and
                 from different samples, is plotted to a PNG file with the same name as
                 the table and a .png extension.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(seq_paired=True, out_file=False, log=False,
                                                        multiproc=True)],
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
    group_contam.add_argument('--distlog', action='store', dest='dist_file', default=None,
                              help='''Specify to write the pairwise Hamming distances of all read
                                   pairs within each UMI and length group to this tab-delimited
                                   file. Only barcodes with read pairs in multiple samples are
                                   included. The file grows with the square of the group sizes.
                                   A histogram of the distances is also plotted to a PNG file
                                   with the same name and a .png extension.''')
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

    removeContaminants(**args_dict)
