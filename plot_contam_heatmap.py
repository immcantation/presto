#!/usr/bin/env python3
"""
Plots a heatmap of the read pairs shared between samples from RemoveContaminants fail files

Also plots a bar plot of the read pairs removed from each sample, and a histogram of the
number of samples each contaminant is found in, next to the heatmap with the same name
followed by _barplot.png and _histogram.png. With --runlog, the heatmap is also plotted
with each cell as a percentage of the read pairs of the row sample, followed by _percent.png.

Each read in a contam-fail file holds the list of the other samples it matched in its shared
field, so the counts are the same as in the RemoveContaminants heatmap. Samples are clustered
by their contamination profiles, using the same clustering and plot as RemoveContaminants.py,
which is imported from the bin directory next to this script.

Takes the read 1 contam-fail file of each sample. The sample label of each file is its file
name without the _contam-fail suffix and extension, which is the default RemoveContaminants
label of its read 1 input file.

An optional samplesheet gives the batch of each sample, shown as a colored strip along the
rows and columns. Each sample label is matched to the sample column, first as is and then
without a trailing _R1, _1 or .R1.

A contaminant in the histogram is a unique combination of UMI barcode and read 1 sequence.
Read pairs matched with mismatches or missing characters form separate contaminants, so some
contaminants can be found in a single sample.

Usage:
  plot_contam_heatmap.py *_R1_contam-fail.fastq [-o heatmap.png] [--runlog run.log] [--sf SHARED]
                         [--bf BARCODE] [--samplesheet samples.tsv] [--samplecol sample]
                         [--batchcol batch]
"""

import os
import shutil
import sys
from argparse import ArgumentParser
from collections import Counter

# RemoveContaminants.py is a script, not a module of the presto package, so its directory is
# added to the import path: the installed script on the PATH, or the presto bin directory
# next to this script
_installed = shutil.which('RemoveContaminants.py')
if _installed is not None:
    # Searched last, so other files in that directory never shadow installed modules
    sys.path.append(os.path.dirname(os.path.realpath(_installed)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin'))

from presto.Annotation import parseAnnotation
from presto.Defaults import default_delimiter
try:
    from RemoveContaminants import batchColors, iterSeqFile, matchBatch, plotSharedHeatmap, \
        readBatches
except ModuleNotFoundError:
    sys.exit('RemoveContaminants.py was not found. Install presto so it is on the PATH, or '
             'keep this script next to the presto bin directory.')


def getFailLabel(fail_file):
    """
    Returns the default sample label of a contam-fail file

    Arguments:
      fail_file : contam-fail file name.

    Returns:
      str : the file name without directory, .gz extension, file extension and
            _contam-fail suffix.
    """
    name = os.path.basename(fail_file)
    if name.endswith('.gz'):  name = name[:-3]
    name = os.path.splitext(name)[0]
    if name.endswith('_contam-fail'):  name = name[:-len('_contam-fail')]

    return name


def countShared(fail_file, shared_field='SHARED', barcode_field='BARCODE',
                delimiter=default_delimiter):
    """
    Counts the reads of one contam-fail file shared with each other sample

    Arguments:
      fail_file : contam-fail file name.
      shared_field : the annotation field listing the other samples a read matched.
      barcode_field : the annotation field holding the UMI barcode.
      delimiter : the (field, value, value list) annotation delimiters.

    Returns:
      tuple : (counts, reads, contaminants), where counts maps each other sample label to its
              number of shared reads, reads is the number of reads in the file, and
              contaminants is the set of hashed (barcode, sequence) pairs of the reads with
              a barcode.
    """
    counts, reads, contaminants = {}, 0, set()
    for title, seq, _ in iterSeqFile(fail_file):
        reads += 1
        ann = parseAnnotation(title, delimiter=delimiter)
        shared = ann.get(shared_field, None)
        if shared is not None:
            for other in shared.split(delimiter[2]):
                counts[other] = counts.get(other, 0) + 1
        barcode = ann.get(barcode_field, None)
        if barcode is not None:
            # Hashed to keep memory low with millions of contaminants
            contaminants.add(hash((barcode, seq.upper())))

    return counts, reads, contaminants


def readTotals(run_log):
    """
    Reads the number of read pairs of each sample from the RemoveContaminants standard output

    The pass and fail counts in the final log record include both reads, so each sum is
    halved to count read pairs, as in the RemoveContaminants heatmap diagonal.

    Arguments:
      run_log : file holding the standard output of RemoveContaminants.py.

    Returns:
      dict : maps each sample label to its number of read pairs.
    """
    counts = {}
    with open(run_log) as handle:
        for line in handle:
            field, sep, value = line.strip().partition('> ')
            value = value.strip()
            if not sep or not value.isdigit():
                continue
            for suffix in ('_PASS', '_FAIL'):
                if field.endswith(suffix):
                    label = field[:-len(suffix)]
                    counts[label] = counts.get(label, 0) + int(value)

    return {label: n // 2 for label, n in counts.items()}


# Surface, text and grid colors of the RemoveContaminants heatmap, and the bar color of the
# validated categorical palette
_surface, _text_1, _text_2, _muted, _grid = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e8e7e2'
_bar = '#2a78d6'


def _styleAxes(ax, value_axis):
    """
    Applies the shared style to the axes of a bar plot

    Arguments:
      ax : matplotlib axes.
      value_axis : 'x' or 'y', the axis of the bar values, which gets light gridlines.
    """
    ax.set_facecolor(_surface)
    for name, spine in ax.spines.items():
        spine.set_visible(name == ('left' if value_axis == 'x' else 'bottom'))
        spine.set_color(_muted)
    ax.grid(axis=value_axis, color=_grid, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=_muted, labelcolor=_text_2, labelsize=9, length=0)


def plotRemovedBars(labels, removed, plot_file, batches=None):
    """
    Plots the number of read pairs removed from each sample as horizontal bars

    Samples are sorted from most to fewest removed read pairs, on a log scale so samples
    with few removed read pairs remain visible. Samples without removed read pairs have no bar.

    Arguments:
      labels : list of sample labels.
      removed : list of the number of read pairs removed from each sample.
      plot_file : output PNG file name.
      batches : optional list of the batch of each sample, used to color the bars.

    Returns:
      str : the plot file name.
    """
    from matplotlib.figure import Figure
    from matplotlib.patches import Patch

    order = sorted(range(len(labels)), key=lambda i: (-removed[i], labels[i]))
    n = len(labels)
    colors, batch_colors = [_bar] * n, None
    if batches is not None:
        batches = [b if b not in (None, '') else 'NA' for b in batches]
        batch_colors = batchColors(batches)
        colors = [batch_colors[batches[i]] for i in order]

    row = 0.32 if n <= 40 else max(0.12, 12.0 / n)
    fig = Figure(figsize=(8.5 + (2.0 if batch_colors else 0), 1.6 + row * n), dpi=150,
                 facecolor=_surface)
    ax = fig.subplots()
    values = [removed[i] for i in order]
    ax.barh(range(n), values, height=0.7, color=colors)
    ax.set_yticks(range(n))
    ax.set_yticklabels([labels[i] for i in order], fontsize=min(9, max(5, row * 30)))
    ax.set_ylim(n - 0.5, -0.5)
    # Bars start at 0 and are clipped at the left edge, just below one read pair
    ax.set_xscale('log')
    ax.set_xlim(0.8, max(values + [1]) * 4)
    _styleAxes(ax, 'x')
    ax.xaxis.set_major_formatter('{x:,.0f}')
    ax.xaxis.set_minor_formatter('')
    ax.set_xlabel('Read pairs removed (log scale)', color=_text_2, fontsize=10)
    ax.set_title('Contaminant read pairs removed per sample', loc='left', color=_text_1,
                 fontsize=14, pad=10)

    # Values at the bar ends, while they fit
    if n <= 40:
        for y, v in enumerate(values):
            ax.text(max(v, 0.8), y, ' {:,}'.format(v), va='center', ha='left', fontsize=8,
                    color=_text_2)

    if batch_colors is not None:
        handles = [Patch(facecolor=c, edgecolor='none', label=b) for b, c in batch_colors.items()]
        legend = ax.legend(handles=handles, title='Batch', loc='upper left',
                           bbox_to_anchor=(1.02, 1), frameon=False, fontsize=9,
                           title_fontsize=10, handlelength=1.0, handleheight=1.0)
        legend.get_title().set_color(_text_2)
        for text in legend.get_texts():
            text.set_color(_text_2)

    fig.tight_layout()
    fig.savefig(plot_file, facecolor=_surface)

    return plot_file


def plotSampleHistogram(sample_counts, plot_file):
    """
    Plots the number of contaminants found in each number of samples

    Arguments:
      sample_counts : dictionary mapping a number of samples to its number of contaminants.
      plot_file : output PNG file name.

    Returns:
      str : the plot file name.
    """
    from matplotlib.figure import Figure
    from matplotlib.ticker import MaxNLocator

    fig = Figure(figsize=(8.5, 4.8), dpi=150, facecolor=_surface)
    ax = fig.subplots()
    xs = list(range(1, max(sample_counts, default=1) + 1))
    ys = [sample_counts.get(x, 0) for x in xs]
    ax.bar(xs, ys, width=0.7, color=_bar)
    _styleAxes(ax, 'y')

    # Log scale when counts span more than two orders of magnitude
    positive = [y for y in ys if y > 0]
    if positive and max(positive) / min(positive) > 100:
        ax.set_yscale('log')
        ax.set_ylim(0.7, max(positive) * 3)
    else:
        ax.set_ylim(0, max(positive + [1]) * 1.15)
    ax.yaxis.set_major_formatter('{x:,.0f}')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlim(0.4, len(xs) + 0.6)

    if len(xs) <= 25:
        for x, y in zip(xs, ys):
            if y > 0:
                ax.text(x, y, '{:,}'.format(y), ha='center', va='bottom', fontsize=8,
                        color=_text_2)

    ax.set_xlabel('Number of samples containing the contaminant', color=_text_2, fontsize=10)
    ax.set_ylabel('Contaminants', color=_text_2, fontsize=10)
    ax.set_title('Samples per contaminant', loc='left', color=_text_1, fontsize=14, pad=22)
    ax.text(0, 1.02, '{:,} unique contaminants (UMI barcode and read 1 sequence)'.format(
                sum(sample_counts.values())),
            transform=ax.transAxes, ha='left', va='bottom', color=_text_2, fontsize=9)

    fig.tight_layout()
    fig.savefig(plot_file, facecolor=_surface)

    return plot_file


def getArgParser():
    parser = ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('fail_files', nargs='+',
                        help='The read 1 contam-fail file of each sample written by '
                             'RemoveContaminants.py, optionally gzip compressed.')
    parser.add_argument('-o', '--out', dest='plot_file', default='contam_heatmap.png',
                        help='Output PNG file.')
    parser.add_argument('--runlog', default=None,
                        help='File holding the standard output of RemoveContaminants.py, used '
                             'to write the number of read pairs of each sample in the diagonal '
                             'and to plot the heatmap as percentages of each sample\'s reads.')
    parser.add_argument('--sf', dest='shared_field', default='SHARED',
                        help='The annotation field listing the other samples a read matched.')
    parser.add_argument('--bf', dest='barcode_field', default='BARCODE',
                        help='The annotation field holding the UMI barcode, used to define '
                             'contaminants in the histogram.')
    parser.add_argument('--samplesheet', default=None,
                        help='CSV (.csv) or tab-delimited samplesheet with the batch of each '
                             'sample, shown as a colored strip on the heatmap and as bar '
                             'colors in the bar plot.')
    parser.add_argument('--samplecol', dest='sample_col', default='sample',
                        help='The samplesheet column holding the sample names.')
    parser.add_argument('--batchcol', dest='batch_col', default='batch',
                        help='The samplesheet column holding the batches.')

    return parser


if __name__ == '__main__':
    args = getArgParser().parse_args()
    args.shared_field = args.shared_field.upper()
    args.barcode_field = args.barcode_field.upper()

    labels = [getFailLabel(f) for f in args.fail_files]
    if len(set(labels)) != len(labels):
        sys.exit('Sample labels inferred from the file names must be unique: %s.' %
                 ','.join(labels))

    # Count the reads of each sample shared with every other sample
    index = {label: i for i, label in enumerate(labels)}
    matrix = [[0] * len(labels) for _ in labels]
    removed, contaminant_samples = [], Counter()
    for fail_file, label in zip(args.fail_files, labels):
        counts, reads, contaminants = countShared(fail_file, shared_field=args.shared_field,
                                                  barcode_field=args.barcode_field)
        removed.append(reads)
        if reads and not contaminants:
            sys.exit('No %s field found in %s. Use --bf to set the barcode field.' %
                     (args.barcode_field, fail_file))
        # Each contaminant is counted once per sample containing it
        contaminant_samples.update(contaminants)
        del contaminants
        if reads and not counts:
            sys.exit('No %s field found in %s. Use --sf to set the shared field.' %
                     (args.shared_field, fail_file))
        unknown = sorted(set(counts) - set(index))
        if unknown:
            sys.exit('Samples %s in the %s field of %s do not match any file name. Give the '
                     'read 1 contam-fail file of every sample.' %
                     (','.join(unknown), args.shared_field, fail_file))
        for other, count in counts.items():
            matrix[index[label]][index[other]] += count

    totals = None
    if args.runlog is not None:
        run_totals = readTotals(args.runlog)
        absent = [l for l in labels if l not in run_totals]
        if absent:
            sys.exit('Sample %s is not in the run log %s.' % (absent[0], args.runlog))
        totals = [run_totals[l] for l in labels]

    batches = None
    if args.samplesheet is not None:
        try:
            sheet = readBatches(args.samplesheet, sample_col=args.sample_col,
                                batch_col=args.batch_col)
        except ValueError as e:
            sys.exit(str(e))
        batches = [matchBatch(l, sheet) for l in labels]
        unmatched = [l for l, b in zip(labels, batches) if b is None]
        if len(unmatched) == len(labels):
            sys.exit('No sample matches the %s column of %s, which holds e.g. %s.' %
                     (args.sample_col, args.samplesheet, ','.join(list(sheet)[:3])))
        if unmatched:
            print('Warning: samples without a batch, shown as NA: %s' % ','.join(unmatched))

        # Check that every batch can get its own color before plotting
        try:
            batchColors([b if b is not None else 'NA' for b in batches])
        except ValueError as e:
            sys.exit(str(e))

    plotSharedHeatmap(matrix, labels, args.plot_file, totals=totals, batches=batches)
    base = os.path.splitext(args.plot_file)[0]
    percent_file = None
    if totals is not None:
        percent_file = plotSharedHeatmap(matrix, labels, base + '_percent.png', totals=totals,
                                         batches=batches, percent=True)
    bar_file = plotRemovedBars(labels, removed, base + '_barplot.png', batches=batches)
    hist = Counter(contaminant_samples.values())
    hist_file = plotSampleHistogram(hist, base + '_histogram.png')

    print('Samples: %s' % ','.join(labels))
    print('Contaminants: %i, in a single sample: %i' %
          (len(contaminant_samples), hist.get(1, 0)))
    print('Heatmap written to %s' % args.plot_file)
    if percent_file is not None:
        print('Percentage heatmap written to %s' % percent_file)
    else:
        print('Percentage heatmap skipped: it needs the read totals from --runlog.')
    print('Bar plot written to %s' % bar_file)
    print('Histogram written to %s' % hist_file)
