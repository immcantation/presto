"""
External application wrappers
"""
# Info
__author__    = 'Jason Anthony Vander Heiden, Namita Gupta'
from presto import __version__, __date__

# Imports
import csv
import sys
import tempfile
import pandas as pd
from itertools import chain
from io import StringIO
from subprocess import CalledProcessError, check_output, PIPE, Popen, STDOUT
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_muscle_exec, default_usearch_exec, default_blastn_exec

# Defaults
choices_usearch_method = ['ublast', 'usearch']
default_usearch_method = 'ublast'
default_uclust_ident = 0.9
default_usearch_ident = 0.5
default_evalue = 1e-5
default_max_hits = 100


def runMuscle(seq_list, muscle_exec=default_muscle_exec):
    """
    Multiple aligns a set of sequences using MUSCLE

    Arguments:
    seq_list = a list of SeqRecord objects to align
    muscle_exec = the MUSCLE executable

    Returns:
    a MultipleSeqAlignment object containing the alignment
    """
    # Return sequence if only one sequence in seq_list
    if len(seq_list) < 2:
        align = MultipleSeqAlignment(seq_list)
        return align

    # Set MUSCLE command
    cmd = [muscle_exec, '-diags', '-maxiters', '2']

    # Convert sequences to FASTA and write to string
    stdin_handle = StringIO()
    SeqIO.write(seq_list, stdin_handle, 'fasta')
    stdin_str = stdin_handle.getvalue()
    stdin_handle.close()

    # Open MUSCLE process
    child = Popen(cmd, bufsize=-1, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                  universal_newlines=True)

    # Send sequences to MUSCLE stdin and retrieve stdout, stderr
    stdout_str, __ = child.communicate(stdin_str)

    # Capture sequences from MUSCLE stdout
    stdout_handle = StringIO(stdout_str)
    align = AlignIO.read(stdout_handle, 'fasta')
    stdout_handle.close()

    return align


def runUClust(seq_list, ident=default_uclust_ident, seq_start=0, seq_end=None,
              usearch_exec=default_usearch_exec):
    """
    Cluster a set of sequences using the UCLUST algorithm from USEARCH

    Arguments:
    seq_list = a list of SeqRecord objects to align
    ident = the sequence identity cutoff to be passed to usearch
    seq_start = the start position to trim sequences at before clustering
    seq_end = the end position to trim sequences at before clustering
    usearch_exec = the path to the usearch executable

    Returns:
    a dictionary object containing {sequence id: cluster id}
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
    cmd = [usearch_exec,
           '-cluster_fast', in_handle.name,
           '-uc', out_handle.name,
           '-id', str(ident),
           '-qmask', 'none',
           '-minseqlength', '1',
           '-threads', '1']

    # Write usearch input fasta file
    SeqIO.write(seq_trimmed, in_handle, 'fasta')
    in_handle.seek(0)

    # Run usearch uclust algorithm
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                                  universal_newlines=True)

        # child = Popen(cmd, bufsize=1, stdout=PIPE, stderr=STDOUT, shell=False,
        #               universal_newlines=True)
        # while child.poll() is None:
        #     out = child.stdout.readline()
        #     sys.stdout.write(out)
        #     sys.stdout.flush()
        # child.wait()
    except CalledProcessError:
        return None

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


def runUSearch(seq, ref_file, method=default_usearch_method,
               ident=default_usearch_ident, evalue=default_evalue,
               max_hits=default_max_hits,
               usearch_exec=default_usearch_exec):
    """
    Aligns a sequence against a reference database using the usearch_local algorithm of USEARCH

    Arguments:
    seq = a SeqRecord object to align
    ref_file = the path to the reference database file
    method = the alignment method to use; one of 'usearch' or 'ublast'
    ident = the identity cut-off for usearch_local
    evalue = the E-value cut-off for usearch_local
    max_hits = the maxhits output limit for usearch_local
    usearch_exec = the path to the usearch executable

    Returns:
    a DataFrame of alignment results
    """
    # Set alignment algorithm
    method_dict = {'ublast': '-ublast', 'usearch': '-usearch_local'}
    try:
        method_cmd = method_dict[method]
    except:
        sys.exit('ERROR: invalid method argument for runUSearch')

    # Open temporary files
    in_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    out_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')

    # Define usearch command
    cmd = [usearch_exec,
           method_cmd, in_handle.name,
           '-db', ref_file,
           '-strand', 'plus',
           '-id', str(ident),
           '-evalue', str(evalue),
           '-maxhits', str(max_hits),
           '-maxaccepts', '0',
           '-maxrejects', '0',
           '-userout', out_handle.name,
           '-userfields', 'query+target+qlo+qhi+tlo+thi+alnlen+evalue+id',
           '-qmask', 'none',
           '-dbmask', 'none',
           '-threads', '1']

    # Write usearch input fasta file
    SeqIO.write(seq, in_handle, 'fasta')
    in_handle.seek(0)

    # Run usearch ublast algorithm
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False, universal_newlines=True)

        # child = Popen(cmd, bufsize=1, stdout=PIPE, stderr=STDOUT, shell=False,
        #              universal_newlines=True)
        # while child.poll() is None:
        #    out = child.stdout.readline()
        #    sys.stdout.write(out)
        #    sys.stdout.flush()
        # child.wait()
    except CalledProcessError:
        return None

    # Parse usearch output
    field_names = ['query', 'target', 'query_start', 'query_end',
                   'target_start', 'target_end',
                   'length', 'evalue', 'identity']
    align_df = pd.read_table(out_handle, header=None, names=field_names, encoding='utf-8')
    # Convert to base-zero indices
    align_df[['query_start', 'query_end', 'target_start', 'target_end']] -= 1

    # Close temp file handles
    in_handle.close()
    out_handle.close()

    return align_df


def runMakeUSearchDb(ref_file, usearch_exec=default_usearch_exec):
    """
    Makes a usearch database file

    Arguments:
    ref_file = the path to the reference database file
    usearch_exec = the path to the usearch executable

    Returns:
    A handle to the named temporary file containing the database file
    """
    # Open temporary file
    db_handle = tempfile.NamedTemporaryFile(suffix='.udb')

    # Define usearch command
    cmd = [usearch_exec,
           '-makeudb_ublast', ref_file,
           '-output', db_handle.name,
           '-dbmask', 'none']

    stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                              universal_newlines=True)

    return db_handle


def runBlastn(seq, ref_file, evalue=default_evalue, max_hits=default_max_hits,
              blastn_exec=default_blastn_exec):
    """
    Aligns a sequence against a reference database using BLASTN

    Arguments:
    seq = a SeqRecord objects to align
    ref_dict = a dictionary of reference SeqRecord objects
    evalue = the E-value cut-off for ublast
    maxhits = the maxhits output limit for ublast
    blastn_exec = the path to the usearch executable

    Returns:
    a DataFrame of alignment results
    """
    seq_fasta = seq.format('fasta')

    # Define blastn command
    cmd = ' '.join([blastn_exec,
                    '-query -',
                    #'-query', str(seq.seq),
                    '-subject', ref_file,
                    '-strand plus',
                    '-evalue', str(evalue),
                    '-max_target_seqs', str(max_hits),
                    '-outfmt "6 qseqid sseqid qstart qend sstart send length evalue pident"',
                    '-num_threads 1'])

    # Run blastn
    child = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                  shell=False, universal_newlines=True)
    stdout_str, stderr_str = child.communicate(seq_fasta)
    out_handle = StringIO(stdout_str)

    # Parse blastn output
    field_names = ['query', 'target', 'query_start', 'query_end', 'target_start', 'target_end',
                   'length', 'evalue', 'identity']
    align_df = pd.read_table(out_handle, header=None, names=field_names)
    # Convert to base-zero indices
    align_df[['query_start', 'query_end', 'target_start', 'target_end']] -= 1

    return align_df