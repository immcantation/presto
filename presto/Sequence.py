"""
Sequence processing functions
"""
# Future
from __future__ import division, absolute_import, print_function
from builtins import zip
from future.moves.itertools import zip_longest
from future.utils import iteritems

# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__, default_missing_chars

# Imports
import re
import sys
from itertools import product
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field, \
                            default_missing_chars, default_missing_residues, \
                            default_min_freq, default_min_qual
from presto.Annotation import parseAnnotation


def compilePrimers(primers):
    """
    Translates IUPAC Ambiguous Nucleotide characters to regular expressions and compiles them

    Arguments:
    key = a dictionary of sequences to translate

    Returns:
    dictionary of compiled regular expressions
    """

    primers_regex = {k: re.compile(re.sub(r'([RYSWKMBDHVN])', translateIUPAC, v))
                     for k, v in iteritems(primers)}

    return primers_regex


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
        return ['[' + k + ']' for k, v in iteritems(IUPAC_trans) if v == key][0]
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
    n_score = a tuple of length two defining scores for all matches against an N
              character for (a, b), with the score for character (a) taking precedence;
              if None score symmetrically according to IUPAC character identity
    gap_score = a tuple of length two defining score for all matches against a [-, .]
                character for (a, b), with the score for character (a) taking precedence;
                if None score symmetrically according to IUPAC character identity

    Returns:
    score for the character pair
    """
    # Define ambiguous character translations
    IUPAC_trans = {'AGWSKMBDHV':'R', 'CTSWKMBDHV':'Y', 'CGKMBDHV':'S', 'ATKMBDHV':'W', 'GTBDHV':'K',
                   'ACBDHV':'M', 'CGTDHV':'B', 'AGTHV':'D', 'ACTV':'H', 'ACG':'V', 'ABCDGHKMRSTVWY':'N',
                   '-.':'.'}
    # Create list of tuples of synonymous character pairs
    IUPAC_matches = [p for k, v in iteritems(IUPAC_trans) for p in list(product(k, v))]

    # Check gap condition
    if gap_score is not None and a in '-.':
        return gap_score[0]
    elif gap_score is not None and b in '-.':
        return gap_score[1]

    # Check N-value condition
    if n_score is not None and a == 'N':
        return n_score[0]
    elif n_score is not None and b == 'N':
        return n_score[1]

    # Return symmetric and reflexive score for IUPAC match conditions
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
    n_score = a tuple of length two defining scores for all matches against an X
              character for (a, b), with the score for character (a) taking precedence;
              if None score symmetrically according to IUPAC character identity
    gap_score = a tuple of length two defining score for all matches against a [-, .]
                character for (a, b), with the score for character (a) taking precedence;
                if None score symmetrically according to IUPAC character identity

    Returns:
    score for the character pair
    """
    # Define ambiguous character translations
    IUPAC_trans = {'RN':'B', 'EQ':'Z', 'LI':'J', 'ABCDEFGHIJKLMNOPQRSTUVWYZ':'X',
                   '-.':'.'}
    # Create list of tuples of synonymous character pairs
    IUPAC_matches = [p for k, v in iteritems(IUPAC_trans) for p in list(product(k, v))]

    # Check gap condition
    if gap_score is not None and a in '-.':
        return gap_score[0]
    elif gap_score is not None and b in '-.':
        return gap_score[1]

    # Check X-value condition
    if x_score is not None and a == 'X':
        return x_score[0]
    elif x_score is not None and b == 'X':
        return x_score[1]

    # Return symmetric and reflexive score for IUPAC match conditions
    if a == b:
        return 1
    elif (a, b) in IUPAC_matches:
        return 1
    elif (b, a) in IUPAC_matches:
        return 1
    else:
        return 0


def getDNAScoreDict(n_score=None, gap_score=None):
    """
    Generates a score dictionary

    Arguments:
    n_score = a tuple of length two defining scores for all matches against an N
              character for (a, b), with the score for character (a) taking precedence;
              if None score symmetrically according to IUPAC character identity
    gap_score = a tuple of length two defining score for all matches against a [-, .]
                character for (a, b), with the score for character (a) taking precedence;
                if None score symmetrically according to IUPAC character identity

    Returns:
    a score dictionary of the form {(char1, char2) : score}
    """
    chars = '-.ACGTRYSWKMBDHVN'
    score_dict = {k:scoreDNA(*k, n_score=n_score, gap_score=gap_score)
                  for k in product(chars, repeat=2)}

    return score_dict


def getAAScoreDict(x_score=None, gap_score=None):
    """
    Generates a score dictionary

    Arguments:
    x_score = a tuple of length two defining scores for all matches against an X
              character for (a, b), with the score for character (a) taking precedence;
              if None score symmetrically according to IUPAC character identity
    gap_score = a tuple of length two defining score for all matches against a [-, .]
                character for (a, b), with the score for character (a) taking precedence;
                if None score symmetrically according to IUPAC character identity

    Returns:
    a score dictionary of the form {(char1, char2) : score}
    """
    chars = '-.*ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    score_dict = {k:scoreAA(*k, x_score=x_score, gap_score=gap_score)
                  for k in product(chars, repeat=2)}

    return score_dict


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
    ignore_chars = a set of characters to ignore

    Returns:
    True if the sequences are equal
    """
    equal = True
    #for a, b in zip(seq1.upper(), seq2.upper()):
    for a, b in zip(seq1, seq2):
        if a != b and a not in ignore_chars and b not in ignore_chars:
            equal = False
            break

    return equal


def weightDNA(seq, ignore_chars=default_missing_chars):
    """
    Returns a score for a single sequence excluding missing positions

    Arguments:
    seq = a SeqRecord or Seq object
    ignore_chars = set of characters to ignore when counting sequence length

    Returns:
    The sum of the character scores for the sequence
    """
    score = sum(1 for x in seq if x not in ignore_chars)
    #score = sum()
    #nuc_score = sum([c in 'ACGTRYSWKMBDHV' for c in seq.upper()])
    #gap_score = 0

    return max(score, 1)


def weightAA(seq, ignore_residues=default_missing_residues):
    """
    Returns a score for a single sequence excluding missing positions

    Arguments:
    seq = a SeqRecord or Seq object
    ignore_residues = set of characters to ignore when counting sequence length

    Returns:
    The sum of the character scores for the sequence
    """
    score = sum(1 for x in seq if x not in ignore_residues)
    #score = sum()
    #nuc_score = sum([c in 'ACGTRYSWKMBDHV' for c in seq.upper()])
    #gap_score = 0

    return max(score, 1)


def scoreDNASeqPair(seq1, seq2, max_error=None, max_weight=None,
                 score_dict=getDNAScoreDict()):
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
    # TODO:  remove upper calls for speed. maybe by extending score dict with lowercase.
    # Determine score
    if max_error is None:
        # Return accurate values when max_error is undefined
        chars = zip(seq1.upper(), seq2.upper())
        score = sum([score_dict[(a, b)] for a, b in chars])
        weight = min(weightDNA(seq1), weightDNA(seq2))
        error = 1.0 - float(score) / weight
    else:
        # If max_error defined return when worst case reach
        score = 0
        if not max_weight:  max_weight = min(len(seq1), len(seq2))
        for i, (a, b) in enumerate(zip(seq1, seq2)):
            score += score_dict[(a, b)]
            if (i - float(score)) / max_weight > max_error:
                score, weight, error = 0, 0, 1.0
                break
        else:
            weight = min(weightDNA(seq1), weightDNA(seq2))
            error = 1.0 - float(score) / weight

    return (score, weight, error)


def calculateDiversity(seq_list, score_dict=getDNAScoreDict(n_score=(1, 1), gap_score=(0, 0))):
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
            scores.append(scoreDNASeqPair(seq1, seq2, score_dict=score_dict)[2])

    return sum(scores) / len(scores)


def calculateSetError(seq_list, ref_seq, ignore_chars=default_missing_chars,
                      score_dict=getDNAScoreDict(n_score=(1, 1), gap_score=(0, 0))):
    """
    Counts the occurrence of nucleotide mismatches from a reference in a set of sequences

    Arguments:
    seq_list = a list of SeqRecord objects with aligned sequences
    ref_seq = a SeqRecord object containing the reference sequence to match against
    ignore_chars = list of characters to exclude from mismatch counts
    score_dict = optional dictionary of alignment scores as {(char1, char2): score}

    Returns:
    a float of the error rate for the set
    """
    # Count informative characters in reference sequence
    ref_bases = sum(1 for b in ref_seq if b not in ignore_chars)

    # Return 0 mismatches for single record case
    if len(seq_list) <= 1:
        return 0.0

    # Iterate over seq_list and count mismatches
    total, score = 0, 0
    for seq in seq_list:
        seq_bases = sum(1 for a in seq if a not in ignore_chars)
        total += min(seq_bases, ref_bases)
        score += sum([score_dict[(a, b)] for a, b in zip(seq, ref_seq)
                      if a not in ignore_chars and b not in ignore_chars])

    return 1.0 - float(score) / total


def deleteSeqPositions(seq, positions):
    """
    Deletes a list of positions from a SeqRecord

    Arguments:
    seq = a SeqRecord objects
    positions = a set of positions (indices) to delete

    Returns:
    a modified SeqRecord with the specified positions removed
    """
    seq_del = ''.join([x for i, x in enumerate(seq.seq) if i not in positions])
    record = SeqRecord(Seq(seq_del, IUPAC.ambiguous_dna),
                       id=seq.id, name=seq.name, description=seq.description)

    if 'phred_quality' in seq.letter_annotations:
        qual_del = [x for i, x in enumerate(seq.letter_annotations['phred_quality']) \
                    if i not in positions]
        record.letter_annotations['phred_quality'] = qual_del

    return record


def findGapPositions(seq_list, max_gap, gap_chars=set(['.', '-'])):
    """
    Finds positions in a set of aligned sequences with a high number of gap characters.

    Arguments:
    seq_list = a list of SeqRecord objects with aligned sequences
    max_gap = a float of the maximum gap frequency to consider a position as non-gapped
    gap_chars = set of characters to consider as gaps

    Returns:
    a list of positions (indices) with gap frequency greater than max_gap
    """
    # Return an empty list in the singleton case
    seq_count = float(len(seq_list))
    if seq_count == 1:
        return []

    # Iterate through positions and count gaps
    gap_positions = []
    seq_str = [str(s.seq) for s in seq_list]
    for i, chars in enumerate(zip_longest(*seq_str, fillvalue='-')):
        gap_count = sum([chars.count(c) for c in gap_chars])
        gap_freq = gap_count / seq_count

        # Update gap position over threshold
        if gap_freq > max_gap:
            gap_positions.append(i)

    return gap_positions


def qualityConsensus(seq_list, min_qual=default_min_qual, min_freq=default_min_freq,
                     dependent=False, ignore_chars=default_missing_chars):
    """
    Builds a consensus sequence from a set of sequences

    Arguments:
    seq_list = a list of SeqRecord objects
    min_qual = the quality cutoff to assign a base
    min_freq = the frequency cutoff to assign a base
    dependent = if False assume sequences are independent for quality calculation
    ignore_chars = a set of characters to exclude when building a consensus sequence

    Returns:
    a consensus SeqRecord object
    """
    # Return a copy of the input SeqRecord upon singleton
    if len(seq_list) == 1:
        seq = seq_list[0]
        # Mask low quality nucleotides
        seq_str = str(seq.seq)
        quals = seq.letter_annotations['phred_quality']
        seq_mask = [seq_str[i] if q >= min_qual else 'N' for i, q in enumerate(quals)]
        seq = SeqRecord(Seq(''.join(seq_mask), IUPAC.ambiguous_dna),
                        id=seq.id,
                        name=seq.name,
                        description=seq.description,
                        letter_annotations=seq.letter_annotations)
        return seq

    # Create sequence and annotation iterators
    # Pad unequal length sequences with character '-' and quality 0
    seq_str = [str(s.seq) for s in seq_list]
    seq_iter = zip_longest(*seq_str, fillvalue='-')
    ann_list = [s.letter_annotations['phred_quality'] for s in seq_list]
    ann_iter = zip_longest(*ann_list, fillvalue=0)

    # Build consensus
    consensus_seq = []
    consensus_qual = []
    for chars, quals in zip(seq_iter, ann_iter):
        # Define set of non-missing characters
        char_set = set(chars).difference(ignore_chars)

        # Assign N if no missing characters and proceed to next position
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

        # TODO: write unit test and verify quality score calculation for sets with missing data
        # Calculate per character consensus quality scores
        if dependent:
            qual_cons = {c:int(max(qual_set[c]) * qual_sum[c] / qual_total) for c in qual_set}
        else:
            qual_cons = {c:int(qual_sum[c] * qual_sum[c] / qual_total) for c in qual_set}

        # Select character with highest consensus quality
        cons = [(c, min(q, 90)) for c, q in iteritems(qual_cons) \
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


def frequencyConsensus(seq_list, min_freq=default_min_freq,
                       ignore_chars=default_missing_chars):
    """
    Builds a consensus sequence from a set of sequences

    Arguments:
    set_seq = a list of SeqRecord objects
    min_freq = the frequency cutoff to assign a base
    ignore_chars = a set of characters to exclude when building a consensus sequence


    Returns:
    a consensus SeqRecord object
    """
    # Return a copy of the input SeqRecord upon singleton
    if len(seq_list) == 1:
        return seq_list[0]

    # Build consensus
    seq_str = [str(s.seq) for s in seq_list]
    consensus_seq = []
    for chars in zip_longest(*seq_str, fillvalue='-'):
        # Define set of non-missing characters
        char_set = set(chars).difference(ignore_chars)

        # Assign N if no non-missing characters and proceed to next position
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
    for key, rec in iteritems(seq_dict):
        tag = parseAnnotation(rec.description, delimiter=delimiter)[field]
        set_dict.setdefault(tag, []).append(key)

    return set_dict


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
