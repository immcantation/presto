"""
Sequence processing functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import re
import sys
from collections import OrderedDict
from itertools import product, zip_longest
from Bio import pairwise2
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field, default_primer_field, \
                            default_gap_chars, default_mask_chars, default_missing_chars, \
                            default_min_freq, default_min_qual, \
                            default_gap_penalty, default_max_error, default_max_len, default_start
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation


class PrimerAlignment:
    """
    A class defining a primer alignment result
    """
    # Instantiation
    def __init__(self, seq=None):
        """
        Initializer

        Arguments:
          seq : Bio.SeqRecord.SeqRecord object contained the sequence primers were aligned against.

        Returns:
          presto.Sequence.PrimerAlignment
        """
        self.seq = seq
        self.primer = None
        self.align_seq = None
        self.align_primer = None
        self.start = None
        self.end = None
        self.gaps = 0
        self.error = 1
        self.rev_primer = False
        self.valid = False

    # Set boolean evaluation to valid value
    def __bool__(self):
        """
        Boolean evaluation of the alignment

        Returns:
          int :  evaluates to the value of the valid attribute
        """
        return self.valid

    # Set length evaluation to length of alignment
    def __len__(self):
        """
        Length of alignment

        Returns:
          int : length of align_seq attribute
        """
        if self.align_seq is None:
            return 0
        else:
            return len(self.align_seq)


def translateAmbigDNA(key):
    """
    Translates IUPAC Ambiguous Nucleotide characters to or from character sets

    Arguments:
      key : String or re.search object containing the character set to translate

    Returns:
      str : Character translation
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
        return ['[' + k + ']' for k, v in IUPAC_trans.items() if v == key][0]
    # Return single ambiguous character for character set
    elif key in IUPAC_trans:
        return IUPAC_trans[key]
    else:
        return 'N'


def scoreDNA(a, b, mask_score=None, gap_score=None):
    """
    Returns the score for a pair of IUPAC Ambiguous Nucleotide characters

    Arguments:
      a : First characters
      b : Second character
      n_score : Tuple of length two defining scores for all matches against an N
                character for (a, b), with the score for character (a) taking precedence;
                if None score symmetrically according to IUPAC character identity
      gap_score : Tuple of length two defining score for all matches against a gap (-, .)
                  character for (a, b), with the score for character (a) taking precedence;
                  if None score symmetrically according to IUPAC character identity

    Returns:
      int : Score for the character pair
    """
    # Define ambiguous character translations
    IUPAC_trans = {'AGWSKMBDHV':'R', 'CTSWKMBDHV':'Y', 'CGKMBDHV':'S', 'ATKMBDHV':'W', 'GTBDHV':'K',
                   'ACBDHV':'M', 'CGTDHV':'B', 'AGTHV':'D', 'ACTV':'H', 'ACG':'V', 'ABCDGHKMRSTVWY':'N',
                   '-.':'.'}

    # Create list of tuples of synonymous character pairs
    IUPAC_matches = [p for k, v in IUPAC_trans.items() for p in list(product(k, v))]

    # Check gap and N-value conditions, prioritizing score for first character
    if gap_score is not None and a in '-.':
        return gap_score[0]
    elif mask_score is not None and a in 'nN':
        return mask_score[0]
    elif gap_score is not None and b in '-.':
        return gap_score[1]
    elif mask_score is not None and b in 'nN':
        return mask_score[1]

    # Return symmetric and reflexive score for IUPAC match conditions
    if a == b:
        return 1
    elif (a, b) in IUPAC_matches:
        return 1
    elif (b, a) in IUPAC_matches:
        return 1
    else:
        return 0


def scoreAA(a, b, mask_score=None, gap_score=None):
    """
    Returns the score for a pair of IUPAC Extended Protein characters

    Arguments:
      a : First character
      b : Second character
      mask_score : Tuple of length two defining scores for all matches against an X
                   character for (a, b), with the score for character (a) taking precedence;
                   if None score symmetrically according to IUPAC character identity
      gap_score : Tuple of length two defining score for all matches against a gap (-, .)
                  character for (a, b), with the score for character (a) taking precedence;
                  if None score symmetrically according to IUPAC character identity

    Returns:
      int : Score for the character pair
    """
    # Define ambiguous character translations
    IUPAC_trans = {'RN':'B', 'EQ':'Z', 'LI':'J', 'ABCDEFGHIJKLMNOPQRSTUVWYZ':'X',
                   '-.':'.'}
    # Create list of tuples of synonymous character pairs
    IUPAC_matches = [p for k, v in IUPAC_trans.items() for p in list(product(k, v))]

    # Check gap and X-value conditions, prioritizing score for first character
    if gap_score is not None and a in '-.':
        return gap_score[0]
    elif mask_score is not None and a in 'xX':
        return mask_score[0]
    elif gap_score is not None and b in '-.':
        return gap_score[1]
    elif mask_score is not None and b in 'xX':
        return mask_score[1]

    # Return symmetric and reflexive score for IUPAC match conditions
    if a == b:
        return 1
    elif (a, b) in IUPAC_matches:
        return 1
    elif (b, a) in IUPAC_matches:
        return 1
    else:
        return 0


def getDNAScoreDict(mask_score=None, gap_score=None):
    """
    Generates a score dictionary

    Arguments:
      mask_score : Tuple of length two defining scores for all matches against an N
                   character for (a, b), with the score for character (a) taking precedence;
                   if None score symmetrically according to IUPAC character identity
      gap_score : Tuple of length two defining score for all matches against a [-, .]
                  character for (a, b), with the score for character (a) taking precedence;
                  if None score symmetrically according to IUPAC character identity

    Returns:
      dict : Score dictionary with keys (char1, char2) mapping to scores
    """
    chars = '-.ACGTRYSWKMBDHVN'
    score_dict = {k:scoreDNA(*k, mask_score=mask_score, gap_score=gap_score)
                  for k in product(chars, repeat=2)}

    return score_dict


def getAAScoreDict(mask_score=None, gap_score=None):
    """
    Generates a score dictionary

    Arguments:
      mask_score : Tuple of length two defining scores for all matches against an X
                   character for (a, b), with the score for character (a) taking precedence;
                   if None score symmetrically according to IUPAC character identity
      gap_score : Tuple of length two defining score for all matches against a [-, .]
                  character for (a, b), with the score for character (a) taking precedence;
                  if None score symmetrically according to IUPAC character identity

    Returns:
      dict : Score dictionary with keys (char1, char2) mapping to scores
    """
    chars = '-.*ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    score_dict = {k:scoreAA(*k, mask_score=mask_score, gap_score=gap_score)
                  for k in product(chars, repeat=2)}

    return score_dict


def reverseComplement(seq):
    """
    Takes the reverse complement of a sequence

    Arguments:
      seq : a SeqRecord object, Seq object or string to reverse complement

    Returns:
      Seq : Object of the same type as the input with the reverse complement sequence
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


def checkSeqEqual(seq1, seq2, ignore_chars=default_missing_chars):
    """
    Determine if two sequences are equal, excluding missing positions

    Arguments:
      seq1 : SeqRecord object
      seq2 : SeqRecord object
      ignore_chars : Set of characters to ignore

    Returns:
      bool : True if the sequences are equal
    """
    equal = True
    #for a, b in zip(seq1.upper(), seq2.upper()):
    for a, b in zip_longest(seq1, seq2):
        if a != b and a not in ignore_chars and b not in ignore_chars:
            equal = False
            break

    return equal


# TODO:  can be removed I think.
def weightSeq(seq, ignore_chars=set()):
    """
    Returns the length of a sequencing excluding ignored characters

    Arguments:
      seq : SeqRecord or Seq object
      ignore_chars : Set of characters to ignore when counting sequence length

    Returns:
      int : Sum of the character scores for the sequence
    """
    return sum(1 for x in seq if x not in ignore_chars)


def scoreSeqPair(seq1, seq2, ignore_chars=set(), score_dict=getDNAScoreDict()):
    """
    Determine the error rate for a pair of sequences

    Arguments:
      seq1 : SeqRecord object
      seq2 : SeqRecord object
      ignore_chars : Set of characters to ignore when scoring and counting the weight
      score_dict : Optional dictionary of alignment scores

    Returns:
      Tuple : Tuple of the (score, minimum weight, error rate) for the pair of sequences
    """
    # TODO:  remove upper calls for speed. maybe by extending score dict with lowercase.
    # Determine score
    chars = list(zip(seq1.upper(), seq2.upper()))
    score_list = [score_dict[(a, b)] for a, b in chars \
                  if a not in ignore_chars and b not in ignore_chars]
    score = sum(score_list)
    weight = len(score_list)
    error = 1.0 - float(score) / weight if weight > 0 else 1.0

    return (score, weight, error)


def calculateDiversity(seq_list, score_dict=getDNAScoreDict()):
    """
    calculateDiversity(seq_list, score_dict=getDNAScoreDict())

    Determine the average pairwise error rate for a list of sequences

    Arguments:
      seq_list : List of SeqRecord objects to score
      score_dict : Optional dictionary of alignment scores as {(char1, char2): score}

    Returns:
      float : Average pairwise error rate for the list of sequences
    """
    # Return 0 if less than 2 sequences
    if len(seq_list) <= 1:
        return 0

    scores = []
    for i, seq1 in enumerate(seq_list):
        for seq2 in seq_list[(i + 1):]:
            scores.append(scoreSeqPair(seq1, seq2, score_dict=score_dict)[2])

    return sum(scores) / len(scores)


def calculateSetError(seq_list, ref_seq, ignore_chars=default_mask_chars,
                      score_dict=getDNAScoreDict()):
    """
    calculateSetError(seq_list, ref_seq, ignore_chars=['n', 'N'], score_dict=getDNAScoreDict())

    Counts the occurrence of nucleotide mismatches from a reference in a set of sequences

    Arguments:
      seq_list : list of SeqRecord objects with aligned sequences.
      ref_seq : SeqRecord object containing the reference sequence to match against.
      ignore_chars : list of characters to exclude from mismatch counts.
      score_dict : optional dictionary of alignment scores as {(char1, char2): score}.

    Returns:
      float : error rate for the set.
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

    # Calculate and return error rate
    try:
        return 1.0 - float(score) / total
    except ZeroDivisionError:
        return 1.0


def deleteSeqPositions(seq, positions):
    """
    Deletes a list of positions from a SeqRecord

    Arguments:
      seq : SeqRecord objects
      positions : Set of positions (indices) to delete

    Returns:
      SeqRecord : Modified SeqRecord with the specified positions removed
    """
    seq_del = ''.join([x for i, x in enumerate(seq.seq) if i not in positions])
    record = SeqRecord(Seq(seq_del, IUPAC.ambiguous_dna),
                       id=seq.id, name=seq.name, description=seq.description)

    if 'phred_quality' in seq.letter_annotations:
        qual_del = [x for i, x in enumerate(seq.letter_annotations['phred_quality']) \
                    if i not in positions]
        record.letter_annotations['phred_quality'] = qual_del

    return record


def findGapPositions(seq_list, max_gap, gap_chars=default_gap_chars):
    """
    Finds positions in a set of aligned sequences with a high number of gap characters.

    Arguments:
      seq_list : List of SeqRecord objects with aligned sequences
      max_gap : Float of the maximum gap frequency to consider a position as non-gapped
      gap_chars : Set of characters to consider as gaps

    Returns:
      list : Positions (indices) with gap frequency greater than max_gap
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
      seq_list : List of SeqRecord objects
      min_qual : Quality cutoff to assign a base
      min_freq : Frequency cutoff to assign a base
      dependent : If False assume sequences are independent for quality calculation
      ignore_chars : Set of characters to exclude when building a consensus sequence

    Returns:
      SeqRecord : Consensus SeqRecord object
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
            qual_cons = {c: int(max(qual_set[c]) * qual_sum[c] / qual_total)
                         for c in qual_set}
        else:
            qual_cons = {c: int(qual_sum[c] * qual_sum[c] / qual_total)
                         for c in qual_set}

        # Select character with highest consensus quality
        cons = [(c, min(q, 90)) for c, q in qual_cons.items() \
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
      set_seq : List of SeqRecord objects
      min_freq : Frequency cutoff to assign a base
      ignore_chars : Set of characters to exclude when building a consensus sequence


    Returns:
      SeqRecord : Consensus SeqRecord object
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
      seq_dict : a dictionary index of sequences returned from SeqIO.index()
      field : the annotation field containing set IDs
      delimiter : a tuple of delimiters for (fields, values, value lists)

    Returns:
      dict : Dictionary mapping set name to a list of record names
    """
    set_dict = {}
    for key, rec in seq_dict.items():
        tag = parseAnnotation(rec.description, delimiter=delimiter)[field]
        set_dict.setdefault(tag, []).append(key)

    return set_dict


def subsetSeqSet(seq_iter, field, values, delimiter=default_delimiter):
    """
    Subsets a sequence set by annotation value

    Arguments:
      seq_iter : Iterator or list of SeqRecord objects
      field : Annotation field to select by
      values : List of annotation values that define the retained sequences
      delimiter : Tuple of delimiters for (annotations, field/values, value lists)

    Returns:
      list : Modified list of SeqRecord objects
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
      seq_dict : Dictionary index of sequences returned from SeqIO.index()
      field : Annotation field to select keys by
      values : List of annotation values that define the retained keys
      delimiter : Tuple of delimiters for (annotations, field/values, value lists)

    Returns:
      list : List of keys
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
      key : Dictionary of sequences to translate

    Returns:
      dict : Dictionary of compiled regular expressions
    """

    primers_regex = {k: re.compile(re.sub(r'([RYSWKMBDHVN])', translateAmbigDNA, v))
                     for k, v in primers.items()}

    return primers_regex


def localAlignment(seq_record, primers, primers_regex=None, max_error=default_max_error,
                   max_len=default_max_len, rev_primer=False, skip_rc=False,
                   gap_penalty=default_gap_penalty,
                   score_dict=getDNAScoreDict(mask_score=(0, 1), gap_score=(0, 0))):
    """
    Performs pairwise local alignment of a list of short sequences against a long sequence

    Arguments:
      seq_record : a SeqRecord object to align primers against
      primers : dictionary of {names: short IUPAC ambiguous sequence strings}
      primers_regex : optional dictionary of {names: compiled primer regular expressions}
      max_error : maximum acceptable error rate before aligning reverse complement
      max_len : maximum length of sample sequence to align
      rev_primer : if True align with the tail end of the sequence
      skip_rc : if True do not check reverse complement sequences
      gap_penalty : a tuple of positive (gap open, gap extend) penalties
      score_dict : optional dictionary of alignment scores as {(char1, char2): score}

    Returns:
      presto.Sequence.PrimerAlignment : primer alignment result object
    """
    # Defined undefined parameters
    if primers_regex is None:  primers_regex = compilePrimers(primers)
    seq_record = seq_record.upper()
    rec_len = len(seq_record)
    max_len = min(rec_len, max_len)

    # Create empty return object
    align = PrimerAlignment(seq_record)
    align.rev_primer = rev_primer

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
        rec.annotations['prorient'] = 'F' if not rev_primer else 'RC'

        # Attempt regular expression match first
    for rec in seq_list:
        scan_seq = str(rec.seq)
        scan_seq = scan_seq[:max_len] if not rev_primer else scan_seq[-max_len:]
        for adpt_id, adpt_regex in primers_regex.items():
            adpt_match = adpt_regex.search(scan_seq)
            # Parse matches
            if adpt_match:
                align.seq = rec
                align.seq.annotations['primer'] = adpt_id
                align.primer = adpt_id
                align.align_seq = scan_seq
                align.align_primer = '-' * adpt_match.start(0) + \
                                     primers[adpt_id] + \
                                     '-' * (max_len - adpt_match.end(0))
                align.gaps = 0
                align.error = 0
                align.valid = True

                # Determine start and end positions
                if not rev_primer:
                    align.start = adpt_match.start(0)
                    align.end = adpt_match.end(0)
                else:
                    rev_pos = rec_len - max_len
                    align.start = adpt_match.start(0) + rev_pos
                    align.end = adpt_match.end(0) + rev_pos

                return align

    # Perform local alignment if regular expression match fails
    best_align, best_rec, best_adpt, best_error = None, None, None, None
    for rec in seq_list:
        this_align = dict()
        scan_seq = str(rec.seq)
        scan_seq = scan_seq[:max_len] if not rev_primer else scan_seq[-max_len:]
        for adpt_id, adpt_seq in primers.items():
            pw2_align = pairwise2.align.localds(scan_seq, adpt_seq, score_dict,
                                                -gap_penalty[0], -gap_penalty[1],
                                                one_alignment_only=True)
            if pw2_align:
                this_align.update({adpt_id: pw2_align[0]})
        if not this_align:  continue

        # Determine alignment with lowest error rate
        for x_adpt, x_align in this_align.items():
            x_error = 1.0 - x_align[2] / len(primers[x_adpt])
            # x_gaps = len(x_align[1]) - max_len
            # x_error = 1.0 - (x_align[2] + x_gaps) / primers[x_adpt])
            if best_error is None or x_error < best_error:
                best_align = this_align
                best_rec = rec
                best_adpt = x_adpt
                best_error = x_error

        # Skip rev_primer complement if forward sequence error within defined threshold
        if best_error <= max_error:  break

    # Set return object to lowest error rate alignment
    if best_align:
        # Define input alignment string and gap count
        align_primer = best_align[best_adpt][1]
        align_len = len(align_primer)
        align_gaps = align_len - max_len

        # Populate return object
        align.seq = best_rec
        align.primer = best_adpt
        align.align_seq = str(best_align[best_adpt][0])
        align.align_primer = align_primer
        align.gaps = align_gaps
        align.error = best_error
        align.valid = True

        # Determine start and end positions
        if not rev_primer:
            # TODO:  need to switch to an aligner that outputs start/end for both sequences in alignment
            align.start = align_len - len(align_primer.lstrip('-'))
            align.end = best_align[best_adpt][4] - align_gaps
        else:
            # Count position from tail and end gaps
            rev_pos = rec_len - align_len
            align.start = rev_pos + best_align[best_adpt][3] + align_gaps
            align.end = rev_pos + len(align_primer.rstrip('-'))

    return align


def scoreAlignment(seq_record, primers, start=default_start, rev_primer=False,
                   score_dict=getDNAScoreDict(mask_score=(0, 1), gap_score=(0, 0))):
    """
    Performs a simple fixed position alignment of primers

    Arguments:
      seq_record : a SeqRecord object to align primers against
      primers : dictionary of {names: short IUPAC ambiguous sequence strings}
      start : position where primer alignment starts
      rev_primer : if True align with the tail end of the sequence
      score_dict : optional dictionary of {(char1, char2): score} alignment scores

    Returns:
      presto.Sequence.PrimerAlignment : primer alignment result object
    """
    # Create empty return dictionary
    seq_record = seq_record.upper()
    align = PrimerAlignment(seq_record)
    align.rev_primer = rev_primer

    # Define orientation variables
    seq_record.annotations['seqorient'] = 'F'
    seq_record.annotations['prorient'] = 'F' if not rev_primer else 'RC'

    # Score primers
    this_align = {}
    rec_len = len(seq_record)
    if rev_primer:  end = rec_len - start
    for adpt_id, adpt_seq in primers.items():
        if rev_primer:
            start = end - len(adpt_seq)
        else:
            end = start + len(adpt_seq)
        chars = zip(seq_record[start:end], adpt_seq)
        score = sum([score_dict[(c1, c2)] for c1, c2 in chars])
        this_align.update({adpt_id: (score, start, end)})

    # Determine primer with lowest error rate
    best_align, best_adpt, best_err = None, None, None
    for adpt, algn in this_align.items():
        # adpt_err = 1.0 - float(algn[0]) / weightSeq(primers[adpt])
        err = 1.0 - float(algn[0]) / len(primers[adpt])
        if best_err is None or err < best_err:
            best_align = algn
            best_adpt = adpt
            best_err = err

    # Set return dictionary to lowest error rate alignment
    if best_align:
        # Populate return object
        align.primer = best_adpt if best_err < 1.0 else None
        align.start = best_align[1]
        align.end = best_align[2]
        align.error = best_err
        align.valid = True

        # Determine alignment sequences
        if not rev_primer:
            align.align_seq = str(seq_record.seq[:best_align[2]])
            align.align_primer = '-' * best_align[1] + primers[best_adpt]
        else:
            align.align_seq = str(seq_record.seq[best_align[1]:])
            align.align_primer = primers[best_adpt] + '-' * (rec_len - best_align[2])

    return align


def extractAlignment(seq_record, start, length, rev_primer=False):
    """
    Extracts a subsequence from sequence

    Arguments:
      data : SeqRecord to process.
      start : position where subsequence starts.
      length : the length of the subsequence to extract.
      rev_primer : if True extract relative to the tail end of the sequence.

    Returns:
      presto.Sequence.PrimerAlignment : extraction results as an alignment object
    """
    # Create empty return dictionary
    seq_record = seq_record.upper()
    align = PrimerAlignment(seq_record)
    align.rev_primer = rev_primer

    # Define orientation variables
    seq_record.annotations['seqorient'] = 'F'
    seq_record.annotations['prorient'] = 'F' if not rev_primer else 'RC'

    rec_len = len(seq_record)
    if rev_primer:
        end = rec_len - start
        start = end - length
    else:
        end = start + length

    # Extract region
    region = str(seq_record.seq[start:end])

    # Populate return object
    align.primer = region
    align.start = start
    align.end = end
    align.error = 0
    align.valid = True

    # Determine alignment sequences
    # Determine alignment sequences
    if not rev_primer:
        align.align_seq = str(seq_record.seq[:end])
        align.align_primer = '-' * start + region
    else:
        align.align_seq = str(seq_record.seq[start:])
        align.align_primer = region + '-' * (rec_len - end)

    return align


def maskSeq(align, mode='mask', barcode=False, barcode_field=default_barcode_field,
            primer_field=default_primer_field, delimiter=default_delimiter):
    """
    Create an output sequence with primers masked or cut

    Arguments:
      align : a PrimerAlignment object.
      mode : defines the action taken; one of 'cut', 'mask', 'tag' or 'trim'.
      barcode : if True add sequence preceding primer to description.
      barcode_field : name of the output barcode annotation.
      primer_field : name of the output primer annotation.
      delimiter : a tuple of delimiters for (annotations, field/values, value lists).

    Returns:
      Bio.SeqRecord.SeqRecord : masked sequence.
    """
    seq = align.seq

    # Build output sequence
    if mode == 'tag' or not align.align_primer:
        # Do not modify sequence
        out_seq = seq
    elif mode == 'trim':
        # Remove region before primer
        if not align.rev_primer:
            out_seq = seq[align.start:]
        else:
            out_seq = seq[:align.end]
    elif mode == 'cut':
        # Remove primer and preceding region
        if not align.rev_primer:
            out_seq = seq[align.end:]
        else:
            out_seq = seq[:align.start]
    elif mode == 'mask':
        # Mask primer with Ns and remove preceding region
        if not align.rev_primer:
            mask_len = align.end - align.start + align.gaps
            out_seq = 'N' * mask_len + seq[align.end:]
            if hasattr(seq, 'letter_annotations') and \
                    'phred_quality' in seq.letter_annotations:
                out_seq.letter_annotations['phred_quality'] = \
                    [0] * mask_len + \
                    seq.letter_annotations['phred_quality'][align.end:]
        else:
            mask_len = min(align.end, len(seq)) - align.start + align.gaps
            out_seq = seq[:align.start] + 'N' * mask_len
            if hasattr(seq, 'letter_annotations') and \
                    'phred_quality' in seq.letter_annotations:
                out_seq.letter_annotations['phred_quality'] = \
                    seq.letter_annotations['phred_quality'][:align.start] + \
                    [0] * mask_len

    # Add alignment annotations to output SeqRecord
    out_seq.annotations = seq.annotations
    out_seq.annotations['primer'] = align.primer
    out_seq.annotations['prstart'] = align.start
    out_seq.annotations['error'] = align.error

    # Parse seq annotation and create output annotation
    seq_ann = parseAnnotation(seq.description, delimiter=delimiter)
    out_ann = OrderedDict([('SEQORIENT', seq.annotations['seqorient']),
                           (primer_field, align.primer)])

    # Add ID sequence to description
    if barcode:
        seq_code = seq[:align.start].seq if not align.rev_primer \
            else seq[align.end:].seq
        out_seq.annotations['barcode'] = seq_code
        out_ann[barcode_field] = seq_code

    out_ann = mergeAnnotation(seq_ann, out_ann, delimiter=delimiter)
    out_seq.id = flattenAnnotation(out_ann, delimiter=delimiter)
    out_seq.description = ''

    return out_seq
