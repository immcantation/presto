.. _Tasks:

Miscellaneous Tasks
================================================================================

.. _Tasks-ImportingData:

Importing data from SRA, ENA or GenBank into pRESTO
--------------------------------------------------------------------------------

If you have download a data set from GenBank, SRA or ENA the format of the
sequences headers are different from the raw Roche 454 and Illumina header
format. As such, they may or may not be compatible with pRESTO, depending on
how the headers have been modified by the sequence archive. The
:ref:`ConvertHeaders` allow you to change incompatible header formats into
the pRESTO format. For example, to convert from SRA or ENA headers the
:program:`sra` subcommand would be used:

.. code-block:: none

    ConvertHeaders.py sra -s reads.fastq

:ref:`ConvertHeaders` provides the following conversion subcommands:

============  =================
Subcommand    Formats Converted
============  =================
generic       Headers with an unknown annotation system
454           Roche 454
genbank       NCBI GenBank and RefSeq
illumina      Illumina HiSeq or MiSeq
imgt          IMGT/GENE-DB
sra           NCBI SRA or EBI ENA
============  =================

Reducing file size for submission to IMGT/HighV-QUEST
--------------------------------------------------------------------------------

`IMGT/HighV-QUEST <http://imgt.org/HighV-QUEST>`__ currently limits the size of
uploaded files to 500,000 sequences. To accomodate this limit, you can use
the :program:`count` subcommand of :ref:`SplitSeq` to divide your files into
small pieces.

.. code-block:: none

    SplitSeq.py count -s reads.fastq -n 500000 --fasta

The :option:`-n 500000 <SplitSeq count -n>` argument sets the maximum number of
sequences in each file and the :option:`--fasta <SplitSeq count --fasta>`
tells the tool to output a FASTA, rather than FASTQ, formatted file.

.. note::

    You can usually avoid the necessity of reducing file sizes by removing
    duplicate sequences first using the :ref:`CollapseSeq` tool.

Subsetting sequence files by annotation
--------------------------------------------------------------------------------

The :program:`group` subcommand of :ref:`SplitSeq` allows you to split one file
into multiple files based on the values in a sequence annotation. For example,
splitting one file with multiple ``SAMPLE`` annotations into separate files
(one for each sample) would be accomplished by::

    SplitSeq.py group -s reads.fastq -f SAMPLE

Which will create a set of files labelled ``SAMPLE-M1`` and ``SAMPLE-M2``, if samples are
named ``M1`` and ``M2``.

If you wanted to split based on a numeric value, rather than a set of categorical values,
then you would add the :option:`--num <SplitSeq group --num>` argument. :ref:`SplitSeq`
would then create two files: one containing sequences with values less than the threshold
specified by the :option:`--num <SplitSeq group --num>` argument and one file containing
sequences with values greater than or equal to the threshold::

    SplitSeq.py group -s reads.fastq -f DUPCOUNT --num 2

Which will create two files with the labels ``atleast-2`` and ``under-2``.

Random sampling from sequence files
--------------------------------------------------------------------------------

The :program:`sample` subcommand of :ref:`SplitSeq` may be used to generate a
random sample from a sequence file or set of pair-end files. The example below
will select a random sample of 1,000 sequences (:option:`-n 1000 <SplitSeq sample -n>`)
which all contain the annotation ``SAMPLE=M1``
(:option:`-f SAMPLE <SplitSeq sample -f>` and :option:`-u M1 <SplitSeq sample -u>`)::

    SplitSeq.py sample -s reads.fastq -f SAMPLE -u M1 -n 1000

Performing an analogous sampling of Illumina paired-end reads would be accomplished using
the :program:`samplepair` subcommand::

    SplitSeq.py samplepair -s reads.fastq -f SAMPLE -u M1 -n 1000 --coord illumina

.. note::

    Both the :option:`-f <SplitSeq sample -f>` and :option:`-n <SplitSeq sample -n>`
    arguments will accept a list of values (eg, ``-n 1000 100 10``), allowing you to
    sample multiple times from multiple files in one command.

Cleaning or removing poor quality sequences
--------------------------------------------------------------------------------

Data sets can be cleaned using one or more invocations of :ref:`FilterSeq`,
which provides multiple sequence quality control operations.  Four subcommands
remove sequences from the data that fail to meet some threshold: including length,
(:program:`length`), number of N or gap characters (:program:`missing`),
homopolymeric tract length (:program:`repeats`), or mean Phred quality score
(:program:`quality`). Two subcommands modify sequences
without removing them: :program:`trimqual` truncates the sequences when the mean
Phred quality scores decays under a threshold, and :program:`maskqual` replaces
positions with low Phred quality scores with N characters.

:ref:`FilterSeq` provides the following quality control subcommands:

============ =================
Subcommand   Operation
============ =================
length       Removes short sequences
missing      Removes sequences with too many Ns or gaps
repeats      Removes sequences with long homopolymeric tracts
quality      Removes sequences with low mean quality scores
trimqual     Truncates sequences where quality scores decay
maskqual     Masks low quality positions
============ =================

Assembling paired-end reads that do not overlap
--------------------------------------------------------------------------------

The typical way to assemble paired-end reads is via *de novo* assembly using
the :program:`align` subcommand of :ref:`AssemblePairs`. However, some sequences
with long CDR3 regions may fail to assemble due to insufficient, or completely
absent, overlap between the mate-pairs. The :program:`reference` subcommand can
be used to assemble mate-pairs that do not overlap using the ungapped V-segment
references sequences as a guide.

First, a normal :program:`align` command would be performed. The
:option:`--failed <AssemblePairs align --failed>` argument is added so that
the reads failing *de novo* alignment are output to separate files::

    AssemblePairs.py align -1 reads-1.fastq -2 reads-2.fastq --rc tail \
        --coord illumina --failed -outname align

Then, the files labeled ``assemble-fail``, along with the ungapped V-segment
reference sequences (:option:`-r vref.fasta <AssemblePairs reference -r>`),
would be input into the :program:`reference` subcommand of :ref:`AssemblePairs`::

    AssemblePairs.py reference -1 align-1_assemble-fail.fastq -2 align-2_assemble-fail.fastq \
        -r vref.fasta --coord illumina --outname ref

Note, we have skipped the argument to reverse complement the tail sequence this
time (:option:`--rc tail <AssemblePairs reference --rc>`) as this was
done in the first invocation of :ref:`AssemblePairs`. You may then process the
two ``assemble-pass`` files separately or concatenate them together into a single file::

    cat align_assemble-pass.fastq ref_assemble-pass.fastq > merged_assemble-pass.fastq

.. note::

    The sequences output by the :program:`reference` subcommand will contain
    an appropriate length spacer of Ns between any mate-pairs that do not overlap.
    The `AssemblePairs reference --fill` argument can be specified to force
    :ref:`AssemblePairs` to insert the germline sequence into the missing positions,
    but this should be used with caution as the inserted sequence may not be
    biologically correct.

Assigning isotype annotations from the constant region sequence
--------------------------------------------------------------------------------

:ref:`MaskPrimers` is usually used to remove primer regions and annotate
sequences with primer identifiers. However, it can be used for any other case
where you need to align a set of short sequences against the reads. One example
alternate use is where you either do not know the C-region primer sequences
or do not trust the primer region to provide an accurate isotype assignment.

If you build a FASTA file containing the reverse-complement of short sequences
from the front of CH-1, then you can annotate the reads with these sequence in the same
way you would C-region specific primers::

    MaskPrimers.py align -s reads.fastq -p IGHC.fasta --maxlen 100 --maxerror 0.3 \
        --mode cut --revpr

Where :option:`--revpr <MaskPrimers align --revpr>` tells :ref:`MaskPrimers` to
reverse-complement the "primer" sequences and look for them at the end of the reads,
:option:`--maxlen 100 <MaskPrimers align --maxlen>` restricts the search to the last
100 bp, :option:`--maxerror 0.3 <MaskPrimers align --maxerror>` allows for up to
30% mismatches, and :option:`-p IGHC.fasta <MaskPrimers align -p>` specifies the file
containing the CH-1 sequences.  An example CH-1 sequence file would look like:

.. literalinclude:: ../workflows/data/IGHC.fasta
   :language: none

:download:`Download IGHC.fasta <../workflows/data/IGHC.fasta>`

.. seealso::

    Constant region reference sequences may be downloaded from
    `IMGT <http://imgt.org/vquest/refseqh.html>`__ and the sequence headers can be
    reformated using the :program:`imgt` subcommand of :ref:`ConvertHeaders`.
    Note, you may need to clean-up the reference sequences a bit
    before running :ref:`ConvertHeaders` if you receive an error about duplicate sequence names
    (eg, remove duplicate allele with different artificial splicing). To cut and
    reverse-complement the constant region sequences use something like
    `seqmagick <http://seqmagick.readthedocs.org>`__.
