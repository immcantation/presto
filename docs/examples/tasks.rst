.. _Tasks:

Miscellaneous Tasks
================================================================================

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

    ConvertHeaders.py sra -s file.fastq

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

Removing junk sequences
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

Reducing file size for submission to IMGT/HighV-QUEST
--------------------------------------------------------------------------------

`IMGT/HighV-QUEST <http://imgt.org/HighV-QUEST>`__ currently limits the size of
uploaded files to 500,000 sequences. To accomodate this limit, you can use
the :program:`count` subcommand of :ref:`SplitSeq` to divide your files into
small pieces.

.. code-block:: none

    SplitSeq.py count -s file.fastq -n 500,000 --fasta

The :option:`-n 500,000 <SplitSeq count -n>` argument sets the maximum number of
sequences in each file and the :option:`--fasta <SplitSeq count --fasta>`
tells the tool to output a FASTA, rather than FASTQ, formatted file.

.. note::

    You can usually avoid the necessity of reducing file sizes by removing
    duplicate sequences first using the :ref:`CollapseSeq` tool.

Sampling and subsetting sequence files
--------------------------------------------------------------------------------

In addition to
`splitting files into smaller pieces <Reducing file size for submission to IMGT/HighV-QUEST>`_,
the :ref:`SplitSeq` tool provides several other methods for subsetting and sampling from sequence
files.

Subsetting by annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :program:`group` subcommand allows you to split one file into multiple files based on
the values in a sequence annotation. For example, splitting one file with multiple ``SAMPLE``
annotations into separate files (one for each sample) would be accomplished by::

    SplitSeq.py group -s file.fastq -f SAMPLE

Which will create a set of files labelled ``SAMPLE-M1`` and ``SAMPLE-M2``, if samples are
named ``M1`` and ``M2``.

If you wanted to split based on a numeric value, rather than a set of categorical values,
then you would add the :option:`--num <SplitSeq group --num>` argument. :ref:`SplitSeq`
would then create two files: one containing sequences with values less than the threshold
specified by the :option:`--num <SplitSeq group --num>` argument and one file containing
sequences with values greater than or equal to the threshold::

    SplitSeq.py group -s file.fastq -f DUPCOUNT --num 2

Which will create two files with the labels ``atleast-2`` and ``under-2``.

Random sampling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To generate a random sample from a sequence file, the :program:`sample` subcommand can be
used. The example below will select a random sample of 1,000 sequences
(:option:`-n 1000 <SplitSeq sample -n>`) which all contain the annotation
``SAMPLE=M1`` (:option:`-f SAMPLE <SplitSeq sample -f>` and :option:`-u M1 <SplitSeq sample -u>`)::

    SplitSeq.py sample -s file.fastq -f SAMPLE -u M1 -n 1000

Performing an analagous sampling of Illumina pair-end reads would be accomplished using
the :program:`samplepair` subcommand::

    SplitSeq.py samplepair -s file.fastq -f SAMPLE -u M1 -n 1000 --coord illumina

.. note::

    Both the :option:`-f <SplitSeq sample -f>` and :option:`-n <SplitSeq sample -n>`
    arguments will accept a list of values (eg, ``-n 1000 100 10``), allowing you to
    sample multiple times from multiple files in one command.

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

    AssemblePairs.py align -1 read1.fastq -2 read1.fastq --rc tail \
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

.. todo::

.. code-block:: bash
    :linenos:

    ConvertHeaders
    MaskPrimers
    ParseHeaders


