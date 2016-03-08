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
:option:`sra` subcommand would be used:

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

.. todo::

.. code-block:: bash
    :linenos:

    FilterSeq missing
    FilterSeq repeats
    FilterSeq quality
    FilterSeq trimqual
    FilterSeq maskqual

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

Sampling and subsetting FASTA and FASTQ files
--------------------------------------------------------------------------------

.. todo::

.. code-block:: bash
    :linenos:

    SplitSeq group
    SplitSeq sample
    SplitSeq samplepair

Dealing with insufficient UMI diversity
--------------------------------------------------------------------------------

.. todo::

.. code-block:: bash
    :linenos:

    ClusterSets
    ParseHeaders
    BuildConsensus

Dealing with misaligned V-segment primers and indels in UMI groups
--------------------------------------------------------------------------------

.. todo::

.. code-block:: bash
    :linenos:

    AlignSets
    BuildConsensus

Assembling paired-end reads that do not overlap
--------------------------------------------------------------------------------

.. todo::

.. code-block:: bash
    :linenos:

    AssemblePairs

Assigning isotype annotations from the constant region sequence
--------------------------------------------------------------------------------

.. todo::

.. code-block:: bash
    :linenos:

    ConvertHeaders
    MaskPrimers
    ParseHeaders

Estimating sequencing and PCR error rates with UMI data
--------------------------------------------------------------------------------

.. todo::

.. code-block:: bash
    :linenos:

    EstimateError
