.. _ImportData:

Importing Data
================================================================================

Importing data from SRA, ENA or GenBank
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
sra           NCBI SRA or EMBL-EBI ENA
============  =================
