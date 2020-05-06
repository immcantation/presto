.. _Assembly:

Fixing Assembly Problems
================================================================================

Assembling paired-end reads that do not overlap
--------------------------------------------------------------------------------

The typical way to assemble paired-end reads is via *de novo* assembly using
the :program:`align` subcommand of :ref:`AssemblePairs`. However, some sequences
with long CDR3 regions may fail to assemble due to insufficient, or completely
absent, overlap between the mate-pairs. The :program:`reference` or
:program:`sequential` subcommands can be used to assemble mate-pairs that do not
overlap using the ungapped V-segment references sequences as a guide.

To handle such sequence in two separate steps, a normal :program:`align` command
would be performed first. The :option:`--failed <AssemblePairs align --failed>`
argument is added so that the reads failing *de novo* alignment are output to
separate files::

    AssemblePairs.py align -1 reads-1.fastq -2 reads-2.fastq --rc tail \
        --coord illumina --failed -outname align

Then, the files labeled ``assemble-fail``, along with the ungapped V-segment
reference sequences (:option:`-r vref.fasta <AssemblePairs reference -r>`),
would be input into the :program:`reference` subcommand of :ref:`AssemblePairs`::

    AssemblePairs.py reference -1 align-1_assemble-fail.fastq -2 align-2_assemble-fail.fastq \
         --rc tail -r vref.fasta --coord illumina --outname ref

This will result in two separate ``assemble-pass`` files - one from each step. You may
process them separately or concatenate them together into a single file::

    cat align_assemble-pass.fastq ref_assemble-pass.fastq > merged_assemble-pass.fastq

However, if you intend to processes them together, you may simplify this by perform both
steps using the :program:`sequential` subcommand, which will attempt *de novo* assembly
followed by reference guided assembly if *de novo* assembly fails::

    AssemblePairs.py sequential -1 reads-1.fastq -2 reads-2.fastq --rc tail \
        --coord illumina -r vref.fasta

.. note::

    The sequences output by the :program:`reference` or :program:`sequential` subcommands
    may contain an appropriate length spacer of Ns between any mate-pairs that do not overlap.
    The :option:`--fill <AssemblePairs reference --fill>` argument may be specified to force
    :ref:`AssemblePairs` to insert the germline sequence into the missing positions,
    but this should be used with caution as the inserted sequence may not be
    biologically correct.
