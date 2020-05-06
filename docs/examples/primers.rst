.. _Primers:

Isotype and Primer Annotations
================================================================================

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
        --mode cut --revpr --pf C_CALL

Where :option:`--revpr <MaskPrimers align --revpr>` tells :ref:`MaskPrimers` to
reverse-complement the "primer" sequences and look for them at the end of the reads,
:option:`--maxlen 100 <MaskPrimers align --maxlen>` restricts the search to the last
100 bp, :option:`--maxerror 0.3 <MaskPrimers align --maxerror>` allows for up to
30% mismatches, and :option:`-p IGHC.fasta <MaskPrimers align -p>` specifies the file
containing the CH-1 sequences. The name of the C-region will be added to the sequence
headers as the ``C_CALL`` annotation, where the field name is specified by the
:option:`--pf <MaskPrimers align --pf>` argument. An example CH-1 sequence file would look like:

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
    `seqmagick <http://seqmagick.readthedocs.io>`__.
