Fixing UMI Problems
================================================================================

.. _UMI-Alignment:

Correcting misaligned V-segment primers and indels in UMI groups
--------------------------------------------------------------------------------

Before generating a consensus for a set of reads sharing a UMI barcode,
the sequences must be properly aligned. Sequences may not be aligned if
more than one PCR primer is identified in a UMI read group - leading to
variations in the the start positions of the reads. Ideally, each set of
reads originating from a single mRNA molecule should be amplified with
the same primer. However, different primers in the multiplex pool may be
incorporated into the same UMI read group during amplification if the
primers are sufficiently similar.

.. _UMI-AlignmentFigure:

.. figure:: ../workflows/figures/Primer_Alignment.svg
    :align: center

    **Correction of misaligned sequences.**
    (A) Discrepancies in the location of primer binding (colored bases,
    with primer name indicated to the left) may cause misalignment of
    sequences sharing a UMI.
    (B) Following multiple alignment of the reads the non-primer regions are
    correctly aligned and suitable for UMI consensus generation.

This type of primer misalignment can be corrected using one of two approaches
using the :ref:`AlignSets` tool. The first approach, which is conceptually
simpler but computationally more expensive, is to perform a full multiple
alignment of reach UMI read group using the :program:`muscle` subcommand of
:ref:`AlignSets`. The :option:`--bf BARCODE <AlignSets muscle --bf>` argument tells
:ref:`AlignSets` to multiple align reads sharing the same ``BARCODE`` annotation.
The :option:`--exec ~/bin/muscle <AlignSets muscle --exec>` is a pointer to
where the `MUSCLE <http://www.drive5.com/muscle>`__ executable is located::

        AlignSets.py muscle -s reads.fastq --bf BARCODE --exec ~/bin/muscle

The above approach will also insert gaps into the sequences where an
insertion/deletion has occured in the reads. As such, you will need to provide
as reasonable gap character threshold to :ref:`BuildConsensus`, such as
:option:`--maxgap 0.5 <BuildConsensus --maxgap>`, defining how you want to handle
positions with gap characters when generating a UMI consensus sequence.

.. note::

    Using the :program:`muscle` subcommand, along with the
    :option:`--maxgap <BuildConsensus --maxgap>` argument to :ref:`BuildConsensus`
    will also address issue with insertions/deletions in UMI read groups.
    Though, in UMI read groups with a sufficient number of reads consensus generation
    will resolve insertions/deletions without the need for multiple alignment,
    as any misaligned reads will simply be washed out by the majority.
    Whether to perform a multiple alignment prior to consensus generation is a
    matter of taste. A multiple alignment may improve consensus quality in
    small UMI read groups (eg, less than 4 sequences), but the extent to which
    small UMI read groups should be trusted is debatable.

The second approach will correct only the primer regions and will not address
insertions/deletions within the sequence, but is much quicker to perform. The first
step involves creation of a primer offset table using the :program:`table` subcommand of
:ref:`AlignSets`::

    AlignSets.py table -p primers.fasta --exec ~/bin/muscle

Which performs a multiple alignment on sequences in ``primers.fasta``
(sequences shown in the :ref:`primer alignment figure <UMI-AlignmentFigure>` above)
to generate a file containing a primer offset table:

.. code-block:: none
    :caption: ``primers_offsets-forward.tab``

    VP1    2
    VP2	   0
    VP3	   1

Then the offset table can be input into the :program:`offset` subcommand of
:ref:`AlignSets` to align the reads::

    AlignSets.py offset -s reads.fastq -d primers_offsets-forward.tab \
        --bf BARCODE --pr VPRIMER --mode pad

In the above command we have specified the field containing the primer annotation
using :option:`--pf VPRIMER <AlignSets offset --pf>` and set the behavior
of the tool to add gap characters to align the reads with the
:option:`--mode pad <AlignSets offset --mode>` argument.
These options will generate the correction shown in **(B)** of the
:ref:`primer alignment figure <UMI-AlignmentFigure>` above.  Alternatively,
we could have deleted unalign positions using the argument
:option:`--mode cut <AlignSets offset --mode>`.

.. note::

    You may need to alter how the offset table is generated if you have used the
    :option:`--mode cut <MaskPrimers align --mode>` argument to :ref:`MaskPrimers`
    rather than :option:`--mode mask <MaskPrimers align --mode>`, as this will
    cause the ends of the primer regions, rather than the front, to be the
    cause of the ragged edges within the UMI read groups. For primers that
    have been cut you would add the :option:`--reverse <AlignSets table --reverse>`
    argument to the :program:`table` operation of :ref:`AlignSets`, which will
    create an offset table that is based on the tail end of the primers.

Dealing with insufficient UMI diversity
--------------------------------------------------------------------------------

Due to errors in the UMI region and/or insufficient UMI length, UMI read groups
are not always homogeneous with respect to the mRNA of origin. This can cause
difficulties in generating a valid UMI consensus sequence. In most cases,
the :option:`--prcons <BuildConsensus --prcons>` and
:option:`--maxerror <BuildConsensus --maxerror>`
(or :option:`--maxdiv <BuildConsensus --maxdiv>`) arguments to :ref:`BuildConsensus` are
sufficient to filter out invalid reads and/or entire invalid UMI groups. However, if
there is significant nucleotide diversity within UMI groups due to insufficient
UMI length or low UMI diversity, the :program:`set` command of the :ref:`ClusterSets`
tool can help correct for this. :ref:`ClusterSets` will cluster sequence by
similarity and add an additional annotation dividing sequences within a UMI read group
into sub-clusters::

    ClusterSets.py set -s reads.fastq -f BARCODE -k CLUSTER --exec ~/bin/usearch

The above command will add an annotation to each sequence named ``CLUSTER``
(:option:`-k CLUSTER <ClusterSets set -k>`) containing a cluster identifier
for each sequence within the UMI barcode group.
The :option:`-f BARCODE <ClusterSets set -f>` argument specifies the UMI annotation and
:option:`--exec ~/bin/usearch <ClusterSets set --exec>` is a pointer to
where the `USEARCH <http://www.drive5.com/usearch>`__ executable is located. After
assigning cluster annotations via :ref:`ClusterSets`, the ``BARCODE`` and ``CLUSTER``
fields can be merged using the :program:`copy` operation of :ref:`ParseHeaders`::

    ParseHeaders.py copy -s reads_cluster-pass.fastq -f BARCODE -k CLUSTER --act cat

Which will copy the UMI annotation (:option:`-f BARCODE <ParseHeaders copy -f>`) into
the cluster annotation (:option:`-k CLUSTER <ParseHeaders copy -k>`) and concatenate
them together (:option:`--act cat <ParseHeaders copy --act>`). Thus converting the
annotations from::

    >SEQ1|BARCODE=ATGTCG|CLUSTER=1
    >SEQ2|BARCODE=ATGTCG|CLUSTER=2

To::

    >SEQ1|BARCODE=ATGTCG|CLUSTER=1ATGTCG
    >SEQ2|BARCODE=ATGTCG|CLUSTER=2ATGTCG

You may then specify :option:`--bf CLUSTER <BuildConsensus --bf>` to
:ref:`BuildConsensus` to tell it to generate UMI consensus sequences by
UMI sub-cluster, rather than by UMI barcode annotation.

Combining split UMIs
--------------------------------------------------------------------------------

Typically, a UMI barcode is attached to only one end of a paired-end mate-pair
and can be copied to other read by a simple invocation of :ref:`PairSeq`.
But in some cases, the UMI may be split such that there are two UMIs, each located on a
different mate-pair. To deal with these sorts of UMIs, you would first employ
:ref:`PairSeq` similarly to how you would in the
:ref:`single UMI case <Stern2014-PairSeq-1>`:

    PairSeq.py -1 reads-1.fastq -2 reads-2.fastq --1f BARCODE --2f BARCODE \
        --coord illumina

The main difference from the single UMI case is that the ``BARCODE`` annotation is
being  simultaneously copied from read 1 to read 2 (:option:`--1f BARCODE <PairSeq --1f>`)
andf rom read 2 to read 1 (:option:`--2f BARCODE <PairSeq --2f>`). This creates
a set of annotations that look like::

    >READ1|BARCODE=ATGTCGTT,GGCTAGTC
    >READ2|BARCODE=ATGTCGTT,GGCTAGTC

Alternatively, these annotations can be combined upon copy using the
:option:`--act cat <PairSeq --act>` argument::

    PairSeq.py -1 reads-1.fastq -2 reads-2.fastq --1f BARCODE --2f BARCODE \
        --coord illumina --act cat

Which concatenates the two values in the ``BARCODE`` field,
yielding UMI annotations suitable for input to :ref:`BuildConsensus`::

    >READ1|BARCODE=ATGTCGTTGGCTAGTC
    >READ2|BARCODE=ATGTCGTTGGCTAGTC


Estimating sequencing and PCR error rates with UMI data
--------------------------------------------------------------------------------

The :ref:`EstimateError` tool provides methods for estimating the combined
PCR and sequencing error rates from large UMI read groups. The assumptions being,
that consensus sequences generated from sufficiently large UMI read groups should
be accurate representations of the true sequences, and that the rate of mismatches
from consensus should therefore be an accurate estimate of the error rate in
the data. However, this is not guaranteed to be true, hence this approach can only
be considered an estimate of a data set's error profile. The following command
generates an error profile from UMI read groups with 50 or more sequences
(:option:`-n 50 <EstimateError -n>`), using a majority rule consensus sequence
(:option:`--mode freq <EstimateError --freq>`), and excluding UMI read groups
with high nucleotide diversity (:option:`--maxdiv 0.1 <EstimateError --maxdiv>`)::

    EstimateError.py -s reads.fastq -n 50 --mode freq --maxdiv 0.1

This generates the following tab-delimited files containing error rates broken
down by various criteria:

============================== ==============================
File                           Error profile
============================== ==============================
reads_error-position.tab       Error rates by read position
reads_error-quality.tab        Error rates by quality score
reads_error-nucleotide.tab     Error rates by nucleotide identity
reads_error-set.tab            Error rates by UMI read group size
============================== ==============================
