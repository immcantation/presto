.. _Annotations:

Manipulating Annotations
================================================================================

The :ref:`ParseHeaders` tool provides a collection of methods for performing
simple manipulations of sequence headers that are format in the
:ref:`pRESTO annotation scheme <AnnotationScheme>`.

For converting sequence headers into the pRESTO format see the
:ref:`Import Task <Tasks-ImportingData>` documentation.

Adding a sample annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Addition of annotation values is accomplished using the :program:`add` subcommand
of :ref:`ParseHeaders`::

    ParseHeaders.py add -s reads.fastq -f SAMPLE -u A1

Which will add the annotation ``SAMPLE=A1`` to each sequence of the input file.

Expanding and renaming annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, pRESTO will not delete annotations. If a sequence header already
contains an annotation that a tool is trying to add it will not overwrite that
annotation. Instead, it will append the annotation value to the values already
present in a comma delimited form. For example, after two interations of
:ref:`MaskPrimers` with the default primer field name (``PRIMER``), you will have
an annotation in the following form, reflecting a match against primer ``VH3`` in
the first iteration and primer ``IGHM`` in the second::

    PRIMER=VH3,IGHM

Separating these annotations into two annotations is accomplished via the
:program:`expand` subcommand of :ref:`ParseHeaders`::

    ParseHeaders.py expand -s reads.fastq -f PRIMER

Resulting in the annotations::

    PRIMER1=VH3|PRIMER2=IGHM

Which may then be renamed via the :program:`rename` subcommand:
:program:`expand` subcommand of :ref:`ParseHeaders`::

    ParseHeaders.py rename -s reads_reheader.fastq -f PRIMER1 PRIMER2 \
        -k VPRIMER CPRIMER

Copying, merging and collapsing annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Nested annotation can be generated using the :program:`copy` or :program:`merge`
subcommands of :ref:`ParseHeaders`. The examples that follow will use the starting
annotation::

    UMI=ATGC|CELL=GGCC|COUNT=10,2

The ``UMI`` and ``CELL`` annotations can be combined into a single ``INDEX``
annotation using the following command::

    ParseHeaders.py merge -s reads.fasta -f UMI CELL -k INDEX --delete
    # result> COUNT=10,2|INDEX=ATGC,GGCC

Without the :option:`--delete <ParseHeaders merge --delete>` argument the
original ``UMI`` and ``CELL`` annotations would be kept in the header.

The nested annotation values can then be combined
using the :program:`collapse` subcommand to create various effects::

    ParseHeaders.py collapse -s reads_reheader.fasta -f INDEX --act cat
    # result> INDEX=ATGCGGCC

    ParseHeaders.py collapse -s reads_reheader.fasta -f INDEX --act first
    # result> INDEX=ATGC

    ParseHeaders.py collapse -s reads_reheader.fasta -f COUNT --act sum
    # result> COUNT=12

    ParseHeaders.py collapse -s reads_reheader.fasta -f COUNT --act min
    # result> COUNT=2

Where the :option:`--act <ParseHeaders collapse --act>` argument specifies
the type of collapse action to perform.

The :program:`copy` subcommand is normally used to create duplicate annotations
with different names, but will have a similar effect to the :program:`merge`
subcommand when the target is an existing field::

    ParseHeaders.py copy -s reads.fasta -f UMI -k CELL
    # result> UMI=ATGC|CELL=GGCC,ATGC|COUNT=10,2

Both the :program:`copy` and :program:`merge` subcommands have an
:option:`--act <ParseHeaders collapse --act>` argument which allows
you to perform an action from the :program:`collapse` subcommand
in the same step as the :program:`copy` or :program:`merge`::

    ParseHeaders.py merge -s reads.fasta -f UMI CELL -k INDEX --delete --act cat
    # result> COUNT=10,2|INDEX=ATGCGGCC

    ParseHeaders.py copy -s reads.fasta -f UMI -k CELL --act cat
    # result> UMI=ATGC|CELL=GGCCATGC|COUNT=10,2

Deleting annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unwanted annotations can be deleted using the :program:`delete` subcommand
of :ref:`ParseHeaders`::

    ParseHeaders.py delete -s reads.fastq -f PRIMER

Which will remove the ``PRIMER`` field from each sequence header.