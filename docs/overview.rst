Overview
================================================================================

Scope and Features
--------------------------------------------------------------------------------

pRESTO performs all stages of raw sequence processing prior to alignment against
reference germline sequences. The toolkit is intended to be easy to use, but some
familiarity with commandline applications is expected. Rather than providing a
fixed solution to a small number of common workflows, we have designed pRESTO to
be as flexible as possible. This design philosophy makes pRESTO suitable for many
existing protocols and adaptable to future technologies, but requires users to
construct a sequence of commands and options specific to their experimental
protocol.

pRESTO is composed of a set of standalone tools to perform specific tasks, often
with a series of subcommands providing different behaviors. A brief description
of each tool is shown in the table below.

.. csv-table::
   :file: tools/tool_summary.tab
   :delim: tab
   :header-rows: 1
   :widths: 15, 10, 75

Input and Output
--------------------------------------------------------------------------------

All tools take as input standard FASTA or FASTQ formatted files and output files
in the same formats. This allows pRESTO to work seamlessly with other sequence
processing tools that use either of these data formats; any steps within a
pRESTO workflow can be exchanged for an alternate tool, if desired.

Each tool appends a specific suffix to its output files describing the step and
output. For example, MaskPrimers will append ``_primers-pass`` to the output
file containing successfully aligned sequences and ``_primers-fail`` to the file
containing unaligned sequences. Details regarding the suffixes added by each
tool can be found in the :ref:`Usage` documentation for the tool.

Annotation Scheme
--------------------------------------------------------------------------------

The majority of pRESTO tools manipulate and add sequences-specific annotations
as part of their processing functions using the scheme shown below. Each
annotation is delimited using a reserved character (``|`` by default), with the
annotation field name and values separated by a second reserved character
(``=`` by default), and each value within a field is separated by a third
reserved character (``,`` by default). These annotations follow the sequence
identifier, which itself immediately follows the ``>`` (FASTA) or ``@`` (FASTQ)
symbol denoting the beginning of a new sequence entry. The sequence identifier
is given the reserved field name ``ID``. To mitigate potential analysis
errors, each tool in pRESTO annotates sequences by appending values to existing
annotation fields when they exist, and will not overwrite or delete annotations
unless explicitly performed using the ParseHeaders tool. All reserved characters
can be redefined using the command line options.

**FASTA Annotation**::

    >SEQUENCE_ID|PRIMER=IgHV-6,IgHC-M|BARCODE=DAY7|DUPCOUNT=8
    NNNNCCACGATTGGTGAAGCCCTCGCAGACCCTCTCACTCACCTGTGCCATCTCCGGGGACAGTGTTTCTACCAAAA

**FASTQ Annotation**::

    @SEQUENCE_ID|PRIMER=IgHV-6,IgHC-M|BARCODE=DAY7|DUPCOUNT=8
    NNNNCCACGATTGGTGAAGCCCTCGCAGACCCTCTCACTCACCTGTGCCATCTCCGGGGACAGTGTTTCTACCAAAA
    +
    !!!!nmoomllmlooj\Xlnngookkikloommononnoonnomnnlomononoojlmmkiklonooooooooomoo

Details regarding the annotations added by each tool can be found in the
:ref:`Usage` documentation for the tool.

