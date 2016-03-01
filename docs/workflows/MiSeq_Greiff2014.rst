Illumina MiSeq 2x250bp B cell receptor mRNA without UMIs
================================================================================

Overview of Experimental Data
--------------------------------------------------------------------------------

The following example uses the publicly available data from:

    | **Quantitative assessment of the robustness of next-generation sequencing
      of antibody variable gene repertoires from immunized mice.**
    | Greiff, V. et al.
    | *BMC Immunol. 15, 40 (2014).*


Which may be downloaded from the EBI European Nucleotide Archive under
accession ID: ERP003950. Primers sequences are available
from the publication. A schematic of the read configuration is shown below.

.. todo::

    Graphic of Grieff et al, 2014 read configuration

Example Data
--------------------------------------------------------------------------------

We have hosted a small subset of the data (Accession: ERR346600) on the
pRESTO website in FASTQ format with accompanying primer files. The sample data
set and workflow script may be downloaded from here:

.. todo::

    `Greiff et al, 2014 example <http://clip.med.yale.edu/presto/examples/Example_Data_Greiff2014.zip>`__

You may also retrieve the example data, the first 10,000 sequences of ERR346600,
using the `SRA Toolkit <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`__:

.. code-block:: bash

    sratoolkit/bin/fastq-dump --split-files -X 10000 ERR346600

Overview of the Workflow
--------------------------------------------------------------------------------

In the following sections, we demonstrate each step of the workflow to move
from raw sequence reads to a fully annotated repertoire of complete V(D)J
sequences. The workflow is divided into four high-level tasks:

    1. Assembly of paired-end reads.
    2. Quality control and annotation of primers.
    3. Removal of duplicate sequences and filtering.

A graphical representation of the workflow along with the corresponding
sequence of pRESTO commands is shown below.

.. figure:: figures/MiSeq_Greiff2014_Flowchart.svg
    :align: center

    **Flowchart of processing steps.**
    Each pRESTO tool is shown as a colored box. The workflow is divided into
    three primary tasks: (1) quality control and annotation of raw reads,
    (2) asssembly of paired-end reads, and (3) removal of duplicate sequences and filtering.
    The intermediate files output by each tool are not shown for the sake of brevity.

.. code-block:: shell
    :linenos:
    :caption: Commands

    AssemblePairs.py align -1 R2*pair-pass.fastq -2 R1*pair-pass.fastq \
        --coord sra --rc tail --outname M1 --log AP.log
    FilterSeq.py quality -s M1_assemble-pass.fastq -q 20
    MaskPrimers.py score -s R1_quality-pass.fastq -p Tipton2015_CPrimers.fasta \
        --start 0 --mode cut --log --outname MP1.log
    MaskPrimers.py score -s R2_quality-pass.fastq -p Tipton2015_VPrimers.fasta \
        --start 0 --mode mask --outname --log MP2.log
    ParseHeaders.py rename -s R1*primers-pass.fastq -f PRIMER -k CPRIMER
    ParseHeaders.py rename -s R2*primers-pass.fastq -f PRIMER -k VPRIMER
    CollapseSeq.py -s Final_assemble-pass.fastq -n 20 --inner --uf CPRIMER \
        --cf VPRIMER --act set

    ParseLog.py -l MP[1-2].log -f ID PRIMER ERROR
    ParseLog.py -l AP.log -f ID OVERLAP ERROR PVALUE
    ParseHeaders.py table -s Final*unique.fastq -f ID DUPCOUNT CPRIMER VPRIMER


Assembly of paired-end reads
--------------------------------------------------------------------------------

.. todo::

Quality control and annotation of raw reads
--------------------------------------------------------------------------------

.. todo::

Removal of duplicate sequences and filtering
--------------------------------------------------------------------------------

.. todo::

Stranded Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. todo::

    What to do when the MiSeq data is stranded

Unstranded Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. todo::

    What to do when the MiSeq data is unstranded



Performance
--------------------------------------------------------------------------------

.. todo::
