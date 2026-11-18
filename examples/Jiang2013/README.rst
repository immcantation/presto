This example uses publicly available data from:

    | **Lineage structure of the human antibody repertoire in response to
      influenza vaccination.**
    | Jiang N, He J, and Weinstein JA, et al.
    | *Sci Transl Med. 2013. 5(171):171ra19. doi:10.1126/scitranslmed.3004794.*

Which may be downloaded from the NCBI Sequence Read Archive under
accession ID: SRX190717. For this example, we use the first 50,000
sequences of sample 43 (accession: SRR765688), provided in the example files
as ``SRR765688.fastq.gz``. It was obtained using fastq-dump from the
`SRA Toolkit <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`__,
with SRA Lite disabled::

    vdb-config --simplified-quality-scores no
    fastq-dump -X 50000 SRR765688

The example script expects the uncompressed file, so decompress it first::

    gunzip SRR765688.fastq.gz

Note: Recent versions of the SRA Toolkit may return SRA Lite data by default,
in which quality scores are replaced by simplified placeholder values. The
workflow begins with quality filtering, so real quality scores are required.
If you download the data yourself, run the ``vdb-config`` command above first,
and do not use the FASTQ download from the NCBI SRA website, which also returns
SRA Lite data without real quality scores.

Primer and sample barcode (referred to as MID in Jiang, He and Weinstein et al, 2013)
sequences are available in the published manuscript. This example assumes that the sample
barcodes, forward primers (V-region), and reverse primers (C-region) have been
extracted from the manuscript and placed into three corresponding FASTA files.
