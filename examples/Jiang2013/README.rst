This example uses publicly available data from:

    | **Lineage structure of the human antibody repertoire in response to
      influenza vaccination.**
    | Jiang N, He J, and Weinstein JA, et al.
    | *Sci Transl Med. 2013. 5(171):171ra19. doi:10.1126/scitranslmed.3004794.*

Which may be downloaded from the NCBI Sequence Read Archive under
accession ID: SRX190717. For this example, we will use the first
50,000 sequences of sample 43 (accession: SRR765688), which may downloaded
downloaded using fastq-dump from the
`SRA Toolkit <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`__::

    fastq-dump -X 50000 SRR765688

Primer and sample barcode (referred to as MID in Jiang, He and Weinstein et al, 2013)
sequences are available in the published manuscript. This example assumes that the sample
barcodes, forward primers (V-region), and reverse primers (C-region) have been
extracted from the manuscript and placed into three corresponding FASTA files.