This example uses publicly available data from:

    | **B cells populating the multiple sclerosis brain mature in the draining
      cervical lymph nodes.**
    | Stern JNH, Yaari G, and Vander Heiden JA, et al.
    | *Sci Transl Med. 2014. 6(248):248ra107. doi:10.1126/scitranslmed.3008879.*

Which may be downloaded from the NCBI Sequence Read Archive under
BioProject accession ID: PRJNA248475. For this example, we use the first
25,000 sequences of sample M12 (accession: SRR1383456), provided in the
example files as ``SRR1383456_1.fastq.gz`` and ``SRR1383456_2.fastq.gz``.
They were obtained using fastq-dump from the
`SRA Toolkit <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`__,
with SRA Lite disabled::

    vdb-config --simplified-quality-scores no
    fastq-dump --split-files -X 25000 SRR1383456

The example script expects the uncompressed files, so decompress them first::

    gunzip SRR1383456_1.fastq.gz SRR1383456_2.fastq.gz

Note: Recent versions of the SRA Toolkit may return SRA Lite data by default,
in which quality scores are replaced by simplified placeholder values. The
workflow begins with quality filtering, so real quality scores are required.
If you download the data yourself, run the ``vdb-config`` command above first,
and do not use the FASTQ download from the NCBI SRA website, which also returns
SRA Lite data without real quality scores.

Primers sequences are available online at the
`supplemental website <http://clip.med.yale.edu/papers/Stern2014STM>`__
for the publication.