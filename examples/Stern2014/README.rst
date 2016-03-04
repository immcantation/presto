This example uses publicly available data from:

    | **B cells populating the multiple sclerosis brain mature in the draining
      cervical lymph nodes.**
    | Stern JNH, Yaari G, and Vander Heiden JA, et al.
    | *Sci Transl Med. 2014. 6(248):248ra107. doi:10.1126/scitranslmed.3008879.*

Which may be downloaded from the NCBI Sequence Read Archive under
BioProject accession ID: PRJNA248475. For this example, we will use the first
25,000 sequences of sample M12 (accession: SRR1383456), which may downloaded
downloaded using fastq-dump from the
`SRA Toolkit <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`__::

    fastq-dump --split-files -X 25000 SRR1383456

Primers sequences are available online at the
`supplemental website <http://clip.med.yale.edu/papers/Stern2014STM>`__
for the publication.