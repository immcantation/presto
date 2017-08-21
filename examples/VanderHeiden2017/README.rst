This example uses publicly available data from:

    | **Dysregulation of B Cell Repertoire Formation in Myasthenia Gravis Patients 
      Revealed through Deep Sequencing.**
    | Vander Heiden JA, et al.
    | *J Immunol. 2017 198(4):1460-1473. doi:10.4049/jimmunol.1601415.*

Which may be downloaded from the NCBI Sequence Read Archive under
BioProject accession ID: PRJNA248475. For this example, we will use the first
25,000 sequences of sample HD09_N_AB8KB (accession: SRR4026043), which may downloaded
downloaded using fastq-dump from the
`SRA Toolkit <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`__::

    fastq-dump --split-files -X 25000 SRR4026043

Primers sequences are available online from the 
`protocols/AbSeq <https://bitbucket.org/kleinstein/immcantation/src/tip/protocols/AbSeq>`__ 
directory in the `Immcantation repository <https://bitbucket.org/kleinstein/immcantation>`__.
