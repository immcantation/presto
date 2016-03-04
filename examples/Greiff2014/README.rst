This example uses publicly available data from:

    | **Quantitative assessment of the robustness of next-generation sequencing
      of antibody variable gene repertoires from immunized mice.**
    | Greiff, V. et al.
    | *BMC Immunol. 2014. 15(1):40. doi:10.1186/s12865-014-0040-5.*

Which may be downloaded from the EBI European Nucleotide Archive under
accession ID: ERP003950. For this example, we will use the first 25,000
sequences of sample Replicate-1-M1 (accession: ERR346600), which may
downloaded using fastq-dump from the
`SRA Toolkit <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`__::

    fastq-dump --split-files -X 25000 ERR346600

Primers sequences are available in additional file 1 of the publication.
