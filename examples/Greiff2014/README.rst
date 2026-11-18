This example uses publicly available data from:

    | **Quantitative assessment of the robustness of next-generation sequencing
      of antibody variable gene repertoires from immunized mice.**
    | Greiff, V. et al.
    | *BMC Immunol. 2014. 15(1):40. doi:10.1186/s12865-014-0040-5.*

Which may be downloaded from the EBI European Nucleotide Archive under
accession ID: ERP003950. For this example, we use the first 25,000
sequences of sample Replicate-1-M1 (accession: ERR346600), provided in the
example files as ``ERR346600_1.fastq`` (read 1, C-region end) and
``ERR346600_2.fastq`` (read 2, V-segment end). They were obtained from ENA
with::

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR346/ERR346600/ERR346600_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR346/ERR346600/ERR346600_2.fastq.gz
    zcat ERR346600_1.fastq.gz | head -n 100000 > ERR346600_1.fastq
    zcat ERR346600_2.fastq.gz | head -n 100000 > ERR346600_2.fastq

Note: NCBI SRA stores the mates of this run in the opposite order, so files
produced with ``fastq-dump --split-files -X 25000 ERR346600`` have read 1 and read 2
swapped relative to ENA and to this example.

Primers sequences are available in additional file 1 of the publication.
