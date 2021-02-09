.. image:: https://img.shields.io/pypi/dm/presto
    :target: https://pypi.org/project/presto
.. image:: https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic
    :target: https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html
	
pRESTO - The REpertoire Sequencing TOolkit
================================================================================

pRESTO is a toolkit for processing raw reads from high-throughput sequencing of
B cell and T cell repertoires.

Dramatic improvements in high-throughput sequencing technologies now enable
large-scale characterization of lymphocyte repertoires, defined as the
collection of trans-membrane antigen-receptor proteins located on the surface of
B cells and T cells. The REpertoire Sequencing TOolkit (pRESTO) is composed of a
suite of utilities to handle all stages of sequence processing prior to germline
segment assignment. pRESTO is designed to handle either single reads or
paired-end reads. It includes features for quality control, primer masking,
annotation of reads with sequence embedded barcodes, generation of
unique molecular identifier (UMI) consensus sequences, assembly of paired-end 
reads and identification of duplicate sequences. Numerous options for sequence 
sorting, sampling and conversion operations are also included.
