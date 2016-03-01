#!/usr/bin/env bash
FilterSeq.py quality -s MS12_R1.fastq -q 20
FilterSeq.py quality -s MS12_R2.fastq -q 20
MaskPrimers.py score -s MS12_R1_quality-pass.fastq -p Stern2014_CPrimers.fasta \
    --start 15 --mode cut --barcode --log MP1.log
MaskPrimers.py score -s MS12_R2_quality-pass.fastq -p Stern2014_VPrimers.fasta \
    --start 0 --mode mask --log MP2.log
ParseLog.py -l MP1.log MP2.log -f ID PRIMER BARCODE ERROR
PairSeq.py -1 MS12_R1*primers-pass.fastq -2 MS12_R2*primers-pass.fastq \
    --coord illumina --1f BARCODE
BuildConsensus.py -s MS12_R1*pair-pass.fastq --bf BARCODE --pf PRIMER \
    --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --log BC1.log
BuildConsensus.py -s MS12_R2*pair-pass.fastq --bf BARCODE --pf PRIMER \
    --maxerror 0.1 --maxgap 0.5 --log BC2.log
ParseLog.py -l BC1.log BC2.log -f BARCODE CONSCOUNT PRCONS ERROR
PairSeq.py -1 MS12_R1*consensus-pass.fastq -2 MS12_R2*consensus-pass.fastq \
    --coord presto
AssemblePairs.py align -1 MS12_R2*consensus-pass_pair-pass.fastq \
    -2 MS12_R1*consensus-pass_pair-pass.fastq --coord presto --rc tail \
    --1f CONSCOUNT --2f CONSCOUNT PRCONS --outname MS12 --log AP.log
ParseLog.py -l AP.log -f ID OVERLAP ERROR PVALUE
ParseHeaders.py collapse -s MS12_assemble-pass.fastq -f CONSCOUNT --act min
CollapseSeq.py -s MS12*reheader.fastq -n 20 --inner --uf PRCONS \
    --cf CONSCOUNT --act sum
SplitSeq.py group -s MS12*collapse-unique.fastq -f CONSCOUNT --num 2
ParseHeaders.py table -s MS12*atleast-2.fastq -f ID PRCONS CONSCOUNT DUPCOUNT
