#!/usr/bin/env bash
AssemblePairs.py align -1 ERR346600_2.fastq -2 ERR346600_1.fastq \
    --coord sra --rc tail --outname M1 --log AP.log
FilterSeq.py quality -s M1_assemble-pass.fastq -q 20 --log FS.log
MaskPrimers.py score -s M1*quality-pass.fastq -p Greiff2014_VPrimers.fasta \
    --start 4 --mode mask --log MP1.log
MaskPrimers.py score -s M1*primers-pass.fastq -p Greiff2015_CPrimers.fasta \
    --start 4 --mode cut --revpr --log MP2.log
ParseHeaders.py expand -s M1*primers-pass_primers-pass.fastq -f PRIMER
ParseHeaders.py rename -s M1*reheader.fastq -f PRIMER1 PRIMER2 \
    -k VPRIMER CPRIMER --outname M1-FINAL
CollapseSeq.py -s M1-FINAL_reheader.fastq -n 20 --inner --uf CPRIMER \
    --cf VPRIMER --act set
ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE
ParseLog.py -l FS.log -f ID QUALITY
ParseLog.py -l MP[1-2].log -f ID PRIMER ERROR
ParseHeaders.py table -s M1-FINAL*unique.fastq -f ID DUPCOUNT CPRIMER VPRIMER