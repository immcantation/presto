#!/usr/bin/env bash
FilterSeq.py length -s SRR765688.fastq -n 300
FilterSeq.py quality -s SRR765688*length-pass.fastq -q 20 --outname STEP1
MaskPrimers.py score -s STEP1*quality-pass.fastq -p SRR765688_MID.fasta \
    --start 0 --mode cut --log MPMID.log
MaskPrimers.py align -s STEP1*primers-pass.fastq -p SRX190717_VPrimers.fasta \
    --maxlen 50 --mode cut --log MPVPrimers.log
MaskPrimers.py align -s STEP1*primers-pass_primers-pass.fastq
    -p SRX190717_CPrimers.fasta --maxlen 50 --revpr --skiprc --mode cut \
    --outname STEP2 --log MPCPrimers.log
ParseLog.py -l *.log -f ID PRSTART PRIMER ERROR
ParseHeaders.py expand -s STEP2_primers-pass.fastq -f PRIMER
ParseHeaders.py rename -s STEP2*reheader.fastq -f PRIMER1 PRIMER2 PRIMER3 \
    -k MID VPRIMER CPRIMER --outname STEP3
CollapseSeq.py -s STEP3_reheader.fastq -n 20 --inner --uf MID CPRIMER \
    --cf VPRIMER --act set
SplitSeq.py group -s STEP3_collapse-unique.fastq -f DUPCOUNT --num 2
ParseHeaders.py table -s STEP3*atleast-2.fastq \
    -f ID DUPCOUNT MID CPRIMER VPRIMER
