#!/usr/bin/env bash
FilterSeq.py length -s SRR765688.fastq -n 300 --outname S43 --log FSL.log
FilterSeq.py quality -s S43_length-pass.fastq -q 20 --outname S43 --log FSQ.log
MaskPrimers.py score -s S43_quality-pass.fastq -p SRR765688_MID.fasta \
    --start 0 --mode cut --outname S43 --log MPM.log
MaskPrimers.py align -s S43_primers-pass.fastq -p SRX190717_VPrimers.fasta \
    --maxlen 50 --mode cut --log MPV.log
MaskPrimers.py align -s S43_primers-pass_primers-pass.fastq
    -p SRX190717_CPrimers.fasta --maxlen 50 --revpr --skiprc --mode cut \
    --outname S43-MP --log MPC.log
ParseHeaders.py expand -s S43-MP_primers-pass.fastq -f PRIMER
ParseHeaders.py rename -s S43-MP*reheader.fastq -f PRIMER1 PRIMER2 PRIMER3 \
    -k MID VPRIMER CPRIMER --outname S43-MP
CollapseSeq.py -s S43-MP_reheader.fastq -n 20 --inner --uf MID CPRIMER \
    --cf VPRIMER --act set --outname S43
SplitSeq.py group -s S43_collapse-unique.fastq -f DUPCOUNT --num 2 --outname S43
ParseHeaders.py table -s S43_atleast-2.fastq -f ID DUPCOUNT MID CPRIMER VPRIMER
ParseLog.py -l FSL.log -f ID LENGTH
ParseLog.py -l FSQ.log -f ID QUALITY
ParseLog.py -l MP[MVC].log -f ID PRSTART PRIMER ERROR
