#!/usr/bin/env bash
FilterSeq.py length -s SRR765688.fastq -n 300 --outname S43 --log FSL.log
FilterSeq.py quality -s S43_length-pass.fastq -q 20 --outname S43 --log FSQ.log
MaskPrimers.py score -s S43_quality-pass.fastq -p SRR765688_MIDs.fasta \
    --start 0 --maxerror 0.1 --mode cut --pf MID \
    --outname S43-MID --log MPM.log
MaskPrimers.py align -s S43-MID_primers-pass.fastq -p SRX190717_VPrimers.fasta \
    --maxlen 50 --maxerror 0.3 --mode mask --pf VPRIMER \
    --outname S43-FWD --log MPV.log
MaskPrimers.py align -s S43-FWD_primers-pass.fastq -p SRX190717_CPrimers.fasta \
    --maxlen 50 --maxerror 0.3 --mode cut --pf CPRIMER --revpr --skiprc \
    --outname S43-REV --log MPC.log
CollapseSeq.py -s S43-REV_primers-pass.fastq -n 20 --inner --uf MID CPRIMER \
    --cf VPRIMER --act set --outname S43
SplitSeq.py group -s S43_collapse-unique.fastq -f DUPCOUNT --num 2 --outname S43
ParseHeaders.py table -s S43_atleast-2.fastq -f ID DUPCOUNT MID CPRIMER VPRIMER
ParseLog.py -l FSL.log -f ID LENGTH
ParseLog.py -l FSQ.log -f ID QUALITY
ParseLog.py -l MPM.log MPV.log MPC.log -f ID PRSTART PRIMER ERROR
