#!/usr/bin/env bash
FilterSeq.py quality -s SRR4026043_1.fastq -q 20 --outname HD09N-R1 --log FS1.log
FilterSeq.py quality -s SRR4026043_2.fastq -q 20 --outname HD09N-R2 --log FS2.log
MaskPrimers.py score -s HD09N-R1_quality-pass.fastq -p AbSeq_R1_Human_IG_Primers.fasta \
    --start 0 --mode cut --outname HD09N-R1 --log MP1.log
MaskPrimers.py score -s HD09N-R2_quality-pass.fastq -p AbSeq_R2_TS.fasta \
    --start 17 --barcode --mode cut --maxerror 0.5 --outname HD09N-R2 --log MP2.log
PairSeq.py -1 HD09N-R1_primers-pass.fastq -2 HD09N-R2_primers-pass.fastq \
    --2f BARCODE --coord sra
BuildConsensus.py -s HD09N-R1_primers-pass_pair-pass.fastq --bf BARCODE --pf PRIMER \
    --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname HD09N-R1 --log BC1.log
BuildConsensus.py -s HD09N-R2_primers-pass_pair-pass.fastq --bf BARCODE \
    --maxerror 0.1 --maxgap 0.5 --outname HD09N-R2 --log BC2.log
PairSeq.py -1 HD09N-R1_consensus-pass.fastq -2 HD09N-R2_consensus-pass.fastq \
    --coord presto
AssemblePairs.py sequential -1 HD09N-R2_consensus-pass_pair-pass.fastq \
    -2 HD09N-R1_consensus-pass_pair-pass.fastq -r IMGT_Human_IG_V.fasta \
    --coord presto --rc tail --scanrev --1f CONSCOUNT --2f CONSCOUNT PRCONS \
    --aligner blastn --outname HD09N-C --log AP.log
MaskPrimers.py align -s HD09N-C_assemble-pass.fastq \
    -p AbSeq_Human_IG_InternalCRegion.fasta --maxlen 100 --maxerror 0.3 \
    --mode tag --revpr --skiprc --log MP3.log --outname HD09N-C
ParseHeaders.py rename -s HD09N-C_primers-pass.fastq -f PRIMER -k CREGION
ParseHeaders.py collapse -s HD09N-C_primers-pass_reheader.fastq -f CONSCOUNT --act min
CollapseSeq.py -s HD09N-C_primers-pass_reheader_reheader.fastq -n 20 --inner \
    --uf CREGION --cf CONSCOUNT --act sum --outname HD09N-C
SplitSeq.py group -s HD09N-C_collapse-unique.fastq \
    -f CONSCOUNT --num 2 --outname HD09N-C
ParseHeaders.py table -s HD09N-C_atleast-2.fastq -f ID CREGION CONSCOUNT DUPCOUNT
ParseLog.py -l FS1.log FS2.log -f ID QUALITY
ParseLog.py -l MP1.log MP2.log MP3.log -f ID PRIMER BARCODE ERROR
ParseLog.py -l BC1.log BC2.log -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR
ParseLog.py -l AP.log -f ID REFID LENGTH OVERLAP GAP ERROR IDENTITY
