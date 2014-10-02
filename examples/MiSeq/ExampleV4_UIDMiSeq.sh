#!/bin/bash
# Example script to run the pRESTO pipeline for MiSeq data with UIDs 
# 
# Author:  Jason Anthony Vander Heiden, Gur Yaari
# Date:    2013.9.27

# Define run parameters
ZIP_FILES=true
PRCONS=0.6
MAXDIV=0.1
MINQUAL=20
R1_MAXERR=0.2 
R2_MAXERR=0.2
AP_MAXERR=0.2 
ALPHA=0.1
FS_MISS=10
CS_MISS=10
NPROC=16
        
# Define input and log files
R1_FILE="MiSeqDemo_R1.fastq"
R2_FILE="MiSeqDemo_R2.fastq"
R1_PRIMER_FILE="R1_MockPrimers.fasta"
R2_PRIMER_FILE="R2_MockPrimers.fasta"
R2_OFFSET_FILE="R2_MockPrimers_offsets-reverse.tab"
RUNLOG="Pipeline.log"

# Start
echo '' > $RUNLOG
echo "INPUT1:" $R1_FILE
echo "INPUT2:" $R2_FILE

# Filter low quality reads
echo "   1: FilterSeq quality      $(date +'%H:%M %D')"
FilterSeq.py quality -s $R1_FILE -q $MINQUAL --nproc $NPROC --outname R1 \
    --log QualityLogR1.log --clean >> $RUNLOG
FilterSeq.py quality -s $R2_FILE -q $MINQUAL --nproc $NPROC --outname R2 \
    --log QualityLogR2.log --clean >> $RUNLOG

# Identify primers and UID 
echo "   2: MaskPrimers score      $(date +'%H:%M %D')"
MaskPrimers.py score -s R1_quality-pass.fastq -p $R1_PRIMER_FILE --mode cut \
    --start 15 --barcode --maxerror $R1_MAXERR --nproc $NPROC \
    --log PrimerLogR1.log --clean >> $RUNLOG
MaskPrimers.py score -s R2_quality-pass.fastq -p $R2_PRIMER_FILE --mode cut \
    --start 0 --maxerror $R1_MAXERR --nproc $NPROC \
    --log PrimerLogR2.log --clean >> $RUNLOG
   
# Assign UIDs to read 1 sequences
echo "   3: PairSeq                $(date +'%H:%M %D')"
PairSeq.py -1 R1*primers-pass.fastq -2 R2*primers-pass.fastq -f BARCODE \
    --coord illumina --clean >> $RUNLOG

# Align sequence start positions by primer alignments
echo "   4: AlignSets              $(date +'%H:%M %D')" 
#AlignSets.py table -p $R2_PRIMER_FILE --reverse --exec /usr/local/bin/muscle3.8.31_i86linux64 > /dev/null
AlignSets.py offset -s R2*pair-pass.fastq -d $R2_OFFSET_FILE --bf BARCODE --pf PRIMER \
    --nproc $NPROC >> $RUNLOG

# Build UID consensus sequences
echo "   5: BuildConsensus         $(date +'%H:%M %D')" 
BuildConsensus.py -s R1*pair-pass.fastq --bf BARCODE --pf PRIMER --prcons $PRCONS \
    -q $MINQUAL --maxdiv $MAXDIV --nproc $NPROC --log ConsensusLogR1.log --clean >> $RUNLOG
BuildConsensus.py -s R2*pair-pass_align-pass.fastq --bf BARCODE --pf PRIMER \
    -q $MINQUAL --maxdiv $MAXDIV --nproc $NPROC --log ConsensusLogR2.log --clean >> $RUNLOG

# Assemble paired ends
echo "   6: AssemblePairs          $(date +'%H:%M %D')" 
AssemblePairs.py align -1 R2*consensus-pass.fastq -2 R1*consensus-pass.fastq \
    --1f CONSCOUNT --2f PRCONS CONSCOUNT --coord presto --rc tail --maxerror $AP_MAXERR --alpha $ALPHA \
    --nproc $NPROC --log AssembleLog.log --outname Assembled --clean >> $RUNLOG

# Remove sequences with many Ns
echo "   7: FilterSeq missing      $(date +'%H:%M %D')" 
FilterSeq.py missing -s Assembled_assemble-pass.fastq -n $FS_MISS --inner \
    --nproc $NPROC --log MissingLog.log >> $RUNLOG

# Rewrite header with minimum of CONSCOUNT
echo "   8: ParseHeaders collapse  $(date +'%H:%M %D')"
ParseHeaders.py collapse -s Assembled*missing-pass.fastq -f CONSCOUNT --act min --fasta > /dev/null

# Remove duplicate sequences
echo "   9: CollapseSeq            $(date +'%H:%M %D')" 
CollapseSeq.py -s Assembled*reheader.fasta -n $CS_MISS --uf PRCONS \
    --cf CONSCOUNT --act sum --outname Assembled --inner >> $RUNLOG

# Filter to sequences with at least 2 supporting sources
echo "  10: SplitSeq group         $(date +'%H:%M %D')" 
SplitSeq.py group -s Assembled_collapse-unique.fasta -f CONSCOUNT --num 2 >> $RUNLOG

# Process log files
echo "  11: ParseLog               $(date +'%H:%M %D')"
ParseLog.py -l QualityLogR[1-2].log -f ID QUALITY > /dev/null &
ParseLog.py -l PrimerLogR[1-2].log -f ID BARCODE PRIMER ERROR > /dev/null &
ParseLog.py -l ConsensusLogR[1-2].log -f BARCODE BCCOUNT CONSCOUNT PRIMER PRCOUNT PRFREQ DIVERSITY > /dev/null &
ParseLog.py -l AssembleLog.log -f ID OVERLAP LENGTH PVAL ERROR > /dev/null &
ParseLog.py -l MissingLog.log -f ID MISSING > /dev/null &
ParseLog.py -l Pipeline.log -f END SEQUENCES PAIRS SETS PASS FAIL UNIQUE DUPLICATE UNDETERMINED PARTS OUTPUT > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    tar -cf LogFiles.tar *.log 
    gzip LogFiles.tar
    rm *.log
    
    tar -cf TempFiles.tar R[1-2]_*.fastq *under* *duplicate* *undetermined* *reheader*
    gzip TempFiles.tar
    rm R[1-2]_*.fastq *under* *duplicate* *undetermined* *reheader*
fi
    
# End
echo -e "DONE\n" 
cd ../

