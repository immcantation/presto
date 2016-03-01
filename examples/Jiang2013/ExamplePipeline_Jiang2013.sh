#!/bin/bash
# Example pRESTO pipeline for Jiang2013 data with MIDs.
# Data from Jiang et al, 2013, Sci Trans Med.
#
# Author:  Jason Anthony Vander Heiden
# Date:    2015.10.10

NAME="SRR765688"
FILE="SRR765688.fastq"
MIDS_PRIMER_FILE="SRR765688_MIDs.fasta"
V_PRIMER_FILE="SRX190717_VPrimers.fasta"
C_PRIMER_FILE="SRX190717_CPrimers.fasta"
RUNLOG="Pipeline.log"

# Define run parameters and input files
READ_FILE=$(readlink -f SRR765688.fastq)
MID_PRIMERS=$(readlink -f SRR765688_MIDs.fasta)
FWD_PRIMERS=$(readlink -f SRX190717_VPrimers.fasta)
REV_PRIMERS=$(readlink -f SRX190717_CPrimers.fasta)
OUTDIR="output"
PIPELINE_LOG="Pipeline.log"
ZIP_FILES=true

# Make output directory and empty log files
mkdir -p $OUTDIR; cd $OUTDIR
echo '' > $PIPELINE_LOG

# Start
echo "DIRECTORY: ${OUTDIR}"
echo -e "START"
STEP=0

# Remove short reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq length"
FilterSeq.py length -s $READ_FILE -n 300 --outname FLTR --outdir . \
	--log LengthLog.log >> $PIPELINE_LOG

# Remove low quality reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
FilterSeq.py quality -s FLTR_length-pass.fastq -q 20 --outname FLTR \
	--log QualityLog.log >> $PIPELINE_LOG

# Identify and remove MIDs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
MaskPrimers.py score -s FLTR_quality-pass.fastq -p $MID_PRIMERS --mode cut \
    --start 0 --maxerror 0.1 --outname MID \
	--log PrimerMIDLog.log >> $PIPELINE_LOG

# Identify and mask forward (V-region) primers
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s MID_primers-pass.fastq -p $FWD_PRIMERS \
    --mode mask --maxlen 50 --maxerror 0.3 --outname PRMR \
    --log PrimerForwardLog.log >> $PIPELINE_LOG

# Identify and remove reverse (C-region) primers
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s PRMR_primers-pass.fastq -p $REV_PRIMERS \
    --mode cut --maxlen 50 --maxerror 0.3 --revpr --skiprc \
	--log PrimerReverseLog.log >> $PIPELINE_LOG

# Expand primer field
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders expand"
ParseHeaders.py expand -s PRMR_primers-pass_primers-pass.fastq \
    -f PRIMER --outname PRMR >> $PIPELINE_LOG

# Rename primer fields
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders expand"
ParseHeaders.py rename -s PRMR_reheader.fastq -f PRIMER1 PRIMER2 PRIMER3 \
	-k MID VPRIMER CPRIMER --outname FIN >> $PIPELINE_LOG

# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
CollapseSeq.py -s FIN_reheader.fastq -n 20 --uf MID CPRIMER --inner \
    --outname FIN >> $PIPELINE_LOG

# Filter to sequences with at least 2 supporting reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s FIN_collapse-unique.fastq -f DUPCOUNT --num 2 \
    >> $PIPELINE_LOG

# Split file by MID
#printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
#SplitSeq.py group -s FIN_collapse-unique.fastq -f MID --outname FIN \
#    >> $PIPELINE_LOG

# Create tables of final repertoire files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s FIN_reheader.fastq \
    -f ID MID CPRIMER VPRIMER --outname Final \
    >> $PIPELINE_LOG
ParseHeaders.py table -s FIN_collapse-unique.fastq \
    -f ID MID CPRIMER VPRIMER DUPCOUNT --outname Final-Unique \
    >> $PIPELINE_LOG
ParseHeaders.py table -s FIN_collapse-unique_atleast-2.fastq \
    -f ID MID CPRIMER VPRIMER DUPCOUNT --outname Final-Unique-Atleast2 \
    >> $PIPELINE_LOG

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l LengthLog.log -f ID LENGTH > /dev/null &
ParseLog.py -l QualityLog.log -f ID QUALITY > /dev/null &
ParseLog.py -l PrimerMIDLog.log PrimerForwardLog.log PrimerReverseLog.log \
    -f ID PRIMER ERROR > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls *Log.log)
    tar -zcf LogFiles.tar $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP

    TEMP_FILES_ZIP=$(ls *.fastq | grep -vP "FIN_reheader.fastq\|FIN_collapse-unique.fastq\|FIN_collapse-unique_atleast-2.fastq")
    tar -zcf TempFiles.tar $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
fi

# End
printf "DONE\n\n"
cd ../
