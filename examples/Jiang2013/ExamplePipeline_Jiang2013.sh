#!/usr/bin/env bash
# Example pRESTO pipeline for Roche 454 data with sample multiplexing
# Data from Jiang, He and Weinstein et al, 2013, Sci Trans Med.
#
# Author:  Jason Anthony Vander Heiden
# Date:    2018.03.14

# Define run parameters and input files
READ_FILE=$(readlink -f SRR765688.fastq)
MID_PRIMERS=$(readlink -f SRR765688_MIDs.fasta)
FWD_PRIMERS=$(readlink -f SRX190717_VPrimers.fasta)
REV_PRIMERS=$(readlink -f SRX190717_CPrimers.fasta)
OUTDIR="output"
OUTNAME="S43"
NPROC=4
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
FilterSeq.py length -s $READ_FILE -n 300 \
	--outname $OUTNAME --log FSL.log --nproc $NPROC --outdir . >> $PIPELINE_LOG

# Remove low quality reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
FilterSeq.py quality -s "${OUTNAME}_length-pass.fastq" -q 20 \
	--outname $OUTNAME --log FSQ.log --nproc $NPROC >> $PIPELINE_LOG

# Identify and remove MIDs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
MaskPrimers.py score -s "${OUTNAME}_quality-pass.fastq" -p $MID_PRIMERS \
    --start 0 --maxerror 0.1 --mode cut --pf MID \
    --outname "${OUTNAME}-MID" --log MPM.log --nproc $NPROC >> $PIPELINE_LOG

# Identify and mask forward (V-region) primers
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s "${OUTNAME}-MID_primers-pass.fastq" -p $FWD_PRIMERS \
    --maxlen 50 --maxerror 0.3 --mode mask --pf VPRIMER \
    --outname "${OUTNAME}-FWD" --log MPV.log --nproc $NPROC >> $PIPELINE_LOG

# Identify and remove reverse (C-region) primers
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s "${OUTNAME}-FWD_primers-pass.fastq" -p $REV_PRIMERS \
    --maxlen 50 --maxerror 0.3 --mode cut --revpr --skiprc --pf CPRIMER \
	--outname "${OUTNAME}-REV" --log MPC.log --nproc $NPROC >> $PIPELINE_LOG

# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
CollapseSeq.py -s "${OUTNAME}-REV_primers-pass.fastq" -n 20 --uf MID CPRIMER \
	--inner --outname $OUTNAME >> $PIPELINE_LOG

# Filter to sequences with at least 2 supporting reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s "${OUTNAME}_collapse-unique.fastq" -f DUPCOUNT --num 2 \
    --outname $OUTNAME >> $PIPELINE_LOG

# Split file by MID
#printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
#SplitSeq.py group -s "${OUTNAME}_collapse-unique.fastq" -f MID >> $PIPELINE_LOG

# Create tables of final repertoire files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s "${OUTNAME}_atleast-2.fastq" \
    -f ID MID CPRIMER VPRIMER DUPCOUNT >> $PIPELINE_LOG

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l FSL.log -f ID LENGTH > /dev/null &
ParseLog.py -l FSQ.log -f ID QUALITY > /dev/null &
ParseLog.py -l MP[MVC].log -f ID PRSTART PRIMER ERROR > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls FS?.log MP?.log)
    tar -zcf LogFiles.tar.gz $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP

    TEMP_FILES_ZIP=$(ls *.fastq | grep -vP "collapse-unique.fastq|atleast-2.fastq")
    tar -zcf TempFiles.tar.gz $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
fi

# End
printf "DONE\n\n"
cd ../
