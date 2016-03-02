#!/bin/bash
# Example pRESTO pipeline for stranded Stern2014 data without UMIs.
# Data from Greiff et al, 2014, BMC Immunol.
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2016.03.01

# Define run parameters and input files
R1_FILE=$(readlink -f ERR346600_1.fastq)
R2_FILE=$(readlink -f ERR346600_2.fastq)
R1_PRIMERS=$(readlink -f Greiff2014_CPrimers.fasta)
R2_PRIMERS=$(readlink -f Greiff2014_VPrimers.fasta)
OUTDIR="output"
PIPELINE_LOG="Pipeline.log"
ZIP_FILES=true

# Make output directory and empty log files
mkdir -p $OUTDIR; cd $OUTDIR
echo '' > $PIPELINE_LOG

# Start
echo "OUTPUT DIRECTORY: ${OUTDIR}"
echo -e "START"
STEP=0

# Assemble paired ends via mate-pair alignment
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
AssemblePairs.py align -1 $R2_FILE -2 $R1_FILE --coord sra --rc tail \
	--minlen 8 --maxerror 0.3 --alpha 1e-5 --scanrev \
	--outname M1 --outdir . --log AP.log >> $PIPELINE_LOG

# Remove low quality reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
FilterSeq.py quality -s M1_assemble-pass.fastq -q 20 \
	--log FS.log >> $PIPELINE_LOG

# Identify primers and UIDs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
MaskPrimers.py score -s M1*quality-pass.fastq -p $R2_PRIMERS \
	--mode mask --start 4 --maxerror 0.2 \
	--log MP1.log >> $PIPELINE_LOG
MaskPrimers.py score -s M1*primers-pass.fastq -p $R1_PRIMERS \
	--mode cut --start 4 --maxerror 0.2 --revpr \
	--log MP2.log >> $PIPELINE_LOG

# Expand primer field
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders expand"
ParseHeaders.py expand -s M1*primers-pass_primers-pass.fastq \
    -f PRIMER >> $PIPELINE_LOG

# Rename primer fields
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders rename"
ParseHeaders.py rename -s M1*reheader.fastq -f PRIMER1 PRIMER2 \
	-k VPRIMER CPRIMER --outname M1-FINAL >> $PIPELINE_LOG
	
# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
CollapseSeq.py -s M1-FINAL_reheader.fastq -n 20 \
	--uf CPRIMER --cf VPRIMER --act set --inner \
	--outname M1-FINAL >> $PIPELINE_LOG

# Subset to sequences observed at least twice
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s M1-FINAL_collapse-unique.fastq -f DUPCOUNT --num 2 \
	--outname M1-FINAL >> $PIPELINE_LOG

# Create annotation table of final unique sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s M1-FINAL_atleast-2.fastq \
    -f ID CPRIMER DUPCOUNT >> $PIPELINE_LOG

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE > /dev/null &
ParseLog.py -l FS.log -f ID QUALITY > /dev/null &
ParseLog.py -l MP[1-2].log -f ID BARCODE PRIMER ERROR > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls AP.log FS.log MP[1-2].log)
    tar -zcf LogFiles.tar $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP

    TEMP_FILES_ZIP=$(ls *.fastq | grep -v "collapse-unique.fastq\|atleast-2.fastq")
    tar -zcf TempFiles.tar $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
fi

# End
printf "DONE\n\n"
cd ../
