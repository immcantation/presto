#!/usr/bin/env bash
# Example pRESTO pipeline for UMI barcoded Illumina Miseq 2x250 data
# Data from Stern, Yaari and Vander Heiden et al, 2014, Sci Trans Med.
#
# Author:  Jason Anthony Vander Heiden, Gur Yaari
# Date:    2016.03.04

# Define run parameters and input files
R1_FILE=$(realpath SRR1383456_1.fastq)
R2_FILE=$(realpath SRR1383456_2.fastq)
R1_PRIMERS=$(realpath Stern2014_CPrimers.fasta)
R2_PRIMERS=$(realpath Stern2014_VPrimers.fasta)
OUTDIR="output"
OUTNAME="MS12"
NPROC=4
PIPELINE_LOG="Pipeline.log"
ZIP_FILES=true

# Make output directory and empty log files
mkdir -p $OUTDIR; cd $OUTDIR
echo '' > $PIPELINE_LOG

# Start
echo "OUTPUT DIRECTORY: ${OUTDIR}"
echo -e "START"
STEP=0

# Remove low quality reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
FilterSeq.py quality -s $R1_FILE -q 20 --outname "${OUTNAME}_R1" --outdir . \
	--log FS1.log --nproc $NPROC >> $PIPELINE_LOG
FilterSeq.py quality -s $R2_FILE -q 20 --outname "${OUTNAME}_R2" --outdir . \
	--log FS2.log --nproc $NPROC >> $PIPELINE_LOG

# Identify primers and UIDs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
MaskPrimers.py score -s "${OUTNAME}_R1_quality-pass.fastq" -p $R1_PRIMERS \
	--mode cut --barcode --start 15 --maxerror 0.2 --outname "${OUTNAME}_R1" \
	--log MP1.log --nproc $NPROC >> $PIPELINE_LOG
MaskPrimers.py score -s "${OUTNAME}_R2_quality-pass.fastq" -p $R2_PRIMERS \
	--mode mask --start 0 --maxerror 0.2 --outname "${OUTNAME}_R2" \
	--log MP2.log --nproc $NPROC >> $PIPELINE_LOG

# Assign UID to read 2 sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 "${OUTNAME}_R1_primers-pass.fastq" -2 "${OUTNAME}_R2_primers-pass.fastq" \
    --1f BARCODE --coord sra >> $PIPELINE_LOG

# Multiple align UID read groups
#printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AlignSets muscle"
#MUSCLE_EXEC=$(which muscle)
#AlignSets.py muscle -s "${OUTNAME}_R1_primers-pass_pair-pass.fastq" --exec $MUSCLE_EXEC \
#	--log AS1.log --outname "${OUTNAME}_R1" >> $PIPELINE_LOG
#AlignSets.py muscle -s "${OUTNAME}_R2_primers-pass_pair-pass.fastq" --exec $MUSCLE_EXEC \
#	--log AS2.log --outname "${OUTNAME}_R2" >> $PIPELINE_LOG

# Build UID consensus sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "BuildConsensus"
BuildConsensus.py -s "${OUTNAME}_R1_primers-pass_pair-pass.fastq" --bf BARCODE --pf PRIMER \
	--prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname "${OUTNAME}_R1" \
	--log BC1.log --nproc $NPROC >> $PIPELINE_LOG
BuildConsensus.py -s "${OUTNAME}_R2_primers-pass_pair-pass.fastq" --bf BARCODE --pf PRIMER \
	--maxerror 0.1 --maxgap 0.5 --outname "${OUTNAME}_R2" \
	--log BC2.log --nproc $NPROC >> $PIPELINE_LOG

# Synchronize consensus sequence files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 "${OUTNAME}_R1_consensus-pass.fastq" -2 "${OUTNAME}_R2_consensus-pass.fastq" \
    --coord presto >> $PIPELINE_LOG

# Assemble paired ends via mate-pair alignment
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
AssemblePairs.py align -1 "${OUTNAME}_R2_consensus-pass_pair-pass.fastq" \
	-2 "${OUTNAME}_R1_consensus-pass_pair-pass.fastq" --1f CONSCOUNT --2f PRCONS CONSCOUNT \
	--coord presto --rc tail --outname $OUTNAME --log AP.log --nproc $NPROC >> $PIPELINE_LOG

# Rewrite header with minimum of CONSCOUNT
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders collapse"
ParseHeaders.py collapse -s "${OUTNAME}_assemble-pass.fastq" -f CONSCOUNT --act min \
    --outname $OUTNAME > /dev/null

# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
CollapseSeq.py -s "${OUTNAME}_reheader.fastq" -n 20 \
	--uf PRCONS --cf CONSCOUNT --act sum --inner \
	--outname $OUTNAME >> $PIPELINE_LOG

# Filter to sequences with at least 2 supporting reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s "${OUTNAME}_collapse-unique.fastq" -f CONSCOUNT --num 2 \
    --outname $OUTNAME >> $PIPELINE_LOG

# Create tables of final repertoire files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s "${OUTNAME}_atleast-2.fastq" -f ID PRCONS CONSCOUNT DUPCOUNT \
    >> $PIPELINE_LOG

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l FS[1-2].log -f ID QUALITY > /dev/null &
ParseLog.py -l MP[1-2].log -f ID BARCODE PRIMER ERROR > /dev/null &
ParseLog.py -l BC[1-2].log -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR \
	> /dev/null &
ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE FIELDS1 FIELDS2 \
    > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls FS[1-2].log MP[1-2].log BC[1-2].log AP.log)
    tar -zcf LogFiles.tar.gz $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP

    TEMP_FILES_ZIP=$(ls *.fastq | grep -vP "collapse-unique.fastq|atleast-2.fastq")
    tar -zcf TempFiles.tar.gz $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
fi

# End
printf "DONE\n\n"
cd ../
