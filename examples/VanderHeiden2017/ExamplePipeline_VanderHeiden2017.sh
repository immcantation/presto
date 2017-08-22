#!/usr/bin/env bash
# Example pRESTO pipeline for UMI barcoded Illumina Miseq 325+275 5'RACE data
# Data from Vander Heiden et al, 2017, J Immunol.
#
# Author:  Jason Anthony Vander Heiden
# Date:    2017.08.21

# Define run parameters and input files
R1_FILE=$(readlink -f SRR4026043_1.fastq)
R2_FILE=$(readlink -f SRR4026043_2.fastq)
R1_PRIMERS=$(readlink -f AbSeq_R1_Human_IG_Primers.fasta)
R2_PRIMERS=$(readlink -f AbSeq_R2_TS.fasta)
CREGION_FILE=$(readlink -f AbSeq_Human_IG_InternalCRegion.fasta)
VREF_FILE=$(readlink -f IMGT_Human_IG_V.fasta)
OUTDIR="output"
OUTNAME="HD09N"
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
FilterSeq.py quality -s $R1_FILE -q 20 --outname "${OUTNAME}-R1" --outdir . \
	--log FS1.log --nproc $NPROC >> $PIPELINE_LOG
FilterSeq.py quality -s $R2_FILE -q 20 --outname "${OUTNAME}-R2" --outdir . \
	--log FS2.log --nproc $NPROC >> $PIPELINE_LOG

# Identify primers and UIDs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
MaskPrimers.py score -s "${OUTNAME}-R1_quality-pass.fastq" -p $R1_PRIMERS \
	--mode cut --start 0 --maxerror 0.2 --outname "${OUTNAME}-R1" \
	--log MP1.log --nproc $NPROC >> $PIPELINE_LOG
MaskPrimers.py score -s "${OUTNAME}-R2_quality-pass.fastq" -p $R2_PRIMERS \
	--mode cut --start 17 --barcode --maxerror 0.5 --outname "${OUTNAME}-R2" \
	--log MP2.log --nproc $NPROC >> $PIPELINE_LOG

# Assign UID to read 2 sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 "${OUTNAME}-R1_primers-pass.fastq" -2 "${OUTNAME}-R2_primers-pass.fastq" \
    --2f BARCODE --coord sra >> $PIPELINE_LOG

# Build UID consensus sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "BuildConsensus"
BuildConsensus.py -s "${OUTNAME}-R1_primers-pass_pair-pass.fastq" --bf BARCODE --pf PRIMER \
	--prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname "${OUTNAME}-R1" \
	--log BC1.log --nproc $NPROC >> $PIPELINE_LOG
BuildConsensus.py -s "${OUTNAME}-R2_primers-pass_pair-pass.fastq" --bf BARCODE \
	--maxerror 0.1 --maxgap 0.5 --outname "${OUTNAME}-R2" \
	--log BC2.log --nproc $NPROC >> $PIPELINE_LOG

# Synchronize consensus sequence files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 "${OUTNAME}-R1_consensus-pass.fastq" -2 "${OUTNAME}-R2_consensus-pass.fastq" \
    --coord presto >> $PIPELINE_LOG

# Assemble paired ends via mate-pair alignment
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs sequential"
AssemblePairs.py sequential -1 "${OUTNAME}-R2_consensus-pass_pair-pass.fastq" \
    -2 "${OUTNAME}-R1_consensus-pass_pair-pass.fastq" -r $VREF_FILE \
    --coord presto --rc tail --scanrev --1f CONSCOUNT --2f PRCONS CONSCOUNT \
	--aligner blastn --outname "${OUTNAME}-C" --log AP.log --nproc $NPROC \
	>> $PIPELINE_LOG

# Annotate with internal C-region
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s "${OUTNAME}-C_assemble-pass.fastq" -p $CREGION_FILE \
    --maxlen 100 --maxerror 0.3 --mode tag --revpr --skiprc \
    --log MP3.log --outname "${OUTNAME}-C" --nproc $NPROC \
    >> $PIPELINE_LOG

# Renamer primer field
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders rename"
ParseHeaders.py rename -s "${OUTNAME}-C_primers-pass.fastq" -f PRIMER -k CREGION \
	> /dev/null

# Rewrite header with minimum of CONSCOUNT
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders collapse"
ParseHeaders.py collapse -s "${OUTNAME}-C_primers-pass_reheader.fastq" \
	-f CONSCOUNT --act min > /dev/null
    
# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
CollapseSeq.py -s "${OUTNAME}-C_primers-pass_reheader_reheader.fastq" -n 20 \
	--uf CREGION --cf CONSCOUNT --act sum --inner \
	--outname "${OUTNAME}-C" >> $PIPELINE_LOG

# Filter to sequences with at least 2 supporting reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s "${OUTNAME}-C_collapse-unique.fastq" -f CONSCOUNT --num 2 \
    --outname "${OUTNAME}-C" >> $PIPELINE_LOG

# Create tables of final repertoire files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s "${OUTNAME}-C_atleast-2.fastq" -f ID PRCONS CONSCOUNT DUPCOUNT \
    >> $PIPELINE_LOG

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l FS[1-2].log -f ID QUALITY > /dev/null &
ParseLog.py -l MP[1-3].log -f ID BARCODE PRIMER ERROR > /dev/null &
ParseLog.py -l BC[1-2].log -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR \
	> /dev/null &
ParseLog.py -l AP.log -f ID REFID LENGTH OVERLAP GAP ERROR PVALUE EVALUE1 EVALUE2 IDENTITY FIELDS1 FIELDS2 \
    > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls FS[1-2].log MP[1-3].log BC[1-2].log AP.log)
    tar -zcf LogFiles.tar.gz $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP

    TEMP_FILES_ZIP=$(ls *.fastq | grep -vP "collapse-unique.fastq|atleast-2.fastq")
    tar -zcf TempFiles.tar.gz $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
fi

# End
printf "DONE\n\n"
cd ../
