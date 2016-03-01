#!/bin/bash
# Example pRESTO pipeline for Stern2014 data with UIDs.
# Data from Stern et al, 2014, Sci Trans Med.
#
# Author:  Jason Anthony Vander Heiden, Gur Yaari
# Date:    2016.03.01

# Define run parameters and input files
R1_FILE=$(readlink -f SRR1383456_1.fastq)
R2_FILE=$(readlink -f SRR1383456_2.fastq)
R1_PRIMERS=$(readlink -f Stern2014_CPrimers.fasta)
R2_PRIMERS=$(readlink -f Stern2014_VPrimers.fasta)
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

# Remove low quality reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
FilterSeq.py quality -s $R1_FILE -q 20 --outname R1 --outdir . \
	--log QualityLogR1.log >> $PIPELINE_LOG
FilterSeq.py quality -s $R2_FILE -q 20 --outname R2 --outdir . \
	--log QualityLogR2.log >> $PIPELINE_LOG

# Identify primers and UIDs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
MaskPrimers.py score -s R1_quality-pass.fastq -p $R1_PRIMERS \
	--mode cut --barcode --start 15 --maxerror 0.2 --outname R1 \
	--log PrimerLogR1.log >> $PIPELINE_LOG
MaskPrimers.py score -s R2_quality-pass.fastq -p $R2_PRIMERS \
	--mode mask --start 0 --maxerror 0.2 --outname R2 \
	--log PrimerLogR2.log >> $PIPELINE_LOG

# Assign UID to read 2 sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 R1_primers-pass.fastq -2 R2_primers-pass.fastq \
    --1f BARCODE --coord sra >> $PIPELINE_LOG

# Multiple align UID read groups
#printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AlignSets muscle"
#MUSCLE_EXEC=$(which muscle)
#AlignSets.py muscle -s R1_primers-pass_pair-pass.fastq --exec $MUSCLE_EXEC \
#	--log AlignLogR1.log --outname R1 >> $PIPELINE_LOG
#AlignSets.py muscle -s R2_primers-pass_pair-pass.fastq --exec $MUSCLE_EXEC \
#	--log AlignLogR2.log --outname R2 >> $PIPELINE_LOG

# Build UID consensus sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "BuildConsensus"
BuildConsensus.py -s R1_primers-pass_pair-pass.fastq --bf BARCODE --pf PRIMER \
	--prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname R1 \
	--log ConsensusLogR1.log >> $PIPELINE_LOG
BuildConsensus.py -s R2_primers-pass_pair-pass.fastq --bf BARCODE --pf PRIMER \
	--maxerror 0.1 --maxgap 0.5 --outname R2 \
	--log ConsensusLogR2.log >> $PIPELINE_LOG

# Synchronize consensus sequence files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 R1_consensus-pass.fastq -2 R2_consensus-pass.fastq \
    --coord presto >> $PIPELINE_LOG

# Assemble paired ends via mate-pair alignment
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
AssemblePairs.py align -1 R2_consensus-pass_pair-pass.fastq \
	-2 R1_consensus-pass_pair-pass.fastq --1f CONSCOUNT --2f PRCONS CONSCOUNT \
	--coord presto --rc tail --minlen 8 --maxerror 0.3 --alpha 1e-5 --scanrev \
	--outname FIN --log AssembleLog.log >> $PIPELINE_LOG

# Rewrite header with minimum of CONSCOUNT
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders collapse"
ParseHeaders.py collapse -s FIN_assemble-pass.fastq -f CONSCOUNT --act min \
    --outname FIN > /dev/null

# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
CollapseSeq.py -s FIN_reheader.fastq -n 20 \
	--uf PRCONS --cf CONSCOUNT --act sum --inner \
	--outname FIN >> $PIPELINE_LOG

# Filter to sequences with at least 2 supporting reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s FIN_collapse-unique.fastq -f CONSCOUNT --num 2 \
    >> $PIPELINE_LOG

# Create tables of final repertoire files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s FIN_reheader.fastq \
    -f ID PRCONS CONSCOUNT --outname Final \
    >> $PIPELINE_LOG
ParseHeaders.py table -s FIN_collapse-unique.fastq \
    -f ID PRCONS CONSCOUNT DUPCOUNT --outname Final-Unique \
    >> $PIPELINE_LOG
ParseHeaders.py table -s FIN_collapse-unique_atleast-2.fastq \
    -f ID PRCONS CONSCOUNT DUPCOUNT --outname Final-Unique-Atleast2 \
    >> $PIPELINE_LOG

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l QualityLogR[1-2].log -f ID QUALITY > /dev/null &
ParseLog.py -l PrimerLogR[1-2].log -f ID BARCODE PRIMER ERROR > /dev/null &
ParseLog.py -l ConsensusLogR[1-2].log -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR \
	> /dev/null &
ParseLog.py -l AssembleLog.log -f ID LENGTH OVERLAP ERROR PVALUE FIELDS1 FIELDS2 \
    > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls *LogR[1-2].log *Log.log)
    tar -zcf LogFiles.tar $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP

    TEMP_FILES_ZIP=$(ls *.fastq | grep -v "FIN_reheader.fastq\|FIN_collapse-unique.fastq\|FIN_collapse-unique_atleast-2.fastq")
    tar -zcf TempFiles.tar $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
fi

# End
printf "DONE\n\n"
cd ../
