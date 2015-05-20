#!/bin/bash
# Wrapper script to run the pRESTO pipeline script on multiple inputs
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2014.11.20

# Define run parameters
DATA_DIR=/scratch2/kleinstein/shlomchik_salmonella
SCRIPT=/scratch2/kleinstein/shlomchik_salmonella/scripts/PrestoPipelineV4.6_AbSeqV3.sh
PRIMERS1=/scratch2/kleinstein/shlomchik_salmonella/primers/AbSeqV3_Mouse_R1CPrimers.fasta
PRIMERS2=/scratch2/kleinstein/shlomchik_salmonella/primers/AbSeqV3_Mouse_R2TSPrimer.fasta
RUN_ID=TG725
LOG_FILE=${RUN_ID}_RunLog.out
FOLDERS=$(ls -d $DATA_DIR/raw/$RUN_ID/TG725*HC| xargs -n 1 basename)
NPROC=20

echo "" > $LOG_FILE 
for F in $FOLDERS
do
    echo "FOLDER: $F" | tee -a $LOG_FILE 
    echo `date` | tee -a $LOG_FILE
    R1=$DATA_DIR/raw/$RUN_ID/$F/*L001_R1_001.fastq
    R2=$DATA_DIR/raw/$RUN_ID/$F/*L001_R2_001.fastq
    OUTDIR=$DATA_DIR/presto/$RUN_ID/$F
    OUTNAME=$F
    $SCRIPT $R1 $R2 $PRIMERS1 $PRIMERS2 $OUTDIR $OUTNAME $NPROC | tee -a $LOG_FILE
done
