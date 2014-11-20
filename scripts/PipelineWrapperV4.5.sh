#!/bin/bash
# Wrapper script to run the pRESTO pipeline script on multiple inputs
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2014.11.20

DATADIR=/scratch2/kleinstein/shlomchik_salmonella
SCRIPT=/scratch2/kleinstein/shlomchik_salmonella/scripts/RunPipelineV4.5_AbSeqV3.sh
PRIMERS1=/scratch2/kleinstein/shlomchik_salmonella/primers/AbSeqV3_Mouse_R1CPrimers.fasta
PRIMERS2=/scratch2/kleinstein/shlomchik_salmonella/primers/AbSeqV3_Mouse_R2TSPrimer.fasta
RUNID=TG725
LOGFILE=${RUNID}_RunLog.out
FOLDERS=$(ls -d $DATADIR/data/$RUNID/TG725*HC| xargs -n 1 basename)
NPROC=20

echo "" > $LOGFILE 
for F in $FOLDERS
do
    echo "FOLDER: $F" | tee -a $LOGFILE 
    echo `date` | tee -a $LOGFILE
    R1=$DATADIR/data/$RUNID/$F/*L001_R1_001.fastq
    R2=$DATADIR/data/$RUNID/$F/*L001_R2_001.fastq
    OUTDIR=$DATADIR/results/$RUNID/$F
    OUTNAME=$F
    $SCRIPT $R1 $R2 $PRIMERS1 $PRIMERS2 $OUTDIR $OUTNAME $NPROC | tee -a $LOGFILE
done
