#!/bin/bash
# Wrapper script to run the pRESTO pipeline script on multiple inputs
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2014.11.18

DATADIR=/scratch2/kleinstein/shlomchik_salmonella
SCRIPT=/scratch2/kleinstein/shlomchik_salmonella/scripts/RunPipelineV4.5_AbSeqV3.sh
RUNID=RQ2410
LOGFILE=${RUNID}_RunLog.out
FOLDERS=$(ls -d $DATADIR/data/$RUNID/Sample_*| xargs -n 1 basename)
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
    $SCRIPT $R1 $R2 $OUTDIR $OUTNAME $NPROC | tee -a $LOGFILE
done
