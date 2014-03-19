#!/bin/bash
# Wrapper script to run the pRESTO pipeline script on multiple inputs
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2014.3.19

LOGFILE=run.out
DATADIR=/scratch2/kleinstein/shlomchik_salmonella
SAMPLE=TG873
FOLDERS=$(ls -d $DATADIR/data/$SAMPLE/Sample_*| xargs -n 1 basename)
NPROC=20

echo "" > $LOGFILE 
for F in $FOLDERS
do
  echo "FOLDER: $F" | tee -a $LOGFILE 
  echo `date` | tee -a $LOGFILE
  R1=$DATADIR/data/$SAMPLE/$F/*L001_R1_001.fastq
  R2=$DATADIR/data/$SAMPLE/$F/*L001_R2_001.fastq
  OUT=$DATADIR/results/$SAMPLE/$F
$DATADIR/scripts/RunPipelineV4_AbVitroV3.0.sh $R1 $R2 $OUT $NPROC | tee -a $LOGFILE
done
