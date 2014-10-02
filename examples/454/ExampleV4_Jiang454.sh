#!/bin/bash
# Example script to run the pRESTO pipeline for 454 data
# 
# Author:  Jason Anthony Vander Heiden, Gur Yaari
# Date:    2014.3.18

# Define run parameters
ZIP_FILES=true
MINLEN=250
MINQUAL=20
NPROC=16

# Define input and log files
NAME="SRR765688"
FILE="SRR765688.fastq"
MIDS_PRIMER_FILE="SRR765688_MIDs.fasta"
V_PRIMER_FILE="SRX190717_VPrimers.fasta"
C_PRIMER_FILE="SRX190717_CPrimers.fasta"
RUNLOG="Pipeline.log"

# Start
echo '' > $RUNLOG
echo "INPUT:" $FILE

# Filter low quality reads 
echo "   1: FilterSeq length     $(date +'%H:%M %D')"
FilterSeq.py length -s $FILE -n $MINLEN --nproc $NPROC \
    --log FilterLen.log --clean --outdir . >> $RUNLOG

echo "   2: FilterSeq quality    $(date +'%H:%M %D')"
FilterSeq.py quality -s $NAME*length-pass.fastq -q $MINQUAL --nproc $NPROC \
    --log FilterQual.log --clean >> $RUNLOG

# Remove sample barcode
echo "   3: MaskPrimers score    $(date +'%H:%M %D')"
MaskPrimers.py score -s $NAME*quality-pass.fastq -p $MIDS_PRIMER_FILE \
    --mode cut --start 0 --nproc $NPROC --log PrimerMID.log --clean >> $RUNLOG

# Remove V-region primers
echo "   4: MaskPrimers align    $(date +'%H:%M %D')"
MaskPrimers.py align -s $NAME*quality-pass_primers-pass.fastq -p $V_PRIMER_FILE \
    --mode cut --maxlen 50 --nproc $NPROC --log PrimerV.log --clean >> $RUNLOG

# Remove C-region primers
echo "   5: MaskPrimers align    $(date +'%H:%M %D')"
MaskPrimers.py align -s $NAME*primers-pass_primers-pass.fastq -p $C_PRIMER_FILE \
    --mode cut --maxlen 50 --revpr --skiprc --nproc $NPROC \
    --outname "${NAME}-ENDMP" --log PrimerC.log --clean >> $RUNLOG

# Split by primers and group by barcode
echo "   6: ParseHeaders expand  $(date +'%H:%M %D')" 
ParseHeaders.py expand -s "${NAME}-ENDMP_primers-pass.fastq" -f PRIMER > /dev/null

echo "   7: ParseHeaders rename  $(date +'%H:%M %D')"
ParseHeaders.py rename -s "${NAME}-ENDMP_primers-pass_reheader.fastq" \
    -f PRIMER1 PRIMER2 PRIMER3 -k BARCODE VPRIMER CPRIMER --outname "${NAME}-ENDMP" > /dev/null

echo "   8: SplitSeq group       $(date +'%H:%M %D')"
SplitSeq.py group -s "${NAME}-ENDMP_reheader.fastq" -f BARCODE --fasta >> $RUNLOG

# Filter duplicates and singletons
echo "   9: CollapseSeq          $(date +'%H:%M %D')"
CollapseSeq.py -s *reheader_MID*.fasta -n 10 --inner --uf CPRIMER \
    --cf BARCODE VPRIMER --act set set >> $RUNLOG

echo "  10: SplitSeq group       $(date +'%H:%M %D')"
SplitSeq.py group -s *unique.fasta -f DUPCOUNT --num 2 >> $RUNLOG

# Process log files
echo "  11: ParseLog             $(date +'%H:%M %D')"
ParseLog.py -l FilterLen.log -f ID LENGTH > /dev/null &
ParseLog.py -l FilterQual.log -f ID QUALITY > /dev/null &
ParseLog.py -l Primer*.log -f ID PRSTART PRIMER ERROR > /dev/null &
ParseLog.py -l Pipeline.log -f END SEQUENCES PAIRS SETS PASS FAIL UNIQUE DUPLICATE UNDETERMINED PARTS OUTPUT > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
	tar -cf LogFiles.tar *.log 
	gzip LogFiles.tar
	rm *.log
	
	tar -cf TempFiles.tar *.fastq *under* *duplicate* *undetermined*
    gzip TempFiles.tar
    rm *.fastq *under* *duplicate* *undetermined*
fi

# End
echo -e "DONE\n" 
cd ../
