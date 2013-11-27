#!/bin/bash
# Super script to run the pRESTO pipeline on AbVitro library v3.0 data 
# 
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2013.11.25
# 
# Required Arguments:
#   $1 = name of the directory containing read files 

# Define run parameters
ZIP_FILES=true
ALIGN_STEP=false
PRCONS=0.6
MAXDIV=0.1
FSQUAL=20
BCQUAL=0
R1_MAXERR=0.2 
R2_MAXERR=0.2
AP_MAXERR=0.2 
ALPHA=0.1
FS_MISS=20
CS_MISS=0
NPROC=30
   
# Define input files
FILEDIR=$1
R1_SUFFIX="R1_001"
R2_SUFFIX="R2_001"
R1_PRIMER_FILE='/scratch2/kleinstein/shlomchik_salmonella/primers/JPrimers_v2.fasta'
R2_PRIMER_FILE='/scratch2/kleinstein/shlomchik_salmonella/primers/VPrimers_shifted_unique.fasta'

# Define script execution command and log files
RUNLOG="/$(pwd)/${FILEDIR}/Pipeline.log"
RUN="nice -19"
#TIMELOG="/$(pwd)/${FILEDIR}/Runtime.log"
#RUN="nice -19 /usr/bin/time -o ${TIMELOG} -a -f '%C\t%E\t%P\t%Mkb'"

# Start
cd $FILEDIR
echo "DIRECTORY:" $FILEDIR

# Filter low quality reads
echo "   1: FilterSeq quality      $(date +'%H:%M %D')"
$RUN FilterSeq.py quality -s *"${R1_SUFFIX}.fastq" -q $FSQUAL --nproc $NPROC --clean --outname R1 >> $RUNLOG
$RUN FilterSeq.py quality -s *"${R2_SUFFIX}.fastq" -q $FSQUAL --nproc $NPROC --clean --outname R2 >> $RUNLOG

# Identify primers and UID 
echo "   2: MaskPrimers score      $(date +'%H:%M %D')"
$RUN MaskPrimers.py score -s R1_quality-pass.fastq -p $R1_PRIMER_FILE --mode cut --start 0 \
    --maxerror $R1_MAXERR --nproc $NPROC --log PrimerLogR1.log --clean >> $RUNLOG
$RUN MaskPrimers.py score -s R2_quality-pass.fastq -p $R2_PRIMER_FILE --mode cut --barcode --start 17 \
    --maxerror $R1_MAXERR --nproc $NPROC --log PrimerLogR2.log --clean >> $RUNLOG
   
# Assign UIDs to read 1 sequences
echo "   3: PairSeq                $(date +'%H:%M %D')"
$RUN PairSeq.py -1 R2*primers-pass.fastq -2 R1*primers-pass.fastq -f BARCODE --coord illumina --clean >> $RUNLOG

if $LEN_STEP; then
	# Multiple align UID read groups
	echo "   4: AlignSets muscle       $(date +'%H:%M %D')"     
	$RUN AlignSets.py muscle -s R1*pair-pass.fastq --exec $MUSCLE_EXEC --nproc $NPROC \
	    --log AlignLogR1.log --clean >> $RUNLOG
	$RUN AlignSets.py muscle -s R2*pair-pass.fastq --exec $MUSCLE_EXEC --nproc $NPROC \
	    --log AlignLogR2.log --clean >> $RUNLOG
	
	# Build UID consensus sequences
	echo "   5: BuildConsensus         $(date +'%H:%M %D')" 
	$RUN BuildConsensus.py -s R1*align-pass.fastq --bf BARCODE --pf PRIMER --prcons $PRCONS \
	    -q $BCQUAL --maxdiv $MAXDIV --nproc $NPROC --log ConsensusLogR1.log --clean >> $RUNLOG
	$RUN BuildConsensus.py -s R2*align-pass.fastq --bf BARCODE --pf PRIMER --prcons $PRCONS \
	    -q $BCQUAL --maxdiv $MAXDIV --nproc $NPROC --log ConsensusLogR2.log --clean >> $RUNLOG
else
	# Build UID consensus sequences
	echo "   5: BuildConsensus         $(date +'%H:%M %D')" 
	$RUN BuildConsensus.py -s R1*pair-pass.fastq --bf BARCODE --pf PRIMER --prcons $PRCONS \
	    -q $BCQUAL --maxdiv $MAXDIV --nproc $NPROC --log ConsensusLogR1.log --clean >> $RUNLOG
	$RUN BuildConsensus.py -s R2*pair-pass.fastq --bf BARCODE --pf PRIMER \
	    -q $BCQUAL --maxdiv $MAXDIV --nproc $NPROC --log ConsensusLogR2.log --clean >> $RUNLOG
fi

# Assemble paired ends
echo "   6: AssemblePairs          $(date +'%H:%M %D')" 
$RUN AssemblePairs.py align -1 R2*consensus-pass.fastq -2 R1*consensus-pass.fastq \
    --1f CONSCOUNT --2f PRCONS CONSCOUNT --coord presto --rc tail --maxerror $AP_MAXERR --alpha $ALPHA \
    --nproc $NPROC --log AssembleLog.log --clean >> $RUNLOG

# Remove sequences with many Ns
echo "   7: FilterSeq missing      $(date +'%H:%M %D')" 
$RUN FilterSeq.py missing -s *assemble-pass.fastq -n $FS_MISS --inner \
    --nproc $NPROC >> $RUNLOG

# Rewrite header with minimum of CONSCOUNT
echo "   8: ParseHeaders collapse  $(date +'%H:%M %D')"
$RUN ParseHeaders.py collapse -s *missing-pass.fastq -f CONSCOUNT --act min \
    --outname Assembled --fasta > /dev/null

# Remove duplicate sequences
echo "   9: CollapseSeq            $(date +'%H:%M %D')" 
$RUN CollapseSeq.py -s Assembled_reheader.fasta -n $CS_MISS --uf PRCONS \
    --cf CONSCOUNT --act sum --outname Assembled --inner >> $RUNLOG

# Filter to sequences with at least 2 supporting sources
echo "  10: SplitSeq group         $(date +'%H:%M %D')" 
$RUN SplitSeq.py group -s Assembled_collapse-unique.fasta -f CONSCOUNT --num 2 >> $RUNLOG

# Create table of final repertoire
echo "  11: ParseHeaders table     $(date +'%H:%M %D')"
$RUN ParseHeaders.py table -s Assembled_collapse-unique_atleast-2.fasta -f ID PRCONS CONSCOUNT DUPCOUNT >> $RUNLOG

# Process log files
echo "  12: ParseLog               $(date +'%H:%M %D')"
$RUN ParseLog.py -l PrimerLogR[1-2].log -f BARCODE PRIMER ERROR > /dev/null &
$RUN ParseLog.py -l ConsensusLogR[1-2].log -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCOUNT PRFREQ DIVERSITY > /dev/null &
$RUN ParseLog.py -l AssembleLog.log -f ID OVERLAP LENGTH PVAL ERROR > /dev/null &
$RUN ParseLog.py -l Pipeline.log -f END SEQUENCES PAIRS SETS PASS FAIL PARTS OUTPUT > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    tar -cf LogFiles.tar *.log 
    gzip LogFiles.tar
    rm *.log
    
    tar -cf TempFiles.tar R[1-2]_*.fastq *under* *duplicate* *undetermined*
    gzip TempFiles.tar
    rm R[1-2]_*.fastq *under* *duplicate* *undetermined*
fi
    
# End
echo -e "DONE\n" 
cd ../

