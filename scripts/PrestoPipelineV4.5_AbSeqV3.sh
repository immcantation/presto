#!/bin/bash
# Super script to run the pRESTO 0.4.5 pipeline on AbVitro AbSeq (V3) data
# 
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2014.12.13
# 
# Required Arguments:
#   $1 = read 1 file (C-region start sequence)
#   $2 = read 2 file (V-region start sequence)
#   $3 = read 1 primer file
#   $4 = read 2 primer file
#   $5 = output directory
#   $6 = output file prefix
#   $7 = number of subprocesses for multiprocessing tools

# Capture command line parameters
R1_FILE=$(readlink -f $1)
R2_FILE=$(readlink -f $2)
R1_PRIMERS=$(readlink -f $3)
R2_PRIMERS=$(readlink -f $4)
OUTDIR=$5
OUTNAME=$6
NPROC=$7

# Define pipeline steps
LOG_RUNTIMES=true
ZIP_FILES=true
FILTER_LOWQUAL=true
ALIGN_UIDSETS=false
REFERENCE_ASSEMBLY=true
MASK_LOWQUAL=false

# FilterSeq run parameters
FS_QUAL=20
FS_MASK=30
FS_MISS=20

# MaskPrimers run parameters
MP_UIDLEN=17
MP_R1_MAXERR=0.2
MP_R2_MAXERR=0.5

# AlignSets run parameters
MUSCLE_EXEC=$HOME/bin/muscle

# BuildConsensus run parameters
BC_PRCONS_FLAG=true
BC_DIV_FLAG=true
BC_MAXDIV=0.1
BC_PRCONS=0.6
BC_QUAL=0

# AssemblePairs-align run parameters
AP_ALN_SCANREV=true
AP_ALN_MAXERR=0.3
AP_ALN_MINLEN=8
AP_ALN_ALPHA=1e-5

# AssemblePairs-reference run parameters
AP_REF_MAXERR=0.5
AP_REF_EVALUE=1e-5
AP_REF_MAXHITS=100
REF_FILE="/scratch2/kleinstein/germlines/IMGT_Human_IGV_2014-08-23.fasta"
#REF_FILE="/scratch2/kleinstein/germlines/IMGT_Mouse_IGV_2014-11-22.fasta"
USEARCH_EXEC=$HOME/bin/usearch

# CollapseSeq run parameters
CS_KEEP=true
CS_MISS=0

# Define log files
RUNLOG="Pipeline.log"
TIMELOG="Time.log"

# Make output directory and empty log files
mkdir -p $OUTDIR; cd $OUTDIR
echo '' > $RUNLOG
if $LOG_RUNTIMES; then
	echo '' > $TIMELOG
	RUN="nice -19 /usr/bin/time -o ${TIMELOG} -a -f %C\t%E\t%P\t%Mkb"
else
	RUN="nice -19"
fi

# Start
echo "DIRECTORY: ${OUTDIR}"
echo "VERSIONS:"
echo "  $(AlignSets.py -v 2>&1)"
echo "  $(AssemblePairs.py -v 2>&1)"
echo "  $(BuildConsensus.py -v 2>&1)"
echo "  $(CollapseSeq.py -v 2>&1)"
echo "  $(FilterSeq.py -v 2>&1)"
echo "  $(MaskPrimers.py -v 2>&1)"
echo "  $(PairSeq.py -v 2>&1)"
echo "  $(ParseHeaders.py -v 2>&1)"
echo "  $(ParseLog.py -v 2>&1)"
echo "  $(SplitSeq.py -v 2>&1)"
echo -e "\nSTART"
STEP=0

# Remove low quality reads
if $FILTER_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
    #OUTPREFIX="$(printf '%02d' $STEP)--${OUTNAME}"
    $RUN FilterSeq.py quality -s $R1_FILE -q $FS_QUAL --nproc $NPROC \
        --outname "${OUTNAME}-R1" --outdir . --clean --log QualityLogR1.log >> $RUNLOG
    $RUN FilterSeq.py quality -s $R2_FILE -q $FS_QUAL --nproc $NPROC \
        --outname "${OUTNAME}-R2" --outdir . --clean --log QualityLogR2.log  >> $RUNLOG
    MPR1_FILE="${OUTNAME}-R1_quality-pass.fastq"
    MPR2_FILE="${OUTNAME}-R2_quality-pass.fastq"
else
    MPR1_FILE=$R1_FILE
    MPR2_FILE=$R2_FILE
fi

# Identify primers and UID 
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
$RUN MaskPrimers.py score -s $MPR1_FILE -p $R1_PRIMERS --mode cut \
    --start 0 --maxerror $MP_R1_MAXERR --nproc $NPROC --log PrimerLogR1.log \
    --outname "${OUTNAME}-R1" --outdir . --clean >> $RUNLOG
$RUN MaskPrimers.py score -s $MPR2_FILE -p $R2_PRIMERS --mode cut \
    --start $MP_UIDLEN --barcode --maxerror $MP_R2_MAXERR --nproc $NPROC --log PrimerLogR2.log \
    --outname "${OUTNAME}-R2" --outdir . --clean >> $RUNLOG

# Assign UIDs to read 1 sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
$RUN PairSeq.py -1 "${OUTNAME}-R2_primers-pass.fastq" -2 "${OUTNAME}-R1_primers-pass.fastq" \
    -f BARCODE --coord illumina --clean >> $RUNLOG

# Multiple align UID read groups
if $ALIGN_UIDSETS; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AlignSets muscle"
	$RUN AlignSets.py muscle -s "${OUTNAME}-R1_primers-pass_pair-pass.fastq" --exec $MUSCLE_EXEC \
	    --nproc $NPROC --log AlignLogR1.log --outname "${OUTNAME}-R1" --clean >> $RUNLOG
	$RUN AlignSets.py muscle -s "${OUTNAME}-R2_primers-pass_pair-pass.fastq" --exec $MUSCLE_EXEC \
	    --nproc $NPROC --log AlignLogR2.log --outname "${OUTNAME}-R2" --clean >> $RUNLOG
	BCR1_FILE="${OUTNAME}-R1_align-pass.fastq"
	BCR2_FILE="${OUTNAME}-R2_align-pass.fastq"
else
	BCR1_FILE="${OUTNAME}-R1_primers-pass_pair-pass.fastq"
	BCR2_FILE="${OUTNAME}-R2_primers-pass_pair-pass.fastq"
fi

# Build UID consensus sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "BuildConsensus"
if $BC_DIV_FLAG; then
    if $BC_PRCONS_FLAG; then
        $RUN BuildConsensus.py -s $BCR1_FILE --bf BARCODE --pf PRIMER --prcons $BC_PRCONS \
            -q $BC_QUAL --maxdiv $BC_MAXDIV --nproc $NPROC --log ConsensusLogR1.log \
            --outname "${OUTNAME}-R1" --clean >> $RUNLOG
    else
        $RUN BuildConsensus.py -s $BCR1_FILE --bf BARCODE --pf PRIMER \
            -q $BC_QUAL --maxdiv $BC_MAXDIV --nproc $NPROC --log ConsensusLogR1.log \
            --outname "${OUTNAME}-R1" --clean >> $RUNLOG
    fi

	$RUN BuildConsensus.py -s $BCR2_FILE --bf BARCODE --pf PRIMER \
	    -q $BC_QUAL --maxdiv $BC_MAXDIV --nproc $NPROC --log ConsensusLogR2.log \
	    --outname "${OUTNAME}-R2" --clean >> $RUNLOG
else
    if $BC_PRCONS_FLAG; then
        $RUN BuildConsensus.py -s $BCR1_FILE --bf BARCODE --pf PRIMER --prcons $BC_PRCONS \
            -q $BC_QUAL --nproc $NPROC --log ConsensusLogR1.log \
            --outname "${OUTNAME}-R1" --clean >> $RUNLOG
    else
        $RUN BuildConsensus.py -s $BCR1_FILE --bf BARCODE --pf PRIMER \
            -q $BC_QUAL --nproc $NPROC --log ConsensusLogR1.log \
            --outname "${OUTNAME}-R1" --clean >> $RUNLOG
    fi

	$RUN BuildConsensus.py -s $BCR2_FILE --bf BARCODE --pf PRIMER \
    	-q $BC_QUAL --nproc $NPROC --log ConsensusLogR2.log \
    	--outname "${OUTNAME}-R2" --clean >> $RUNLOG
fi

# Assemble paired ends via mate-pair alignment
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"

if $BC_PRCONS_FLAG; then
    PRFIELD="PRCONS"
else
    PRFIELD="PRIMER"
fi

if $AP_ALN_SCANREV; then
    $RUN AssemblePairs.py align -1 "${OUTNAME}-R2_consensus-pass.fastq" \
        -2 "${OUTNAME}-R1_consensus-pass.fastq" --1f CONSCOUNT --2f $PRFIELD CONSCOUNT \
        --coord presto --rc tail --minlen $AP_ALN_MINLEN --maxerror $AP_ALN_MAXERR \
        --alpha $AP_ALN_ALPHA --nproc $NPROC --log AssembleAlignLog.log \
        --outname "${OUTNAME}-ALN" --scanrev >> $RUNLOG
else
    $RUN AssemblePairs.py align -1 "${OUTNAME}-R2_consensus-pass.fastq" \
        -2 "${OUTNAME}-R1_consensus-pass.fastq" --1f CONSCOUNT --2f $PRFIELD CONSCOUNT \
        --coord presto --rc tail --minlen $AP_ALN_MINLEN --maxerror $AP_ALN_MAXERR \
        --alpha $AP_ALN_ALPHA --nproc $NPROC --log AssembleAlignLog.log \
        --outname "${OUTNAME}-ALN" >> $RUNLOG
fi

# Assemble paired ends via alignment against V-region reference database
if $REFERENCE_ASSEMBLY; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs reference"
    $RUN AssemblePairs.py reference -1 "${OUTNAME}-ALN-1_assemble-fail.fastq" \
        -2 "${OUTNAME}-ALN-2_assemble-fail.fastq" -r $REF_FILE \
        --1f CONSCOUNT --2f $PRFIELD CONSCOUNT --coord presto \
        --maxerror $AP_REF_MAXERR --evalue $AP_REF_EVALUE --maxhits $AP_REF_MAXHITS \
        --nproc $NPROC --log AssembleReferenceLog.log --outname "${OUTNAME}-REF" \
        --exec $USEARCH_EXEC >> $RUNLOG
    cat "${OUTNAME}-ALN_assemble-pass.fastq" "${OUTNAME}-REF_assemble-pass.fastq" > \
        "${OUTNAME}-CAT_assemble-pass.fastq"
    PH_FILE="${OUTNAME}-CAT_assemble-pass.fastq"
else
    PH_FILE="${OUTNAME}-ALN_assemble-pass.fastq"
fi

# Rewrite header with minimum of CONSCOUNT
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders collapse"
$RUN ParseHeaders.py collapse -s $PH_FILE -f CONSCOUNT --act min \
    --outname "${OUTNAME}-FIN" > /dev/null

# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
if $CS_KEEP; then
    $RUN CollapseSeq.py -s "${OUTNAME}-FIN_reheader.fastq" -n $CS_MISS --uf PRCONS \
        --cf CONSCOUNT --act sum --outname "${OUTNAME}-FIN" --inner --keepmiss >> $RUNLOG
else
    $RUN CollapseSeq.py -s "${OUTNAME}-FIN_reheader.fastq" -n $CS_MISS --uf PRCONS \
        --cf CONSCOUNT --act sum --outname "${OUTNAME}-FIN" --inner >> $RUNLOG
fi

# Filter to sequences with at least 2 supporting sources
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
$RUN SplitSeq.py group -s "${OUTNAME}-FIN_collapse-unique.fastq" -f CONSCOUNT --num 2 >> $RUNLOG

# Mask low quality positions
if $MASK_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq maskqual"
    $RUN FilterSeq.py maskqual -s "${OUTNAME}-FIN_collapse-unique_atleast-2.fastq" -q $FS_MASK \
        --nproc $NPROC --outname "${OUTNAME}-FIN" --clean >> $RUNLOG
    FSMISS_FILE="${OUTNAME}-FIN_maskqual-pass.fastq"
else
    FSMISS_FILE="${OUTNAME}-FIN_collapse-unique_atleast-2.fastq"
fi

# Remove sequences with many Ns
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq missing"
$RUN FilterSeq.py missing -s $FSMISS_FILE -n $FS_MISS --inner --nproc $NPROC \
    --log MissingLog.log --outname "${OUTNAME}-FIN" --clean >> $RUNLOG

# Create table of final repertoire
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
mv "${OUTNAME}-FIN_missing-pass.fastq" "${OUTNAME}_high-fidelity.fastq"
$RUN ParseHeaders.py table -s "${OUTNAME}_high-fidelity.fastq" \
    -f ID PRIMER PRCONS CONSCOUNT DUPCOUNT >> $RUNLOG

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
if $FILTER_LOWQUAL; then
    $RUN ParseLog.py -l QualityLogR[1-2].log -f ID QUALITY > /dev/null &
fi
$RUN ParseLog.py -l PrimerLogR[1-2].log -f ID BARCODE PRIMER ERROR > /dev/null &
$RUN ParseLog.py -l ConsensusLogR[1-2].log -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ DIVERSITY > /dev/null &
$RUN ParseLog.py -l AssembleAlignLog.log -f ID OVERLAP LENGTH ERROR PVALUE HEADPOS TAILPOS HEADFIELDS TAILFIELDS > /dev/null &
if $REFERENCE_ASSEMBLY; then
    $RUN ParseLog.py -l AssembleReferenceLog.log -f ID REFID OVERLAP LENGTH GAP ERROR HEADEVALUE TAILEVALUE HEADPOS TAILPOS HEADFIELDS TAILFIELDS > /dev/null &
fi
$RUN ParseLog.py -l MissingLog.log -f ID MISSING > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    tar -cf LogFiles.tar *LogR[1-2].log *Log.log
    rm *LogR[1-2].log *Log.log
    gzip LogFiles.tar

    if $CS_KEEP; then
        tar -cf TempFiles.tar *R[1-2]_*.fastq *assemble-pass* *under* *duplicate* *reheader* *fail*
        rm *R[1-2]_*.fastq *assemble-pass* *under* *duplicate* *reheader* *fail*
    else
        tar -cf TempFiles.tar *R[1-2]_*.fastq *assemble-pass* *under* *duplicate* *undetermined* *reheader* *fail*
        rm *R[1-2]_*.fastq *assemble-pass* *under* *duplicate* *undetermined* *reheader* *fail*
    fi
    gzip TempFiles.tar

fi

# End
printf "DONE\n"
cd ../

