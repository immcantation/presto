Illumina MiSeq 2x250bp B cell receptor mRNA without UMIs
================================================================================

Stranded Data
--------------------------------------------------------------------------------

.. todo::

    Example workflow for MiSeq stranded data without UMIs


.. code-block:: bash
    :linenos:

    #!/bin/bash
    # Super script to run the pRESTO 0.4.7 pipeline on MiSeq data without UMIs
    #
    # Author:  Jason Anthony Vander Heiden
    # Date:    2015.05.19
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
    ZIP_FILES=true
    FILTER_LOWQUAL=true
    REFERENCE_ASSEMBLY=true
    MASK_LOWQUAL=false

    # FilterSeq run parameters
    FS_QUAL=20
    FS_MASK=30

    # MaskPrimers run parameters
    MP_R1_MODE="cut"
    MP_R2_MODE="mask"
    MP_R1_START=17
    MP_R2_START=0
    MP_R1_MAXERR=0.2
    MP_R2_MAXERR=0.2

    # AssemblePairs-align run parameters
    AP_ALN_SCANREV=true
    AP_ALN_MAXERR=0.3
    AP_ALN_MINLEN=8
    AP_ALN_ALPHA=1e-5

    # AssemblePairs-reference run parameters
    AP_REF_MINIDENT=0.5
    AP_REF_EVALUE=1e-5
    AP_REF_MAXHITS=100
    REF_FILE="/scratch2/kleinstein/germlines/IMGT_Human_IGV_ungapped_2014-08-23.fasta"
    #REF_FILE="/scratch2/kleinstein/germlines/IMGT_Mouse_IGV_ungapped_2014-11-22.fasta"
    USEARCH_EXEC=$HOME/bin/usearch

    # CollapseSeq run parameters
    CS_KEEP=true
    CS_MISS=0

    # Define log files
    PIPELINE_LOG="Pipeline.log"
    ERROR_LOG="Pipeline.err"

    # Make output directory and empty log files
    mkdir -p $OUTDIR; cd $OUTDIR
    echo '' > $PIPELINE_LOG
    echo '' > $ERROR_LOG

    # Start
    echo "DIRECTORY: ${OUTDIR}"
    echo "VERSIONS:"
    echo "  $(AlignSets.py --version 2>&1)"
    echo "  $(AssemblePairs.py --version 2>&1)"
    echo "  $(BuildConsensus.py --version 2>&1)"
    echo "  $(ClusterSets.py --version 2>&1)"
    echo "  $(CollapseSeq.py --version 2>&1)"
    echo "  $(ConvertHeaders.py --version 2>&1)"
    echo "  $(FilterSeq.py --version 2>&1)"
    echo "  $(MaskPrimers.py --version 2>&1)"
    echo "  $(PairSeq.py --version 2>&1)"
    echo "  $(ParseHeaders.py --version 2>&1)"
    echo "  $(ParseLog.py --version 2>&1)"
    echo "  $(SplitSeq.py --version 2>&1)"
    echo -e "\nSTART"
    STEP=0

    # Remove low quality reads
    if $FILTER_LOWQUAL; then
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        #OUTPREFIX="$(printf '%02d' $STEP)--${OUTNAME}"
        FilterSeq.py quality -s $R1_FILE -q $FS_QUAL --nproc $NPROC \
            --outname "${OUTNAME}-R1" --outdir . --log QualityLogR1.log \
            >> $PIPELINE_LOG  2> $ERROR_LOG
        FilterSeq.py quality -s $R2_FILE -q $FS_QUAL --nproc $NPROC \
            --outname "${OUTNAME}-R2" --outdir . --log QualityLogR2.log  \
            >> $PIPELINE_LOG  2> $ERROR_LOG
        MPR1_FILE="${OUTNAME}-R1_quality-pass.fastq"
        MPR2_FILE="${OUTNAME}-R2_quality-pass.fastq"
    else
        MPR1_FILE=$R1_FILE
        MPR2_FILE=$R2_FILE
    fi

    # Identify primers and UID
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
    MaskPrimers.py score -s $MPR1_FILE -p $R1_PRIMERS --mode $MP_R1_MODE \
        --start $MP_R1_START --maxerror $MP_R1_MAXERR --nproc $NPROC --log PrimerLogR1.log \
        --outname "${OUTNAME}-R1" --outdir . >> $PIPELINE_LOG 2> $ERROR_LOG
    MaskPrimers.py score -s $MPR2_FILE -p $R2_PRIMERS --mode $MP_R2_MODE \
        --start $MP_R2_START --maxerror $MP_R2_MAXERR --nproc $NPROC --log PrimerLogR2.log \
        --outname "${OUTNAME}-R2" --outdir . >> $PIPELINE_LOG 2> $ERROR_LOG

    # Assign UIDs to read 1 sequences
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
    PairSeq.py -1 "${OUTNAME}-R2_primers-pass.fastq" -2 "${OUTNAME}-R1_primers-pass.fastq" \
        --coord illumina >> $PIPELINE_LOG 2> $ERROR_LOG

    # Assemble paired ends via mate-pair alignment
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"

    if $AP_ALN_SCANREV; then
        AssemblePairs.py align -1 "${OUTNAME}-R2_primers-pass_pair-pass.fastq" \
            -2 "${OUTNAME}-R1_primers-pass_pair-pass.fastq" --2f PRIMER \
            --coord illumina --rc tail --minlen $AP_ALN_MINLEN --maxerror $AP_ALN_MAXERR \
            --alpha $AP_ALN_ALPHA --nproc $NPROC --log AssembleAlignLog.log \
            --outname "${OUTNAME}-ALN" --scanrev --failed >> $PIPELINE_LOG 2> $ERROR_LOG
    else
        AssemblePairs.py align -1 "${OUTNAME}-R2_primers-pass_pair-pass.fastq" \
            -2 "${OUTNAME}-R1_primers-pass_pair-pass.fastq" --2f PRIMER \
            --coord illumina --rc tail --minlen $AP_ALN_MINLEN --maxerror $AP_ALN_MAXERR \
            --alpha $AP_ALN_ALPHA --nproc $NPROC --log AssembleAlignLog.log \
            --outname "${OUTNAME}-ALN" --failed >> $PIPELINE_LOG 2> $ERROR_LOG
    fi

    # Assemble paired ends via alignment against V-region reference database
    if $REFERENCE_ASSEMBLY; then
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs reference"
        AssemblePairs.py reference -1 "${OUTNAME}-ALN-1_assemble-fail.fastq" \
            -2 "${OUTNAME}-ALN-2_assemble-fail.fastq" -r $REF_FILE --2f PRIMER --coord illumina \
            --minident $AP_REF_MINIDENT --evalue $AP_REF_EVALUE --maxhits $AP_REF_MAXHITS \
            --nproc $NPROC --log AssembleReferenceLog.log --outname "${OUTNAME}-REF" \
            --exec $USEARCH_EXEC --failed >> $PIPELINE_LOG 2> $ERROR_LOG
        cat "${OUTNAME}-ALN_assemble-pass.fastq" "${OUTNAME}-REF_assemble-pass.fastq" > \
            "${OUTNAME}-CAT_assemble-pass.fastq"
        FS_FILE="${OUTNAME}-CAT_assemble-pass.fastq"
    else
        FS_FILE="${OUTNAME}-ALN_assemble-pass.fastq"
    fi

    # Mask low quality positions
    if $MASK_LOWQUAL; then
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq maskqual"
        FilterSeq.py maskqual -s $FS_FILE -q $FS_MASK --nproc $NPROC \
            --outname "${OUTNAME}-FIN" --log MaskqualLog.log >> $PIPELINE_LOG 2> $ERROR_LOG
        CS_FILE="${OUTNAME}-FIN_maskqual-pass.fastq"
    else
        CS_FILE=$FS_FILE
    fi

    # Remove duplicate sequences
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
    if $CS_KEEP; then
        CollapseSeq.py -s $CS_FILE -n $CS_MISS --uf PRIMER --inner --keepmiss \
        --outname "${OUTNAME}-FIN" >> $PIPELINE_LOG 2> $ERROR_LOG
    else
        CollapseSeq.py -s $CS_FILE -n $CS_MISS --uf PRIMER --inner \
        --outname "${OUTNAME}-FIN" >> $PIPELINE_LOG 2> $ERROR_LOG
    fi

    # Create table of final repertoire
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
    ParseHeaders.py table -s "${OUTNAME}-FIN_collapse-unique.fastq" \
        -f ID PRIMER DUPCOUNT --outname "Unique" >> $PIPELINE_LOG 2> $ERROR_LOG

    # Process log files
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
    if $FILTER_LOWQUAL; then
        ParseLog.py -l QualityLogR[1-2].log -f ID QUALITY > /dev/null &
    fi
    ParseLog.py -l PrimerLogR[1-2].log -f ID BARCODE PRIMER ERROR \
        > /dev/null  2> $ERROR_LOG &
    ParseLog.py -l AssembleAlignLog.log -f ID LENGTH OVERLAP ERROR PVALUE FIELDS1 FIELDS2 \
        > /dev/null  2> $ERROR_LOG &
    if $REFERENCE_ASSEMBLY; then
        ParseLog.py -l AssembleReferenceLog.log -f ID REFID LENGTH OVERLAP GAP EVALUE1 EVALUE2 IDENTITY FIELDS1 FIELDS2 \
        > /dev/null  2> $ERROR_LOG &
    fi
    if $MASK_LOWQUAL; then
        ParseLog.py -l MaskqualLog.log -f ID MASKED > /dev/null  2> $ERROR_LOG &
    fi
    wait

    # Zip intermediate and log files
    if $ZIP_FILES; then
        LOG_FILES_ZIP=$(ls *LogR[1-2].log *Log.log)
        tar -zcf LogFiles.tar $LOG_FILES_ZIP
        rm $LOG_FILES_ZIP

        TEMP_FILES_ZIP=$(ls *.fastq | grep -v "collapse-unique.fastq")
        tar -zcf TempFiles.tar $TEMP_FILES_ZIP
        rm $TEMP_FILES_ZIP
    fi

    # End
    printf "DONE\n\n"
    cd ../


Unstranded Data
--------------------------------------------------------------------------------

.. todo::

    Example workflow for MiSeq unstranded data without UMIs


.. code-block:: bash

    #!/bin/bash
    # Super script to run the pRESTO 0.4.5 pipeline on unstranded MiSeq data without UIDs
    #
    # Author:  Jason Anthony Vander Heiden
    # Date:    2014.11.26
    #
    # Required Arguments:
    #   $1 = read 1 file
    #   $2 = read 2 file
    #   $3 = read 1 primer file (V-region)
    #   $4 = read 2 primer file (C-region)
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
    QUAL_STEP=true
    MASK_STEP=false
    MISS_STEP=false

    # Define pRESTO run parameters
    FS_QUAL=20
    FS_MASK=20
    FS_MISS=20
    AP_SCANREV=true
    AP_MAXERR=0.3
    AP_ALPHA=0.01
    MP1_MAXLEN=100
    MP2_MAXLEN=100
    MP1_MAXERR=0.2
    MP2_MAXERR=0.2
    CS_MISS=20
    MUSCLE_EXEC=$HOME/bin/muscle


    # Define script execution command and log files
    mkdir -p $OUTDIR; cd $OUTDIR
    RUNLOG="Pipeline.log"
    echo '' > $RUNLOG
    if $LOG_RUNTIMES; then
        TIMELOG="Time.log"
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

    # Filter low quality reads
    echo -e "\nSTART"
    STEP=0

    # Remove low quality reads
    if $QUAL_STEP; then
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        #OUTPREFIX="$(printf '%02d' $STEP)--${OUTNAME}"
        $RUN FilterSeq.py quality -s $R1_FILE -q $FS_QUAL --nproc $NPROC \
            --outname "${OUTNAME}-R1" --outdir . >> $RUNLOG
        $RUN FilterSeq.py quality -s $R2_FILE -q $FS_QUAL --nproc $NPROC \
            --outname "${OUTNAME}-R2" --outdir . >> $RUNLOG
        APR1_FILE="${OUTNAME}-R1_quality-pass.fastq"
        APR2_FILE="${OUTNAME}-R2_quality-pass.fastq"
    else
        APR1_FILE=$R1_FILE
        APR2_FILE=$R2_FILE
    fi

    # Assemble paired ends
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
    if $AP_SCANREV; then
        $RUN AssemblePairs.py align -1 $APR2_FILE -2 $APR1_FILE --coord illumina --rc tail \
            --maxerror $AP_MAXERR --alpha $AP_ALPHA --nproc $NPROC \
            --log AssembleLog.log --outname "${OUTNAME}" --outdir . --scanrev >> $RUNLOG
    else
        $RUN AssemblePairs.py align -1 $APR2_FILE -2 $APR1_FILE --coord illumina --rc tail \
            --maxerror $AP_MAXERR --alpha $AP_ALPHA --nproc $NPROC \
            --log AssembleLog.log --outname "${OUTNAME}" --outdir . >> $RUNLOG
    fi

    # Identify and mask V-region primers
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align (PR1)"
    $RUN MaskPrimers.py align -s "${OUTNAME}_assemble-pass.fastq" -p $R1_PRIMERS --mode mask \
        --maxerror $MP1_MAXERR --maxlen $MP1_MAXLEN --nproc $NPROC \
        --log PrimerPR1Log.log --outname "${OUTNAME}-PR1" >> $RUNLOG

    # Rename V-region primer field
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders rename (PR1)"
    $RUN ParseHeaders.py rename -s "${OUTNAME}-PR1_primers-pass.fastq" \
        -f PRIMER -k VPRIMER --outname "${OUTNAME}-PR1" > /dev/null

    # Identify and mask C-region primers
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align (PR2)"
    $RUN MaskPrimers.py align -s "${OUTNAME}-PR1_reheader.fastq" -p $R2_PRIMERS --mode mask \
        --maxerror $MP2_MAXERR --maxlen $MP2_MAXLEN --nproc $NPROC --revpr --skiprc \
        --log PrimerPR2Log.log --outname "${OUTNAME}-PR2" >> $RUNLOG

    # Rename C-region primer field
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders rename (PR2)"
    $RUN ParseHeaders.py rename -s "${OUTNAME}-PR2_primers-pass.fastq" \
        -f PRIMER -k CPRIMER --outname "${OUTNAME}-PR2" > /dev/null

    # Remove duplicate sequences
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
    $RUN CollapseSeq.py -s "${OUTNAME}-PR2_reheader.fastq" -n $CS_MISS --uf CPRIMER \
        --cf VPRIMER --act set --outname "${OUTNAME}" --inner >> $RUNLOG

    # Mask low quality positions
    if $MASK_STEP; then
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq maskqual"
        $RUN FilterSeq.py maskqual -s "${OUTNAME}-collapse-unique.fastq" -q $FS_MASK \
            --nproc $NPROC --outname "${OUTNAME}-unique" >> $RUNLOG
        FSMISS_FILE="${OUTNAME}-unique_maskqual-pass.fastq"
    else
        FSMISS_FILE="${OUTNAME}_collapse-unique.fastq"
    fi

    # Remove sequences with many Ns
    if $MISS_STEP; then
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq missing"
        $RUN FilterSeq.py missing -s $FSMISS_FILE -n $FS_MISS --inner --nproc $NPROC \
            --log MissingLog.log --outname "${OUTNAME}" --fasta >> $RUNLOG
        FINAL_FILE="${OUTNAME}_missing-pass.fastq"
    else
        FINAL_FILE=$FSMISS_FILE
    fi

    # Rename final file
    mv $FINAL_FILE "${OUTNAME}_final.fastq"

    # Create table of final repertoire
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
    $RUN ParseHeaders.py table -s "${OUTNAME}_final.fastq" \
        -f ID CPRIMER VPRIMER DUPCOUNT >> $RUNLOG

    # Split final file into sets of singletons and sequences with at least 2 reads
    #printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
    #$RUN SplitSeq.py group -s "${OUTNAME}_final.fastq" -f DUPCOUNT --num 2 >> $RUNLOG

    # Process log files
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
    $RUN ParseLog.py -l AssembleLog.log -f ID REFID OVERLAP LENGTH GAP PVALUE ERROR HEADFIELDS TAILFIELDS > /dev/null &
    $RUN ParseLog.py -l PrimerPR[1-2]Log.log -f ID SEQORIENT PRORIENT PRSTART PRIMER ERROR > /dev/null &
    if $MISS_STEP; then
        $RUN ParseLog.py -l MissingLog.log -f ID MISSING > /dev/null &
    fi
    wait

    if $ZIP_FILES; then
        tar -cf LogFiles.tar *Log.log
        gzip LogFiles.tar
        rm *Log.log

        tar -cf TempFiles.tar *quality* *duplicate* *undetermined* *reheader* *fail*
        gzip TempFiles.tar
        rm *quality* *duplicate* *undetermined* *reheader* *fail*
    fi

    # End
    printf "DONE\n"
    cd ../


