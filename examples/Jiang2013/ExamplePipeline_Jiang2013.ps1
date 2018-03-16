# Set path to pRESTO scripts
$anacondapath = 'C:\Anaconda3\Scripts'
 
# Define run parameters and output files
$NAME = 'SRR765688'
$READ_FILE = (Resolve-Path SRR765688.fastq).Path
$MID_PRIMERS = (Resolve-Path SRR765688_MIDs.fasta).Path
$FWD_PRIMERS= (Resolve-Path SRX190717_VPrimers.fasta).Path
$REV_PRIMERS= (Resolve-Path SRX190717_CPrimers.fasta).Path

# Make output directory and empty log file
$OUTDIR = (mkdir -Force output).FullName
echo '' > $OUTDIR\Pipeline.log
$PIPELINE_LOG = (Resolve-Path $OUTDIR\Pipeline.log).Path

# Start
cd $OUTDIR
echo ('DIRECTORY: ' + (pwd).Path)
  
# Remove short reads
Write-Output "STEP1 Remove short reads"
python $anacondapath\FilterSeq.py length -s $READ_FILE -n 300 --outname FLTR --outdir $OUTDIR >> $PIPELINE_LOG
# Remove low quality reads
Write-Output "STEP2 Remove low quality reads"
python $anacondapath\FilterSeq.py quality -s FLTR_length-pass.fastq  -q 20 --outname FLTR >> $PIPELINE_LOG
# Identify and remove MIDs
Write-Output "STEP3 Identify and remove MIDs"
python $anacondapath\MaskPrimers.py score -s FLTR_quality-pass.fastq -p $MID_PRIMERS --mode cut --start 0 --maxerror 0.1 --outname MID  >> $PIPELINE_LOG
# Identify and mask forward (V-region) primers
Write-Output "STEP4 Identify and mask forward (V-region) primers"
python $anacondapath\MaskPrimers.py align -s MID_primers-pass.fastq -p $FWD_PRIMERS --mode mask --maxlen 50 --maxerror 0.3 --outname PRMR  >> $PIPELINE_LOG
# Identify and remove reverse (C-region) primers
Write-Output "STEP5 Identify and remove reverse (C-region) primers"
python $anacondapath\MaskPrimers.py align -s PRMR_primers-pass.fastq -p $REV_PRIMERS --mode cut --maxlen 50 --maxerror 0.3 --revpr --skiprc  --outname FIN >> $PIPELINE_LOG
# Remove duplicate sequences
Write-Output "STEP8 Remove duplicate sequences"
python $anacondapath\CollapseSeq.py -s FIN_primers-pass.fastq -n 20 --uf MID CPRIMER --inner --outname FIN  >> $PIPELINE_LOG
# Filter to sequences with at least 2 supporting reads
Write-Output "STEP9 Filter to sequences with at least 2 supporting reads"
python $anacondapath\SplitSeq.py group -s FIN_collapse-unique.fastq -f DUPCOUNT --num 2  >> $PIPELINE_LOG

# Create tables of final repertoire files
Write-Output "STEP10 Create tables of final repertoire files"
python $anacondapath\ParseHeaders.py table -s FIN_reheader.fastq -f ID MID CPRIMER VPRIMER --outname Final
python $anacondapath\ParseHeaders.py table -s FIN_collapse-unique.fastq -f ID MID CPRIMER VPRIMER DUPCOUNT --outname Final-Unique
python $anacondapath\ParseHeaders.py table -s FIN_collapse-unique_atleast-2.fastq -f ID MID CPRIMER VPRIMER DUPCOUNT --outname Final-Unique-Atleast2
