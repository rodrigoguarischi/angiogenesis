#!/bin/sh

### align.sh ######################################################################################
# Performs alignment of raw RNA-seq reads (fastq) using STAR 2-pass with most recent ENSEMBL 
# assembly and respective annotation. The aim is to use the same set of alignment files for 
# expression and editing analysis.
# 

### HISTORY #######################################################################################
# Version Date
# 0.1 15-Feb-2016
#

### MAIN ##########################################################################################

# Set base directory
baseDir=/work/angiogenesis/rodrigoguarischi;

scriptDir=$baseDir/git/scripts;
alignmentDir=$baseDir/alignment;
fastqDir=$baseDir/alignment/raw;

# Make links to raw FASTQ files (Run script create_links.pl on $fastqDir)
mkdir -p $fastqDir;
cd $fastqDir;
$scriptDir/create_links.pl;

# Build genome index file. Expect one file named ensembl_genome.fa on it
starExe=$baseDir/tools/STAR-2.5.1b/bin/Linux_x86_64_static/STAR;
genomeDir=$baseDir/mm10/genome_database;
cd $genomeDir;
$starExe --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeDir/ensembl_genome.fa --runThreadN 20;

# Run 1 pass alignment
genomeAnnotation=$baseDir/mm10/annotation/ensembl_genes.gtf;
runDir1=$baseDir/alignment/alg_1pass;
mkdir -p $runDir1;
cd $runDir1;
for file in $(ls $fastqDir/*.gz); do file_bn=$(basename $file .fastq.gz); nice $starExe --genomeDir $genomeDir --sjdbGTFfile $genomeAnnotation --readFilesCommand zcat --readFilesIn $file --runThreadN 40 --outFileNamePrefix ./$file_bn. --outSAMtype BAM SortedByCoordinate; done;

# Run 2 pass alignment
runDir2=$baseDir/alignment/alg_2pass;
mkdir -p $runDir2;
cd $runDir2;
SJ_1pass=$(ls $runDir1/*.tab);
for file in $(ls $fastqDir/*.gz); do file_bn=$(basename $file .fastq.gz); nice $starExe --genomeDir $genomeDir --sjdbGTFfile $genomeAnnotation --sjdbFileChrStartEnd $SJ_1pass --readFilesCommand zcat --readFilesIn $file --runThreadN 40 --outFileNamePrefix ./$file_bn. --outSAMtype BAM SortedByCoordinate; done;

# Add read groups and sort (picard-tools)
resultsDir=$baseDir/alignment/results;
logFile=$baseDir/alignment/results/Add_read_groups_and_sort.log;
mkdir -p $resultsDir;
cd $resultsDir;

echo -e "## Adding Read Groups ##\n" > $logFile;

for file in $(ls $runDir2/*.bam); do file_bn=$(basename $file .Aligned.sortedByCoord.out.bam); run_id=$(echo $file_bn | perl -pe "s/_[ACGT]{6}$//"); sample=$(echo $file_bn | cut -d'-' -f 1); barcode=$(echo $file_bn | perl -pe "s/.*([ACGT]{6}$)/\1/"); lane_id=$(samtools view $file | head -n 1 | cut -f 1 | perl -pe "s/:\d+:\d+:\d+$//"); nice picard-tools AddOrReplaceReadGroups I=$file O=$resultsDir/$file_bn.rg_added_sorted.bam SO=coordinate RGID=$run_id RGLB=$sample RGPL='illumina' RGPU=$barcode.$lane_id RGSM=$sample & done 2>> $logFile;

wait;

echo -e "\n\n" >> $logFile;

