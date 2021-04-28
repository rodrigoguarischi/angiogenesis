#!/bin/sh

### summarize_by_genes.sh #########################################################################
# Performs counting on alignment files (BAM) using HTSeq-count.
# 

### HISTORY #######################################################################################
# Version Date
# 0.1 17-Feb-2016
#

### MAIN ##########################################################################################

# Set base directory
baseDir=/work/angiogenesis/rodrigoguarischi;

alignmentsDir=$baseDir/alignment/results;
genomeAnnotation=$baseDir/mm10/annotation/ensembl_genes.gtf;
countsDir=$baseDir/mRNA/counts_by_gene;

mkdir -p $countsDir;
cd $countsDir;

logFile=$baseDir/mRNA/count_tables/counting.log;

echo -e "## Summarizing counts by gene ##\n" > $logFile;

for file in $(ls $alignmentsDir/*.rg_added_sorted.bam); do file_bn=$(basename $file .rg_added_sorted.bam); run_id=$(echo $file_bn | perl -pe "s/_[ACGT]{6}$//"); nice htseq-count --quiet --format bam --stranded reverse $file $genomeAnnotation > $countsDir/$run_id.count & done 2>> $logFile;

wait;

echo -e "\n\n" >> $logFile;

