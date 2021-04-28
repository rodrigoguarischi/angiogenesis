#!/bin/sh

### summarize_by_exon.sh ######################################################################
# This script is intended to summarize alignment files (BAM) by exons for further use by
# DEXSeq package. First step is to "collapse" reference annotation for later use on counting
# 

### HISTORY #######################################################################################
# Version Date
# 0.1 26-Feb-2016
#

### MAIN ##########################################################################################

# Set base directory
baseDir=/work/angiogenesis/rodrigoguarischi;

# Collapse reference annotation
dexseqScriptsDir=$baseDir/tools/DEXSeq/inst/python_scripts;
annotationDir=$baseDir/mm10/annotation;

cd $annotationDir;

python $dexseqScriptsDir/dexseq_prepare_annotation.py $annotationDir/ensembl_genes.gtf $annotationDir/ensembl_genes.dexseq.gff

# Count number of reads by exon on collapsed reference generated
dexseqCountsDir=$baseDir/mRNA/counts_by_exon;
bamDir=$baseDir/alignment/results;
mkdir -p $dexseqCountsDir;

cd $dexseqCountsDir;

for file in $(ls $bamDir/*.bam); do sample_id=$(basename $file .rg_added_sorted.bam | perl -pe "s/_[ACGT]{6}//i" ); samtools view -h $file | python $dexseqScriptsDir/dexseq_count.py -s reverse $annotationDir/ensembl_genes.dexseq.gff /dev/stdin $dexseqCountsDir/$sample_id.txt 2>> $dexseqCountsDir/counting.log & done;


