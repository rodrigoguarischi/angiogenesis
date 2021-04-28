#!/bin/sh

### transcripts_discovery.sh ######################################################################
# Starting from same alignment files (BAM) used on expression analysis perform transcriptome
# reconstruction to identify genes outside ENSEMBL annotation
# 

### HISTORY #######################################################################################
# Version Date
# 0.1 22-Feb-2016
#

### MAIN ##########################################################################################

# Set base directory
baseDir=/work/angiogenesis/rodrigoguarischi;

# Create directories and link files
transcriptomeDir=$baseDir/mRNA/transcripts_discovert;
mkdir -p $transcriptomeDir;
cd $transcriptomeDir;

# Merge BAM files by samples using Picard-tools
bamDir=$baseDir/alignment/results;
inputDir=$baseDir/mRNA/transcripts_discovert/alignments_merged_by_sample;
mkdir -p $inputDir;
cd $inputDir;

samples=$(for file in $(ls $bamDir/*.bam); do sample_id=$(basename $file .rg_added_sorted.bam | perl -pe "s/-\d+.*$//"); echo -e "$sample_id"; done | uniq ); for sample_i in $samples; do bam_files=$(for bam_i in $(ls $bamDir/$sample_i*); do echo -n "I=$bam_i "; done); nice picard-tools MergeSamFiles $bam_files O=$inputDir/$sample_i.merged.bam & done 2> merge.log;

wait;

# Run StringTie for each sample
stringTieExe=$baseDir/tools/stringtie-1.2.1.Linux_x86_64/stringtie;
genomeAnnotation=$baseDir/mm10/annotation/ensembl_genes.gtf;

reconstructedDir=$baseDir/mRNA/transcripts_discovert/results;
mkdir -p $reconstructedDir;
cd $reconstructedDir;

for file in $(ls $inputDir/*.bam); do sample_id=$(basename $file .merged.bam); $stringTieExe $file -G $genomeAnnotation -o $sample_id.gtf -p 5 & done;

wait;

# Merge all reconstruction files into one
gtf_files=$(ls *.gtf); $stringTieExe --merge $gtf_files -G $genomeAnnotation -o merged.gtf;

# bedtools intersect -a merged.gtf -b ../../../mm10/annotation/ensembl_genes.gtf -v | awk 'BEGIN{FS="\t"; OFS="\t"} $3=="exon" {print $0, $5-$4, "chr"$1":"$4"-"$5}' | sort -k 10nr -t$'\t' | less -S

# Evaluating Transcript Discovery
bedtools intersect -a merged.gtf -b ../../../mm10/annotation/ensembl_genes.gtf -v | awk 'BEGIN{FS="\t"; OFS="\t"} $3=="exon" {print $0, $5-$4+1, "chr"$1":"$4"-"$5}' | sort -k 10nr -t$'\t' > intergenic_exons.txt
bedtools intersect -a merged.gtf -b ../../../mm10/annotation/ensembl_genes.gtf -v -s | awk 'BEGIN{FS="\t"; OFS="\t"} $3=="exon" {print $0, $5-$4+1, "chr"$1":"$4"-"$5}' | sort -k 10nr -t$'\t' > intergenic_exons_stranded.txt

