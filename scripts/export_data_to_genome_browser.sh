#!/bin/sh

### export_data_to_genome_browser.sh ##############################################################
# Prepare data (alignments and de novo transcripts) to be imported on UCSC genome Browser
# 

### HISTORY #######################################################################################
# Version Date
# 0.1 24-Feb-2016
#

### MAIN ##########################################################################################

# Set base directory
baseDir=/work/angiogenesis/rodrigoguarischi;
gbDir=$baseDir/genome_browser;
toolsDir=$baseDir/tools/ucscTools;

mkdir -p $gbDir;
cd $gbDir;

# Pipe to create bed from StringTie output. Steps: 
# 	1) Convert UCSC chromosome IDS to ENSEMBL IDS
#	2) Conver GTF file into genePred
# 	3) Convert genePred into Bed
#	4) Sort file by chr name and position (bigBed requires sorted entry)
#	5) Change colours of new/old transcripts. Light green for new transcripts and dark green for known transcripts.
ens2ucsc -o mm10 -f gtf ../mRNA/transcripts_dicovery/results/merged.gtf | $toolsDir/gtfToGenePred /dev/stdin /dev/stdout | $toolsDir/genePredToBed /dev/stdin /dev/stdout | sort -k1,1 -k2,2n | awk 'BEGIN{FS="\t";OFS="\t"} { if ($4 ~ "ENSMUST.*") {print $1,$2,$3,$4,$5,$6,$7,$8,"0,68,27",$10,$11,$12} else {print $1,$2,$3,$4,$5,$6,$7,$8,"102,194,164",$10,$11,$12} }' > angio_transcripts.bed 

# Pipe to incorporate sequence on bed files
ucsc2ens -o mm10 -f bed angio_transcripts.bed | bedtools getfasta -s -name -split -fi ../mm10/genome_database/ensembl_genome.fa -bed /dev/stdin -fo /dev/stdout | grep -v ">" | perl -pe "s/([ACGT]{60})/\1<\/br>/gi" > angio_transcripts.fasta; 
paste angio_transcripts.bed angio_transcripts.fasta > angio_transcripts_seq.bed

# Convert bed to bigBed
$toolsDir/fetchChromSizes mm10 > mm10.chrom.sizes
$toolsDir/bedToBigBed -tab -type=bed12+1 -as=bed12wSeq.as angio_transcripts_seq.bed -extraIndex=name mm10.chrom.sizes angio_transcripts_seq.bb

# Create track with new transcripts only
grep -P "MSTRG.\d+" angio_transcripts_seq.bed > angio_transcripts_new_transcripts.bed
$toolsDir/bedToBigBed -tab -type=bed12+1 -as=bed12wSeq.as angio_transcripts_new_transcripts.bed -extraIndex=name mm10.chrom.sizes angio_transcripts_seq_new_transcripts.bb

# Pipe to export BAM files to UCSC
# 1) Filter out low quality mapping
# 2) Convert Ens IDs to UCSC
# 3) Compact output back to BAM format
mergedBySampleDir=$baseDir/mRNA/transcripts_dicovery/alignments_merged_by_sample;
for file in $(ls $mergedBySampleDir/*.bam); do file_bn=$(basename $file .bam); samtools view -h -q 255 $file | ens2ucsc -o mm10 -f sam | samtools view -S -b /dev/stdin > $file_bn.hq.bam & done;

wait;

# Create index for each file
for file in $(ls *.hq.bam); do samtools index $file & done;

wait;

