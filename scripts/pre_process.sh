# Following GATK Best-pratices, pre-processing step consists of:
# 1) Mark Duplicates (Picard tool)
# 2) Split'N'Trim (GATK)
# 3) Recalibrate Bases (GATK)

# >>>>>> REMOVER <<<<<
# # Prepare reference genome - create dictionary and index
# cd $genomeDir
# picard-tools CreateSequenceDictionary R=$genomeDir/ensembl_genome.fa O=$genomeDir/ensembl_genome.dict
# samtools faidx $genomeDir/ensembl_genome.fa
# >>>>>> REMOVER <<<<<<

# Mark duplicates, and create index
preProcessingDir=$baseDir/preProcessing

echo -e "## Mark Duplicates ##\n" >> $logFile;

for file in $(ls $preProcessingDir/*.rg_added_sorted.bam); do file_bn=$(basename $file .rg_added_sorted.bam); picard-tools MarkDuplicates I=$file O=$preProcessingDir/$file_bn.rg_added_sorted.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$preProcessingDir/$file_bn.metrics & done 2>> $logFile;

wait;

echo -e "\n\n" >> $logFile;


# Split'N'Trim and reassign mapping qualities
cd $preProcessingDir
gatkJar=$baseDir/tools/GATK-3.4.0/GenomeAnalysisTK.jar
reference_file=$genomeDir/ensembl_genome.fa

echo -e "## Split'N'Trim ##\n" >> $logFile;
for file in $(ls $preProcessingDir/*.rg_added_sorted.dedupped.bam); do file_bn=$(basename $file .rg_added_sorted.dedupped.bam); java -jar $gatkJar -T SplitNCigarReads -R $reference_file -I $file -o $preProcessingDir/$file_bn.rg_added_sorted.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS & done 2>> $logFile;

wait;

echo -e "\n\n" >> $logFile;


# Base Recalibration
mm10_vcf=$genomeDir/Mus_musculus.vcf
echo -e "## Base Recalibrator - Modeling ##\n" >> $logFile;

for file in $(ls $preProcessingDir/*.rg_added_sorted.dedupped.split.bam); do file_bn=$(basename $file .rg_added_sorted.dedupped.split.bam); java -jar $gatkJar -T BaseRecalibrator -R $reference_file -I $file -knownSites $mm10_vcf -o $preProcessingDir/$file_bn.recal_data.grp & done 2>> $logFile;

wait;

echo -e "\n\n" >> $logFile;

echo -e "## Base Recalibrator - Write BAM ##\n" >> $logFile;

for file in $(ls $preProcessingDir/*.rg_added_sorted.dedupped.split.bam); do file_bn=$(basename $file .rg_added_sorted.dedupped.split.bam); java -jar $gatkJar -T PrintReads -R $reference_file -I $file -BQSR $preProcessingDir/$file_bn.recal_data.grp -o $preProcessingDir/$file_bn.rg_added_sorted.dedupped.split.recal.bam & done 2>> $logFile;

wait;

echo -e "\n\n" >> $logFile;


# Join files by sample and index BAMs
variantDiscoveryDir=$baseDir/variantDiscovery
logFile=$variantDiscoveryDir/run.log
mkdir $variantDiscoveryDir
cd $variantDiscoveryDir
echo -e "## Variant Discovery - Merge BAMs ##\n" > $logFile;

for sample in 'P12_1' 'P12_2' 'P15_1' 'P15_2' 'P17_1' 'P17_2' 'R12_1' 'R12_2' 'R12_5_1' 'R12_5_2' 'R15_1' 'R15_2' 'R17_1' 'R17_2'; do sample_files=$(for file in $(ls $preProcessingDir/$sample-*.dedupped.split.recal.bam); do echo -n "I=$file "; done); (picard-tools MergeSamFiles $sample_files O=$variantDiscoveryDir/$sample.merged.bam; picard-tools BuildBamIndex I=$file) & done 2>> $logFile;

wait;

echo -e "\n\n" >> $logFile;


# Run variant calling
echo -e "## Variant Discovery - Variant calling ##\n" >> $logFile;

for file in $(ls $variantDiscoveryDir/*.merged.bam); do file_bn=$(basename $file .merged.bam); java -jar $gatkJar -T HaplotypeCaller -R $reference_file -I $file -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --dbsnp $mm10_vcf -o $variantDiscoveryDir/$file_bn.vcf & done 2>> $logFile;


wait;

echo -e "\n\n" >> $logFile;


# Variant Filtering
echo -e "## Variant Discovery - Variant filtering ##\n" >> $logFile;

for file in $(ls $variantDiscoveryDir/*.vcf); do file_bn=$(basename $file .vcf); java -jar $gatkJar -T VariantFiltration -R $reference_file -V $file -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $variantDiscoveryDir/$file_bn.filtered.vcf & done 2>> $logFile;

wait;

echo -e "\n\n" >> $logFile;



# Run variant calling - Joint Genotyping
echo -e "## Variant Discovery - Variant calling (Joint Genotyping) ##\n" >> $logFile;

sample_files=$(for file in $(ls $variantDiscoveryDir/*.merged.bam); do echo -n "-I $file "; done); java -jar $gatkJar -T HaplotypeCaller -R $reference_file $sample_files -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --dbsnp $mm10_vcf -o $variantDiscoveryDir/joint_genotyping.vcf 2>> $logFile;

echo -e "\n\n" >> $logFile;


# Variant Filtering
echo -e "## Variant Discovery - Variant filtering (Joint Genotyping) ##\n" >> $logFile;

java -jar $gatkJar -T VariantFiltration -R $reference_file -V $variantDiscoveryDir/joint_genotyping.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $variantDiscoveryDir/joint_genotyping.filtered.vcf 2>> $logFile;

echo -e "\n\n" >> $logFile;


