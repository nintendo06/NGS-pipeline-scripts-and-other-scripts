#!/bin/bash

# Specify input and output files/directories
reference_genome="reference.fasta"
input_fastq="sample.fastq"
output_directory="output"

# Step 1: Index the reference genome
bwa index $reference_genome

# Step 2: Align reads to the reference genome
bwa mem -t 4 $reference_genome $input_fastq > aligned.sam

# Step 3: Convert SAM to BAM format and sort
samtools view -S -b aligned.sam > aligned.bam
samtools sort -@ 4 aligned.bam > sorted.bam

# Step 4: Index the sorted BAM file
samtools index sorted.bam

# Step 5: Mark duplicates (optional, if using GATK)
gatk MarkDuplicates \
  -I sorted.bam \
  -O marked_duplicates.bam \
  -M marked_dup_metrics.txt

# Step 6: Base quality score recalibration (optional, if using GATK)
gatk BaseRecalibrator \
  -R $reference_genome \
  -I marked_duplicates.bam \
  --known-sites known_sites.vcf.gz \
  -O recalibrated_data.table

gatk ApplyBQSR \
  -R $reference_genome \
  -I marked_duplicates.bam \
  --bqsr-recal-file recalibrated_data.table \
  -O recalibrated.bam

# Step 7: Variant calling
gatk HaplotypeCaller \
  -R $reference_genome \
  -I recalibrated.bam \
  -O variants.vcf

# Step 8: Filter variants (optional)
gatk VariantFiltration \
  -R $reference_genome \
  -V variants.vcf \
  -O filtered_variants.vcf \
  --filter-expression "QUAL < 30.0" \
  --filter-name "LowQual"

# Move output files to the desired output directory
mkdir -p $output_directory
mv sorted.bam sorted.bam.bai marked_duplicates.bam marked_dup_metrics.txt recalibrated.bam recalibrated_data.table variants.vcf filtered_variants.vcf $output_directory

# Cleanup intermediate files
rm aligned.sam aligned.bam aligned.bam.bai


# Manuplation of vcf file using bcf tools:

#!/bin/bash

# View the contents of a VCF file
bcftools view input.vcf

# Filter variants based on a minimum quality score of 30
bcftools filter -i 'QUAL > 30' input.vcf > filtered.vcf

# Merge multiple VCF files
bcftools merge input1.vcf input2.vcf > merged.vcf

# Subset specific samples from a VCF file
bcftools view -s sample1,sample2 input.vcf > subset.vcf

# Convert VCF to BCF format
bcftools convert input.vcf -o output.bcf

# Convert VCF to plaintext format
bcftools convert input.vcf -o output.txt

# Calculate basic statistics for a VCF file
bcftools stats input.vcf > stats.txt
