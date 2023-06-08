#!/bin/bash

# Set the path to the GATK executable
GATK="/path/to/gatk/gatk"

# Set the reference genome file path
REF_GENOME="/path/to/reference_genome.fa"

# Set the output directory
OUTPUT_DIR="/path/to/output_directory"

# List all input BAM files
BAM_DIR="/path/to/input_directory"

# Iterate over each BAM file in the input directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Get the sample name from the BAM file name
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

    # Set the output VCF file path
    OUTPUT_VCF="$OUTPUT_DIR/$SAMPLE_NAME.vcf"

    # Perform variant calling using GATK's HaplotypeCaller
    $GATK HaplotypeCaller \
        -R $REF_GENOME \
        -I $BAM_FILE \
        -O $OUTPUT_VCF
