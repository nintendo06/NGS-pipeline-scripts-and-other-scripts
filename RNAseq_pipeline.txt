#!/bin/bash

# RNA-seq pipeline script

# Set the necessary software paths
FASTQ_DIR=/path/to/fastq/files
GENOME_INDEX=/path/to/genome/index
GTF_FILE=/path/to/annotation.gtf
OUTPUT_DIR=/path/to/output

# Step 1: Quality control with FastQC
fastqc -o $OUTPUT_DIR $FASTQ_DIR/*.fastq

# Step 2: Adapter and quality trimming with Trim Galore
for file in $FASTQ_DIR/*.fastq; do
    trim_galore --quality 20 --illumina --output_dir $OUTPUT_DIR $file
done

# Step 3: Alignment with HISAT2
for file in $OUTPUT_DIR/*trimmed.fq; do
    base=$(basename $file "_trimmed.fq")
    hisat2 -x $GENOME_INDEX -U $file -S $OUTPUT_DIR/$base.sam
done

# Step 4: Convert SAM to BAM and sort with Samtools
for file in $OUTPUT_DIR/*.sam; do
    base=$(basename $file ".sam")
    samtools view -S -b $file > $OUTPUT_DIR/$base.bam
    samtools sort -o $OUTPUT_DIR/$base.sorted.bam $OUTPUT_DIR/$base.bam
done

# Step 5: Index the sorted BAM files with Samtools
for file in $OUTPUT_DIR/*.sorted.bam; do
    samtools index $file
done

# Step 6: Quantify gene expression with featureCounts
featureCounts -T 4 -a $GTF_FILE -o $OUTPUT_DIR/counts.txt $OUTPUT_DIR/*.sorted.bam

# Step 7: Perform differential expression analysis with DESeq2 or other tools
Rscript differential_expression_analysis.R --input $OUTPUT_DIR/counts.txt --output $OUTPUT_DIR/differential_expression_results.txt


# METHOD TWO
#!/bin/bash

# Specify input and output files/directories
input_directory="input"
output_directory="output"
reference_genome="reference.fasta"
annotation_file="annotation.gtf"

# Step 1: Quality control with FastQC
mkdir -p $output_directory/fastqc
fastqc -o $output_directory/fastqc $input_directory/*.fastq


# Step 2: Alignment with HISAT2
mkdir -p $output_directory/alignment
hisat2-build $reference_genome $output_directory/alignment/reference_index

for file in $input_directory/*.fastq; do
    sample_name=$(basename "$file" .fastq)
    hisat2 -p 4 --dta-cufflinks -x $output_directory/alignment/reference_index -U $file -S $output_directory/alignment/$sample_name.sam
done

# Step 3: Assembly and quantification with StringTie
mkdir -p $output_directory/quantification
for file in $output_directory/alignment/*.sam; do
    sample_name=$(basename "$file" .sam)
    stringtie $file -p 4 -G $annotation_file -o $output_directory/quantification/$sample_name.gtf -l $sample_name
done

# Step 4: Merge transcript assemblies
stringtie --merge -p 4 -G $annotation_file -o $output_directory/quantification/merged.gtf $output_directory/quantification/*.gtf

# Step 5: Estimate transcript abundances
for file in $output_directory/alignment/*.sam; do
    sample_name=$(basename "$file" .sam)
    stringtie -e -B -p 4 -G $output_directory/quantification/merged.gtf -o $output_directory/quantification/$sample_name.gtf $file
done

# Step 6: Run differential expression analysis with DESeq2
Rscript <<EOF
library(DESeq2)

# Read count matrix
count_data <- as.matrix(read.table("$output_directory/quantification/count_matrix.txt", header=TRUE, row.names=1))

# Sample information
sample_info <- read.table("$output_directory/quantification/sample_info.txt", header=TRUE, row.names=1)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=count_data, colData=sample_info, design=~condition)
dds <- DESeq(dds)

# Differential expression analysis
res <- results(dds)
res <- res[order(res$padj), ]
write.table(res, file="$output_directory/quantification/differential_expression_results.txt", sep="\t", quote=FALSE)
