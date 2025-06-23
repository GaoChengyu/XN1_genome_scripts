#!/bin/bash
set -e  # Exit immediately if any command fails

# Configuration
REFERENCE="../../ChrInfo/XN1.Chr.fasta"
THREADS=16
SAMPLES=("DC1" "YL1")  # Add more samples as needed

# 1. Alignment and SNP calling for each sample
for SAMPLE in "${SAMPLES[@]}"; do
    # Align using minimap2 (asm5 preset for closely related genomes)
    minimap2 -a -x asm5 --cs -r2k -t ${THREADS} ${REFERENCE} ../DC1_JKI.genome.fa > ${SAMPLE}.sam
    
    # Convert SAM to sorted BAM
    samtools sort -@ ${THREADS} -O bam -o ${SAMPLE}.bam ${SAMPLE}.sam
    
    # Index BAM file
    samtools index ${SAMPLE}.bam
    
    # Generate raw variants using bcftools
    bcftools mpileup -Ou ${SAMPLE}.bam -f ${REFERENCE} | \
        bcftools call -mv -o ${SAMPLE}.vcf
    
    # Compress and index VCF
    bcftools sort ${SAMPLE}.vcf -o ${SAMPLE}.vccf.gz -O z
    bcftools index ${SAMPLE}.vcf.gz
done

# 2. Merge all samples' VCF files
bcftools merge -O z -o merged.vcf.gz *.vcf.gz

# 3. Identify sample-specific SNPs (using first sample as reference)
TARGET_SAMPLE="DC1"
BACKGROUND_SAMPLE="YL1"

# Create directory for common sites analysis
mkdir -p common_sites_analysis
cd common_sites_analysis

# Extract common sites
bcftools isec -p . ../${TARGET_SAMPLE}.vcf.gz ../${BACKGROUND_SAMPLE}.vcf.gz

# Filter sites present only in target sample (flag 10 in sites.txt)
awk 'BEGIN{FS=OFS="\t"} $5=="10"' sites.txt > only_${TARGET_SAMPLE}_sites.txt

# Optional: Clean up intermediate files
# find .. -name "*.sam" -delete
# find .. -name "*.bam" -delete