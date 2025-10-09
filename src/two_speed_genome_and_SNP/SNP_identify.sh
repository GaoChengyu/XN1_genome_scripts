#!/bin/bash
set -e  # Exit immediately if any command fails to prevent partial execution

# SNP Identification Pipeline
# Description: Performs alignment, variant calling, and identifies sample-specific SNPs
# Usage: ./SNP_identify.sh <reference_fasta> <output_directory> <threads> <sample1_fasta> [sample2_fasta ...]

# Check if minimum number of arguments are provided
if [ $# -lt 4 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: $0 <reference_fasta> <output_directory> <threads> <sample1_fasta> [sample2_fasta ...]"
    echo "  <reference_fasta>: Path to reference genome in FASTA format"
    echo "  <output_directory>: Directory for all output files"
    echo "  <threads>: Number of CPU threads for parallel processing"
    echo "  <sample1_fasta>: Path to first sample genome in FASTA format"
    echo "  [sample2_fasta ...]: Additional sample genomes (optional)"
    exit 1
fi

# Assign input parameters to variables
REFERENCE="$1"
OUTPUT_DIR="$2"
THREADS="$3"
shift 3  # Remove first three arguments, remaining are sample files

# Store sample files in array
SAMPLE_FILES=("$@")

# Validate input files exist
echo "Validating input files..."
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference genome file not found: $REFERENCE"
    exit 1
fi

for sample_file in "${SAMPLE_FILES[@]}"; do
    if [ ! -f "$sample_file" ]; then
        echo "Error: Sample file not found: $sample_file"
        exit 1
    fi
done

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "Error: Cannot access output directory '$OUTPUT_DIR'"; exit 1; }

echo "Starting SNP Identification Pipeline"
echo "====================================="
echo "Reference genome: $REFERENCE"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "Samples to process: ${#SAMPLE_FILES[@]}"
for i in "${!SAMPLE_FILES[@]}"; do
    echo "  Sample $((i+1)): ${SAMPLE_FILES[i]}"
done
echo ""

# Extract sample names from file paths for use in output files
declare -a SAMPLE_NAMES
for sample_file in "${SAMPLE_FILES[@]}"; do
    # Extract base filename without extension
    base_name=$(basename "$sample_file")
    sample_name="${base_name%.*}"
    sample_name="${sample_name%.*}"  # Remove double extension if present
    SAMPLE_NAMES+=("$sample_name")
done

###############################################################################
### 1. ALIGNMENT AND SNP CALLING FOR EACH SAMPLE
###############################################################################

echo "Step 1: Performing alignment and SNP calling for each sample..."
echo "This step aligns each sample genome to the reference and calls variants."

for i in "${!SAMPLE_FILES[@]}"; do
    SAMPLE_FILE="${SAMPLE_FILES[i]}"
    SAMPLE_NAME="${SAMPLE_NAMES[i]}"
    
    echo "Processing sample: $SAMPLE_NAME"
    
    # Step 1.1: Alignment using minimap2
    # -a: Output in SAM format (required for downstream processing)
    # -x asm5: Optimized preset for closely related genomes (~5% divergence)
    # --cs: Output cigar string for precise alignment information
    # -r2k: Set chaining and alignment bandwidth (2k bases)
    # -t: Number of threads for parallel alignment
    echo "  Aligning sample to reference with minimap2..."
    minimap2 -a -x asm5 --cs -r2k -t "$THREADS" "$REFERENCE" "$SAMPLE_FILE" > "${SAMPLE_NAME}.sam"
    
    # Check if alignment was successful
    if [ ! -s "${SAMPLE_NAME}.sam" ]; then
        echo "Error: Alignment failed for sample $SAMPLE_NAME"
        exit 1
    fi
    
    # Step 1.2: Convert SAM to sorted BAM format
    # SAMtools sort: Converts SAM to BAM and sorts by coordinate
    # -@: Number of threads for sorting
    # -O bam: Output format is BAM
    # Sorting is required for efficient downstream processing and indexing
    echo "  Converting SAM to sorted BAM format..."
    samtools sort -@ "$THREADS" -O bam -o "${SAMPLE_NAME}.bam" "${SAMPLE_NAME}.sam"
    
    # Step 1.3: Index BAM file for efficient access
    # BAM index allows random access to specific genomic regions
    # Required for variant calling and visualization tools
    echo "  Indexing BAM file..."
    samtools index "${SAMPLE_NAME}.bam"
    
    # Step 1.4: Generate raw variants using bcftools
    # bcftools mpileup: Generates genotype likelihoods from BAM alignment
    # -Ou: Output uncompressed BCF for piping
    # -f: Reference genome for base alignment
    # bcftools call: Calls variants from the pileup data
    # -mv: Use multi-allelic caller and output variant sites only
    # -o: Output VCF file with raw variant calls
    echo "  Calling variants with bcftools..."
    bcftools mpileup -Ou "${SAMPLE_NAME}.bam" -f "$REFERENCE" | \
        bcftools call -mv -o "${SAMPLE_NAME}.vcf"
    
    # Check if variant calling was successful
    if [ ! -s "${SAMPLE_NAME}.vcf" ]; then
        echo "Warning: No variants called for sample $SAMPLE_NAME"
        # Create empty VCF to allow pipeline to continue
        touch "${SAMPLE_NAME}.vcf"
    fi
    
    # Step 1.5: Compress and index VCF file
    # Compression reduces file size and enables indexing
    # Indexing allows efficient querying of specific genomic regions
    echo "  Compressing and indexing VCF file..."
    bcftools sort "${SAMPLE_NAME}.vcf" -o "${SAMPLE_NAME}.vcf.gz" -O z
    bcftools index "${SAMPLE_NAME}.vcf.gz"
    
    echo "  Completed processing for sample: $SAMPLE_NAME"
    echo ""
done

###############################################################################
### 2. MERGE ALL SAMPLES' VCF FILES
###############################################################################

echo "Step 2: Merging VCF files from all samples..."
echo "This step combines variant calls from all samples into a unified VCF file."

# Check if we have multiple samples to merge
if [ ${#SAMPLE_FILES[@]} -gt 1 ]; then
    # bcftools merge: Combines multiple VCF/BCF files
    # -O z: Output compressed VCF format
    # -o: Output file for merged variants
    # This creates a multi-sample VCF with genotypes for all samples
    echo "Merging VCF files from ${#SAMPLE_FILES[@]} samples..."
    bcftools merge -O z -o "merged.vcf.gz" *.vcf.gz
    
    # Check if merge was successful
    if [ ! -s "merged.vcf.gz" ]; then
        echo "Warning: Merged VCF file is empty"
    else
        echo "Merged VCF created: merged.vcf.gz"
    fi
else
    echo "Only one sample provided, skipping merge step"
    # For single sample, use its VCF as the "merged" file for downstream compatibility
    cp "${SAMPLE_NAMES[0]}.vcf.gz" "merged.vcf.gz"
    cp "${SAMPLE_NAMES[0]}.vcf.gz.csi" "merged.vcf.gz.csi"
fi

###############################################################################
### 3. IDENTIFY SAMPLE-SPECIFIC SNPS
###############################################################################

echo "Step 3: Identifying sample-specific SNPs..."
echo "This step identifies SNPs unique to specific samples by comparing VCF files."

# Only proceed with sample-specific analysis if we have at least 2 samples
if [ ${#SAMPLE_FILES[@]} -ge 2 ]; then
    # Use first sample as target and second as background for comparison
    # Users can modify this logic based on their experimental design
    TARGET_SAMPLE="${SAMPLE_NAMES[0]}"
    BACKGROUND_SAMPLE="${SAMPLE_NAMES[1]}"
    
    echo "Comparing $TARGET_SAMPLE (target) vs $BACKGROUND_SAMPLE (background)"
    
    # Create directory for common sites analysis
    mkdir -p "common_sites_analysis"
    cd "common_sites_analysis" || { echo "Error: Cannot access common_sites_analysis directory"; exit 1; }
    
    # Step 3.1: Extract common and unique sites using bcftools isec
    # bcftools isec: Intersects multiple VCF files to find common/unique variants
    # -p .: Output prefix (current directory)
    # This creates multiple files categorizing variants by their presence/absence
    echo "  Intersecting VCF files to identify unique variants..."
    bcftools isec -p . "../${TARGET_SAMPLE}.vcf.gz" "../${BACKGROUND_SAMPLE}.vcf.gz"
    
    # Step 3.2: Filter sites present only in target sample
    # The 'sites.txt' file contains information about variant presence:
    # - Flag values indicate which files contain each variant
    # - Flag "10" indicates variant is only in the first file (target sample)
    # This identifies sample-specific SNPs not found in the background
    echo "  Extracting variants unique to target sample..."
    if [ -f "sites.txt" ]; then
        awk 'BEGIN{FS=OFS="\t"} $5=="10"' "sites.txt" > "only_${TARGET_SAMPLE}_sites.txt"
        
        # Count unique variants
        unique_count=$(wc -l < "only_${TARGET_SAMPLE}_sites.txt" 2>/dev/null || echo 0)
        echo "  Found $unique_count variants unique to $TARGET_SAMPLE"
    else
        echo "Warning: sites.txt not found - intersection may have failed"
        touch "only_${TARGET_SAMPLE}_sites.txt"
    fi
    
    cd ..  # Return to main output directory
else
    echo "Only one sample provided, skipping sample-specific SNP analysis"
    mkdir -p "common_sites_analysis"
    touch "common_sites_analysis/only_${SAMPLE_NAMES[0]}_sites.txt"
fi

###############################################################################
### 4. GENERATE SUMMARY STATISTICS
###############################################################################

echo "Step 4: Generating summary statistics..."
echo "This step provides an overview of the SNP calling results."

# Create summary report
SUMMARY_FILE="snp_analysis_summary.txt"
{
    echo "SNP Identification Pipeline Summary"
    echo "==================================="
    echo "Reference genome: $REFERENCE"
    echo "Output directory: $OUTPUT_DIR"
    echo "Number of samples: ${#SAMPLE_FILES[@]}"
    echo "Date: $(date)"
    echo ""
    echo "Sample Details:"
    for i in "${!SAMPLE_NAMES[@]}"; do
        echo "  ${SAMPLE_NAMES[i]}: ${SAMPLE_FILES[i]}"
    done
    echo ""
    echo "Variant Statistics:"
    
    # Count variants for each sample
    for sample_name in "${SAMPLE_NAMES[@]}"; do
        if [ -f "${sample_name}.vcf.gz" ]; then
            variant_count=$(bcftools view "${sample_name}.vcf.gz" | grep -v '^#' | wc -l 2>/dev/null || echo 0)
            echo "  ${sample_name}: $variant_count variants"
        fi
    done
    
    # Report on sample-specific variants if available
    if [ -f "common_sites_analysis/only_${TARGET_SAMPLE}_sites.txt" ]; then
        unique_count=$(wc -l < "common_sites_analysis/only_${TARGET_SAMPLE}_sites.txt" 2>/dev/null || echo 0)
        echo ""
        echo "Sample-specific variants:"
        echo "  ${TARGET_SAMPLE}-specific: $unique_count variants"
    fi
    
    echo ""
    echo "Output Files:"
    echo "  Individual sample VCFs: ${SAMPLE_NAMES[*]}.vcf.gz"
    echo "  Merged VCF: merged.vcf.gz"
    echo "  Sample-specific SNPs: common_sites_analysis/only_${TARGET_SAMPLE}_sites.txt"
    echo "  This summary: $SUMMARY_FILE"
    
} > "$SUMMARY_FILE"

echo "Summary report generated: $SUMMARY_FILE"

###############################################################################
### 5. OPTIONAL CLEANUP OF INTERMEDIATE FILES
###############################################################################

# Uncomment the following lines to remove intermediate files and save disk space
# echo "Step 5: Cleaning up intermediate files..."
# find . -name "*.sam" -delete
# find . -name "*.bam" -delete
# find . -name "*.bai" -delete
# find . -name "*.vcf" -delete  # Keep only compressed VCFs

echo ""
echo "SNP Identification Pipeline Completed Successfully!"
echo "==================================================="
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Generated files:"
echo "  - Individual sample VCFs: ${SAMPLE_NAMES[*]}.vcf.gz"
echo "  - Merged variants: merged.vcf.gz"
echo "  - Sample-specific SNPs: common_sites_analysis/only_${TARGET_SAMPLE}_sites.txt"
echo "  - Summary report: $SUMMARY_FILE"
echo ""
echo "Next steps:"
echo "  - Validate SNP calls using IGV or other visualization tools"
echo "  - Annotate variants with functional consequences"
echo "  - Perform population genetics or association analysis"