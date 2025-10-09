#!/bin/bash
# SV Validation Pipeline Using Pre-aligned BAM Files v1.2
# Input: Reference genome, SV annotation BED, directory of aligned BAM files
# Output: Validated SV classification (TP/FP), coverage matrix, visualization

###############################################################################
### PARAMETER CONFIGURATION
###############################################################################
REF="XN1.Chr.fasta"                      # Reference genome in FASTA format
SV_BED="XN1_DEL.bed"          # BED file of insertion SVs (4+ columns: chr, start, end, SV_ID)
MIN_COVERAGE=5                       # Minimum read depth to validate presence
MIN_SAMPLES=2                        # Minimum samples supporting coverage for FP classification
BAM_DIR="bam_dir"      # Directory containing sorted, indexed BAM files
THREADS=8                            # CPU threads for parallel processing

###############################################################################
### 1. INPUT VALIDATION
###############################################################################
# Check if BAM directory exists
if [ ! -d "$BAM_DIR" ]; then
    echo "Error: BAM directory $BAM_DIR not found"
    exit 1
fi

# Check if BAM files are present
BAM_FILES=($BAM_DIR/*.bam)
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "Error: No BAM files found in $BAM_DIR"
    exit 1
fi

# Verify BAM indexes exist or create them
echo "Verifying BAM indexes..."
for BAM in "${BAM_FILES[@]}"; do
    if [ ! -f "${BAM}.bai" ]; then
        echo "Indexing ${BAM##*/}..."
        samtools index -@ $THREADS "$BAM"
    fi
done

###############################################################################
### 2. SV REGION COVERAGE ANALYSIS
###############################################################################
echo "Calculating SV coverage..."

# Generate coverage matrix across all SV regions
bedtools multicov -bams "${BAM_FILES[@]}" \
    -bed "$SV_BED" -f 0.9 -r > sv_coverage_matrix.tsv

# Check if coverage matrix was generated
if [ ! -s "sv_coverage_matrix.tsv" ]; then
    echo "Error: Failed to generate coverage matrix"
    exit 1
fi

###############################################################################
### 3. FALSE POSITIVE SV IDENTIFICATION
###############################################################################
echo "Identifying false positive SVs..."

awk -v min_cov=$MIN_COVERAGE -v min_sam=$MIN_SAMPLES '
BEGIN {
    OFS="\t";
    # Print header for output file
    print "chr", "start", "end", "sv_id", "covered_samples", "validation_status"
}
# Skip header row from bedtools output
NR>1 {
    sample_count = 0;  # Reset counter for each SV region
    
    # Iterate coverage columns (columns 5 to NF)
    for(i=5; i<=NF; i++) {
        if($i >= min_cov) sample_count++  # Count samples meeting coverage threshold
    }
    
    # Classification logic:
    # FP = False Positive (present in reference)
    # TP = True Positive (genuine insertion)
    status = (sample_count >= min_sam) ? "FP" : "TP";
    
    # Output results with original coordinates and ID
    print $1, $2, $3, $4, sample_count, status
}' sv_coverage_matrix.tsv > sv_validation_results.tsv

###############################################################################
### 4. RESULT VISUALIZATION (REQUIRES R ENVIRONMENT)
###############################################################################
echo "Generating visualization..."

Rscript -e '
library(ggplot2)
# Import validation results
data <- read.delim("sv_validation_results.tsv", header=TRUE)

# Create bar plot showing validation distribution
pdf("SV_Validation_Report.pdf", width=10, height=6)
ggplot(data, aes(x=covered_samples, fill=validation_status)) +
  geom_bar(alpha=0.85) +
  geom_vline(xintercept='$MIN_SAMPLES'-0.5, linetype="dashed", color="red", linewidth=0.8) +
  scale_fill_manual(values=c("FP"="#FF6B6B", "TP"="#4ECDC4")) +
  labs(title="Insertion SV Validation by Long-Read Coverage",
       subtitle=paste("Threshold:", '$MIN_SAMPLES', "samples with >=", '$MIN_COVERAGE', "x coverage"),
       x="Number of Supporting Samples", 
       y="SV Count",
       fill="Validation Status") +
  theme_minimal(base_size=12) +
  theme(legend.position="top")
dev.off()
'

###############################################################################
### OUTPUT FILES SUMMARY
###############################################################################
echo "### PROCESS COMPLETED ###"
echo "Coverage matrix: sv_coverage_matrix.tsv"
echo "Validation results: sv_validation_results.tsv"
echo "Visualization report: SV_Validation_Report.pdf"