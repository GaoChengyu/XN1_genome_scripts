#!/bin/bash
# Description: Aligns transcripts to a reference genome, processes alignments, and collapses redundant transcripts using TAMA
# Usage: ./part04_tama_del_redundancy.sh <input_fasta> <output_directory> [mode]

# Check if correct number of arguments are provided
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: ./part04_tama_del_redundancy.sh <input_fasta> <output_directory> [mode]"
    echo "  <input_fasta>: Path to the input transcript FASTA file"
    echo "  <output_directory>: Path to the output directory for all generated files"
    echo "  [mode]: Optional processing mode: 'strict' or 'lenient' (default: strict)"
    exit 1
fi

# Assign input parameters to variables
INPUT_FASTA="$1"
OUTPUT_DIR="$2"

# Set processing mode (default to strict if not provided)
if [ $# -eq 3 ]; then
    MODE="$3"
else
    MODE="strict"
fi

# Validate mode parameter
if [ "$MODE" != "strict" ] && [ "$MODE" != "lenient" ]; then
    echo "Error: Invalid mode '$MODE'. Must be 'strict' or 'lenient'."
    exit 1
fi

# Extract sample name from input filename (remove extension)
SAMPLE=$(basename "$INPUT_FASTA" | sed 's/\.[^.]*$//')

# Configuration Section - Modify these variables as needed
THREADS=16                         # Number of threads for parallel processing
REF_GENOME="XN1.genome.chr.fasta"  # Path to reference genome
TAMA_PATH="/mnt/d/sf/tama-master"   # Path to TAMA installation

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Change to output directory to ensure all files are created there
cd "$OUTPUT_DIR" || { echo "Error: Cannot access output directory '$OUTPUT_DIR'"; exit 1; }

# Copy input file to output directory for processing
cp "$INPUT_FASTA" "${SAMPLE}.fa"

echo "Starting transcript processing pipeline"
echo "Input file: $INPUT_FASTA"
echo "Output directory: $OUTPUT_DIR"
echo "Sample: $SAMPLE"
echo "Mode: $MODE"
echo "Threads: $THREADS"

# Step 1: Alignment using Minimap2
# Aligns transcript sequences to reference genome using splice-aware alignment
# -ax splice: Optimized for splice alignment of RNA-seq data
# -uf: Consider only sense strand transcripts (use forward strand)
# --secondary=no: Suppress secondary alignments to reduce redundancy
# -C5: Increased sensitivity for splice site detection
# Output: SAM format alignment file
echo "Step 1: Aligning transcripts to reference genome with Minimap2..."
minimap2 -ax splice -t $THREADS -uf --secondary=no -C5 \
    $REF_GENOME \
    "${SAMPLE}.fa" > "${SAMPLE}_aligned.sam"

# Step 2: SAM Processing and Sorting
# Converts SAM to BAM format, sorts by coordinate, and converts back to SAM
# Sorting is required for many downstream analysis tools
# BAM format is more efficient for storage and processing
echo "Step 2: Processing and sorting alignment files..."

# Convert SAM to sorted BAM for efficient storage
samtools sort -@ $THREADS \
    -o "${SAMPLE}_aligned_sorted.bam" \
    "${SAMPLE}_aligned.sam"

# Convert back to SAM format for TAMA compatibility
# Some tools like TAMA require SAM format rather than BAM
samtools view -@ $THREADS -h \
    "${SAMPLE}_aligned_sorted.bam" > "${SAMPLE}_aligned_sorted.sam"

# Step 3: TAMA Collapse - Remove Redundant Transcripts
# Collapses redundant transcript isoforms based on splice junctions and ends
# Identifies and merges transcripts that likely represent the same isoform
echo "Step 3: Collapsing redundant transcripts with TAMA..."

if [ "$MODE" == "strict" ]; then
    # Strict mode parameters - for high-confidence, polycistronic-filtered transcripts
    # Uses more permissive parameters to capture true biological variation
    echo "Using strict collapse parameters for high-confidence transcripts..."
    # -s: Input SAM file
    # -f: Reference genome FASTA file
    # -p: Output prefix
    # -x no_cap: Don't require Cap analysis (for Iso-Seq data)
    # -e common_ends: Collapse transcripts with common ends
    # -c 0: No 5' end threshold
    # -i 0: No 3' end threshold
    # -m 20: Exon/splice junction tolerance (20bp)
    # -a 99999999: Large 5' end tolerance to capture all variations
    python ${TAMA_PATH}/tama_collapse.py \
        -s "${SAMPLE}_aligned_sorted.sam" \
        -f $REF_GENOME \
        -p "${SAMPLE}_tama" \
        -x no_cap \
        -e common_ends \
        -c 0 \
        -i 0 \
        -m 20 \
        -a 99999999
else
    # Lenient mode parameters - for non-polycistronic-filtered transcripts
    # Uses more stringent parameters to reduce false positives
    echo "Using lenient collapse parameters for standard transcripts..."
    # -m 0: No exon/splice junction tolerance (exact match required)
    # -z 0: No 3' end tolerance (exact match required)
    # -a 1000: 5' end tolerance of 1000bp
    python ${TAMA_PATH}/tama_collapse.py \
        -s "${SAMPLE}_aligned_sorted.sam" \
        -f $REF_GENOME \
        -p "${SAMPLE}_tama" \
        -x no_cap \
        -e common_ends \
        -m 0 \
        -z 0 \
        -a 1000
fi

# Step 4: TAMA Merge - Integrate Transcripts from Multiple Samples
# Merges transcript sets from multiple samples/conditions
# Creates a unified transcriptome annotation
echo "Step 4: Merging transcript sets from multiple samples..."

# Create file list for TAMA merge
# Format: <BED_file> <cap_status> <color_RGB> <sample_name>
cat > filelist.txt <<EOF
${SAMPLE}_tama.bed    no_cap    1,1,1    ${SAMPLE}_sample
# Add other samples here following the same format
# Example: condi_tama.bed no_cap    2,1,1    condi_sample
EOF

echo "File list created for TAMA merge:"
cat filelist.txt

# Merge parameters for creating unified transcriptome:
# -f filelist.txt: Input file containing list of samples to merge
# -p isoseq_tama: Output prefix for merged files
# -e common_ends: Merge transcripts with common ends
# -d merge_dup: Merge duplicate transcripts
# -a 9999: 5' end tolerance for merging
# -m 0: No exon/splice junction tolerance (strict matching)
# -z 0: No 3' end tolerance (strict matching)
echo "Running TAMA merge with the following parameters..."
python ${TAMA_PATH}/tama_merge.py \
    -f filelist.txt \
    -p "isoseq_tama" \
    -e common_ends \
    -d merge_dup \
    -a 9999 \
    -m 0 \
    -z 0

# Cleanup temporary files (optional)
# Uncomment the following lines to remove intermediate files and save disk space
# echo "Cleaning up temporary files..."
# rm "${SAMPLE}_aligned.sam" "${SAMPLE}_aligned_sorted.bam" filelist.txt

echo "Processing complete!"
echo "Input file: $INPUT_FASTA"
echo "Output directory: $OUTPUT_DIR"
echo "Sample: $SAMPLE"
echo "Processing mode: $MODE"
echo "Generated files:"
echo "  - Alignment: ${SAMPLE}_aligned_sorted.sam"
echo "  - Collapsed transcripts: ${SAMPLE}_tama.bed"
echo "  - Merged transcripts: isoseq_tama.bed"
echo "  - Additional TAMA output files with various suffixes"
