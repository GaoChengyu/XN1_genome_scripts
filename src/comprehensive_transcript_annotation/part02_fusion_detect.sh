#!/bin/bash
# Description: Processes polished transcripts through alignment, fusion detection, and filtering
# Usage: ./part02_fusion_detect.sh <input_fasta> <output_directory> <ref_genome> <cluster_report>

# Check if correct number of arguments are provided
if [ $# -ne 4 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: ./part02_fusion_detect.sh <input_fasta> <output_directory>"
    echo "  <input_fasta>: Path to the input polished transcripts FASTA file (e.g., sample_clean.fa)"
    echo "  <output_directory>: Path to the output directory for all generated files"
    echo "  <ref_genome>: Reference genome path"
    echo "  <cluster_report>: Cluster report from part01_data_filter.sh"
    exit 1
fi

# Assign input parameters to variables
INPUT_FASTA="$1"
OUTPUT_DIR="$2"
REF_GENOME="$3" 
CLUSTER_REPORT="$4"

# Extract sample name from input filename for use in output files
SAMPLE=$(basename "$INPUT_FASTA" | sed 's/_clean.fa//')

# Configuration Section - Modify these as needed
THREADS=30                       # Number of threads for parallel processing

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Change to output directory to ensure all files are created there
cd "$OUTPUT_DIR" || { echo "Error: Cannot access output directory '$OUTPUT_DIR'"; exit 1; }

# Step 1: Alignment using Minimap2
# Aligns polished transcripts to reference genome using splice-aware alignment
# -ax splice: Optimize for splice alignment
# -uf: Use forward strand for spliced alignment
# --secondary=no: Suppress secondary alignments to reduce redundancy
# -C5: Output CIGAR with >5 operations for detailed alignment information
echo "Step 1: Aligning transcripts to reference genome with Minimap2"
minimap2 -ax splice -t $THREADS -uf --secondary=no -C5 \
    $REF_GENOME \
    "$INPUT_FASTA" > "${SAMPLE}_clean.sam"

# Step 2: SAM Processing and Format Conversion
# Converts SAM to sorted BAM and back to SAM for compatibility with downstream tools
echo "Step 2: Processing and converting SAM file formats"

# Convert SAM to sorted BAM for efficient storage and processing
samtools sort -@ $THREADS \
    -o "${SAMPLE}_clean.sort.bam" \
    "${SAMPLE}_clean.sam"

# Convert back to SAM format for fusion_finder compatibility
# Some tools like fusion_finder require SAM format rather than BAM
samtools view -@ $THREADS -h \
    "${SAMPLE}_clean.sort.bam" > "${SAMPLE}_clean.sort.sam"

# Step 3: Fusion Transcript Detection
# Identifies potential fusion transcripts using cluster information
# Requires cds_Cupcake package installed from GitHub
# --input: Polished transcripts FASTA file
# -s: Sorted SAM alignment file
# --cluster_report_csv: Cluster report from initial Iso-Seq processing
# -o: Output prefix for fusion detection results
echo "Step 3: Detecting fusion transcripts"
fusion_finder.py \
    --input "$INPUT_FASTA" \
    -s "${SAMPLE}_clean.sort.sam" \
    -o "${SAMPLE}_fusion" \
    --cluster_report_csv $CLUSTER_REPORT

# Step 4: Filter Fusion Transcripts
# Removes detected fusion transcripts from the final transcript set
echo "Step 4: Filtering out fusion transcripts"

# Extract fusion transcript headers from the group file
# The group file contains information about potential fusion events
awk 'BEGIN{FS="\t";OFS="\t"}{print $2}' "${SAMPLE}_fusion.group.txt" | 
    awk '{gsub(",","\n",$0);print $0}' > "${SAMPLE}_del_fusion_header.txt"

# Remove fusion transcripts from the original FASTA using pyfasta
# --header: Operate on sequence headers
# --exclude: Exclude sequences listed in the header file
# --file: File containing headers of sequences to exclude
pyfasta extract --header \
    --fasta "$INPUT_FASTA" \
    --exclude --file "${SAMPLE}_del_fusion_header.txt" \
    > "${SAMPLE}_clean_nofusion.fa"

# Step 5: Cleanup Temporary Files (Optional)
# Remove intermediate SAM files to save disk space
# Comment out these lines if you need to keep intermediate files for debugging
echo "Step 5: Cleaning up temporary files"
rm "${SAMPLE}_clean.sam" "${SAMPLE}_clean.sort.sam"

echo "Processing complete for sample: $SAMPLE"
echo "Input file: $INPUT_FASTA"
echo "Output directory: $OUTPUT_DIR"
echo "Final filtered transcripts: ${SAMPLE}_clean_nofusion.fa"
echo "Fusion detection results: ${SAMPLE}_fusion.*"


