#!/bin/bash
# Description: Automated pipeline for Iso-Seq data processing and alignment using SMRT Tools and GMAP
# Usage: ./part01_data_filter.sh <input_file> <output_directory>

# Check if correct number of arguments are provided
if [ $# -ne 3 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: ./part01_data_filter.sh <input_file> <output_directory> <reference_genome>"
    echo "  <input_file>: Path to the input subreads BAM file"
    echo "  <output_directory>: Path to the output directory for all generated files"
    echo "  <reference_genome>: Path to reference genome
    exit 1
fi

# Assign input parameters to variables
INPUT_FILE="$1"
OUTPUT_DIR="$2"
REF_PATH="$3"

# Configuration Section - Modify these variables as needed
THREADS=30                       # Number of threads for parallel processing
SAMPLE="condi"                      # Sample identifier
CONDITION="condi"                # Experimental condition identifier
REF_GENOME="XN1.fasta"           # Reference genome file name
PRIMER="primer.fasta"            # Primer sequences file

# Note: The second subreads file path is now constructed differently since we're only taking one main input
# For the polish step, we'll use the same input file or modify as needed
SECOND_SUBREADS_FILE="$1" # Second subreads file for polishing

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Change to output directory to ensure all files are created there
cd "$OUTPUT_DIR" || { echo "Error: Cannot access output directory '$OUTPUT_DIR'"; exit 1; }

# Step 1: CCS Processing - Generate Circular Consensus Sequences
# This step processes raw subreads to generate high-quality consensus sequences
echo "Step 1: Generating Circular Consensus Sequences (CCS)"
ccs "$INPUT_FILE" \
    "${SAMPLE}.ccs.1.bam" \
    --min-rq=0.70 \
    --min-passes 1 \
    --noPolish \
    -j $THREADS

# Step 2: Lima Processing - Remove primers and orient sequences
# Demultiplexes and removes primers from CCS sequences
echo "Step 2: Removing primers and orienting sequences with Lima"
lima "${SAMPLE}.ccs.1.bam" \
    "$PRIMER" \
    "${SAMPLE}.fl.bam" \
    -j $THREADS \
    --isoseq \
    --peek-guess

# Step 3: Refine Processing - Remove noise and polyA tails
# Filters out low-quality reads and removes polyA tails
echo "Step 3: Refining sequences - removing noise and polyA tails"
isoseq3 refine "${SAMPLE}.fl.primer_5p--primer_3p.bam" \
    "$PRIMER" \
    "${SAMPLE}.flnc.bam" \
    -j $THREADS \
    --require-polya

# Step 4: Cluster Processing - Cluster FLNC reads
# Clusters full-length non-chimeric reads to collapse redundancies
echo "Step 4: Clustering FLNC reads"
isoseq3 cluster "${SAMPLE}.flnc.bam" \
    "${SAMPLE}.clustered.bam" \
    --verbose \
    -j $THREADS

# Step 5: Polish Processing - Polish transcripts using subreads
# Polishes clustered transcripts using original subreads for higher accuracy
echo "Step 5: Polishing transcripts using subreads"
isoseq3 polish -j $THREADS \
    "${CONDITION}.clustered.bam" \
    "$SECOND_SUBREADS_FILE" \
    "${CONDITION}.polished.bam"

# Step 6: GMAP Index Preparation
# Builds GMAP index for the reference genome
echo "Step 6: Building GMAP index"
gmap_build -D gmap_index \
    -d $REF_GENOME \
    $REF_PATH

# Step 7: Combine High and Low Quality Transcripts
# Combines both high-quality and low-quality polished transcripts for alignment
echo "Step 7: Combining high and low quality transcripts"
gunzip "${CONDITION}.polished.hq.fasta.gz" 2>/dev/null || true
gunzip "${CONDITION}.polished.lq.fasta.gz" 2>/dev/null || true
cat "${CONDITION}.polished.hq.fasta" "${CONDITION}.polished.lq.fasta" > "${CONDITION}.polished.hlq.fasta"

# Step 8: GMAP Alignment
# Aligns transcripts to reference genome using GMAP
echo "Step 8: Aligning transcripts to reference genome with GMAP"
gmap -D "gmap_index/${REF_GENOME}" \
    -d $REF_GENOME \
    -f samse \
    -t $THREADS \
    -n 1 \
    --no-chimeras \
    --max-intronlength-middle=20000 \
    --max-intronlength-ends=20000 \
    --min-intronlength=20 \
    --split-large-introns \
    -z sense_force \
    "${CONDITION}.polished.hlq.fasta" > "${CONDITION}.aligned.sam"

# Step 9: SAM Processing and Sorting
# Converts SAM to BAM and sorts the alignment file
echo "Step 9: Processing and sorting SAM file"
samtools sort -@ $THREADS \
    -o "${CONDITION}.aligned.sorted.bam" \
    "${CONDITION}.aligned.sam"

samtools view -@ $THREADS \
    -h "${CONDITION}.aligned.sorted.bam" > "${CONDITION}.aligned.sorted.sam"

echo "Pipeline completed successfully!"
echo "Output files are located in: $OUTPUT_DIR"
