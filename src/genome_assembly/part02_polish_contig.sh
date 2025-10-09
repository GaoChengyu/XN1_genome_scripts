#!/bin/bash

### PacBio/Nanopore Long-Read Polishing Pipeline ###
# Description: Performs iterative polishing of genome assemblies using long-reads (Racon) and short-reads (Pilon)
# Usage: ./part02_polish_contig.sh <output_dir> <initial_assembly> <long_reads> <short_read1> <short_read2> [threads] [pilon_mem]

# =================================================================
# Parameter Validation
# =================================================================

# Check if correct number of arguments are provided
if [ $# -lt 5 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: $0 <output_dir> <initial_assembly> <long_reads> <short_read1> <short_read2> [threads] [pilon_mem]"
    echo "  <output_dir>: Directory for all output files"
    echo "  <initial_assembly>: Initial assembly FASTA file to polish"
    echo "  <long_reads>: Raw long-reads file (PacBio/Nanopore) in FASTA/FASTQ format"
    echo "  <short_read1>: Illumina R1 reads (FASTQ, can be gzipped)"
    echo "  <short_read2>: Illumina R2 reads (FASTQ, can be gzipped)"
    echo "  [threads]: Optional - Number of CPU threads (default: 16)"
    echo "  [pilon_mem]: Optional - Memory for Pilon (default: 52G)"
    exit 1
fi

# Assign input parameters to variables
OUTPUT_DIR="$1"
INITIAL_ASSEMBLY="$2"
LONG_READS="$3"
SHORT_READ1="$4"
SHORT_READ2="$5"

# Set optional parameters with defaults
THREADS="${6:-16}"              # Number of CPU threads (default: 16)
PILON_MEM="${7:-52G}"           # Memory allocation for Pilon (default: 52G)
OUTPUT_PREFIX="polished"        # Prefix for output files

# =================================================================
# Input Validation
# =================================================================

# Validate input files exist
echo "Validating input files..."
for file in "$INITIAL_ASSEMBLY" "$LONG_READS" "$SHORT_READ1" "$SHORT_READ2"; do
    if [ ! -f "$file" ]; then
        echo "Error: Input file not found: $file"
        exit 1
    fi
done

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "Error: Cannot access output directory '$OUTPUT_DIR'"; exit 1; }

echo "Starting polishing pipeline with parameters:"
echo "  Output directory: $OUTPUT_DIR"
echo "  Initial assembly: $INITIAL_ASSEMBLY"
echo "  Long reads: $LONG_READS"
echo "  Short reads: $SHORT_READ1, $SHORT_READ2"
echo "  Threads: $THREADS"
echo "  Pilon memory: $PILON_MEM"

# =================================================================
# Long-read Polishing with Racon (3 iterations)
# =================================================================
# Purpose: Corrects large-scale errors and improves consensus using raw long-read data
# Racon uses partial-order alignment to generate consensus from read overlaps
# Multiple iterations progressively improve consensus quality

echo "=== STARTING LONG-READ POLISHING WITH RACON ==="

# Copy initial assembly to working directory
cp "$INITIAL_ASSEMBLY" "${OUTPUT_PREFIX}_run0.cns.fa"
current_assembly="${OUTPUT_PREFIX}_run0.cns.fa"

for i in {1..3}; do
    echo "Long-read polishing iteration $i/3..."
    
    # Step 1: Create reference index for minimap2
    # Minimap2 index enables fast mapping of reads to the assembly
    echo "  Creating reference index..."
    minimap2 -d "ref_long_$i.mmi" "$current_assembly"
    
    # Step 2: Map long reads to current assembly
    # -ax map-pb: Optimized preset for PacBio reads (use -ax map-ont for Nanopore)
    # This generates SAM format alignments for consensus generation
    echo "  Mapping long reads to assembly..."
    minimap2 -ax map-pb -t "$THREADS" "ref_long_$i.mmi" "$LONG_READS" > "align_long_$i.sam"
    
    # Step 3: Polish assembly with Racon
    # Racon generates consensus sequence using the alignment and raw reads
    # It corrects indels and substitutions based on read evidence
    echo "  Generating consensus with Racon..."
    racon -t "$THREADS" "$LONG_READS" "align_long_$i.sam" "$current_assembly" > "${OUTPUT_PREFIX}_run$i.cns.fa"
    
    # Update assembly for next iteration
    current_assembly="${OUTPUT_PREFIX}_run$i.cns.fa"
    echo "  Completed iteration $i. Output: $current_assembly"
    
    # Cleanup intermediate files from this iteration
    rm "ref_long_$i.mmi" "align_long_$i.sam"
done

LONG_POLISHED="$current_assembly"
echo "Long-read polishing completed. Final output: $LONG_POLISHED"

# =================================================================
# Short-read Polishing with Pilon (3 iterations)
# =================================================================
# Purpose: Corrects small-scale errors and improves base-level accuracy using Illumina data
# Pilon uses high-quality short reads to fix SNPs, small indels, and gaps
# Multiple iterations ensure comprehensive error correction

echo "=== STARTING SHORT-READ POLISHING WITH PILON ==="

current_assembly="$LONG_POLISHED"

for i in {1..3}; do
    echo "Short-read polishing iteration $i/3..."
    
    # Step 1: Map short reads to current assembly using minimap2
    # -ax sr: Optimized preset for short reads
    # Output is piped directly to samtools for BAM conversion
    echo "  Mapping short reads to assembly..."
    minimap2 -ax sr -t "$THREADS" "$current_assembly" "$SHORT_READ1" "$SHORT_READ2" \
        | samtools view -@ "$THREADS" -bS - > "ngs_align_$i.bam"
    
    # Step 2: Sort BAM file by coordinate
    # Sorting is required for duplicate marking and efficient access
    echo "  Sorting BAM file..."
    samtools sort -@ "$THREADS" "ngs_align_$i.bam" -o "ngs_sorted_$i.bam"
    
    # Step 3: Mark PCR duplicates
    # Sambamba identifies and marks duplicate reads from PCR amplification
    # This prevents over-counting of duplicated fragments
    echo "  Marking PCR duplicates..."
    if command -v sambamba &> /dev/null; then
        sambamba markdup -t "$THREADS" "ngs_sorted_$i.bam" "ngs_markdup_$i.bam"
    else
        echo "  Warning: sambamba not found, using samtools markdup instead..."
        samtools markdup -@ "$THREADS" "ngs_sorted_$i.bam" "ngs_markdup_$i.bam"
    fi
    
    # Step 4: Filter high-quality alignments
    # -q 30: Keep only alignments with MAPQ ≥ 30 (high confidence)
    # This ensures only reliable alignments are used for polishing
    echo "  Filtering high-quality alignments..."
    samtools view -@ "$THREADS" -b -q 30 "ngs_markdup_$i.bam" > "ngs_filtered_$i.bam"
    samtools index -@ "$THREADS" "ngs_filtered_$i.bam"
    
    # Step 5: Polish assembly with Pilon
    # Pilon uses the aligned short reads to correct various error types:
    # --fix all: Correct SNPs, indels, gaps, and local misassemblies
    # --changes: Output list of changes made
    # --vcf: Generate VCF file with variant calls
    echo "  Polishing with Pilon..."
    pilon_out="${OUTPUT_PREFIX}_pilon$i"
    
    # Check if pilon.jar is available
    if [ ! -f "pilon.jar" ]; then
        # Try to find pilon in PATH
        if command -v pilon &> /dev/null; then
            pilon_cmd="pilon"
        else
            echo "Error: pilon.jar not found in current directory and 'pilon' not in PATH"
            echo "Please ensure Pilon is installed and accessible"
            exit 1
        fi
    else
        pilon_cmd="java -Xmx$PILON_MEM -jar pilon.jar"
    fi
    
    $pilon_cmd \
        --genome "$current_assembly" \
        --frags "ngs_filtered_$i.bam" \
        --output "$pilon_out" \
        --fix all \
        --changes \
        --vcf \
        --threads "$THREADS" > "pilon_iter$i.log" 2>&1
    
    # Check if Pilon produced output
    if [ ! -f "${pilon_out}.fasta" ]; then
        echo "Error: Pilon iteration $i failed to produce output file"
        echo "Check pilon_iter$i.log for details"
        exit 1
    fi
    
    # Update assembly for next iteration
    current_assembly="${pilon_out}.fasta"
    echo "  Completed iteration $i. Output: $current_assembly"
    
    # Cleanup intermediate files from this iteration
    rm "ngs_align_$i.bam" "ngs_sorted_$i.bam" "ngs_markdup_$i.bam" "ngs_filtered_$i.bam" "ngs_filtered_$i.bam.bai"
done

FINAL_ASSEMBLY="$current_assembly"

# =================================================================
# Final Summary and Output
# =================================================================

echo "=== POLISHING PIPELINE COMPLETED SUCCESSFULLY ==="
echo "Final polished assembly: $FINAL_ASSEMBLY"
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Generated files:"
echo "  - Final assembly: $FINAL_ASSEMBLY"
echo "  - Intermediate polished assemblies: ${OUTPUT_PREFIX}_run*.cns.fa, ${OUTPUT_PREFIX}_pilon*.fasta"
echo "  - Pilon logs: pilon_iter*.log"
echo "  - Pilon change reports: ${OUTPUT_PREFIX}_pilon*.changes"
echo "  - Pilon VCF files: ${OUTPUT_PREFIX}_pilon*.vcf"

# Create a symbolic link to the final assembly for easy access
ln -sf "$FINAL_ASSEMBLY" "final_polished_assembly.fasta"
echo ""
echo "Symbolic link created: final_polished_assembly.fasta -> $FINAL_ASSEMBLY"