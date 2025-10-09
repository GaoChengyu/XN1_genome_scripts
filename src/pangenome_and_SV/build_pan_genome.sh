#!/bin/bash
# Script: Pan-genome Analysis Pipeline
# Description: Constructs a genome graph, identifies core genomic regions, and performs pan-genome analysis
# Usage: ./build_pan_genome.sh <reference_genome> <output_directory> <script_directory> <threads> <genome1> <genome2> ... <genomeN>

# Check if minimum number of arguments are provided
if [ $# -lt 4 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: $0 <reference_genome> <output_directory> <script_directory> <threads> <genome1> <genome2> ... <genomeN>"
    echo "  <reference_genome>: Path to the reference genome FASTA file"
    echo "  <output_directory>: Directory for all output files"
    echo "  <script_directory>: Directory containing helper Python scripts"
    echo "  <threads>: Number of CPU threads for parallel processing"
    echo "  <genome1> <genome2> ...: Paths to additional genome FASTA files"
    exit 1
fi

# Assign input parameters to variables
REF_GENOME="$1"
OUTPUT_DIR="$2"
SCRIPT_DIR="$3"
THREADS="$4"

# Extract additional genome files from remaining arguments
shift 4  # Remove the first 4 arguments
GENOME_LIST=("$@")  # Remaining arguments are genome files

# Validate input files exist
echo "Validating input files..."
if [ ! -f "$REF_GENOME" ]; then
    echo "Error: Reference genome file not found: $REF_GENOME"
    exit 1
fi

for genome in "${GENOME_LIST[@]}"; do
    if [ ! -f "$genome" ]; then
        echo "Error: Genome file not found: $genome"
        exit 1
    fi
done

if [ ! -d "$SCRIPT_DIR" ]; then
    echo "Error: Script directory not found: $SCRIPT_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Change to output directory to ensure all files are created there
cd "$OUTPUT_DIR" || { echo "Error: Cannot access output directory '$OUTPUT_DIR'"; exit 1; }

# Generate output prefix from reference genome name
REF_PREFIX=$(basename "$REF_GENOME" | cut -d. -f1)
OUTPUT_PREFIX="${REF_PREFIX}_pangenome"

echo "Starting Pan-genome Analysis Pipeline"
echo "======================================"
echo "Reference genome: $REF_GENOME"
echo "Output directory: $OUTPUT_DIR"
echo "Script directory: $SCRIPT_DIR"
echo "Threads: $THREADS"
echo "Additional genomes: ${#GENOME_LIST[@]}"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

# Step 1: Build Genome Graph with Minigraph
# Constructs a pangenome graph using the reference genome and additional genomes
# -cxggs: Enable complex graph generation with cigar strings and GFA output
# -t: Number of threads for parallel processing during graph construction
# The graph represents genomic variations and structural differences between genomes
echo "Step 1: Building genome graph with Minigraph..."
echo "This step constructs a pangenome graph representing genomic variations across all input genomes."
minigraph -cxggs -t "$THREADS" \
    "$REF_GENOME" \
    "${GENOME_LIST[@]}" \
    > "${OUTPUT_PREFIX}.gfa"

# Check if graph building was successful
if [ ! -s "${OUTPUT_PREFIX}.gfa" ]; then
    echo "Error: Genome graph construction failed - output file is empty"
    exit 1
fi
echo "Genome graph successfully created: ${OUTPUT_PREFIX}.gfa"
echo ""

# Step 2: Align Genomes to Graph and Generate Coverage Files
# Aligns each genome to the pangenome graph to determine presence/absence of graph segments
# --cov: Generate coverage information for each graph segment
# -x asm: Use assembly-to-graph alignment mode
echo "Step 2: Aligning genomes to pangenome graph and generating coverage files..."
echo "This step aligns each genome to the graph to determine which graph segments are present in each genome."

# Process reference genome first
echo "Processing reference genome: $REF_GENOME"
ref_prefix=$(basename "$REF_GENOME" .fasta | cut -d. -f1)

# Align reference genome to graph
minigraph -t "$THREADS" --cov -x asm \
    "${OUTPUT_PREFIX}.gfa" \
    "$REF_GENOME" \
    > "${ref_prefix}.gaf"

# Convert GAF to TSV coverage format using Python helper script
# The coverage file indicates which graph nodes are covered by this genome
python "$SCRIPT_DIR/comb_coverage01.py" \
    -g "${ref_prefix}.gaf" \
    -a "$ref_prefix" \
    -o "${ref_prefix}Cov.tsv" \
    -r Y

# Process each additional genome
for genome in "${GENOME_LIST[@]}"; do
    genome_prefix=$(basename "$genome" .fasta | cut -d. -f1)
    echo "Processing genome: $genome_prefix"
    
    # Align genome to graph
    minigraph -t "$THREADS" --cov -x asm \
        "${OUTPUT_PREFIX}.gfa" \
        "$genome" \
        > "${genome_prefix}.gaf"
    
    # Check if alignment was successful
    if [ ! -s "${genome_prefix}.gaf" ]; then
        echo "Warning: Alignment failed for genome $genome_prefix - skipping coverage conversion"
        continue
    fi
    
    # Convert GAF to TSV coverage format
    python "$SCRIPT_DIR/comb_coverage01.py" \
        -g "${genome_prefix}.gaf" \
        -a "$genome_prefix" \
        -o "${genome_prefix}Cov.tsv" \
        -r Y
done
echo "Coverage files generated for all genomes"
echo ""

# Step 3: Generate Graph Segment Information (if needed)
# This step may require additional processing to extract graph segment lengths and positions
# The exact method depends on the minigraph version and available utilities
echo "Step 3: Preparing graph segment information..."
echo "Extracting graph segment lengths and positions for core genome analysis."

# Note: This step may need adjustment based on your specific minigraph version
# Some versions may require different commands to extract graph segment information
# The following is a placeholder - adjust based on your minigraph installation
if command -v gfatools &> /dev/null; then
    # Example using gfatools to extract graph information
    gfatools stats "${OUTPUT_PREFIX}.gfa" > graph_stats.txt
    echo "Graph statistics extracted using gfatools"
else
    # Alternative method: parse GFA file directly
    # This extracts node IDs and segment lengths from the GFA file
    grep '^S' "${OUTPUT_PREFIX}.gfa" | awk '{print $2 "\t" length($3)}' > graph_len_prelim.tsv
    echo "Graph segment information extracted from GFA file"
fi
echo ""

# Step 4: Run R Analysis for Core Genome Identification
# Identifies core genomic regions present in all genomes and dispensable regions in subsets
# Uses coverage data to create a binary presence/absence matrix
echo "Step 4: Identifying core and dispensable genomic regions using R..."
echo "This step analyzes coverage data to identify:"
echo "  - Core genome: Regions present in all genomes"
echo "  - Dispensable genome: Regions present in some but not all genomes"

Rscript --vanilla - <<EOF
# Pan-genome Analysis Script
# Identifies core genomic regions present in all genomes and analyzes pan-genome structure

# Load required libraries
library(tidyverse)

# Set working directory to output directory
setwd("$OUTPUT_DIR")

cat("Reading coverage data for all genomes...\n")

# Generate list of coverage files based on processed genomes
# Include reference genome and all successfully processed additional genomes
coverage_files <- c("${ref_prefix}Cov.tsv")

# Add coverage files for additional genomes that were successfully processed
additional_genomes <- c($(printf "\"%s\" " "${GENOME_LIST[@]}" | sed 's/\.fasta/Cov.tsv/g'))
coverage_files <- c(coverage_files, additional_genomes)

# Remove any coverage files that don't exist (in case some genomes failed)
coverage_files <- coverage_files[file.exists(coverage_files)]

cat("Found", length(coverage_files), "coverage files to process\n")

# Read and merge coverage data from all genomes
# Create a binary matrix where 1 = segment present, 0 = segment absent
cat("Creating presence-absence matrix...\n")
datmat <- reduce(
    map(coverage_files, ~ {
        # Read coverage file
        cov_data <- read_tsv(.x, show_col_types = FALSE)
        # Extract genome name from filename
        genome_name <- tools::file_path_sans_ext(basename(.x))
        genome_name <- sub("Cov", "", genome_name)
        
        # Create binary presence/absence data
        cov_data %>% 
            select(nodeid, coverage) %>%
            mutate(coverage = if_else(coverage > 0, 1, 0)) %>%
            rename(!!genome_name := coverage)
    }),
    full_join, 
    by = "nodeid"
)

# Write full node matrix (presence/absence for all segments in all genomes)
write_tsv(datmat, "nodemat.tsv")
cat("Full node matrix written to: nodemat.tsv\n")

# Identify core nodes (present in all genomes)
core_nodeid <- datmat %>%
    filter(if_all(-nodeid, ~ .x == 1)) %>%
    select(nodeid)

write_tsv(core_nodeid, "core_nodemat.tsv")
cat("Core nodes (present in all genomes) written to: core_nodemat.tsv\n")
cat("Number of core nodes:", nrow(core_nodeid), "\n")

# Calculate pan-genome statistics
total_nodes <- nrow(datmat)
core_nodes <- nrow(core_nodeid)
dispensable_nodes <- total_nodes - core_nodes

cat("\nPan-genome Statistics:\n")
cat("Total graph nodes:", total_nodes, "\n")
cat("Core genome nodes:", core_nodes, "\n")
cat("Dispensable genome nodes:", dispensable_nodes, "\n")
cat("Core genome proportion:", round(core_nodes/total_nodes * 100, 2), "%\n")

# Note: The following section requires graph segment length information
# This may need to be adjusted based on how graph_len.tsv is generated
if (file.exists("graph_len.tsv")) {
    cat("Processing graph segment lengths...\n")
    
    # Extract graph segment information (lengths, positions)
    graph_len <- read_delim(
        "graph_len.tsv", 
        col_names = c("nodeid", "length", "chromo", "pos", "rrank"), 
        delim = " "
    )
    
    # Create core region BED file (regions present in all genomes)
    core_graph_len <- graph_len %>%
        inner_join(core_nodeid, by = "nodeid") %>%
        mutate(
            start = as.numeric(pos),
            end = start + length - 1
        ) %>%
        select(chromo, start, end)
    
    write_tsv(core_graph_len, "core_graph_len.bed", col_names = FALSE)
    cat("Core genome regions (BED format) written to: core_graph_len.bed\n")
    
    # Create dispensable region BED file (regions absent in at least one genome)
    dispensable_graph_len <- graph_len %>%
        anti_join(core_nodeid, by = "nodeid") %>%
        mutate(
            start = as.numeric(pos),
            end = start + length - 1
        ) %>%
        select(chromo, start, end)
    
    write_tsv(dispensable_graph_len, "Dispensable_graph_len.bed", col_names = FALSE)
    cat("Dispensable genome regions (BED format) written to: Dispensable_graph_len.bed\n")
    
    # Calculate length-based statistics
    if (nrow(core_graph_len) > 0) {
        core_bp <- sum(graph_len %>% inner_join(core_nodeid, by = "nodeid") %>% pull(length))
        dispensable_bp <- sum(graph_len %>% anti_join(core_nodeid, by = "nodeid") %>% pull(length))
        total_bp <- core_bp + dispensable_bp
        
        cat("\nLength-based Statistics:\n")
        cat("Core genome size:", core_bp, "bp\n")
        cat("Dispensable genome size:", dispensable_bp, "bp\n")
        cat("Total pan-genome size:", total_bp, "bp\n")
        cat("Core genome proportion:", round(core_bp/total_bp * 100, 2), "%\n")
    }
} else {
    cat("Warning: graph_len.tsv not found - skipping BED file generation\n")
    cat("Please ensure graph segment information is available for complete analysis\n")
}

# Generate additional pan-genome statistics by genome
genome_stats <- datmat %>%
    select(-nodeid) %>%
    summarise(across(everything(), ~ sum(.x))) %>%
    pivot_longer(everything(), names_to = "genome", values_to = "nodes_present") %>%
    mutate(total_nodes = total_nodes)

write_tsv(genome_stats, "genome_specific_stats.tsv")
cat("Genome-specific statistics written to: genome_specific_stats.tsv\n")

cat("\nR analysis completed successfully!\n")
EOF

# Check if R analysis completed successfully
if [ $? -ne 0 ]; then
    echo "Warning: R analysis encountered errors, but proceeding with available results"
fi

echo ""
echo "Pan-genome Analysis Pipeline Completed Successfully!"
echo "===================================================="
echo "Output files generated in: $OUTPUT_DIR"
echo ""
echo "Key output files:"
echo "  - ${OUTPUT_PREFIX}.gfa: Pangenome graph in GFA format"
echo "  - nodemat.tsv: Presence/absence matrix of graph segments"
echo "  - core_nodemat.tsv: List of core genome segments"
echo "  - core_graph_len.bed: Core genome regions in BED format"
echo "  - Dispensable_graph_len.bed: Dispensable genome regions in BED format"
echo "  - genome_specific_stats.tsv: Statistics for each genome"
echo ""
echo "Summary:"
echo "  Reference genome: $(basename $REF_GENOME)"
echo "  Additional genomes processed: ${#GENOME_LIST[@]}"
echo "  Output directory: $OUTPUT_DIR"