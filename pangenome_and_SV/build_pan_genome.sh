#!/bin/bash
# Script: Pan-genome Analysis Pipeline
# Description: Constructs a genome graph, identifies core genomic regions, and performs pan-genome analysis
# Author: [Your Name]
# Date: [Date]

# Configuration Section - Modify these as needed
THREADS=48                         # Number of threads for parallel processing
REF_GENOME="XN1.Chr.fasta"         # Reference genome file
GENOME_LIST=("DC1_JKI.genome.fasta" "MC2.fasta" "NBRC30405.genome.fasta" 
             "NL1.genome.fa" "PGHB.genome.fa" "YL1.genome.fa")  # Genome files to include
OUTPUT_PREFIX="seven_genome"        # Output file prefix
SCRIPT_DIR="/mnt/e/python脚本"      # Directory containing helper scripts
WORK_DIR="$(pwd)"                   # Working directory

# Step 1: Build Genome Graph with Minigraph
# -cxggs: Enable cigar generation and gfa output
# -t: Thread count for parallel processing
echo "Building genome graph with minigraph..."
minigraph -cxggs -t$THREADS \
    $REF_GENOME \
    ${GENOME_LIST[@]} \
    > ${OUTPUT_PREFIX}.gfa

# Step 2: Align Genomes to Graph and Generate Coverage Files
for genome in $REF_GENOME ${GENOME_LIST[@]}; do
    genome_prefix=$(basename $genome .fasta | cut -d. -f1)
    
    echo "Processing $genome_prefix genome..."
    
    # Align genome to graph
    minigraph -t $THREADS --cov -x asm \
        ${OUTPUT_PREFIX}.gfa \
        $genome \
        > ${genome_prefix}.gaf
    
    # Convert GAF to TSV coverage format
    python $SCRIPT_DIR/comb_coverage01.py \
        -g ${genome_prefix}.gaf \
        -a $genome_prefix \
        -o ${genome_prefix}Cov.tsv \
        -r Y
done

# Step 3: Run R Analysis for Core Genome Identification
echo "Identifying core genomic regions..."
Rscript --vanilla - <<EOF
# Pan-genome Analysis Script
# Identifies core genomic regions present in all genomes

# Load required libraries
library(tidyverse)

# Set working directory
setwd("$WORK_DIR")

# Read coverage data for all genomes
coverage_files <- c(
    "XN1Cov.tsv", "DC1Cov.tsv", "MC2Cov.tsv", 
    "NBRC2Cov.tsv", "NL1Cov.tsv", "YL1Cov.tsv", "PGHBCov.tsv"
)

# Create node matrix (1 = present, 0 = absent)
datmat <- reduce(
    map(coverage_files, ~ read_tsv(.x) %>% select(nodeid, coverage)),
    full_join, by = "nodeid"
) %>%
    mutate(across(-nodeid, ~ if_else(.x > 0, 1, 0))) %>%
    set_names(c("nodeid", str_replace(coverage_files, "Cov.tsv", "")))

# Write full node matrix
write_tsv(datmat, "nodemat.tsv")

# Identify core nodes (present in all genomes)
core_nodeid <- datmat %>%
    filter(if_all(-nodeid, ~ .x == 1)) %>%
    select(nodeid)

write_tsv(core_nodeid, "core_nodemat.tsv")

# Extract graph segment information
graph_len <- read_delim(
    "graph_len.tsv", 
    col_names = c("nodeid", "length", "chromo", "pos", "rrank"), 
    delim = " "
)

# Create core region BED file
core_graph_len <- graph_len %>%
    inner_join(core_nodeid, by = "nodeid") %>%
    mutate(
        start = as.numeric(pos),
        end = start + length - 1
    ) %>%
    select(chromo, start, end)

write_tsv(core_graph_len, "core_graph_len.bed", col_names = FALSE)

# Create dispensable region BED file
dispensable_graph_len <- graph_len %>%
    anti_join(core_nodeid, by = "nodeid") %>%
    mutate(
        start = as.numeric(pos),
        end = start + length - 1
    ) %>%
    select(chromo, start, end)

write_tsv(dispensable_graph_len, "Dispensable_graph_len.bed", col_names = FALSE)
EOF

# Save results
write_tsv(datpan, "pan_genome_results.tsv")
EOF