#!/bin/bash

# LTR Prediction Pipeline
# Integrates LTR_Finder, LTR_Harvest, and LTR_retriever for comprehensive LTR retrotransposon discovery
# Usage: ./ltr_pipeline.sh <genome.fasta>

# Check input parameter
if [ $# -ne 1 ]; then
    echo "Error: Please specify genome file"
    echo "Usage: ./ltr_pipeline.sh genome.fasta"
    exit 1
fi

# Initialize variables
GENOME="$1"
PREFIX=$(basename "$GENOME" .fasta)  # Extract base name for output files

####################################
# Step 1: Run LTR_Finder prediction
####################################
echo "Running LTR_Finder..."
ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.85 "$GENOME" > "${PREFIX}.finder.scn"

####################################
# Step 2: Run LTR_Harvest prediction
# Part A: Create enhanced suffix array
####################################
echo "Building suffix array..."
gt suffixerator -db "$GENOME" -indexname "$PREFIX" -tis -suf -lcp -des -ssp -sds -dna

# Part B: Run LTR_Harvest
####################################
echo "Running LTR_Harvest..."
nohup gt ltrharvest -index "$PREFIX" -similar 85 -vic 10 -seed 20 -seqids yes \
      -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 \
      -motif TGCA -motifmis 1 > "${PREFIX}.harvest.TGCA.scn" 1>ltrharvest.log 2>&1 &

wait  # Ensure harvest completes before next step

####################################
# Step 3: Integrate predictions with LTR_retriever
# Outputs to ltr_retriever.log
####################################
echo "Integrating results with LTR_retriever..."
nohup LTR_retriever -genome "$GENOME" -inharvest "${PREFIX}.harvest.TGCA.scn" \
      -infinder "${PREFIX}.finder.scn" -threads 10 1>ltr_retriever.log 2>&1 &

echo "Pipeline completed successfully!"