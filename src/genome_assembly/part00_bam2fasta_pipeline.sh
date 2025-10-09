#!/usr/bin/zsh
# bam2fasta Pipeline
# Description: Integrated workflow for PacBio long-read file format transformation
# Usage: ./part00_bam2fasta_pipeline.sh <output_fasta_file> <input_bam_file>
# Requirements: bam2fastx v3.0.0

# Verify input arguments
if [[ $# -ne 2 ]]; then
    echo "Error: Invalid arguments"
    echo "Usage: $0 <output_fasta_file> <input_bam_file>"
    exit 1
fi

# Check input file existence
if [[ ! -f $2 ]]; then
    echo "Error: Input file $2 not found"
    exit 1
fi


bam2fasta -o $1 $2