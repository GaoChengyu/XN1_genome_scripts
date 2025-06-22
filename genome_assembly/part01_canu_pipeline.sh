#!/usr/bin/zsh
# Canu Pipeline: Correction, Trimming and Assembly
# Description: Integrated workflow for PacBio long-read processing using Canu
# Usage: ./part01_canu_pipeline.sh <output_directory> <input_reads.fastq>
# Requirements: Canu v2.0+

###############################################################################
### PARAMETER VALIDATION
###############################################################################
# Verify input arguments
if [[ $# -ne 2 ]]; then
    echo "Error: Invalid arguments"
    echo "Usage: $0 <output_directory> <input_reads.fastq>"
    exit 1
fi

# Check input file existence
if [[ ! -f $2 ]]; then
    echo "Error: Input file $2 not found"
    exit 1
fi

###############################################################################
### CONFIGURATION
###############################################################################
# Global parameters
REF_GENOME_SIZE="52m"         # Estimated genome size (correction/trimming steps)
FINAL_GENOME_SIZE="54m"       # Adjusted genome size (assembly step)
MIN_READ_LENGTH=2000          # Minimum read length to consider (bp)
MIN_OVERLAP=500               # Minimum overlap length (bp)

# Resource allocation
CORRECTION_THREADS=4          # CPUs for correction step
CORRECTION_MEM="14"           # GB RAM for correction
ASSEMBLY_THREADS=8            # CPUs for trim/assembly
ASSEMBLY_MEM="32"             # GB RAM for trim/assembly

# Disable memory-intensive features
export CANU_DISABLE_GNUPLOT=1 # Bypass gnuplot dependency

###############################################################################
### 1. READ CORRECTION STAGE
###############################################################################
# Objective: Improve raw read accuracy using overlapping read consensus
# Input: Raw PacBio reads
# Output: Error-corrected reads in ${1}/correction directory
echo "\n=== STARTING READ CORRECTION ==="
canu -correct \
    -p "canu_corrected" \          # File prefix for outputs
    -d "${1}/correction" \          # Output subdirectory
    corThreads=${CORRECTION_THREADS} \
    corMemory=${CORRECTION_MEM} \
    genomeSize=${REF_GENOME_SIZE} \
    minReadLength=${MIN_READ_LENGTH} \
    minOverlapLength=${MIN_OVERLAP} \
    mhapPipe=false \               # Disable MHAP for smaller genomes
    purgeOverlaps=false \          # Retain all overlaps
    saveOverlaps=true \            # Preserve overlap data for downstream
    corOutCoverage=200 \           # Target correction coverage depth
    corMinCoverage=2 \             # Minimum coverage for correction
    -pacbio "${2}"                 # Input raw PacBio reads

# Check correction success
if [[ ! -f "${1}/correction/canu_corrected.correctedReads.fasta.gz" ]]; then
    echo "Error: Correction step failed"
    exit 2
fi

###############################################################################
### 2. READ TRIMMING STAGE
###############################################################################
# Objective: Remove low-quality regions from corrected reads
# Input: Corrected reads from stage 1
# Output: Trimmed reads in ${1}/trimming directory
echo "\n=== STARTING READ TRIMMING ==="
canu -trim \
    -p "canu_trimmed" \            # File prefix for outputs
    -d "${1}/trimming" \            # Output subdirectory
    corThreads=${ASSEMBLY_THREADS} \
    corMemory=${ASSEMBLY_MEM} \
    genomeSize=${REF_GENOME_SIZE} \
    minReadLength=${MIN_READ_LENGTH} \
    minOverlapLength=${MIN_OVERLAP} \
    -corrected \                   # Indicate input is corrected reads
    -pacbio "${1}/correction/canu_corrected.correctedReads.fasta.gz"

# Check trimming success
if [[ ! -f "${1}/trimming/canu_trimmed.trimmedReads.fasta.gz" ]]; then
    echo "Error: Trimming step failed"
    exit 3
fi

###############################################################################
### 3. ASSEMBLY STAGE
###############################################################################
# Objective: Construct final genome assembly
# Input: Trimmed corrected reads from stage 2
# Output: Assembled contigs in ${1}/assembly directory
echo "\n=== STARTING GENOME ASSEMBLY ==="
canu -assemble \
    -p "canu_assembly" \           # File prefix for outputs
    -d "${1}/assembly" \            # Output subdirectory
    corThreads=${ASSEMBLY_THREADS} \
    corMemory=${ASSEMBLY_MEM} \
    genomeSize=${FINAL_GENOME_SIZE} \
    correctedErrorRate=0.055 \     # Adjusted error rate for corrected reads
    -corrected \                   # Indicate input is corrected reads
    -pacbio "${1}/trimming/canu_trimmed.trimmedReads.fasta.gz"

# Verify assembly output
if [[ ! -f "${1}/assembly/canu_assembly.contigs.fasta" ]]; then
    echo "Error: Assembly step failed"
    exit 4
fi

echo "\n=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "Final assembly: ${1}/assembly/canu_assembly.contigs.fasta"