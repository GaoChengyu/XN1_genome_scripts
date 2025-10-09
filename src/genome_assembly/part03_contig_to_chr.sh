#!/bin/bash
set -e  # Exit immediately on any error

# Check required arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <contigs.fasta> <R1_reads> <R2_reads>"
    exit 1
fi

# Configuration parameters
CONTIGS_FASTA="$1"
R1_READS="$2"
R2_READS="$3"
THREADS=50  # Adjust according to your system resources
PREFIX="aligned"

# 1. Indexing phase
echo "[$(date)] Building indexes..."
# Create FASTA index for downstream tools
samtools faidx "$CONTIGS_FASTA"

# Build chromap index for faster alignment
chromap -i -r "$CONTIGS_FASTA" -o contigs.index

# 2. Alignment using chromap (Hi-C optimized)
echo "[$(date)] Performing Hi-C alignment..."
chromap \
    --preset hic \
    -r "$CONTIGS_FASTA" \
    -x contigs.index \
    --remove-pcr-duplicates \
    -1 "$R1_READS" \
    -2 "$R2_READS" \
    --SAM \
    -o "${PREFIX}.sam" \
    -t "$THREADS"

# 3. SAM to sorted BAM conversion
echo "[$(date)] Converting SAM to sorted BAM..."
samtools view -bh "${PREFIX}.sam" | \
    samtools sort -@ "$THREADS" -n - > "${PREFIX}.bam"

# Clean up intermediate files
rm "${PREFIX}.sam"

# 4. Scaffolding with YAHS
echo "[$(date)] Performing scaffolding with YAHS..."
yahs "$CONTIGS_FASTA" "${PREFIX}.bam"

# 5. Prepare files for Juicer
echo "[$(date)] Preparing Juicer input..."
juicer pre -a -o out_JBAT \
    yahs.out.bin \
    yahs.out_scaffolds_final.agp \
    "${CONTIGS_FASTA}.fai" > out_JBAT.log 2>&1

# 6. Generate .hic file using juicer_tools
echo "[$(date)] Generating .hic file..."
# Extract chromosome sizes from log for hic conversion
CHROM_SIZES=$(grep 'PRE_C_SIZE' out_JBAT.log | awk '{print $2 " " $3}')

# Create hic file
java -Xmx32G -jar /mnt/d/sf/juicer_tools_1.19.02.jar pre \
    out_JBAT.txt out_JBAT.hic.part \
    <(echo "$CHROM_SIZES") && \
    mv out_JBAT.hic.part out_JBAT.hic

echo "[$(date)] Pipeline completed successfully!"