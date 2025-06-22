#!/bin/bash
# Description: Aligns transcripts to a reference genome, processes alignments, and collapses redundant transcripts using TAMA

# Configuration Section - Modify these variables as needed
THREADS=16                         # Number of threads for parallel processing
REF_GENOME="XN1.genome.chr.fasta"  # Path to reference genome
SAMPLE="js"                         # Sample identifier
INPUT_FASTA="js_clean_nofusion_noptmrna.fa"  # Input transcript file
TAMA_PATH="/mnt/d/sf/tama-master"   # Path to TAMA installation
MODE="strict"                       # Processing mode: "strict" or "lenient"

# Step 1: Alignment using Minimap2
# -ax splice: Optimized for splice alignment
# -uf: Consider only sense strand transcripts
# --secondary=no: Suppress secondary alignments
# -C5: Increased sensitivity for splice site detection
echo "Aligning transcripts to reference genome..."
minimap2 -ax splice -t $THREADS -uf --secondary=no -C5 \
    $REF_GENOME \
    $INPUT_FASTA > ${SAMPLE}_aligned.sam

# Step 2: SAM Processing and Sorting
echo "Processing and sorting alignment files..."
samtools sort -@ $THREADS \
    -o ${SAMPLE}_aligned_sorted.bam \
    ${SAMPLE}_aligned.sam
samtools view -@ $THREADS -h \
    ${SAMPLE}_aligned_sorted.bam > ${SAMPLE}_aligned_sorted.sam

# Step 3: TAMA Collapse - Remove Redundant Transcripts
echo "Collapsing redundant transcripts with TAMA..."
if [ "$MODE" == "strict" ]; then
    # Strict mode parameters (for polycistronic-filtered transcripts)
    # -e common_ends: Collapse transcripts with common ends
    # -m 20: Exon/splice junction tolerance
    # -a 99999999: 5' end tolerance
    python ${TAMA_PATH}/tama_collapse.py \
        -s ${SAMPLE}_aligned_sorted.sam \
        -f $REF_GENOME \
        -p ${SAMPLE}_tama \
        -x no_cap \
        -e common_ends \
        -c 0 \
        -i 0 \
        -m 20 \
        -a 99999999
else
    # Lenient mode parameters (for non-polycistronic-filtered transcripts)
    # -m 0: No exon/splice junction tolerance
    # -z 0: No 3' end tolerance
    # -a 1000: 5' end tolerance
    python ${TAMA_PATH}/tama_collapse.py \
        -s ${SAMPLE}_aligned_sorted.sam \
        -f $REF_GENOME \
        -p ${SAMPLE}_tama \
        -x no_cap \
        -e common_ends \
        -m 0 \
        -z 0 \
        -a 1000
fi

# Step 4: TAMA Merge - Integrate Transcripts from Multiple Samples
# (Example for merging two samples: js and condi)
echo "Merging transcript sets from multiple samples..."
cat > filelist.txt <<EOF
js_tama.bed    no_cap    1,1,1    js_sample
condi_tama.bed no_cap    2,1,1    condi_sample
EOF

# Merge parameters:
# -e common_ends: Merge transcripts with common ends
# -d merge_dup: Merge duplicate transcripts
# -a 9999: 5' end tolerance
# -m 0: No exon/splice junction tolerance
# -z 0: No 3' end tolerance
python ${TAMA_PATH}/tama_merge.py \
    -f filelist.txt \
    -p isoseq_tama \
    -e common_ends \
    -d merge_dup \
    -a 9999 \
    -m 0 \
    -z 0

echo "Processing complete!"
echo "Final collapsed transcripts: ${SAMPLE}_tama.bed"
echo "Merged transcripts: isoseq_tama.bed"