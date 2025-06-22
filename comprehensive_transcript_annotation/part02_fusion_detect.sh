#!/bin/bash
# Description: Processes polished transcripts through alignment, fusion detection, and filtering

# Configuration Section - Modify these as needed
THREADS=30                       # Number of threads for parallel processing
REF_GENOME="XN1.genome.chr.fasta"  # Reference genome path
SAMPLE="condi"                   # Sample identifier (change per sample)
CLUSTER_REPORT="condi.clustered.cluster_report.csv"  # Cluster report from Step 1

# Step 1: Alignment using Minimap2
# -ax splice: Optimize for splice alignment
# -uf: Use forward strand for spliced alignment
# --secondary=no: Suppress secondary alignments
# -C5: Output CIGAR with >5 operations
minimap2 -ax splice -t $THREADS -uf --secondary=no -C5 \
    $REF_GENOME \
    ${SAMPLE}_clean.fa > ${SAMPLE}_clean.sam

# Step 2: SAM Processing
# Sort SAM file and convert to BAM
samtools sort -@ $THREADS \
    -o ${SAMPLE}_clean.sort.bam \
    ${SAMPLE}_clean.sam

# Convert back to SAM for fusion_finder compatibility
samtools view -@ $THREADS -h \
    ${SAMPLE}_clean.sort.bam > ${SAMPLE}_clean.sort.sam

# Step 3: Fusion Transcript Detection
# Requires cds_Cupcake installed from GitHub
# --input: Polished transcripts
# -s: Sorted SAM alignment file
# --cluster_report_csv: Cluster report from initial processing
fusion_finder.py \
    --input ${SAMPLE}_clean.fa \
    -s ${SAMPLE}_clean.sort.sam \
    -o ${SAMPLE}_fusion \
    --cluster_report_csv $CLUSTER_REPORT

# Step 4: Filter Fusion Transcripts
# Extract fusion transcript headers
awk 'BEGIN{FS="\t";OFS="\t"}{print $2}' ${SAMPLE}_fusion.group.txt | 
    awk '{gsub(",","\n",$0);print $0}' > ${SAMPLE}_del_fusion_header.txt

# Remove fusion transcripts using pyfasta
pyfasta extract --header \
    --fasta ${SAMPLE}_clean.fa \
    --exclude --file ${SAMPLE}_del_fusion_header.txt \
    > ${SAMPLE}_clean_nofusion.fa

# Cleanup temporary files (optional)
rm ${SAMPLE}_clean.sam ${SAMPLE}_clean.sort.sam

echo "Processing complete for sample: $SAMPLE"
echo "Final filtered transcripts: ${SAMPLE}_clean_nofusion.fa"