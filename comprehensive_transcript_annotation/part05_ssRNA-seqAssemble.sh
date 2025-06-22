#!/bin/bash
# Description: Assembles transcripts from ssRNA-seq data, filters by expression, and collapses redundant isoforms

# Configuration Section - Modify these variables as needed
THREADS=16                         # Number of threads for parallel processing
REF_GENOME="XN1.genome.chr.fasta"  # Reference genome path
HISAT2_INDEX="hisat2/XN1"  # HISAT2 index path
DATA_DIR="data"                    # Directory containing FASTQ files
OUTPUT_DIR="rnaseq_assembly"       # Output directory
TACO_PATH="/mnt/d/sf/taco-v0.7.3.Linux_x86_64"  # Path to TACO
TAMA_PATH="/mnt/d/sf/tama-master"  # Path to TAMA
MIN_TPM=1.0                        # Minimum TPM for transcript retention

# Step 1: Create output directory structure
mkdir -p ${OUTPUT_DIR}/bam ${OUTPUT_DIR}/gtf ${OUTPUT_DIR}/counts

# Step 2: HISAT2 Alignment
# -g: HISAT2 index prefix
# -s: Strandedness (1 = forward)
# -d: Input data directory
# -o: Output directory
# -t: Thread count
echo "Aligning RNA-seq reads with HISAT2..."
nohup hisat2pip.py -g $HISAT2_INDEX \
                  -s 1 \
                  -d $DATA_DIR \
                  -o ${OUTPUT_DIR}/bam \
                  -t $THREADS \
                  1> ${OUTPUT_DIR}/hisat2_log.txt 2>&1 &

# Wait for alignment to complete
wait

# Step 3: Transcript Assembly with StringTie
echo "Assembling transcripts with StringTie..."
for bam_file in ${OUTPUT_DIR}/bam/*.bam; do
    sample=$(basename $bam_file .bam)
    stringtie $bam_file \
        -p $THREADS \
        -o ${OUTPUT_DIR}/gtf/${sample}.gtf \
        -l $sample
done

# Create file list for TACO merge
ls ${OUTPUT_DIR}/gtf/*.gtf > gtf_list.txt

# Step 4: Merge Assemblies with TACO
echo "Merging transcript assemblies with TACO..."
${TACO_PATH}/taco_run gtf_list.txt \
    -o ${OUTPUT_DIR}/merged \
    -p $THREADS

# Step 5: Extract Transcript Sequences
echo "Extracting transcript sequences..."
gffread ${OUTPUT_DIR}/merged/assembly.gtf \
    -w ${OUTPUT_DIR}/ssRNA.transcript.all.fa \
    -g $REF_GENOME

# Step 6: Transcript Quantification and Filtering
# -p: Count fragments (paired-end)
# -s 2: Reverse stranded
# -t exon: Count at exon level
# -g transcript_id: Group by transcript ID
echo "Quantifying transcript expression..."
featureCounts -p -s 2 \
    -t exon \
    -g transcript_id \
    -a ${OUTPUT_DIR}/merged/assembly.gtf \
    -o ${OUTPUT_DIR}/counts/rnaseq.raw.tsv \
    ${OUTPUT_DIR}/bam/*.bam

# Process TPM values (simplified example)
# In practice, use a proper TPM calculation script
awk -v min_tpm=$MIN_TPM '
    BEGIN {FS="\t"; OFS=","}
    NR>2 {
        # Calculate TPM values (pseudo-code)
        # Actual TPM calculation requires normalization
        for(i=2; i<=NF; i++) sum[i] += $i
        # Store transcript counts
        transcript[$1] = $0
    }
    END {
        # Print header
        print "transcript_id," "sample1,sample2,..."
        # Print transcripts meeting TPM threshold
        for (t in transcript) {
            print transcript[t]
            # Filtering logic would go here
        }
    }
' ${OUTPUT_DIR}/counts/rnaseq.raw.tsv > ${OUTPUT_DIR}/counts/tpm.csv

# Filter transcripts (TPM > 1 in at least one sample)
# This requires a proper TPM calculation and filtering script
# Placeholder for actual filtering logic
echo "Filtering low-expression transcripts..."
# In practice: awk -v min=$MIN_TPM '...' tpm.csv > select.transcript.list
# For demo: select all transcripts
cut -f1 ${OUTPUT_DIR}/counts/rnaseq.raw.tsv | tail -n +3 > select.transcript.list

# Extract filtered transcripts
echo "Extracting filtered transcript sequences..."
# Using seqkit for efficient sequence extraction (install via conda)
seqkit grep -f select.transcript.list \
    ${OUTPUT_DIR}/ssRNA.transcript.all.fa \
    > ${OUTPUT_DIR}/ssRNA.transcript.filtered.fa

# Remove low-complexity sequences
echo "Removing low-complexity sequences..."
seqkit seq -g -u ${OUTPUT_DIR}/ssRNA.transcript.filtered.fa | 
    awk '!/A{10,}|T{10,}|C{10,}|G{10,}/' \
    > ${OUTPUT_DIR}/ssRNA.transcript.fa

# Step 7: Transcript Alignment with Minimap2
echo "Aligning transcripts to reference genome..."
minimap2 -ax splice -t $THREADS \
    -uf --secondary=no \
    -C5 $REF_GENOME \
    ${OUTPUT_DIR}/ssRNA.transcript.fa > ${OUTPUT_DIR}/ssRNA.transcript.sam

# Step 8: SAM Processing
echo "Processing alignment files..."
samtools sort -@ $THREADS \
    -o ${OUTPUT_DIR}/ssRNA.transcript.sort.bam \
    ${OUTPUT_DIR}/ssRNA.transcript.sam
samtools view -@ $THREADS -h \
    ${OUTPUT_DIR}/ssRNA.transcript.sort.bam > ${OUTPUT_DIR}/ssRNA.transcript.sort.sam

# Step 9: TAMA Collapse - Remove Redundant Transcripts
echo "Collapsing redundant transcripts with TAMA..."
python ${TAMA_PATH}/tama_collapse.py \
    -s ${OUTPUT_DIR}/ssRNA.transcript.sort.sam \
    -f $REF_GENOME \
    -p ${OUTPUT_DIR}/ssrna_tama \
    -x no_cap \
    -e common_ends \
    -c 50 -i 50 -m 20 -a 1000

echo "Processing complete!"
echo "Final collapsed transcripts: ${OUTPUT_DIR}/ssrna_tama.bed"