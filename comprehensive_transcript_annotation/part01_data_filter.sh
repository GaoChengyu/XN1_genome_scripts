#!/bin/bash
# Description: Automated pipeline for Iso-Seq data processing and alignment using SMRT Tools and GMAP

# Configuration Section - Modify these variables as needed
THREADS=30                       # Number of threads for parallel processing
SAMPLE="js"                       # Sample identifier
CONDITION="condi"                # Experimental condition identifier
REF_GENOME="XN1.fasta"         # Reference genome file name
REF_PATH="XN1.genome.chr.fasta"  # Path to reference genome
SUBSAMPLES=("YLFSJSthird" "YLFScond/all/FISO22H001113_1A")  # Subdirectories for different conditions
SUBREADS_FILES=("m64164_220319_063828.subreads.bam" "m64270e_220523_002302.subreads.bam")  # Subreads files
PRIMER="primer.fasta"             # Primer sequences file

# Step 1: CCS Processing - Generate Circular Consensus Sequences
ccs "../../RNA_seq_data/${SAMPLE}_thirdIso_Seq/X101SC21070545-Z01-J002_PacBio_Rawdata_XXXX/${SUBSAMPLES[0]}/${SUBREADS_FILES[0]}" \
    "${SAMPLE}.ccs.1.bam" \
    --min-rq=0.70 \
    --min-passes 1 \
    --noPolish \
    -j $THREADS

# Step 2: Lima Processing - Remove primers and orient sequences
lima "${SAMPLE}.ccs.1.bam" \
    $PRIMER \
    "${SAMPLE}.fl.bam" \
    -j $THREADS \
    --isoseq \
    --peek-guess

# Step 3: Refine Processing - Remove noise and polyA tails
isoseq3 refine "${SAMPLE}.fl.primer_5p--primer_3p.bam" \
    $PRIMER \
    "${SAMPLE}.flnc.bam" \
    -j $THREADS \
    --require-polya

# Step 4: Cluster Processing - Cluster FLNC reads
isoseq3 cluster "${SAMPLE}.flnc.bam" \
    "${SAMPLE}.clustered.bam" \
    --verbose \
    -j $THREADS

# Step 5: Polish Processing - Polish transcripts using subreads
isoseq3 polish -j $THREADS \
    "${CONDITION}.clustered.bam" \
    "../../RNA_seq_data/${CONDITION}_thirdIso_Seq/X101SC21070545-Z01-J003/X101SC21070545-Z01-J003_PacBio_Rawdata_XXXX/${SUBSAMPLES[1]}/${SUBREADS_FILES[1]}" \
    "${CONDITION}.polished.bam"

# Step 6: GMAP Index Preparation
gmap_build -D gmap_index \
    -d $REF_GENOME \
    $REF_PATH

# Step 7: Combine High and Low Quality Transcripts
gunzip "${CONDITION}.polished.hq.fasta.gz"
gunzip "${CONDITION}.polished.lq.fasta.gz"
cat "${CONDITION}.polished.hq.fasta" "${CONDITION}.polished.lq.fasta" > "${CONDITION}.polished.hlq.fasta"

# Step 8: GMAP Alignment
gmap -D "gmap_index/${REF_GENOME}" \
    -d $REF_GENOME \
    -f samse \
    -t $THREADS \
    -n 1 \
    --no-chimeras \
    --max-intronlength-middle=20000 \
    --max-intronlength-ends=20000 \
    --min-intronlength=20 \
    --split-large-introns \
    -z sense_force \
    "${CONDITION}.polished.hlq.fasta" > "${CONDITION}.aligned.sam"

# Step 9: SAM Processing and Sorting
samtools sort -@ $THREADS \
    -o "${CONDITION}.aligned.sorted.bam" \
    "${CONDITION}.aligned.sam"

samtools view -@ $THREADS \
    -h "${CONDITION}.aligned.sorted.bam" > "${CONDITION}.aligned.sorted.sam"