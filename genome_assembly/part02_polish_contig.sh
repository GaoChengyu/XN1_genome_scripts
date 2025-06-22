#!/bin/bash

### PacBio/Nanopore Long-Read Polishing Pipeline ###
# This script performs iterative polishing of genome assemblies using:
# 1. Long-read (PacBio/Nanopore) based polishing with Racon
# 2. Short-read (Illumina) based polishing with Pilon
# Requires: minimap2, racon, samtools, sambamba, pilon

# =================================================================
# Configuration Section - Modify parameters as needed
# =================================================================
THREADS=16                      # Number of CPU threads to use
LONG_READS="long_reads.fasta"    # Raw long-reads file (FASTA/FASTQ)
SHORT_READ1="short_1_clean.fq.gz"  # Illumina R1 reads
SHORT_READ2="short_2_clean.fq.gz"  # Illumina R2 reads
INITIAL_ASSEMBLY="contig.fasta"  # Initial assembly to polish
OUTPUT_PREFIX="sample"          # Prefix for output files
PILON_MEM="52G"                 # Memory allocation for Pilon

# =================================================================
# Long-read Polishing (Racon - 3 iterations)
# =================================================================
# Note: Iterative polishing improves consensus accuracy using raw signals
current_assembly="$INITIAL_ASSEMBLY"

for i in {1..3}; do
  echo "Starting long-read polishing iteration $i"
  
  # 1. Create reference index
  minimap2 -d "ref_$i.mmi" "$current_assembly"
  
  # 2. Map long reads to assembly
  minimap2 -ax map-pb -t "$THREADS" "ref_$i.mmi" "$LONG_READS" > "align_$i.sam"
  
  # 3. Polish assembly with Racon
  racon -t "$THREADS" "$LONG_READS" "align_$i.sam" "$current_assembly" > "${OUTPUT_PREFIX}_run$i.cns.fa"
  
  # Update assembly for next iteration
  current_assembly="${OUTPUT_PREFIX}_run$i.cns.fa"
  echo "Completed iteration $i. Output: $current_assembly"
done

LONG_POLISHED="$current_assembly"  # Final long-read polished assembly

# =================================================================
# Short-read Polishing (Pilon - 3 iterations)
# =================================================================
# Note: Short-read polishing corrects small errors and improves base accuracy
current_assembly="$LONG_POLISHED"

for i in {1..3}; do
  echo "Starting short-read polishing iteration $i"
  
  # 1. Map short reads to assembly
  minimap2 -ax sr -t "$THREADS" "$current_assembly" "$SHORT_READ1" "$SHORT_READ2" \
    | samtools view -@ "$THREADS" -bS - > "ngs.polish$i.bam"
  
  # 2. Sort BAM file
  samtools sort -@ "$THREADS" "ngs.polish$i.bam" -o "ngs.polish$i.sorted.bam"
  
  # 3. Mark PCR duplicates (for non-PCR-free libraries)
  sambamba markdup -t "$THREADS" "ngs.polish$i.sorted.bam" "ngs.polish${i}_markdup.bam"
  
  # 4. Filter high-quality alignments (MAPQ ≥ 30)
  samtools view -@ "$THREADS" -b -q 30 "ngs.polish${i}_markdup.bam" > "ngs.polish${i}_filtered.bam"
  samtools index -@ "$THREADS" "ngs.polish${i}_filtered.bam"
  
  # 5. Polish with Pilon
  pilon_out="pilon_polished$i"
  java -Xmx"$PILON_MEM" -jar pilon.jar \
    --genome "$current_assembly" \
    --frags "ngs.polish${i}_filtered.bam" \
    --output "$pilon_out" \
    --fix all \
    --changes \
    --vcf \
    --threads "$THREADS" > "pilon_iter$i.log" 2>&1
  
  # Update assembly for next iteration
  current_assembly="${pilon_out}.fasta"
  echo "Completed iteration $i. Output: $current_assembly"
done

FINAL_ASSEMBLY="$current_assembly"
echo "Polishing complete! Final assembly: $FINAL_ASSEMBLY"