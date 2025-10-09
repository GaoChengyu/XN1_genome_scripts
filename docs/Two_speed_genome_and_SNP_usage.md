# two_speed_genome_and_SNP
## SNP_identify.sh

A comprehensive pipeline for identifying single nucleotide polymorphisms (SNPs) through comparative genomic analysis. This workflow performs reference-based alignment, variant calling, and identifies sample-specific SNPs using high-quality genomic data.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- Minimap2 - Sequence alignment
- SAMtools - BAM file processing and manipulation
- bcftools - Variant calling and VCF file operations

### Quick Start

```{}
./SNP_identify.sh <reference_fasta> <output_directory> <threads> <sample1_fasta> [sample2_fasta ...]
#<reference_fasta>: Reference genome in FASTA format
#<threads>: Number of CPU threads for parallel processing
#<sample1_fasta>: First sample genome in FASTA format
```

## two_speed_genome_identify.r

An R script for identifying and analyzing two-speed genomes using mixture models and hidden Markov models. This toolkit detects regions of differential evolutionary rates in genomes, distinguishing between slow-evolving (conserved) and fast-evolving (variable) genomic regions.

### Prerequisites

```{}
# Install required packages
install.packages(c("mixtools", "depmixS4", "fBasics", "stats4", "shape"))
```