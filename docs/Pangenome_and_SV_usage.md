# pangenome_and_SV
## build_pan_genome.sh

A comprehensive pipeline for constructing pangenome graphs and analyzing core and dispensable genomic regions across multiple genomes. This workflow uses Minigraph for graph construction and R for statistical analysis to identify conserved and variable genomic regions.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- Minigraph - Genome graph construction and alignment
- Python 3 - For running helper scripts
- R with tidyverse - For statistical analysis and visualization
- Helper Scripts - comb_coverage01.py (must be in script directory)

### Quick Start

```{}
./build_pan_genome.sh <reference_genome> <output_directory> <script_directory> <threads> <genome1> <genome2> ... <genomeN>
#<reference_genome>: Reference genome FASTA file
#<output_directory>: Directory for all output files
#<script_directory>: Directory containing helper Python scripts
#<threads>: Number of CPU threads for parallel processing
#<genome1> <genome2> ...: Additional genome FASTA files
```

## SV_distrabution.py

A comprehensive Python tool for analyzing the distribution and impact of structural variants (SVs) across genomic regions. This tool provides statistical analysis of SV enrichment in different functional regions and identifies genes potentially affected by regulatory or coding sequence changes.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- Python
- Required packages: pandas, pybedtools, numpy

### Quick Start

```{}
python SV_distrabution.py --gff3 <annotation.gff3> --bed <svs.bed> --genome_size <chrom_sizes.tsv>
#--gff3: Gene annotation file in GFF3 format
#--bed: Structural variants file in BED format
#--genome_size: Genome size file in TSV format
```

## SV_Validation.sh

A pipeline for validating structural variants (SVs) using pre-aligned BAM files from multiple samples. This tool classifies SVs as true positives (TP) or false positives (FP) based on read coverage evidence across multiple samples and generates comprehensive validation reports.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- bedtools - Genome arithmetic operations
- samtools - BAM file processing and indexing
- R with ggplot2 - Statistical computing and visualization

### Quick Start

```{}
./SV_Validation.sh -r <reference.fasta> -s <svs.bed> -b <bam_directory>
#-r, --reference: Reference genome FASTA file
#-s, --sv-bed: Structural variants BED file
#-b, --bam-dir: Directory containing BAM files
```

## comb_coverage01.py

The script was downloaded from https://github.com/AnimalGenomicsETH/bovine-graphs and utilized to implement the `build_pan_genome.sh` pipeline.
