# repeat_identified
## part01_earlGrey_code.sh

Using earlGrey to detect and classify TEs

### Prerequisites

- earlGrey
### Quick Start

```{}
earlGrey -g <REF_GENOME> -s <Prefix> -o <output_directory> -d yes
```
## part02_LTR_Prediction.sh

A comprehensive pipeline for discovering and annotating Long Terminal Repeat (LTR) retrotransposons in genomic sequences. This workflow integrates three specialized tools - LTR_Finder, LTR_Harvest, and LTR_retriever - to provide robust and complementary LTR element predictions.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- LTR_Finder - De novo LTR retrotransposon finder
- GenomeTools - Contains LTR_Harvest and suffix array utilities
- LTR_retriever - Integration and quality filtering tool

### Quick Start

```{}
./part02_LTR_Prediction.sh <genome.fasta>
#<genome.fasta>: Input genome sequence in FASTA format
```
## part03_TE_Annotation_Refinement.sh

A comprehensive pipeline for refining Transposable Element (TE) annotations using PacBio long-read data. This tool validates TE predictions, resolves overlapping annotations, and generates high-confidence TE annotations with comprehensive statistical summaries.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- bedtools - Genome arithmetic operations
- samtools - BAM file processing and indexing
- pandas - Data processing and analysis

### Quick Start

```{}
./part03_TE_Annotation_Refinement.sh <TE_GFF> <REF_GENOME> <PACBIO_BAM>
#<TE_GFF>: Input TE annotation file in GFF format
#<REF_GENOME>: Reference genome FASTA file
#<PACBIO_BAM>: PacBio aligned BAM file
```
