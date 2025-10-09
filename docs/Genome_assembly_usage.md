# genome_assembly
## part00_bam2fasta_pipeline.sh

Integrated workflow for PacBio long-read file format transformation: bam to fasta.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- bam2fastx - transform raw data format

### Quick Start

```{}
./part00_bam2fasta_pipeline.sh <output_fasta_file> <input_bam_file>
#<output_fasta_file>: PacBio long-read FASTA file
#<input_bam_file>: PacBio long-read BAM file
```
## part01_canu_pipeline.sh

A comprehensive and automated pipeline for processing PacBio long-read sequencing data using Canu. This workflow performs error correction, quality trimming, and de novo assembly in three sequential stages to produce high-quality genome assemblies from raw PacBio reads.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- Canu (v2.0 or higher) - Long-read assembler

### Quick Start

```{}
./part01_canu_pipeline.sh <output_directory> <input_reads.fastq>
#<output_directory>: Directory where all output files will be written
#<input_reads.fastq>: Raw PacBio sequencing reads
```
## part02_polish_contig.sh

A comprehensive pipeline for iterative polishing of genome assemblies using both long-read (PacBio/Nanopore) and short-read (Illumina) sequencing data. This workflow combines Racon for long-read-based consensus improvement and Pilon for short-read-based error correction to produce high-quality polished genome assemblies.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- minimap2 - Sequence alignment
- racon - Long-read consensus polishing
- samtools - BAM/SAM processing
- sambamba (optional) - Duplicate marking (samtools markdup used as fallback)
- Pilon - Short-read based polishing
- Java - Required for Pilon execution

### Quick Start

```{}
./part02_polish_contig.sh <output_dir> <initial_assembly> <long_reads> <short_read1> <short_read2>
#<output_directory>: Directory where all output files will be written
#<initial_assembly>: Initial genome assembly FASTA file to polish
#<long_reads>: Raw long-reads file (PacBio/Nanopore)
#<short_read1>: Illumina forward reads (R1)
#<short_read2>: Illumina reverse reads (R2)
```

## part03_contig_to_chr.sh

A high-performance pipeline for chromosome-scale genome scaffolding using Hi-C data. This workflow leverages Chromap for fast Hi-C read alignment and YAHS for scaffolding, producing high-quality chromosome-length scaffolds and visualization files for downstream analysis.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- Chromap - Hi-C optimized read aligner
- SAMtools - BAM/SAM processing and indexing
- YAHS - Hi-C scaffolding tool
- Juicer Tools - Hi-C data processing and .hic file generation
- Java - Required for Juicer Tools

### Quick Start

```{}
./hic_scaffolding_pipeline.sh <contigs.fasta> <R1_reads> <R2_reads>
#<contigs.fasta>: Initial genome assembly contigs in FASTA format
#<R1_reads>: Hi-C forward reads (R1)
#<R2_reads>: Hi-C reverse reads (R2)
```
