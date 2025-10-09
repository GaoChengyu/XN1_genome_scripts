# XN1_genome_scripts

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)

> A Detailed Genomic Analysis Workflow for Haploid Fungi: A Case Study of the Apple Leaf Blotch Pathogen

## 📋 Table of Contents

- [Overview](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#-overview)
- [Features](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#-features)
- [Installation](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#-installation)
- [Usage](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#-usage)
- [Data Format](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#-data-format)
- [Project Structure](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#-project-structure)
- [Contributing](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#-contributing)
- [Citation](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#-citation)
- [License](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#-license)
- [Contact](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#%E2%80%8D%EF%B8%8F-contact)
- [Acknowledgments](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/README.md#-acknowledgments)

## 🎯 Overview
The plant pathogenic fungus *Diplocarpon coronariae* is a major causative agent of apple Marssonina leaf blotch, affecting the global apple industry safety. We isolated a purified XN1 strain of *D. coronariae* via single spore technology, assembled its genome to T2T levels, and presented a precise full-length transcriptome using PacBio sequencing and manual curation. Based on high-quality assembly and gene annotation, we discovered that LTR-RTs drove structural variations and gene evolution within and between species of *D. coronariae*. Furthermore, we assembled a young chromosome specific to *D. coronariae*, designated as Chr15.

Here, we have uploaded all the codes related to this study to this repository after sorting, aiming at knowledge sharing and improving the reproducibility of our research.

## ✨ Features
All executable scripts in this code repository are stored in the ./src directory, which is divided into multiple subdirectories with the following features:

- **genome_assembly**： Codes related to the T2T-level assembly of the XN1 genome.
- **comprehensive_transcript_annotation**： Codes related to the construction of the full-length transcriptome of the XN1 genome.
- **repeat_identified**： Codes related to the identification, classification, and filtering of repetitive sequences.
- **pangenome_and_SV**： Codes related to pan-genome construction and structural variation (SV) identification.
- **two_speed_genome_and_SNP**:  Codes related to two-speed genome and SNP identification.
- **plot**： Codes related to fugures plotting.

## 🚀 Installation

```{}
# Clone this repository
git clone https://github.com/GaoChengyu/XN1_genome_scripts.git
cd XN1_genome_scripts/src
```

## 📖 Usage
### genome_assembly
**part00_bam2fasta_pipeline.sh**

Integrated workflow for PacBio long-read file format transformation: bam to fasta.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- bam2fastx - transform raw data format

Quick Start:

```{}
./part00_bam2fasta_pipeline.sh <output_fasta_file> <input_bam_file>
#<output_fasta_file>: PacBio long-read FASTA file
#<input_bam_file>: PacBio long-read BAM file
```
**part01_canu_pipeline.sh**

A comprehensive and automated pipeline for processing PacBio long-read sequencing data using Canu. This workflow performs error correction, quality trimming, and de novo assembly in three sequential stages to produce high-quality genome assemblies from raw PacBio reads.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- Canu (v2.0 or higher) - Long-read assembler

Quick Start:

```{}
./part01_canu_pipeline.sh <output_directory> <input_reads.fastq>
#<output_directory>: Directory where all output files will be written
#<input_reads.fastq>: Raw PacBio sequencing reads
```
**part02_polish_contig.sh**

A comprehensive pipeline for iterative polishing of genome assemblies using both long-read (PacBio/Nanopore) and short-read (Illumina) sequencing data. This workflow combines Racon for long-read-based consensus improvement and Pilon for short-read-based error correction to produce high-quality polished genome assemblies.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- minimap2 - Sequence alignment
- racon - Long-read consensus polishing
- samtools - BAM/SAM processing
- sambamba (optional) - Duplicate marking (samtools markdup used as fallback)
- Pilon - Short-read based polishing
- Java - Required for Pilon execution

Quick Start:

```{}
./part02_polish_contig.sh <output_dir> <initial_assembly> <long_reads> <short_read1> <short_read2>
#<output_directory>: Directory where all output files will be written
#<initial_assembly>: Initial genome assembly FASTA file to polish
#<long_reads>: Raw long-reads file (PacBio/Nanopore)
#<short_read1>: Illumina forward reads (R1)
#<short_read2>: Illumina reverse reads (R2)
```

**part03_contig_to_chr.sh**

A high-performance pipeline for chromosome-scale genome scaffolding using Hi-C data. This workflow leverages Chromap for fast Hi-C read alignment and YAHS for scaffolding, producing high-quality chromosome-length scaffolds and visualization files for downstream analysis.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- Chromap - Hi-C optimized read aligner
- SAMtools - BAM/SAM processing and indexing
- YAHS - Hi-C scaffolding tool
- Juicer Tools - Hi-C data processing and .hic file generation
- Java - Required for Juicer Tools

Quick Start:

```{}
./part03_contig_to_chr.sh <contigs.fasta> <R1_reads> <R2_reads>
#<contigs.fasta>: Initial genome assembly contigs in FASTA format
#<R1_reads>: Hi-C forward reads (R1)
#<R2_reads>: Hi-C reverse reads (R2)
```

### comprehensive_transcript_annotation

**part01_data_filter.sh**

A comprehensive pipeline for processing PacBio Iso-Seq data from raw subreads to aligned transcripts. This workflow performs full-length transcript identification, clustering, polishing, and genome alignment using SMRT Tools and GMAP.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- SMRT Link Tools (v9.0+)
  - ccs - Circular Consensus Sequencing
  - lima - Primer removal and demultiplexing
  - isoseq3 - Iso-Seq processing pipeline
- GMAP/GSNAP - Splice-aware genome alignment
- SAMtools - BAM/SAM processing

Quick Start:

```{}
./part01_data_filter.sh <input_file> <output_directory> <reference_genome>
#<input_file>: Input PacBio subreads BAM file
#<output_directory>: Directory for all output files
#<reference_genome>: Reference genome FASTA file
```


**part02_fusion_detect.sh**

A specialized pipeline for detecting and filtering fusion transcripts from Iso-Seq data. This workflow identifies potential fusion events using alignment patterns and cluster information, producing a clean set of transcripts with fusion artifacts removed.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- Minimap2 - Sequence alignment
- SAMtools - BAM/SAM processing
- cds_Cupcake - Fusion detection tools
- pyfasta - FASTA file manipulation

Quick Start:

```{}
./part02_fusion_detect.sh <input_fasta> <output_directory> <ref_genome> <cluster_report>
#<input_fasta>: Polished transcripts FASTA file from Iso-Seq processing
#<output_directory>: Directory for all output files
#<ref_genome>: Reference genome FASTA file
#<cluster_report>: Cluster report CSV from Iso-Seq processing
```

**part03_filter_PolycistronicmRNA.sh**

A specialized pipeline for identifying and filtering polycistronic mRNAs from transcriptome data. This workflow detects transcripts containing multiple non-overlapping open reading frames (ORFs) that encode distinct proteins, a common feature in prokaryotic and some eukaryotic transcriptomes.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- TransDecoder - ORF prediction and annotation
- BLAST+ - Sequence alignment and database tools
- Python 3 with pyfasta - Sequence manipulation

Quick Start:

```{}
./part03_filter_PolycistronicmRNA.sh <input_fasta> <output_directory>
#<input_fasta>: Filtered transcripts FASTA file from previous pipeline steps
#<output_directory>: Directory for all output files
```

**part04_tama_del_redundancy.sh**

A specialized pipeline for reducing redundancy in transcriptome assemblies using TAMA. This workflow aligns transcripts to a reference genome and collapses redundant isoforms based on splice junctions and transcript ends, producing a non-redundant set of high-confidence transcripts.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- Minimap2 - Sequence alignment
- SAMtools - BAM/SAM processing
- TAMA - Transcript collapse and merge
- Python - Required for running TAMA scripts

Quick Start:

```{}
./part04_tama_del_redundancy.sh <input_fasta> <output_directory> [mode]
#<input_fasta>: Input transcript FASTA file
#<output_directory>: Directory for all output files
#[mode]: Optional processing mode (default: strict)
#- strict: For high-confidence transcripts 
#- lenient: For standard transcripts 
```

**part05_ssRNA-seqAssemble.sh**

This pipeline performs comprehensive transcriptome assembly from strand-specific RNA-seq data. It processes raw sequencing reads through alignment, transcript assembly, expression quantification, and redundancy removal to produce a high-confidence set of transcripts.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- HISAT2 - Spliced read aligner
- StringTie - Transcript assembly
- TACO - Transcript assembly merging
- featureCounts - Read counting
- SAMtools - BAM/SAM processing
- Minimap2 - Sequence alignment
- TAMA - Transcript collapsing
- seqkit - Sequence processing

Quick Start:

```{}
./part05_ssRNA-seqAssemble.sh <input_data_directory> <output_directory>
#<input_data_directory>: Directory containing RNA-seq FASTQ files
#<output_directory>: Directory where all output files will be written
```

**part06_all_transcript_merge.PY**

A Python-based pipeline for identifying and restoring original Iso-Seq transcripts that were inappropriately merged with ssRNA-seq transcripts during TAMA processing. This tool ensures high-quality Iso-Seq transcript structures are preserved while maintaining the benefits of transcriptome merging.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- Python

Quick Start:

```{}
python part06_all_transcript_merge.py <merged_bed> <original_isoseq_bed> <output_directory>
#<merged_bed>: TAMA merged transcript BED file
#<original_isoseq_bed>: Original Iso-Seq transcripts BED file
#<output_directory>: Directory for output files
```

### repeat_identified
**part01_earlGrey_code.sh**

Using earlGrey to detect and classify TEs

Prerequisites:

- earlGrey
Quick Start:

```{}
earlGrey -g <REF_GENOME> -s <Prefix> -o <output_directory> -d yes
```
**part02_LTR_Prediction.sh**

A comprehensive pipeline for discovering and annotating Long Terminal Repeat (LTR) retrotransposons in genomic sequences. This workflow integrates three specialized tools - LTR_Finder, LTR_Harvest, and LTR_retriever - to provide robust and complementary LTR element predictions.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- LTR_Finder - De novo LTR retrotransposon finder
- GenomeTools - Contains LTR_Harvest and suffix array utilities
- LTR_retriever - Integration and quality filtering tool

Quick Start:

```{}
./part02_LTR_Prediction.sh <genome.fasta>
#<genome.fasta>: Input genome sequence in FASTA format
```
**part03_TE_Annotation_Refinement.sh**

A comprehensive pipeline for refining Transposable Element (TE) annotations using PacBio long-read data. This tool validates TE predictions, resolves overlapping annotations, and generates high-confidence TE annotations with comprehensive statistical summaries.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- bedtools - Genome arithmetic operations
- samtools - BAM file processing and indexing
- pandas - Data processing and analysis

Quick Start:

```{}
./part03_TE_Annotation_Refinement.sh <TE_GFF> <REF_GENOME> <PACBIO_BAM>
#<TE_GFF>: Input TE annotation file in GFF format
#<REF_GENOME>: Reference genome FASTA file
#<PACBIO_BAM>: PacBio aligned BAM file
```

### pangenome_and_SV
**build_pan_genome.sh**

A comprehensive pipeline for constructing pangenome graphs and analyzing core and dispensable genomic regions across multiple genomes. This workflow uses Minigraph for graph construction and R for statistical analysis to identify conserved and variable genomic regions.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- Minigraph - Genome graph construction and alignment
- Python 3 - For running helper scripts
- R with tidyverse - For statistical analysis and visualization
- Helper Scripts - comb_coverage01.py (must be in script directory)

Quick Start:

```{}
./build_pan_genome.sh <reference_genome> <output_directory> <script_directory> <threads> <genome1> <genome2> ... <genomeN>
#<reference_genome>: Reference genome FASTA file
#<output_directory>: Directory for all output files
#<script_directory>: Directory containing helper Python scripts
#<threads>: Number of CPU threads for parallel processing
#<genome1> <genome2> ...: Additional genome FASTA files
```

**SV_distrabution.py**

A comprehensive Python tool for analyzing the distribution and impact of structural variants (SVs) across genomic regions. This tool provides statistical analysis of SV enrichment in different functional regions and identifies genes potentially affected by regulatory or coding sequence changes.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- Python
- Required packages: pandas, pybedtools, numpy

Quick Start:

```{}
python SV_distrabution.py --gff3 <annotation.gff3> --bed <svs.bed> --genome_size <chrom_sizes.tsv>
#--gff3: Gene annotation file in GFF3 format
#--bed: Structural variants file in BED format
#--genome_size: Genome size file in TSV format
```

**SV_Validation.sh**

A pipeline for validating structural variants (SVs) using pre-aligned BAM files from multiple samples. This tool classifies SVs as true positives (TP) or false positives (FP) based on read coverage evidence across multiple samples and generates comprehensive validation reports.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- bedtools - Genome arithmetic operations
- samtools - BAM file processing and indexing
- R with ggplot2 - Statistical computing and visualization

Quick Start:

```{}
./SV_Validation.sh -r <reference.fasta> -s <svs.bed> -b <bam_directory>
#-r, --reference: Reference genome FASTA file
#-s, --sv-bed: Structural variants BED file
#-b, --bam-dir: Directory containing BAM files
```

**comb_coverage01.py**

The script was downloaded from https://github.com/AnimalGenomicsETH/bovine-graphs and utilized to implement the `build_pan_genome.sh` pipeline.

### two_speed_genome_and_SNP
**SNP_identify.sh**

A comprehensive pipeline for identifying single nucleotide polymorphisms (SNPs) through comparative genomic analysis. This workflow performs reference-based alignment, variant calling, and identifies sample-specific SNPs using high-quality genomic data.

Prerequisites:

Ensure the following tools are installed and available in your PATH

- Minimap2 - Sequence alignment
- SAMtools - BAM file processing and manipulation
- bcftools - Variant calling and VCF file operations

Quick Start:

```{}
./SNP_identify.sh <reference_fasta> <output_directory> <threads> <sample1_fasta> [sample2_fasta ...]
#<reference_fasta>: Reference genome in FASTA format
#<threads>: Number of CPU threads for parallel processing
#<sample1_fasta>: First sample genome in FASTA format
```

**two_speed_genome_identify.r**

An R script for identifying and analyzing two-speed genomes using mixture models and hidden Markov models. This toolkit detects regions of differential evolutionary rates in genomes, distinguishing between slow-evolving (conserved) and fast-evolving (variable) genomic regions.

Prerequisites:

```{}
# Install required packages
install.packages(c("mixtools", "depmixS4", "fBasics", "stats4", "shape"))
```

### plot
For visualization scripts, we provide accompanying example datasets stored in the /example_data directory to facilitate reproducibility of the reported results.


## 📊 Data Format

For detailed descriptions of data formats and related example data, please refer to the files in the /example_data and /docs.

## 📁 Project Structure

```{}
XN1_genome_scripts/
├── 📁 example_data/           # Example file for visualization
├── 📁 src/                    # Source code
│   ├── genome_assembly/    
│   │   ├── part00_bam2fasta_pipeline.sh 
│   │   ├── part01_canu_pipeline.sh  
│   │   ├── part02_polish_contig.sh   
│   │   └── part03_contig_to_chr.sh 
│   ├── comprehensive_transcript_annotation/   
│   │   ├── part01_data_filter.sh
│   │   ├── part02_fusion_detect.sh    
│   │   ├── part03_filter_PolycistronicmRNA.sh    
│   │   ├── part04_tama_del_redundancy.sh 
│   │   ├── part05_ssRNA-seqAssemble.sh 
│   │   └── part06_all_transcript_merge.PY     
│   ├── repeat_identified/  
│   │   ├── part01_earlGrey_code.sh    
│   │   ├── part02_LTR_Prediction.sh   
│   │   ├── part03_TE_Annotation_Refinement.sh  
│   │   └── part04_Calculation_of_LTR_RT_insertion.r 
│   ├── two_speed_genome_and_SNP/  
│   │   ├── SNP_identify.sh    
│   │   └── two_speed_genome_identify.r
│   ├── pangenome_and_SV/   
│   │   ├── build_pan_genome.sh    
│   │   ├── comb_coverage01.py   
│   │   ├── SV_distrabution.py  
│   │   └── SV_Validation.sh 
│   ├── plot/  
│   │   ├── fig1c.r    
│   │   ├── fig2bc.py   
│   │   ├── fig3d.r  
│   │   ├── fig3fg.r  
│   │   ├── fig3h.r  
│   │   ├── fig4a.r
│   │   ├── fig4d.r
│   │   ├── fig5d.r  
│   │   └── fig5e.sh 
├── 📁 docs/                   # Documentation
│   ├── Genome_assembly_usage.md 
│   ├── Comprehensive_transcript_annotation_usage.md 
│   ├── Repeat_identified_usage.md 
│   ├── Pangenome_and_SV_usage.md 
│   ├── Two_speed_genome_and_SNP_usage.md 
│   └── Plot_usage.md
├── 📄 LICENSE        # MIT License
└── 📄 README.md              # This file
```

## 🤝 Contributing
We welcome contributions from the community! 

## 📝 Citation
If you use this code in your research, please cite our paper:

```{bibtex}
@inproceedings{Gao2025,
  title={High-quality Genome Assembly of Diplocarpon coronariae Unveils LTR Retrotransposon-Driven Structural Dynamics in Fungi Evolution},
  author={CY Gao, X Liu, BS Zhao, H Feng and LL Huang.},
  year={2025},
  organization={NWAFU}
}
```

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/GaoChengyu/XN1_genome_scripts/blob/main/LICENSE) file for details.

## 🙋‍♂️ Contact
For questions, issues, or suggestions:

Issues: [GitHub Issues](https://github.com/GaoChengyu/XN1_genome_scripts/issues)

Email: gaocy9611@nwafu.edu.cn


## 🙏 Acknowledgments
We thank:

The developers of all software using in this work for their excellent tools

Reviewers for their valuable feedback

Contributors who helped improve this project



<div align="center">
  
If this work is helpful, please give it a ⭐!

</div> 
