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
All executable scripts in this code repository are stored in the ./scr directory, which is divided into multiple subdirectories with the following features:

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
cd XN1_genome_scripts/scr
```

## 📖 Usage
### genome_assembly
**part00_bam2fasta_pipeline.sh**

Integrated workflow for PacBio long-read file format transformation: bam to fasta.

Prerequisites:
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
- Canu (v2.0 or higher) - Long-read assembler

Quick Start:

```{}
./part01_canu_pipeline.sh <output_directory> <input_reads.fastq>
#<output_directory>: Directory where all output files will be written
#<input_reads.fastq>: Raw PacBio sequencing reads
```

### comprehensive_transcript_annotation
**part05_ssRNA-seqAssemble.sh**

This pipeline performs comprehensive transcriptome assembly from strand-specific RNA-seq data. It processes raw sequencing reads through alignment, transcript assembly, expression quantification, and redundancy removal to produce a high-confidence set of transcripts.

**part05_ssRNA-seqAssemble.sh**

This pipeline performs comprehensive transcriptome assembly from strand-specific RNA-seq data. It processes raw sequencing reads through alignment, transcript assembly, expression quantification, and redundancy removal to produce a high-confidence set of transcripts.

Prerequisites:
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
**part05_ssRNA-seqAssemble.sh**

This pipeline performs comprehensive transcriptome assembly from strand-specific RNA-seq data. It processes raw sequencing reads through alignment, transcript assembly, expression quantification, and redundancy removal to produce a high-confidence set of transcripts.

## 📊 Data Format

For detailed descriptions of data formats and related example data, please refer to the files in the example_data directory.

## 📁 Project Structure

```{}
XN1_genome_scripts/
├── 📁 configs/                 # Configuration files
│   ├── main.yaml              # Main configuration
│   ├── train_baseline.yaml    # Training config
│   └── hparam_space.yaml      # Hyperparameter search space
├── 📁 example_data/           # Example Data directory
│   ├── raw/                   # Raw data (immutable)
│   ├── processed/             # Processed data
│   └── external/              # External data sources
├── 📁 src/                    # Source code
│   ├── genome_assembly/        # Data processing modules
│   │   ├── __init__.py
│   │   ├── loaders.py         # Data loading utilities
│   │   ├── preprocessors.py   # Data preprocessing
│   │   └── transformers.py    # Feature transformations
│   ├── models/                # Model architectures
│   │   ├── __init__.py
│   │   ├── base_model.py      # Base model class
│   │   ├── your_model.py      # Your main model
│   │   └── layers.py          # Custom layers
│   ├── training/              # Training utilities
│   │   ├── trainers.py        # Training loops
│   │   ├── callbacks.py       # Training callbacks
│   │   └── optimizers.py      # Optimizer configurations
│   ├── evaluation/            # Evaluation metrics
│   │   ├── metrics.py         # Custom metrics
│   │   ├── analyzers.py       # Result analysis
│   │   └── reporters.py       # Report generation
│   └── utils/                 # Utility functions
│       ├── logging.py         # Logging configuration
│       ├── visualization.py   # Plotting utilities
│       └── helpers.py         # Helper functions
├── 📁 notebooks/              # Jupyter notebooks
│   ├── 01_data_exploration.ipynb
│   ├── 02_model_prototyping.ipynb
│   └── 03_result_analysis.ipynb
├── 📁 scripts/                # Utility scripts
│   ├── download_data.py       # Data download script
│   ├── setup_environment.py   # Environment setup
│   └── run_experiments.py     # Batch experiment runner
├── 📁 tests/                  # Test suite
│   ├── test_data.py
│   ├── test_models.py
│   └── test_integration.py
├── 📁 docs/                   # Documentation
│   ├── api.md                 # API reference
│   ├── tutorials.md           # Tutorials and examples
│   └── troubleshooting.md     # Common issues and solutions
├── 📄 environment.yml         # Conda environment
├── 📄 requirements.txt        # Pip requirements
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
