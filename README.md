# XN1_genome_scripts

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

> A Detailed Genomic Analysis Workflow for Haploid Fungi: A Case Study of the Apple Leaf Blotch Pathogen

## 📋 Table of Contents

- [Overview](#🎯 Overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Data Format](#data-format)
- [Configuration](#configuration)
- [Project Structure](#project-structure)
- [Reproducing Results](#reproducing-results)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## 🎯 Overview
The plant pathogenic fungus *Diplocarpon coronariae* is a major causative agent of apple Marssonina leaf blotch, affecting the global apple industry safety. We isolated a purified XN1 strain of *D. coronariae* via single spore technology, assembled its genome to T2T levels, and presented a precise full-length transcriptome using PacBio sequencing and manual curation. Based on high-quality assembly and gene annotation, we discovered that LTR-RTs drove structural variations and gene evolution within and between species of *D. coronariae*. Furthermore, we assembled a young chromosome specific to *D. coronariae*, designated as Chr15.

Here, we have uploaded all the codes related to this study to this repository after sorting, aiming at knowledge sharing and improving the reproducibility of our research.

## Content
- **genome_assembly**： Codes related to the T2T-level assembly of the XN1 genome.
- **comprehensive_transcript_annotation**： Codes related to the construction of the full-length transcriptome of the XN1 genome.
- **repeat_identified**： Codes related to the identification, classification, and filtering of repetitive sequences.
- **pangenome_and_SV**： Codes related to pan-genome construction and structural variation (SV) identification.
- **two_speed_genome_and_SNP**:  Codes related to two-speed genome and SNP identification.
- **plot**： Codes related to fugures plotting.

## 🚀 Installation

## 📖 Usage

## 📊 Data Format

## 📁 Project Structure

```{}
XN1_genome_scripts/
├── 📁 configs/                 # Configuration files
│   ├── main.yaml              # Main configuration
│   ├── train_baseline.yaml    # Training config
│   └── hparam_space.yaml      # Hyperparameter search space
├── 📁 data/                   # Data directory
│   ├── raw/                   # Raw data (immutable)
│   ├── processed/             # Processed data
│   └── external/              # External data sources
├── 📁 src/                    # Source code
│   ├── data/                  # Data processing modules
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
├── 📄 setup.py               # Package installation
├── 📄 pyproject.toml         # Modern package config
├── 📄 Makefile               # Automation commands
├── 📄 Dockerfile             # Container configuration
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

Issues: GitHub Issues

Email: gaocy9611@nwafu.edu.cn


## 🙏 Acknowledgments
We thank:

The developers of [Library/Framework] for their excellent tools

[Dataset providers] for making their data available

Reviewers for their valuable feedback

Contributors who helped improve this project



<div align="center">
If this work is helpful, please give it a ⭐!

</div> 
