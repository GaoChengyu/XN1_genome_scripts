# comprehensive_transcript_annotation

## part01_data_filter.sh

A comprehensive pipeline for processing PacBio Iso-Seq data from raw subreads to aligned transcripts. This workflow performs full-length transcript identification, clustering, polishing, and genome alignment using SMRT Tools and GMAP.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- SMRT Link Tools (v9.0+)
  - ccs - Circular Consensus Sequencing
  - lima - Primer removal and demultiplexing
  - isoseq3 - Iso-Seq processing pipeline
- GMAP/GSNAP - Splice-aware genome alignment
- SAMtools - BAM/SAM processing

### Quick Start

```{}
./part01_data_filter.sh <input_file> <output_directory> <reference_genome>
#<input_file>: Input PacBio subreads BAM file
#<output_directory>: Directory for all output files
#<reference_genome>: Reference genome FASTA file
```


## part02_fusion_detect.sh

A specialized pipeline for detecting and filtering fusion transcripts from Iso-Seq data. This workflow identifies potential fusion events using alignment patterns and cluster information, producing a clean set of transcripts with fusion artifacts removed.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- Minimap2 - Sequence alignment
- SAMtools - BAM/SAM processing
- cds_Cupcake - Fusion detection tools
- pyfasta - FASTA file manipulation

### Quick Start:

```{}
./part02_fusion_detect.sh <input_fasta> <output_directory> <ref_genome> <cluster_report>
#<input_fasta>: Polished transcripts FASTA file from Iso-Seq processing
#<output_directory>: Directory for all output files
#<ref_genome>: Reference genome FASTA file
#<cluster_report>: Cluster report CSV from Iso-Seq processing
```

## part03_filter_PolycistronicmRNA.sh

A specialized pipeline for identifying and filtering polycistronic mRNAs from transcriptome data. This workflow detects transcripts containing multiple non-overlapping open reading frames (ORFs) that encode distinct proteins, a common feature in prokaryotic and some eukaryotic transcriptomes.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- TransDecoder - ORF prediction and annotation
- BLAST+ - Sequence alignment and database tools
- Python 3 with pyfasta - Sequence manipulation

### Quick Start

```{}
./part03_filter_PolycistronicmRNA.sh <input_fasta> <output_directory>
#<input_fasta>: Filtered transcripts FASTA file from previous pipeline steps
#<output_directory>: Directory for all output files
```

## part04_tama_del_redundancy.sh

A specialized pipeline for reducing redundancy in transcriptome assemblies using TAMA. This workflow aligns transcripts to a reference genome and collapses redundant isoforms based on splice junctions and transcript ends, producing a non-redundant set of high-confidence transcripts.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- Minimap2 - Sequence alignment
- SAMtools - BAM/SAM processing
- TAMA - Transcript collapse and merge
- Python - Required for running TAMA scripts

### Quick Start

```{}
./part04_tama_del_redundancy.sh <input_fasta> <output_directory> [mode]
#<input_fasta>: Input transcript FASTA file
#<output_directory>: Directory for all output files
#[mode]: Optional processing mode (default: strict)
#- strict: For high-confidence transcripts 
#- lenient: For standard transcripts 
```

## part05_ssRNA-seqAssemble.sh

This pipeline performs comprehensive transcriptome assembly from strand-specific RNA-seq data. It processes raw sequencing reads through alignment, transcript assembly, expression quantification, and redundancy removal to produce a high-confidence set of transcripts.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- HISAT2 - Spliced read aligner
- StringTie - Transcript assembly
- TACO - Transcript assembly merging
- featureCounts - Read counting
- SAMtools - BAM/SAM processing
- Minimap2 - Sequence alignment
- TAMA - Transcript collapsing
- seqkit - Sequence processing

### Quick Start

```{}
./part05_ssRNA-seqAssemble.sh <input_data_directory> <output_directory>
#<input_data_directory>: Directory containing RNA-seq FASTQ files
#<output_directory>: Directory where all output files will be written
```

## part06_all_transcript_merge.PY

A Python-based pipeline for identifying and restoring original Iso-Seq transcripts that were inappropriately merged with ssRNA-seq transcripts during TAMA processing. This tool ensures high-quality Iso-Seq transcript structures are preserved while maintaining the benefits of transcriptome merging.

### Prerequisites

Ensure the following tools are installed and available in your PATH

- Python

### Quick Start

```{}
python part06_all_transcript_merge.py <merged_bed> <original_isoseq_bed> <output_directory>
#<merged_bed>: TAMA merged transcript BED file
#<original_isoseq_bed>: Original Iso-Seq transcripts BED file
#<output_directory>: Directory for output files
```
