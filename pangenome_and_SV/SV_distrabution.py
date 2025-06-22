import pandas as pd
from pybedtools import BedTool
import logging
import argparse
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Load gene annotation file (GFF3 format)
def load_gene_annotation(gff3_file):
    try:
        gff3 = pd.read_csv(gff3_file, sep='\t', comment='#', header=None,
                           names=["seqid", "source", "type", "start", "end", 
                                  "score", "strand", "phase", "attributes"])
    except pd.errors.ParserError:
        raise ValueError(f"{gff3_file} is not TAB-delimited. Please check the file format.")
    
    # Ensure start and end positions are integers
    gff3["start"] = gff3["start"].astype(int)
    gff3["end"] = gff3["end"].astype(int)
    
    # Extract gene, exon, and CDS regions - use copy() to avoid SettingWithCopyWarning
    genes = gff3[gff3["type"] == "gene"].copy()
    exons = gff3[gff3["type"] == "exon"].copy()
    cds = gff3[gff3["type"] == "CDS"].copy()
    
    # Extract promoter regions (2000 bp upstream of genes)
    # Fix SettingWithCopyWarning: use loc for assignment
    genes.loc[:, "promoter_start"] = genes.apply(
        lambda row: max(0, row["start"] - 2000) if row["strand"] == "+" else row["end"], axis=1
    )
    genes.loc[:, "promoter_end"] = genes.apply(
        lambda row: row["start"] if row["strand"] == "+" else row["end"] + 2000, axis=1
    )
    promoters = genes[["seqid", "promoter_start", "promoter_end", "strand", "attributes"]].copy()
    
    return genes, exons, cds, promoters

# Load SV file (BED format)
def load_sv_bed(sv_file):
    try:
        # Only read first 3 columns (chromosome, start, end)
        sv = pd.read_csv(sv_file, sep='\t', header=None, usecols=[0, 1, 2],
                         names=["chrom", "start", "end"])
    except pd.errors.ParserError:
        raise ValueError(f"{sv_file} is not TAB-delimited. Please check the file format.")
    
    # Ensure start and end positions are integers
    sv["start"] = sv["start"].astype(int)
    sv["end"] = sv["end"].astype(int)
    return sv

# Calculate total region length (merge overlapping regions) - add sorting functionality
def calculate_region_length(bed_df, chrom_order=None):
    """Calculate total region length (after merging overlapping regions)"""
    if bed_df.empty:
        return 0
    
    # Sort by chromosome and start position
    if chrom_order:
        # Create custom sorter
        chrom_sorter = {chrom: idx for idx, chrom in enumerate(chrom_order)}
        bed_df["chrom_order"] = bed_df["chrom"].map(chrom_sorter)
        bed_df = bed_df.sort_values(by=["chrom_order", "start"]).drop(columns=["chrom_order"])
    else:
        bed_df = bed_df.sort_values(by=["chrom", "start"])
    
    bed = BedTool.from_dataframe(bed_df)
    merged = bed.merge()
    total = 0
    for interval in merged:
        total += interval.end - interval.start
    return total

# Map SVs to gene regions (add normalization)
def map_sv_to_gene_regions(sv, genes, exons, cds, promoters, genome_size):
    # Get chromosome order from genome size file
    chrom_order = list(genome_size.keys())
    
    # Ensure all DataFrames have 'chrom' column (standardize chromosome column name)
    genes = genes.rename(columns={"seqid": "chrom"})
    exons = exons.rename(columns={"seqid": "chrom"})
    cds = cds.rename(columns={"seqid": "chrom"})
    promoters = promoters.rename(columns={"seqid": "chrom"})
    sv = sv.rename(columns={sv.columns[0]: "chrom"})
    
    # Sort SV data
    sv