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
    sv_sorted = sv.sort_values(by=["chrom", "start"])
    sv_bed = BedTool.from_dataframe(sv_sorted[["chrom", "start", "end"]])
    
    # Calculate region lengths (merge overlapping regions) with sorting
    promoters_length = calculate_region_length(
        promoters[["chrom", "promoter_start", "promoter_end"]].rename(
            columns={"promoter_start": "start", "promoter_end": "end"}
        ), chrom_order
    )
    
    exons_length = calculate_region_length(
        exons[["chrom", "start", "end"]], chrom_order
    )
    
    cds_length = calculate_region_length(
        cds[["chrom", "start", "end"]], chrom_order
    )
    
    genes_length = calculate_region_length(
        genes[["chrom", "start", "end"]], chrom_order
    )
    
    introns_length = genes_length - exons_length  # Intron length = Gene length - Exon length
    
    # Calculate total covered length (union of genes and promoters)
    covered_regions = pd.concat([
        promoters[["chrom", "promoter_start", "promoter_end"]].rename(
            columns={"promoter_start": "start", "promoter_end": "end"}
        ),
        genes[["chrom", "start", "end"]]
    ])
    covered_length = calculate_region_length(covered_regions, chrom_order)
    
    # Calculate total genome length
    total_genome_length = sum(genome_size.values())
    
    # Calculate length outside gene regions
    outside_length = total_genome_length - covered_length
    
    # Map SVs to gene regions (using sorted SV data)
    genes_bed = BedTool.from_dataframe(genes[["chrom", "start", "end"]])
    exons_bed = BedTool.from_dataframe(exons[["chrom", "start", "end"]])
    cds_bed = BedTool.from_dataframe(cds[["chrom", "start", "end"]])
    introns_bed = genes_bed.subtract(exons_bed)  # Calculate intron regions
    promoters_bed = BedTool.from_dataframe(
        promoters[["chrom", "promoter_start", "promoter_end"]].rename(
            columns={"promoter_start": "start", "promoter_end": "end"}
        )
    )
    
    # Count SVs in gene regions
    sv_in_introns = sv_bed.intersect(introns_bed, u=True).count()
    sv_in_promoters = sv_bed.intersect(promoters_bed, u=True).count()
    sv_in_exons = sv_bed.intersect(exons_bed, u=True).count()
    sv_in_cds = sv_bed.intersect(cds_bed, u=True).count()
    
    # Count SVs outside gene/promoter regions
    sv_outside = len(sv) - sv_in_introns - sv_in_promoters - sv_in_exons - sv_in_cds
    
    # Normalization: Calculate SVs per million base pairs
    def normalize(count, length):
        return (count / length) * 1_000_000 if length > 0 else 0.0
    
    # Return raw counts and normalized values
    return {
        "raw_counts": {
            "outside_genes_and_promoters": sv_outside,
            "in_promoters": sv_in_promoters,
            "in_introns": sv_in_introns,
            "in_exons": sv_in_exons,
            "in_cds": sv_in_cds
        },
        "region_lengths": {
            "outside_genes_and_promoters": outside_length,
            "in_promoters": promoters_length,
            "in_introns": introns_length,
            "in_exons": exons_length,
            "in_cds": cds_length
        },
        "normalized_counts": {
            "outside_genes_and_promoters": normalize(sv_outside, outside_length),
            "in_promoters": normalize(sv_in_promoters, promoters_length),
            "in_introns": normalize(sv_in_introns, introns_length),
            "in_exons": normalize(sv_in_exons, exons_length),
            "in_cds": normalize(sv_in_cds, cds_length)
        }
    }

# Identify genes affected by nonsynonymous mutations from SVs
def find_nonsynonymous_mutations(sv, cds):
    # Standardize chromosome column names
    sv = sv.rename(columns={sv.columns[0]: "chrom"})
    cds = cds.rename(columns={"seqid": "chrom"})
    
    # Convert SVs and CDS regions to BedTool objects
    sv_bed = BedTool.from_dataframe(sv[["chrom", "start", "end"]])
    cds_bed = BedTool.from_dataframe(cds[["chrom", "start", "end", "attributes"]])
    
    # Find SVs overlapping CDS regions
    overlapping_sv = sv_bed.intersect(cds_bed, wa=True, wb=True)
    
    # Parse overlap results
    nonsynonymous_genes = set()
    for overlap in overlapping_sv:
        sv_chrom, sv_start, sv_end = overlap.fields[:3]
        cds_chrom, cds_start, cds_end, attributes = overlap.fields[3:7]
        
        # Extract gene name (in GFF3 format, gene name is in attributes field)
        gene_name = None
        for attr in attributes.split(";"):
            if attr.startswith("Parent="):
                gene_name = attr.split("=")[1]
                break
        if gene_name:
            nonsynonymous_genes.add(gene_name)
    
    return nonsynonymous_genes

# Find gene IDs with SVs in promoter regions
def find_promoter_genes(sv, promoters):
    # Standardize chromosome column names
    sv = sv.rename(columns={sv.columns[0]: "chrom"})
    promoters = promoters.rename(columns={"seqid": "chrom"})
    
    # Convert SVs and promoter regions to BedTool objects
    sv_bed = BedTool.from_dataframe(sv[["chrom", "start", "end"]])
    promoters_bed = BedTool.from_dataframe(promoters[["chrom", "promoter_start", "promoter_end", "attributes"]])
    
    # Find SVs overlapping promoter regions
    overlapping_sv = sv_bed.intersect(promoters_bed, wa=True, wb=True)
    
    # Parse results and extract gene IDs
    promoter_genes = set()
    for overlap in overlapping_sv:
        sv_chrom, sv_start, sv_end = overlap.fields[:3]
        promoter_chrom, promoter_start, promoter_end, attributes = overlap.fields[3:7]
        
        # Extract gene name (adjust based on GFF3 attributes)
        gene_name = None
        for attr in attributes.split(";"):
            if attr.startswith("ID=") or attr.startswith("Parent="):
                gene_name = attr.split("=")[1]
                break
        if gene_name:
            promoter_genes.add(gene_name)
    
    return promoter_genes

# Save results to files (updated to include normalized results)
def save_results(sv_stats, nonsynonymous_genes, promoter_genes, output_stats_file, output_genes_file, output_promoter_genes_file):
    # Save SV statistics (raw counts and normalized values)
    with open(output_stats_file, "w") as f:
        # Write raw counts
        f.write("=== Raw SV Counts ===\n")
        for region, count in sv_stats["raw_counts"].items():
            f.write(f"{region}: {count}\n")
        
        # Write region lengths
        f.write("\n=== Region Lengths (bp) ===\n")
        for region, length in sv_stats["region_lengths"].items():
            f.write(f"{region}: {length}\n")
        
        # Write normalized results (per million base pairs)
        f.write("\n=== Normalized SV Density (per Mb) ===\n")
        for region, density in sv_stats["normalized_counts"].items():
            f.write(f"{region}: {density:.4f}\n")  # Increased decimal precision
    
    # Save list of affected nonsynonymous mutation genes
    with open(output_genes_file, "w") as f:
        for gene in nonsynonymous_genes:
            f.write(f"{gene}\n")
    
    # Save list of genes with SVs in promoter regions
    with open(output_promoter_genes_file, "w") as f:
        for gene in promoter_genes:
            f.write(f"{gene}\n")

# Load genome size file
def load_genome_size(genome_size_file):
    genome_size = {}
    with open(genome_size_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                chrom = parts[0]
                try:
                    size = int(parts[1])
                    genome_size[chrom] = size
                except ValueError:
                    logging.warning(f"Skipping invalid line: {line.strip()}")
    return genome_size

# Main function (add genome size loading)
def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Analyze SV distribution in gene regions and their impact")
    parser.add_argument("--gff3", required=True, help="Gene annotation file (GFF3 format)")
    parser.add_argument("--bed", required=True, help="SV file (BED format)")
    parser.add_argument("--genome_size", required=True, help="Genome size file (TSV: chromosome<tab>size)")
    parser.add_argument("--output_stats", default="sv_stats.txt", help="SV statistics output file")
    parser.add_argument("--output_genes", default="nonsynonymous_genes.txt", help="Output file for nonsynonymous mutation genes")
    parser.add_argument("--output_promoter_genes", default="promoter_genes.txt", help="Output file for promoter-affected genes")
    args = parser.parse_args()
    
    # Load genome size
    logging.info("Loading genome size file: %s", args.genome_size)
    genome_size = load_genome_size(args.genome_size)
    
    # Load data
    logging.info("Loading gene annotation file: %s", args.gff3)
    genes, exons, cds, promoters = load_gene_annotation(args.gff3)
    
    logging.info("Loading SV file: %s", args.bed)
    sv = load_sv_bed(args.bed)
    
    # Map SVs to gene regions (with normalization)
    logging.info("Analyzing SV distribution in gene regions (with normalization)")
    sv_stats = map_sv_to_gene_regions(sv, genes, exons, cds, promoters, genome_size)
    
    # Identify genes with nonsynonymous mutations
    logging.info("Identifying genes affected by nonsynonymous mutations")
    nonsynonymous_genes = find_nonsynonymous_mutations(sv, cds)
    
    # Identify genes with SVs in promoters
    logging.info("Identifying genes with SVs in promoter regions")
    promoter_genes = find_promoter_genes(sv, promoters)
    
    # Save results
    save_results(sv_stats, nonsynonymous_genes, promoter_genes, args.output_stats, args.output_genes, args.output_promoter_genes)
    logging.info("Results saved to %s, %s, and %s", args.output_stats, args.output_genes, args.output_promoter_genes)

if __name__ == "__main__":
    main()