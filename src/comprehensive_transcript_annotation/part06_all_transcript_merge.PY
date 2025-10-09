#!/usr/bin/env python3
"""
Description: Identifies and restores original Iso-Seq transcripts that were merged with ssRNA-seq transcripts
This script compares TAMA-merged transcripts with original Iso-Seq transcripts and identifies cases where
Iso-Seq transcripts were inappropriately merged with ssRNA-seq transcripts during the TAMA merge process.
It generates a list of transcripts to be restored to their original Iso-Seq structure.

Usage: 
    python part06_all_transcript_merge.py <merged_bed> <original_isoseq_bed> <output_directory> [--length_threshold LENGTH] [--exon_threshold EXON]

Arguments:
    merged_bed: TAMA merged transcript BED file (contains both Iso-Seq and ssRNA-seq transcripts)
    original_isoseq_bed: Original Iso-Seq transcripts BED file (before merging)
    output_directory: Directory where output files will be written
    --length_threshold: Maximum allowed length difference before restoration (default: 50)
    --exon_threshold: Maximum allowed exon count difference before restoration (default: 0)
"""

import pandas as pd
import argparse
import os
import sys

def parse_arguments():
    """Parse command line arguments and validate input files."""
    parser = argparse.ArgumentParser(
        description="Identify and restore original Iso-Seq transcripts that were merged with ssRNA-seq transcripts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python part06_all_transcript_merge.py merged_transcripts.bed original_isoseq.bed /path/to/output
    python part06_all_transcript_merge.py merged.bed isoseq.bed ./results --length_threshold 100 --exon_threshold 1
        """
    )
    
    # Required arguments
    parser.add_argument(
        "merged_bed",
        help="TAMA merged transcript BED file (contains transcripts from both Iso-Seq and ssRNA-seq)"
    )
    parser.add_argument(
        "original_isoseq_bed", 
        help="Original Iso-Seq transcripts BED file (before TAMA merging)"
    )
    parser.add_argument(
        "output_directory",
        help="Directory where output files will be written"
    )
    
    # Optional arguments with default values
    parser.add_argument(
        "--length_threshold",
        type=int,
        default=50,
        help="Maximum allowed length difference (bp) between merged and original transcript before restoration (default: 50)"
    )
    parser.add_argument(
        "--exon_threshold", 
        type=int,
        default=0,
        help="Maximum allowed exon count difference between merged and original transcript before restoration (default: 0)"
    )
    
    return parser.parse_args()

def validate_input_files(merged_bed, original_isoseq_bed):
    """Validate that input files exist and are accessible."""
    if not os.path.exists(merged_bed):
        raise FileNotFoundError(f"Merged BED file not found: {merged_bed}")
    if not os.path.exists(original_isoseq_bed):
        raise FileNotFoundError(f"Original Iso-Seq BED file not found: {original_isoseq_bed}")
    
    # Check if files are readable
    try:
        with open(merged_bed, 'r') as f:
            pass
    except IOError:
        raise IOError(f"Cannot read merged BED file: {merged_bed}")
    
    try:
        with open(original_isoseq_bed, 'r') as f:
            pass
    except IOError:
        raise IOError(f"Cannot read original Iso-Seq BED file: {original_isoseq_bed}")

def create_output_directory(output_dir):
    """Create output directory if it doesn't exist."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

def load_bed_files(merged_bed_path, original_isoseq_bed_path):
    """
    Load BED files into pandas DataFrames.
    
    Parameters:
    -----------
    merged_bed_path : str
        Path to TAMA merged transcript BED file
    original_isoseq_bed_path : str
        Path to original Iso-Seq transcripts BED file
        
    Returns:
    --------
    tuple
        (merged_df, orig_df) - DataFrames containing the BED file data
    """
    print("Loading BED files...")
    
    # Define column names for BED format
    # BED format typically has 12 columns: 
    # chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
    bed_columns = [
        'chrom', 'start', 'end', 'name', 'score', 'strand',
        'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
        'blockSizes', 'blockStarts'
    ]
    
    # Load merged BED file - may contain an additional 'source' column
    try:
        merged_df = pd.read_csv(
            merged_bed_path, 
            sep='\t', 
            header=None,
            comment='#',  # Skip comment lines
            names=bed_columns + ['source']  # Add source column for TAMA merged files
        )
        print(f"  Loaded {len(merged_df)} transcripts from merged BED file")
    except pd.errors.ParserError:
        # If the source column is missing, load without it
        merged_df = pd.read_csv(
            merged_bed_path,
            sep='\t',
            header=None,
            comment='#',
            names=bed_columns
        )
        merged_df['source'] = 'unknown'  # Add default source
        print(f"  Loaded {len(merged_df)} transcripts from merged BED file (no source column found)")
    
    # Load original Iso-Seq BED file
    orig_df = pd.read_csv(
        original_isoseq_bed_path,
        sep='\t',
        header=None,
        comment='#',
        names=bed_columns
    )
    print(f"  Loaded {len(orig_df)} transcripts from original Iso-Seq BED file")
    
    return merged_df, orig_df

def identify_transcripts_for_restoration(merged_df, orig_df, length_threshold, exon_threshold):
    """
    Identify transcripts that need restoration by comparing merged and original structures.
    
    Criteria for restoration:
    1. Transcript must be from Iso-Seq source in the merged file
    2. Either:
       - Length difference > length_threshold, OR
       - Exon count difference > exon_threshold
       
    Parameters:
    -----------
    merged_df : DataFrame
        Merged transcripts from TAMA
    orig_df : DataFrame
        Original Iso-Seq transcripts
    length_threshold : int
        Maximum allowed length difference (bp)
    exon_threshold : int
        Maximum allowed exon count difference
        
    Returns:
    --------
    list
        List of transcript IDs that need restoration
    """
    print("Identifying transcripts for restoration...")
    
    # Filter for Iso-Seq transcripts in the merged set
    # TAMA typically marks the source of transcripts in the 'source' column
    merged_iso = merged_df[merged_df['source'].str.contains('isoseq', case=False, na=False)]
    print(f"  Found {len(merged_iso)} Iso-Seq transcripts in merged dataset")
    
    # Create a set of original transcript names for fast lookup
    orig_names = set(orig_df['name'])
    
    # Find transcripts that need replacement by comparing structures
    replace_ids = []
    
    for idx, row in merged_iso.iterrows():
        # Find the corresponding transcript in the original Iso-Seq data
        orig_matches = orig_df[orig_df['name'] == row['name']]
        
        if not orig_matches.empty:
            orig_row = orig_matches.iloc[0]
            
            # Calculate differences between merged and original transcript
            merged_length = row['end'] - row['start']
            orig_length = orig_row['end'] - orig_row['start']
            length_diff = abs(merged_length - orig_length)
            
            # Get exon count difference
            exon_diff = abs(row['blockCount'] - orig_row['blockCount'])
            
            # Check if transcript structure changed significantly
            # If length difference exceeds threshold OR exon count differs, mark for restoration
            if length_diff > length_threshold or exon_diff > exon_threshold:
                replace_ids.append(row['name'])
                
                if len(replace_ids) <= 5:  # Print first few for debugging
                    print(f"    Marking for restoration: {row['name']} "
                          f"(length_diff: {length_diff}bp, exon_diff: {exon_diff})")
    
    print(f"  Identified {len(replace_ids)} transcripts for restoration")
    return replace_ids

def write_output_files(replace_ids, orig_df, output_dir):
    """
    Write output files for transcript restoration.
    
    Generates:
    1. replace_old.list - List of merged transcript IDs to be replaced
    2. replace.bed - Original Iso-Seq transcripts to restore
    
    Parameters:
    -----------
    replace_ids : list
        List of transcript IDs to restore
    orig_df : DataFrame
        Original Iso-Seq transcripts DataFrame
    output_dir : str
        Output directory path
    """
    print("Writing output files...")
    
    # Define output file paths
    output_replace_list = os.path.join(output_dir, "replace_old.list")
    output_replace_bed = os.path.join(output_dir, "replace.bed")
    
    # Write list of transcripts to replace
    with open(output_replace_list, 'w') as f:
        for transcript_id in replace_ids:
            f.write(f"{transcript_id}\n")
    print(f"  Wrote {len(replace_ids)} transcript IDs to: {output_replace_list}")
    
    # Extract and write the original Iso-Seq transcripts that need restoration
    replace_df = orig_df[orig_df['name'].isin(replace_ids)]
    replace_df.to_csv(output_replace_bed, sep='\t', header=False, index=False)
    print(f"  Wrote {len(replace_df)} original transcripts to: {output_replace_bed}")
    
    return output_replace_list, output_replace_bed

def main():
    """Main function to execute the transcript restoration pipeline."""
    
    # Parse command line arguments
    args = parse_arguments()
    
    try:
        # Validate input files
        validate_input_files(args.merged_bed, args.original_isoseq_bed)
        
        # Create output directory
        create_output_directory(args.output_directory)
        
        # Load BED files
        merged_df, orig_df = load_bed_files(args.merged_bed, args.original_isoseq_bed)
        
        # Identify transcripts that need restoration
        replace_ids = identify_transcripts_for_restoration(
            merged_df, orig_df, args.length_threshold, args.exon_threshold
        )
        
        # Write output files
        output_files = write_output_files(replace_ids, orig_df, args.output_directory)
        
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()