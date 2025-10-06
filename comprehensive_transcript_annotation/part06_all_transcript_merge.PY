#!/usr/bin/env python3
# Description: Identifies and restores original Iso-Seq transcripts that were merged with ssRNA-seq transcripts

import pandas as pd

# Configuration
MERGED_BED = "tama.bed"          # TAMA merged transcript file
ORIGINAL_ISOSEQ_BED = "isoseq_tama.bed"  # Original Iso-Seq transcripts
OUTPUT_REPLACE_LIST = "replace_old.list"  # List of merged IDs to replace
OUTPUT_REPLACE_BED = "replace.bed"  # Original Iso-Seq transcripts to restore

def main():
    # Load transcript sets
    merged_df = pd.read_csv(MERGED_BED, sep='\t', header=None, 
                          names=['chrom', 'start', 'end', 'name', 'score', 'strand',
                                 'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
                                 'blockSizes', 'blockStarts', 'source'])
    orig_df = pd.read_csv(ORIGINAL_ISOSEQ_BED, sep='\t', header=None,
                         names=['chrom', 'start', 'end', 'name', 'score', 'strand',
                                'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
                                'blockSizes', 'blockStarts'])

    # Identify merged Iso-Seq transcripts
    # Criteria: Transcripts from Iso-Seq source that were modified during merging
    merged_iso = merged_df[merged_df['source'].str.contains('isoseq')]
    orig_names = set(orig_df['name'])
    
    # Find transcripts that need replacement
    replace_ids = []
    for _, row in merged_iso.iterrows():
        orig_row = orig_df[orig_df['name'] == row['name']]
        if not orig_row.empty:
            # Check if structure changed significantly
            orig = orig_row.iloc[0]
            length_diff = abs((row['end'] - row['start']) - (orig['end'] - orig['start']))
            exon_diff = abs(row['blockCount'] - orig['blockCount'])
            
            # Restore if significant changes detected
            if length_diff > 50 or exon_diff > 0:
                replace_ids.append(row['name'])
    
    # Output replacement list
    with open(OUTPUT_REPLACE_LIST, 'w') as f:
        for tid in replace_ids:
            f.write(f"{tid}\n")
    
    # Output original Iso-Seq transcripts for restoration
    replace_df = orig_df[orig_df['name'].isin(replace_ids)]
    replace_df.to_csv(OUTPUT_REPLACE_BED, sep='\t', header=False, index=False)
    
    print(f"Identified {len(replace_ids)} transcripts for restoration")

if __name__ == "__main__":
    main()
