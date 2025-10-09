#!/bin/bash
# TE_Annotation_Refinement
# Input: TE annotation GFF, PacBio aligned BAM, reference genome
# Output: High-confidence TE annotations (with class info), statistical summary
# Dependencies: bedtools, samtools

###############################################################################
### 0. Environment Setup
###############################################################################
OUT_DIR="TE_PacBio_Refinement"
mkdir -p $OUT_DIR

# Verify required tools are installed
for tool in bedtools samtools; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool is not installed. Please install it first."
        exit 1
    fi
done

###############################################################################
### 1. Parameter Configuration
###############################################################################
TE_GFF="$1"                    # Input TE annotation file (GFF format)
REF_GENOME="$2"                # Reference genome FASTA
PACBIO_BAM="$3"                # PacBio aligned BAM file
MIN_COVERAGE=10                # Minimum long-read coverage depth
MIN_TE_LENGTH=300              # Minimum TE length threshold
MAX_OVERLAP=0.3                # Maximum allowed overlap ratio

###############################################################################
### 2. Data Preprocessing (Extract class info from GFF)
###############################################################################
echo "### STEP 2: Data Preprocessing ###"

# Convert GFF to BED format while preserving class information
awk -v OFS='\t' '
    {
        start = $4 - 1;  # Convert to 0-based coordinate system
        end = $5;
        class = $3;      # Use the third column as TE class
        print $1, start, end, "TE_" NR, class;  # Generate unique TE IDs
    }
' "$TE_GFF" > "$OUT_DIR/TE_annotations.bed"

# Validate conversion success
if [ ! -s "$OUT_DIR/TE_annotations.bed" ]; then
    echo "Error: Failed to extract TE annotations"
    exit 1
fi

echo "Number of TEs extracted: $(wc -l < "$OUT_DIR/TE_annotations.bed")"

###############################################################################
### 3. TE Validation Using PacBio Data
###############################################################################
echo "### STEP 3: PacBio Validation ###"

# Verify BAM file exists
if [ ! -f "$PACBIO_BAM" ]; then
    echo "Error: PacBio BAM file not found"
    exit 1
fi

# Create BAM index if missing
if [ ! -f "${PACBIO_BAM}.bai" ]; then
    echo "Indexing BAM file..."
    samtools index "$PACBIO_BAM"
fi

# Calculate coverage depth for each TE region
echo "Calculating TE coverage depth..."
bedtools coverage -a "$OUT_DIR/TE_annotations.bed" -b "$PACBIO_BAM" -mean > "$OUT_DIR/TE_coverage.bed"

# Identify low-coverage and short TEs
awk -v min_cov=$MIN_COVERAGE -v min_len=$MIN_TE_LENGTH '
    BEGIN {OFS = "\t"} {
        te_length = $3 - $2;
        cov = $(NF);  # Last column is mean coverage
        status = "PASS";
        if (cov < min_cov) status = "LOW_COV";
        if (te_length < min_len) status = "SHORT";
        print $1, $2, $3, $4, $5, te_length, cov, status;
    }
' "$OUT_DIR/TE_coverage.bed" > "$OUT_DIR/TE_coverage_with_status.tsv"

###############################################################################
### 4. Overlapping TE Detection and Resolution
###############################################################################
echo "### STEP 4: Overlapping TE Processing ###"

# Sort TE annotations for consistent processing
bedtools sort -i "$OUT_DIR/TE_annotations.bed" > "$OUT_DIR/sorted_TEs.bed"

# Detect overlapping TE regions
echo "Detecting overlapping TEs..."
bedtools intersect -a "$OUT_DIR/sorted_TEs.bed" -b "$OUT_DIR/sorted_TEs.bed" -wo |
awk -v max_ov=$MAX_OVERLAP '
    ($1 == $6 && $4 != $9) {  # Same chromosome, different TE
        overlap = $12;
        len1 = $3 - $2;
        len2 = $8 - $7;
        ov_ratio1 = overlap / len1;
        ov_ratio2 = overlap / len2;
        if (ov_ratio1 > max_ov || ov_ratio2 > max_ov) {
            print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, overlap, ov_ratio1, ov_ratio2;
        }
    }
' > "$OUT_DIR/TE_overlaps.tsv"

# Check if any overlapping TEs were detected
if [ ! -s "$OUT_DIR/TE_overlaps.tsv" ]; then
    echo "No overlapping TEs detected - using original annotations"
    cp "$OUT_DIR/TE_annotations.bed" "$OUT_DIR/final_TE_base.bed"
else
    echo "Merging overlapping TEs..."
    python -c "
import pandas as pd
import os
from collections import defaultdict

# Load TE annotations
te_bed = pd.read_csv('$OUT_DIR/sorted_TEs.bed', sep='\t', header=None, 
                     names=['chr', 'start', 'end', 'id', 'class'])
te_bed['length'] = te_bed['end'] - te_bed['start']

# Load overlap results
df_ov = pd.read_csv('$OUT_DIR/TE_overlaps.tsv', sep='\t', header=None, 
                    names=['chr1', 'start1', 'end1', 'id1', 'class1',
                           'chr2', 'start2', 'end2', 'id2', 'class2',
                           'overlap', 'ratio1', 'ratio2'])

# Build overlap graph
graph = defaultdict(set)
for _, row in df_ov.iterrows():
    graph[row['id1']].add(row['id2'])
    graph[row['id2']].add(row['id1'])

# Find connected components
visited = set()
clusters = []
for te_id in te_bed['id']:
    if te_id not in visited:
        cluster = set()
        stack = [te_id]
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                cluster.add(node)
                stack.extend(graph.get(node, []))
        clusters.append(cluster)

# Merge each cluster
merged_regions = []
for cluster_idx, cluster in enumerate(clusters):
    cluster_df = te_bed[te_bed['id'].isin(cluster)]
    if cluster_df.empty:
        continue
    
    # Calculate merged region coordinates
    merged_start = cluster_df['start'].min()
    merged_end = cluster_df['end'].max()
    merged_length = merged_end - merged_start
    
    # Determine dominant class (length-weighted)
    class_votes = cluster_df.groupby('class')['length'].sum()
    main_class = class_votes.idxmax() if not class_votes.empty else 'Merged'
    
    # Generate unique ID
    merged_id = f'Merged_{cluster_idx + 1}'
    
    merged_regions.append({
        'chr': cluster_df['chr'].iloc[0],
        'start': merged_start,
        'end': merged_end,
        'id': merged_id,
        'class': main_class,
        'length': merged_length,
        'components': len(cluster),
        'original_ids': ','.join(cluster_df['id'].tolist()),
        'original_classes': ','.join(cluster_df['class'].unique())
    })

# Save merged results
merged_df = pd.DataFrame(merged_regions)
merged_df[['chr', 'start', 'end', 'id', 'class']].to_csv(
    '$OUT_DIR/merged_TEs.bed', sep='\t', index=False
)

# Create final TE base: unmerged TEs + merged TEs
all_merged_ids = set()
for _, row in merged_df.iterrows():
    all_merged_ids.update(row['original_ids'].split(','))

# Filter unmerged TEs
unmerged_tes = te_bed[~te_bed['id'].isin(all_merged_ids)]

# Combine unmerged and merged TEs
final_base = pd.concat([unmerged_tes[['chr', 'start', 'end', 'id', 'class']], 
                       merged_df[['chr', 'start', 'end', 'id', 'class']]])

# Save final TE base set
final_base.to_csv('$OUT_DIR/final_TE_base.bed', sep='\t', index=False, header=False)
"
fi

###############################################################################
### 5. Comprehensive Filtering and Result Generation
###############################################################################
echo "### STEP 5: Comprehensive Filtering ###"

python -c "
import pandas as pd
import os

# Load base TE annotations (post-overlap processing)
if os.path.exists('$OUT_DIR/final_TE_base.bed'):
    te_base = pd.read_csv('$OUT_DIR/final_TE_base.bed', sep='\t', 
                          header=None, names=['chr', 'start', 'end', 'id', 'class'])
else:
    te_base = pd.read_csv('$OUT_DIR/TE_annotations.bed', sep='\t', 
                          header=None, names=['chr', 'start', 'end', 'id', 'class'])

# Load coverage data
cov_data = pd.read_csv('$OUT_DIR/TE_coverage_with_status.tsv', sep='\t', 
                      names=['chr', 'start', 'end', 'id', 'class', 'te_length', 'coverage', 'cov_status'])

# Keep only essential columns for merging
cov_data = cov_data[['id', 'te_length', 'coverage', 'cov_status']]

# Merge coverage information using TE ID as key
final_df = pd.merge(te_base, cov_data, on='id', how='left')

# Handle new TEs created during merging
final_df['te_length'] = final_df['te_length'].fillna(final_df['end'] - final_df['start'])
final_df['coverage'] = final_df['coverage'].fillna(0)
final_df['cov_status'] = final_df['cov_status'].fillna('PASS')

# Apply filtering rules
final_df['status'] = 'KEEP'
final_df.loc[final_df['cov_status'] == 'LOW_COV', 'status'] = 'FILTERED'
final_df.loc[final_df['cov_status'] == 'SHORT', 'status'] = 'FILTERED'

# Save final annotations
final_df.to_csv('$OUT_DIR/final_TE_annotations.tsv', sep='\t', index=False)

# Generate high-confidence TE set
high_conf = final_df[final_df['status'] == 'KEEP']
high_conf[['chr', 'start', 'end', 'id', 'class']].to_csv(
    '$OUT_DIR/high_confidence_TEs.bed', sep='\t', index=False, header=False
)
"

###############################################################################
### 6. Statistical Summary Generation
###############################################################################
echo "### STEP 6: Generating Statistical Summary ###"

python <<EOF
import pandas as pd
import numpy as np

# Calculate reference genome size
print('Calculating genome size...')
genome_size = 0
with open('$REF_GENOME', 'r') as f:
    for line in f:
        if not line.startswith('>'):
            genome_size += len(line.strip())

# Load final TE annotations
print('Loading TE annotations...')
te_df = pd.read_csv('$OUT_DIR/final_TE_annotations.tsv', sep='\t')

# Calculate class-specific statistics
print('Calculating class statistics...')
class_stats = te_df.groupby('class').agg(
    count=('id', 'size'),
    total_length=('te_length', 'sum'),
    avg_length=('te_length', 'mean'),
    min_length=('te_length', 'min'),
    max_length=('te_length', 'max')
).reset_index()

# Calculate genome percentage
class_stats['genome_percentage'] = (class_stats['total_length'] / genome_size) * 100

# Calculate status-specific statistics
print('Calculating status statistics...')
status_stats = te_df.groupby('status').agg(
    count=('id', 'size'),
    total_length=('te_length', 'sum')
).reset_index()
status_stats['genome_percentage'] = (status_stats['total_length'] / genome_size) * 100

# Generate comprehensive report
print('Generating report...')
with open('$OUT_DIR/TE_statistics_summary.txt', 'w') as f:
    f.write('TE Annotation Statistics Report\n')
    f.write('=' * 50 + '\n\n')
    
    f.write(f'Reference Genome Size: {genome_size:,} bp\n')
    f.write(f'Total TE Elements: {len(te_df):,}\n')
    f.write(f'Total TE Sequence Length: {te_df["te_length"].sum():,} bp\n')
    f.write(f'Genome Coverage by TEs: {te_df["te_length"].sum() / genome_size * 100:.2f}%\n\n')
    
    f.write('By TE Class:\n')
    f.write('Class\tCount\tTotal Length\tAvg Length\tMin Length\tMax Length\tGenome %\n')
    for _, row in class_stats.iterrows():
        f.write(f"{row['class']}\t{row['count']}\t{row['total_length']:,}\t{row['avg_length']:.1f}\t")
        f.write(f"{row['min_length']}\t{row['max_length']}\t{row['genome_percentage']:.2f}%\n")
    
    f.write('\nBy Filtering Status:\n')
    f.write('Status\tCount\tTotal Length\tGenome %\n')
    for _, row in status_stats.iterrows():
        f.write(f"{row['status']}\t{row['count']}\t{row['total_length']:,}\t{row['genome_percentage']:.2f}%\n")
    
    # Add class distribution data
    f.write('\nClass Distribution:\n')
    for _, row in class_stats.iterrows():
        pct = row['count'] / len(te_df) * 100
        f.write(f"{row['class']}: {pct:.1f}%\n")
    
    # Add length distribution data
    f.write('\nLength Distribution:\n')
    bins = [0, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 1000000]
    te_df['length_bin'] = pd.cut(te_df['te_length'], bins=bins)
    bin_counts = te_df['length_bin'].value_counts().sort_index()
    for bin_name, count in bin_counts.items():
        f.write(f"{bin_name}: {count} elements\n")
EOF

echo "### PROCESS COMPLETED SUCCESSFULLY! ###"
echo "High-confidence TE annotations: $OUT_DIR/high_confidence_TEs.bed"
echo "Comprehensive statistical summary: $OUT_DIR/TE_statistics_summary.txt"