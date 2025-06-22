#!/bin/bash
# Description: Predicts ORFs, identifies non-overlapping transcripts, detects polycistronic mRNAs, and filters transcripts

# Configuration Section - Modify these as needed
THREADS=16                       # Number of threads for parallel processing
SAMPLE="js"                       # Sample identifier
REF_CDS="cds.fasta"         # Reference CDS sequences
BLAST_DB="cdsdb/cds"             # BLAST database name/path

# Step 1: ORF Prediction with TransDecoder
# -t: Input transcript file
TransDecoder.LongOrfs -t ${SAMPLE}_clean_nofusion.fa

# Step 2: Filter Non-overlapping ORFs (Placeholder for custom script)
# This step requires the custom getorf.ipynb script
# Input: ${SAMPLE}_clean_nofusion.fa and TransDecoder output
# Output: uniq.orf.fa
echo "Run getorf.ipynb manually:"
echo "Input files: ${SAMPLE}_clean_nofusion.fa and transdecoder_dir/longest_orfs.gff3"
echo "Output file: uniq.orf.fa"

# Step 3: Calculate CDS Lengths
# Processes CDS fasta to calculate sequence lengths
awk 'BEGIN{FS="\t";OFS="\t"} {
    if ($0~/^>/) {
        id = $0
        next
    } else {
        len[id] += length($0)
    }
} END {
    for (i in len) {
        print substr(i,2), len[i]
    }
}' $REF_CDS > cds.length.txt

# Step 4: BLAST Database Creation
# -dbtype nucl: Nucleotide database
# -in: Input CDS sequences
makeblastdb -dbtype nucl -in $REF_CDS -out $BLAST_DB

# Step 5: BLAST Alignment
# -query: Filtered ORF sequences
# -outfmt 6: Tabular format
# -evalue: Significance threshold
blastn -query uniq.orf.fa \
       -out blast.txt \
       -db $BLAST_DB \
       -outfmt 6 \
       -evalue 1e-5 \
       -num_threads $THREADS

# Step 6: Polycistronic mRNA Detection (Placeholder for custom script)
# This step requires the custom judgeptmRNA.ipynb script
# Input: blast.txt and cds.length.txt
# Output: PolycistronicmRNA.list
echo "Run judgeptmRNA.ipynb manually:"
echo "Input files: blast.txt and cds.length.txt"
echo "Output file: PolycistronicmRNA.list"

# Step 7: Filter Polycistronic Transcripts
# Remove polycistronic mRNAs from main set
pyfasta extract --header \
                --fasta ${SAMPLE}_clean_nofusion.fa \
                --exclude --file PolycistronicmRNA.list \
                > ${SAMPLE}_clean_nofusion_noptmrna.fa

# Extract polycistronic mRNAs to separate file
pyfasta extract --header \
                --fasta ${SAMPLE}_clean_nofusion.fa \
                --file PolycistronicmRNA.list \
                > ${SAMPLE}_ptmrna.fa

echo "Processing complete for sample: $SAMPLE"
echo "Final filtered transcripts: ${SAMPLE}_clean_nofusion_noptmrna.fa"
echo "Polycistronic transcripts: ${SAMPLE}_ptmrna.fa"