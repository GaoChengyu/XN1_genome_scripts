#!/bin/bash
# Description: Predicts ORFs, identifies non-overlapping transcripts, detects polycistronic mRNAs, and filters transcripts
# Usage: ./part03_filter_PolycistronicmRNA.sh <input_fasta> <output_directory>

# Check if correct number of arguments are provided
if [ $# -ne 2 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: ./part03_filter_PolycistronicmRNA.sh <input_fasta> <output_directory>"
    echo "  <input_fasta>: Path to the input filtered transcripts FASTA file (e.g., sample_clean_nofusion.fa)"
    echo "  <output_directory>: Path to the output directory for all generated files"
    exit 1
fi

# Assign input parameters to variables
INPUT_FASTA="$1"
OUTPUT_DIR="$2"

# Extract sample name from input filename (removes '_clean_nofusion.fa' suffix)
SAMPLE=$(basename "$INPUT_FASTA" | sed 's/_clean_nofusion.fa//')

# Configuration Section - Modify these as needed
THREADS=16                       # Number of threads for parallel processing
REF_CDS="cds.fasta"              # Reference CDS sequences file
BLAST_DB="cdsdb/cds"             # BLAST database name/path

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Change to output directory to ensure all files are created there
cd "$OUTPUT_DIR" || { echo "Error: Cannot access output directory '$OUTPUT_DIR'"; exit 1; }

# Copy input file to output directory for processing
cp "$INPUT_FASTA" "${SAMPLE}_clean_nofusion.fa"

# Step 1: ORF Prediction with TransDecoder
# Identifies long open reading frames (ORFs) in transcript sequences
# -t: Input transcript file
echo "Step 1: Predicting Open Reading Frames with TransDecoder"
TransDecoder.LongOrfs -t "${SAMPLE}_clean_nofusion.fa"

# Step 2: Filter Non-overlapping ORFs
# Identifies transcripts with multiple non-overlapping ORFs using custom Python script
# Criteria: Transcript must contain at least two ORFs that don't overlap with each other
echo "Step 2: Filtering for non-overlapping ORFs"
python <<EOF
from pyfasta import Fasta

def getfasta(fafile):
    """Load FASTA file and return Fasta object"""
    fa = Fasta(fafile)
    return fa

# Define file paths - using output directory for all files
gff_file = "${SAMPLE}_clean_nofusion.fa.transdecoder_dir/longest_orfs.gff3"
fasta_file = "${SAMPLE}_clean_nofusion.fa"
output_file = "uniq.orf.fa"

def getORF(orfgff):
    """Generator function to read ORF GFF file line by line"""
    with open(orfgff, 'r') as f:
        for line in f:
            yield line

def getorfinfo(orfgff):
    """Extract CDS information from GFF file and organize by transcript ID"""
    CDSdic = {}
    for line in getORF(orfgff):
        if line != "\n":
            lin = line.strip().split("\t")
            ID_, source, type_, start, end, im1, strand, im2, feature = lin
            if strand == "+" and type_ == "CDS":
                range_ = (start, end)
                CDSdic.setdefault(ID_, []).append(range_)
    return CDSdic

def getintersect(orfgff):
    """Identify non-overlapping ORFs within transcripts"""
    cdsdic = getorfinfo(orfgff)
    posdic = {}
    for key in cdsdic.keys():
        sort_range = sorted(cdsdic[key], key=lambda x: x[0])
        
        for i in range(0, len(sort_range)):
            i_match_num = 0
            for j in range(0, len(sort_range)):
                # Check if ORF i doesn't overlap with ORF j
                if int(sort_range[i][0]) < int(sort_range[j][0]) and int(sort_range[i][1]) < int(sort_range[j][0]):
                    i_match_num += 1
                elif int(sort_range[i][0]) > int(sort_range[j][0]) and int(sort_range[j][1]) < int(sort_range[i][0]):
                    i_match_num += 1
            # If ORF doesn't overlap with any other ORFs in the transcript
            if i_match_num == len(sort_range) - 1:
                posdic.setdefault(key, []).append(sort_range[i])
    # Filter for transcripts with at least 2 non-overlapping ORFs
    posdic_filter = dict(filter(lambda x: len(x[1]) > 1, posdic.items()))
    return posdic_filter

def getorffa(fafile, orf):
    """Extract ORF sequences from FASTA file based on GFF coordinates"""
    orfseq = {}
    fa = getfasta(fafile)
    posdic = getintersect(orf)
    for k1 in posdic.keys():
        for pos in posdic[k1]:
            start = int(pos[0]) - 1  # Convert to 0-based indexing
            end = int(pos[1])
            seq = fa[k1][start:end]
            orfseq.setdefault(k1, []).append(seq)
    return orfseq

def writeorfseq(fafile, orf, outputfile):
    """Write non-overlapping ORF sequences to output FASTA file"""
    orfseq = getorffa(fafile, orf)
    with open(outputfile, 'w') as w:
        for key in orfseq.keys():
            for i in range(0, len(orfseq[key])):
                w.write(f">{key}.p{i}\n")
                w.write(orfseq[key][i] + "\n")

writeorfseq(fasta_file, gff_file, output_file)
EOF

# Step 3: Calculate CDS Lengths
# Calculate length of each CDS in reference for coverage calculations
echo "Step 3: Calculating reference CDS lengths"
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
# Create BLAST database from reference CDS sequences
# -dbtype nucl: Nucleotide database
# -in: Input CDS sequences
echo "Step 4: Creating BLAST database from reference CDS"
makeblastdb -dbtype nucl -in $REF_CDS -out $BLAST_DB

# Step 5: BLAST Alignment
# Align filtered ORFs against reference CDS database
# -query: Filtered ORF sequences
# -outfmt 6: Tabular format for easy parsing
# -evalue: E-value threshold for significant matches
echo "Step 5: Aligning ORFs to reference CDS using BLAST"
blastn -query uniq.orf.fa \
       -out blast.txt \
       -db $BLAST_DB \
       -outfmt 6 \
       -evalue 1e-5 \
       -num_threads $THREADS

# Step 6: Polycistronic mRNA Detection
# Identify polycistronic mRNAs using custom criteria:
# - Transcript must have ≥2 non-overlapping ORFs
# - Each ORF must align to annotated protein-coding genes with >50% coverage
echo "Step 6: Detecting polycistronic mRNAs"
python <<EOF
def readblast(blastfile):
    """Generator function to read BLAST results line by line"""
    with open(blastfile, 'r') as blast:
        for line in blast:
            yield line

def readlen(lenfile):
    """Read CDS lengths into dictionary"""
    with open(lenfile, "r") as cdslen:
        cdslendic = {}
        for line in cdslen:
            line = line.strip().split()
            cdslendic[line[0]] = line[1]
    return cdslendic

def filter_cds(blastfile, lenfile):
    """Filter BLAST results by coverage (>50% of CDS length)"""
    cdslendic = readlen(lenfile)
    filter_judge = []
    for line in readblast(blastfile):
        line = line.strip().split()
        name = line[1]
        alignlen = int(line[9]) - int(line[8])
        coverage = alignlen / int(cdslendic[name])
        if coverage > 0.5:
            filter_judge.append("yes")
        else:
            filter_judge.append("no")
    return filter_judge

def del_50cds(blastfile, lenfile):
    """Filter out BLAST hits with less than 50% coverage"""
    filter_judge = filter_cds(blastfile, lenfile)
    i = 0
    for line in readblast(blastfile):
        if filter_judge[i] == "yes":
            yield line
        i += 1

def judgeorf(blastfile, lenfile):
    """Identify polycistronic transcripts based on alignment patterns"""
    dic1 = {}
    for line in del_50cds(blastfile, lenfile):
        geneid = line.strip().split("\t")[0].split(".")[0]  # Transcript ID
        transid = line.strip().split("\t")[0].split(".")[1]  # ORF number
        cdsid = line.strip().split("\t")[1]  # CDS ID
        dic1.setdefault(geneid, []).append((transid, cdsid))
    
    for key in dic1.keys():
        if len(dic1[key]) > 1:  # Only consider transcripts with multiple ORFs
            dic2 = {}
            for val in dic1[key]:
                dic2.setdefault(val[0], []).append(val[1])
            
            alllist = []
            maxlen = []
            for value in dic2.values():
                alllist += value
                maxlen.append(len(value))
            
            nowlen = len(set(alllist))  # Number of unique CDS IDs
            maxnum = max(maxlen)  # Maximum number of CDS IDs per ORF
            
            # If different ORFs align to different CDS sequences, likely polycistronic
            if nowlen > maxnum:
                yield key

def output(blastfile, lenfile, path):
    """Write polycistronic transcript IDs to output file"""
    with open(f'{path}/PolycistronicmRNA.list', 'w') as w1:
        for key in judgeorf(blastfile, lenfile):
            w1.write(key + "\n")

blast_file = "blast.txt"
len_file = "cds.length.txt"
output(blast_file, len_file, ".")
EOF

# Step 7: Filter Polycistronic Transcripts
# Separate polycistronic mRNAs from moncistronic transcripts
echo "Step 7: Filtering polycistronic transcripts"

# Remove polycistronic mRNAs from main transcript set
pyfasta extract --header \
                --fasta "${SAMPLE}_clean_nofusion.fa" \
                --exclude --file PolycistronicmRNA.list \
                > "${SAMPLE}_clean_nofusion_noptmrna.fa"

# Extract polycistronic mRNAs to separate file for further analysis
pyfasta extract --header \
                --fasta "${SAMPLE}_clean_nofusion.fa" \
                --file PolycistronicmRNA.list \
                > "${SAMPLE}_ptmrna.fa"

echo "Processing complete for sample: $SAMPLE"
echo "Input file: $INPUT_FASTA"
echo "Output directory: $OUTPUT_DIR"
echo "Final filtered transcripts: ${SAMPLE}_clean_nofusion_noptmrna.fa"
echo "Polycistronic transcripts: ${SAMPLE}_ptmrna.fa"
echo "Non-overlapping ORFs: uniq.orf.fa"
echo "Polycistronic transcript list: PolycistronicmRNA.list"
