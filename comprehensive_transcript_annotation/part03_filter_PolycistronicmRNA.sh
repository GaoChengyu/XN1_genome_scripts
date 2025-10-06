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

# Step 2: Filter Non-overlapping ORFs
# Input: ${SAMPLE}_clean_nofusion.fa and TransDecoder output
# Output: uniq.orf.fa
echo "Input files: ${SAMPLE}_clean_nofusion.fa and transdecoder_dir/longest_orfs.gff3"
echo "Output file: uniq.orf.fa"

python <<EOF
from pyfasta import Fasta


def getfasta(fafile):
    fa=Fasta(fafile)
    return fa

gff="/mnt/d/data/XN21/isoseq/ptORF/condi/condi_clean_nofusion.fa.transdecoder_dir/longest_orfs.gff3"
fas="/mnt/d/data/XN21/isoseq/ptORF/condi/condi_clean_nofusion.fa"
out="/mnt/d/data/XN21/isoseq/ptORF/condi/uniq.orf.fa"
def getORF(orfgff):
    with open(orfgff,'r') as f:
        for line in f:
            yield line

def getorfinfo(orfgff):
    CDSdic={}
    for line in getORF(orfgff):
        if line!="\n":
            lin=line.strip().split("\t")
            ID_,source,type_,start,end,im1,strand,im2,feature=lin
            if strand=="+" and type_=="CDS":
                range_=(start,end)
                CDSdic.setdefault(ID_,[]).append(range_)
    return CDSdic

#筛选不重叠ORF标准 在一个转录本全部ORF中，至少存在两个转录本，与其他所有转录本都不重叠
def getintersect(orfgff):
    cdsdic=getorfinfo(orfgff)
    posdic={}
    for key in cdsdic.keys():
        sort_range=sorted(cdsdic[key],key=lambda x:x[0])

        for i in range(0,len(sort_range)):
            i_match_num=0
            for j in range(0,len(sort_range)):
                if int(sort_range[i][0])<int(sort_range[j][0]) and int(sort_range[i][1])<int(sort_range[j][0]):
                    i_match_num+=1
                elif int(sort_range[i][0])>int(sort_range[j][0]) and int(sort_range[j][1])<int(sort_range[i][0]):
                    i_match_num+=1
            if i_match_num==len(sort_range)-1:
                posdic.setdefault(key,[]).append(sort_range[i])
    posdic_filter=dict(filter(lambda x:len(x[1])>1,posdic.items()))
    return posdic_filter

def getorffa(fafile,orf):
    orfseq={}
    fa=getfasta(fafile)
    posdic=getintersect(orf)
    for k1 in posdic.keys():
        for pos in posdic[k1]:
            start=int(pos[0])-1
            end=int(pos[1])
            seq=fa[k1][start:end]
            orfseq.setdefault(k1,[]).append(seq)
    return orfseq

def writeorfseq(fafile,orf,outputfile):
    orfseq=getorffa(fafile,orf)
    with open(outputfile,'w') as w:
        for key in orfseq.keys():
            for i in range(0,len(orfseq[key])):
                w.write(f">{key}.p{i}\n")
                w.write(orfseq[key][i]+"\n")

writeorfseq(fas,gff,out)
EOF
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
