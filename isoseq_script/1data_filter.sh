
####步骤1：使用官方软件SMRT Link(5.1.0) 中的pbtranscript工具对测序下机初始数据subreads进行加工。
#CCS：对原始数据进行过滤
ccs ../../RNA_seq_data/js_thirdIso_Seq/X101SC21070545-Z01-J002_PacBio_Rawdata_XXXX/YLFSJSthird/m64164_220319_063828.subreads.bam js.ccs.1.bam --min-rq=0.70 --min-passes 1  --noPolish -j 30
#lima：去除引物，并对序列进行orient（5'→3'），demultiplex
lima js.ccs.1.bam primer.fasta js.fl.bam -j 16 --isoseq --peek-guess  #需要一个引物fasta文件
#refine：去除full length 的噪音。从FL reads中除去polyA尾(可选)和嵌合体（concatemers），产生FLNC转录本。即FL到FLNC
isoseq3 refine condi.fl.primer_5p--primer_3p.bam primer.fasta condi.flnc.bam -j 16 --require-polya
#cluster：聚类FLNC reads产生unpolished transcripts (FLNC to UNPOLISHED)
isoseq3 cluster js.flnc.bam js.clustered.bam --verbose -j 16
#polish：使用subreads来polish转录本 (UNPOLISHED to POLISHED)   ###v3.7.0 没有polish 搞不懂为什么 作者也没写为啥去除了 只能用3.4.0
isoseq3 polish -j 16 condi.clustered.bam /mnt/d/data/XN21/RNA_seq_data/condi_thirdIso_Seq/X101SC21070545-Z01-J003/X101SC21070545-Z01-J003_PacBio_Rawdata_XXXX/YLFScond/all/FISO22H001113_1A/m64270e_220523_002302.subreads.bam condi.polished.bam


## gmap建立基因组索引
gmap_build -D gmap_index -d McXN21.fasta ../../pacbio/MC.XN21.genome.chr.fasta
#解压高、低质量的一致性序列
gunzip condi.polished.hq.fasta.gz
gunzip condi.polished.lq.fasta.gz
cat condi.polished.hq.fasta condi.polished.lq.fasta >condi.polished.hlq.fasta
#gmap比对  使用--split-output=all参数，获得所有比对信息（单一，多比对）
gmap -D ../gmap_index/McXN21.fasta -d McXN21.fasta -f samse -t 30 -n 1  --no-chimeras --max-intronlength-middle=20000 --max-intronlength-ends=20000 --min-intronlength=20 --split-large-introns  -z sense_force ../condi.polished.hlq.fasta >condi.aligned.sam
samtools sort -@ 30 condi.aligned.sam >condi.aligned.sorted.sam
samtools view -@ 20 -h condi.aligned.sorted.bam >condi.aligned.sorted.sam



