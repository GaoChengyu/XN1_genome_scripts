###minimap2比对转录本到参考基因组上   转录本为去除多顺反子转录本


# -ax 不同的数据类型指定不同的模式，如基因组的subreads，可以用map-pb\map-ont，Iso-seq的reads用splice,二代数据reads用sr等
# -uf 对于iso-seq数据，只考虑正义链转录本
# -C5 对于识别不同物种的剪切位点，灵敏度更高
minimap2 -ax splice -t 16 -uf --secondary=no -C5 ../../MC.XN21.genome.chr.fasta  js_clean_nofusion_noptmrna.fa > js_clean_nofusion_noptmrna.sam
samtools sort -@ 30 js_clean_nofusion_noptmrna.sam >js_clean_nofusion_noptmrna.sort.bam
samtools view -@ 30 -h js_clean_nofusion_noptmrna.sort.bam >js_clean_nofusion_noptmrna.sort.sam



###tama 从github下载 https://github.com/GenomeRIK/tama/wiki/Tama-Collapse   python2环境运行 去除冗余
python /mnt/d/sf/tama-master/tama_collapse.py -s condi_clean_nofusion_noptmrna.sort.sam -f ../../MC.XN21.genome.chr.fasta -p condi_tama -x no_cap -e common_ends -c 0 -i 0 -m 20 -a 99999999
#得到一个终极文件bed


###将不同来源的bed文件整合 例如来自孢子和菌丝

python /mnt/d/sf/tama-master/tama_merge.py -f filelist.txt -p MC_tama -e common_ends -d merge_dup





####下一步打算通过ssRNAseq数据预测转录本，与isoseq转录本整合，补充





####################分界线
#####还有一种思路 因为做的是无参全长转录本 可以完全抛弃已有参考基因 所以不做多顺反子鉴定 直接用过滤了融合基因的序列进行下一步
minimap2 -ax splice -t 16 -uf --secondary=no -C5 ../../../MC.XN21.genome.chr.fasta  condi_clean_nofusion.fa > condi_clean_nofusion.sam
samtools sort -@ 30 condi_clean_nofusion.sam >condi_clean_nofusion.sort.bam
samtools view -@ 30 -h condi_clean_nofusion.sort.bam >condi_clean_nofusion.sort.sam

-c C Coverage (default 99)
Coverage is defined by the percentage of the original read that is mapped to the genome. In some cases regions at the beginnning or end of the read will be clipped from the mapping representation. The amount of this clipping is what defines the coverage value.
-i I Identity (default 85)
Identity is defined as the number of bases of the read that match the bases of the genome as per the genomic mapping. Basically it is how close the read sequence is to the genomic sequence.
-a A 5 prime threshold (default 10)
The 5 prime threshold is the amount of tolerance at the 5' end of the transcript for grouping reads to be collapsed.
-m M Exon/Splice junction threshold (default 10)
The Exon/Splice junction threshold is the amount of tolerance for the splice junctions of the transcript for grouping reads to be collapsed.

#参数选择
###保持3’和内含子剪接的零容忍！
python /mnt/d/sf/tama-master/tama_collapse.py -s condi_clean_nofusion.sort.sam -f MC.XN21.genome.chr.fasta -p condi_tama -x no_cap -e common_ends -m 0 -z 0 -a 1000



######tama_merge.py


-f filelist.txt

The filelist file contains the name of the files you want to merge as well as some additional information. The format for the file should be like this (tab separated, do not include header):

  file_name    cap_flag    merge_priority(start,junctions,end)    source_name
  annotation_capped.bed        capped  1,1,1   caplib
  annotation_nocap.bed        no_cap  2,1,1   nocaplib



python /mnt/d/sf/tama-master/tama_merge.py -f filelist.txt -p isoseq_tama -e common_ends -d merge_dup -a 9999 -m 0 -z 0
#增加-a参数 增加5'的容忍度 不增加-z和-m，3'和剪接的容忍度  防治错误折叠更多转录本  特别针对全长转录组！！！一定保留更多3‘不同的转录本！
