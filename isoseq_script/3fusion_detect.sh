#不同样本的二代polish后的转录本
###minimap2比对转录本到参考基因组上
minimap2 -ax splice -t 16 -uf --secondary=no -C5 ../../MC.XN21.genome.chr.fasta  condi_clean.fa > condi_clean.sam
samtools sort -@ 30 condi_clean.sam > condi_clean.sort.bam
samtools view -@ 30 -h condi_clean.sort.bam >condi_clean.sort.sam

#GitHub安装cds_Cupcake  （genome环境下）
fusion_finder.py --input condi_clean.fa -s condi_clean.sort.sam -o condi_fusion --cluster_report_csv condi.clustered.cluster_report.csv
#--cluster_report_csv 后接第一大步骤cluster生成的数据

#删除融合基因转录本
awk 'BEGIN{FS="\t";OFS="\t"}{print $2}' js_fusion.group.txt|awk '{gsub(",","\n",$0);print $0}' >del_fusion_header.txt #生成融合基因的头部文件
#利用pyfasta删除不需要的序列
pyfasta extract --header --fasta js_clean.fa --exclude --file del_fusion_header.txt > js_clean_nofusion.fa


#合并二代polish后的不同样本的转录本##这步不做
#cat condi/resultT/condi_clean.fa js/result/js_clean.fa |awk 'BEGIN{i=0}{if ($0~/>/){i+=1;$0=">transcript/"i}}1' > all.isoseq.polish.fa