#pacbio下机数据bam
#利用bam2fastx
bam2fasta -o subreads.fasta subreads.bam

#canu 2.2组装rawdata
#correct
#!/usr/bin/zsh
canu -correct \
    -p ath -d $1 corThreads=4 corMemory=14 gnuplot=undef \
    genomeSize=42m minReadLength=2000 minOverlapLength=500 \
    mhapPipe=false purgeOverlaps=false saveOverlaps=true \
    corOutCoverage=200 corMinCoverage=2 \
    -pacbio $2

#trim
canu -trim \
	-p fp -d $1 corThreads=8 corMemory=32 gnuplot=undef \
	genomeSize=42m minReadLength=2000 minOverlapLength=500 \
	-corrected -pacbio $2

#assemble
#!/usr/bin/zsh
canu -assemble \
	-p assemble -d $1 corThreads=8 corMemory=32 gnuplot=undef \
	genomeSize=38m correctedErrorRate=0.055 \
	-corrected -pacbio $2

#剔除组装后基因组序列中的bubbles
seqkit grep -nrp "suggestBubble=no" assemble.contigs.fasta >no_bubble.contig.fasta

#用minimap2和racon 用自身数据进行polish

#先建立索引，在比对。
x ：非常中要的一个选项，软件预测的一些值，针对不同的数据选择不同的值
map-pb/map-ont: pb或者ont数据与参考序列比对:
ava-pb/ava-ont: 寻找pd数据或者ont数据之间的overlap关系；
asm5/asm10/asm20: 拼接结果与参考序列进行比对，适合~0.1/1/5% 序列分歧度；
splice: 长reads的切割比对
sr: 短reads比对
-d :创建索引文件名
-a ：指定输出格式为sa格式，默认为PAF
-Q ：sam文件中不输出碱基质量
-R ：reads Group信息，与bwa比对中的-R一致
-t：线程数，默认为3
####

#run1
minimap2 -d contig.mmi improved3.fasta       # indexing
minimap2 -ax map-pb -t 8  contig.mmi raw_reads.fasta > alin_1.sam
racon -t 8 raw_reads.fasta alin_1.sam improved3.fasta > sample_run1.cns.fa

# run2
minimap2 -d ref_2.mmi sample_run1.cns.fa
minimap2 -ax map-pb -t 8  ref_2.mmi  raw_reads.fasta > alin_2.sam
racon -t 8  raw_reads.fasta alin_2.sam  sample_run1.cns.fa  >  sample_run2.cns.fa

# run3
minimap2 -d ref_3.mmi sample_run2.cns.fa
minimap2 -ax map-pb -t 8  ref_3.mmi  raw_reads.fasta > alin_3.sam
racon -t 8 raw_reads.fasta alin_3.sam sample_run2.cns.fa  >  sample_run3.cns.fa

###用pbmm gcpp三代数据进行polish （python2）
1. samtools faidx assembly.fa
2. pbmm2 align pacbio.subreads.bam assembly.fa contigs.fasta.bam
3. samtools sort -@ 32 contigs.fasta.bam -o contigs.fasta.sort.bam
4. samtools index -@ 16 contigs.fasta.sort.bam
5. gcpp -j32 --algorithm=arrow contigs.fasta.sorted.bam -r contig.fasta -o myConsensus.fasta

#用二代数据进行polish
#第一步：比对
minimap2 -ax sr -t 4 contig_polish.fa ../hph121_1_clean.fq.gz ../hph121_2_clean.fq.gz |samtools view -@ 8 -bS> ngs.polish1.bam

samtools sort -@ 8 ngs.polish1.bam -o ngs.polish1.sorted.bam
#samtools index -@ 8 ngs.polish1.sorted.bam

#第二步：标记重复（非PCR-free建库)

sambamba去除pcr重复，快速去除PCR重复，较samtools，picard更快


sambamba markdup -t 4 ngs.polish1.sorted.bam ngs.polish1_markdup.bam

#第三步 高质量序列过滤
samtools view -h -@ 10 -q 30 ngs.polish1_markdup.bam |samtools view -b >ngs.polish1_filter.bam
samtools index -@ 8 ngs.polish1_filter.bam
#第四步  pilon纠错（返回第一步迭代三次）
nohup java -Xmx52G -jar /home/gaocy/miniconda3/envs/genome/share/pilon-1.24-0/pilon.jar --genome contig_polish.fa --frags ngs.polish1_filter.bam --fix all --output pilon_polished --vcf 1> pilon.log 2>&1 &

####按从大到小对contig排序
使用自编脚本  切记！！！！！！！！！ ################输出一个基因组完全组装后的fasta文件


###############
##考虑全手动安装
先采用 RepeatModeler v. 1.0.11（http://www.repeatmasker.org）对基因组中的重复序列进行从头预测
###Create a Database for RepeatModeler
BuildDatabase -name NC_repeat Nc_2polish.fasta

###Run RepeatModeler
nohup RepeatModeler -database NC_repeat -threads 8 1>repeatmodeler.log 2>&1 &
#选择consensi.fa.classified进行下一步的合并
###合并预测重复序列与重复序列数据库
cat secies.repeat.fasta RMRBSeqs.fasta > all.repeat.fasta

###利用RepeatMasker v. 4.0.9 (Tarailo‐Graovac and Chen 2009)软件通过与重复系列数据库进行比对鉴定各基因组中的重复序列
nohup RepeatMasker -e rmblast -lib all.repeat.fasta -pa 16 -xsmall -dir dirname(##输出文件储存路径) genome.fasta  1 >repetmasker.log 2>&1 &       #-xsmall 软屏蔽 软屏蔽有利于之后的基因预测  降低假阳性
###这几部出现报错可能是rmblast版本冲突 删除conda版本 手动安装

############
ltr转座子鉴定

#######LTRfinder v. 1.5.9 (Ellinghaus et al. 2008)软件进行预测，https://github.com/xzhub/LTR_Finder下载，命令行分别为：
ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.85 ${GENOME} > ${i}.finder.scn
###ltr_finder速度很慢 可以用多线程版本LTR_FINDER_parallel   github下载
perl /mnt/d/sf/LTR_FINDER_parallel-1.1/LTR_FINDER_parallel -seq MC.XN21.genome.chr.gapless.fasta -threads 16 -harvest_out

##
gt suffixerator -db ${genome.fasta} -indexname ${genomename} -tis -suf -lcp -des -ssp -sds -dna
nohup gt ltrharvest -index ${genomename} -similar 85 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 > mc.harvest.TGCA.scn 1>ltrharvest.log 2>&1 &

#####ltr_retriever整合两种ltr预测结果
nohup LTR_retriever -genome Marssonina_coronaria.genome.fasta -inharvest mc.harvest.TGCA.scn -infinder Mc.finder.scn -threads 10 1>ltr_retriever.log 2>&1 &



#####转录组分析
###自编脚本hisat2pip.py对转录组数据回帖


#####后续基因模型预测时需要拼接好的转录本和内含子位置信息，因此，采用 stringTie(Kovaka et al. 2019)软件对 RNA-seq 数据进行转录本拼接

自编脚本stringtie.sh #在组装的转录本中，也会给出定量的结果，对于组装的新转录本和基因，默认采用STRG加数字编号进行区分

stringtie --merge -o merge.assembly.gtf -p 20 *.gtf  #单个样本组装完成后，会合并所有样本的转录本组装结果，得到一个非冗余的转录本集合


####利用STAR鉴定内含子

###构建索引

#--runMode genomeGenerate option directs STAR to run genome indices generation job  
#default alignReads
STAR --runMode genomeGenerate --runThreadN 2 --genomeDir STAR --genomeFastaFiles Marssonina_coronaria.genome.fasta --sjdbGTFfile merge.assembly.gtf

#比对 预测 
自编脚本 STAR.sh  跑完以后可以得到一个SJ.out.tab, 该文件为剪切位点
#####获得introns.gff   
方法1  用awk
awk 'BEGIN{OFS="\t"} $7 >= 2{if($4==1){st="+"}else{st="-"} print $1,"STAR","intron",$2,$3,$7,st,".","."}' SJ.out.tab > STAR.gff
方法2
star_to_gff.pl --star  SJ.out.tab --gff SJ.gff --label introns



#######组装质量评估
#1 计算平均subreads深度
minimap2 -ax map-pb -t 16 ../../MC.XN21.genome.chr.fasta ../../XN21.subreads.fasta.fasta >depth.sam
samtools sort -@ 30 depth.sam >depth.sort.bam
samtools depth depth.sorted.bam -a > depth.txt
cat depth.txt|wc -l
awk -F"\t" '{(total+=$3)};END{print total}' depth.txt |head
awk '$3 != 0 {count++} END {print count}' depth.txt
#2 单行perl确定组装基因组gap
perl -ne 'chomp;if( />(.*)/){$head = $1; $i=0; next};@a=split("",$_); foreach(@a){$i++;if($_ eq "N" && $s ==0 ){print "$head\t$i"; $s =1}elsif($s==1 && $_ ne "N"){print "\t$i\n";$s=0}}' ../../MC.XN21.genome.chr.fasta >gap.bed


######利用 Merqury v. 1.3 (Rhie et al. 2020)软件对拼接的基因组进行精确性评估
1.确定合适的kmer 大小
best_k.sh <genome_size> [tolerable_collision_rate=0.001]
2.构建db 
meryl threads=4 k=18 count ${subreads.fasta} output meryl_db
3
merqury.sh meryl_db genome.fasta merqury_res


#######使用 MUMmer v. 3.1 (Marçais et al. 2018)软件包中的 nucmer 程序 共线性分析
nucmer --maxmatch --maxgap=500 --mincluster=100 -t 6 --prefix=prefix  genome.fasta1 genome.fasta2
show-coords -d -l -L 5000-r -c -T strain1_stran2.delta > strain1_strain2.tab.txt


########组装线粒体基因组
#第一次使用需要 get_organelle_config.py --add embplant_pt,embplant_mt 下载库
#依赖spades不要用conda安装 去github安装 不然报错
get_organelle_from_reads.py -1 1.fq.clean.gz -2 2.fq.clean.gz --max-reads 7.5E7 -w 85 -R 15 -k 21,45,65,85,105,125,145 -F fungus_mt -o dirname 
# if you fail with the default database, use your own seed database and label database with "-s" and "--genes" 

##注释
https://chlorobox.mpimp-golm.mpg.de/geseq.html  参数默认 在ncbi下载酵母（或其他）genebank格式的线粒体序列
在reference Sequences本地上传
submit 会生成gb gff注释
##绘图https://irscope.shinyapps.io/Chloroplot/



#####利用二代数据填补基因组gap
GapCloser -b config.txt -a MC.XN21.genome.chr.fasta -o nogap.fasta -t 16

#config.txt文件内容

#maximal read length
max_rd_len=149
[LIB]
#average insert size
avg_ins=408
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=4
#use only first 100 bps of each read
rd_len_cutoff=100
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#a pair of fastq file, read 1 file should always be followed by read 2 file
q1=/mnt/d/data/XN21/ngs/Result/01.Cleandata/YLFS/YLFS.L350_FDSW21H000542-1r_1.fq.clean.gz
q2=/mnt/d/data/XN21/ngs/Result/01.Cleandata/YLFS/YLFS.L350_FDSW21H000542-1r_2.fq.clean.gz


####另一种挂载策略，对单倍体物种更友好
chromap+yahs （conda安装）
使用脚本yahs.sh


利用minimap2将三代数据比对到基因组上
再利用bwa将二代数据比对到基因组上
samtools计算深度 合并深度
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1$2]=$3;next}{if(a[$1$2]<$3){print $1,$2,$3}else{print $1,$2,a[$1$2]}}' pacbio.depth.txt ngs.depth.txt >combine.depth.txt


#将未锚定的序列合并为Chr00
先将没锚定的contig复制出来chrun.fa
grep "^>" -v chrun.fa| awk '{ ORS = ""; $1 = $1; print $0}' > chrUn.fasta
fold -w 100 chrUn.fasta > chrUnplaced.fasta
sed -i '1 i\>Chr00' chrUnplaced.fasta

#最终完成终极基因组序列

###busco
nohup busco -i genome.fa -c 10 -o genomebusco -m genome -l /mnt/d/db/fungi_odb10/fungi_odb10 --offline --force --augustus --augustus_species magnaporthe_grisea &