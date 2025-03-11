########步骤2：利用参考基因组和链特异性RNA-Seq数据，通过软件TranscriptClean (Wyman and Mortazavi 2019) 对一致性序列进行修复，修复内容包括插入和缺失 (indels)，非匹配位点 (mismatches) 和剪接位点 (splice junction)。
#首先通过软件STAR (Dobin et al. 2013) 将与PacBio测序相同样品测序的不同时期链特异性RNA-Seq数据比对到参考基因组上，通过添加参数 --twopassMode Basic 进行两次比对从而获得更加精确的剪接位点信息文件SJ.out.tab，随后将参考基因组序序列和剪接位点信息文件对一致性序列进行修复。


###STAR   用拼接好的参考基因组序列构建索引
STAR --runMode genomeGenerate --runThreadN 2 --genomeDir starIndex --genomeFastaFiles ../MC.XN21.genome.chr.fasta
#使用自编脚本 STAR.sh  跑完以后可以得到一个SJ.out.tab, 该文件为剪切位点
##如果报错可以解压缩之后重试
nohup STAR.sh star_data starresult starIndex &
#跑完以后可以得到一个SJ.out.tab, 该文件为剪切位点

合并多个文件
cat *.tab |sort|uniq > JS_SJ.out.tab

#####安装TranscriptClean v2.0.3 from github
#依赖
#pyfasta (v0.5.2): https://pypi.python.org/pypi/pyfasta/
#Samtools (v1.9): https://github.com/samtools/samtools/releases/
#PyRanges: https://pyranges.readthedocs.io
#pybedtools (optional for variants) (v0.7.8): https://daler.github.io/pybedtools/
#Bedtools (optional for variants) (v2.25.0): http://bedtools.readthedocs.io/en/latest/content/installation.html

#genome环境下运行

""" Make sure that every chromosome in the SAM file also exists in the
        reference genome. This is a common source of crashes. Also, make sure
        that the VCF chromosome naming convention matches the SAM file if the
        user has provided a VCF. """
#纠错需要的sam文件 是三代polish后产生的一致性序列 gmap比对到参考基因组后 sort的sam文件
nohup python /mnt/d/sf/TranscriptClean-2.0.3/TranscriptClean.py --threads 16 --sam js.aligned.sorted.sam --genome MC.XN21.genome.chr.fasta --spliceJns JS_SJ.out.tab --maxLenIndel=5000 --maxSJOffset=80000 --outprefix result/js &
#--maxLenIndel=5000 --maxSJOffset=80000 要设置很大 避免因为缺失deletion太大无法纠错

##base环境下R安装不报错 小环境下报错
Rscript /mnt/d/sf/TranscriptClean-2.0.3/generate_report.R result/js