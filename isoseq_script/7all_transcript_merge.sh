#用tama将isoseq和ssrnaseq拼接去冗余后的转录本merge，（可以增加5端容忍度，避免降解造成的多转录本）得到最终版转录本
python /mnt/d/sf/tama-master/tama_merge.py -f filelist.txt -p MC_tama -e common_ends -d merge_dup -a 9999 -z 0 -m 0

isoseq和ssrna转录本merge的时候会将一部分isoseq转录本和ssrna合并，这会丢失长读转录本的优势，因此我们用以下方式将被merge的isoseq转录本换回来

运行python脚本 replace_merge2isoseq.ipynb
grep -F -v -w -f replace_old.list MC_tama.bed >no_need_replace.bed
cat no_need_replace.bed replace.bed >all.merge.v2.bed 


这里存在一个问题 将isoseq bed和ssrnaseq bed merge后 由于ssrnaseq读长短的弊端 会将一些原本为两个基因的错误组装成长转录本 并在merge时 将本来isoseq为两个基因的转录本合并成一个基因的转录本 因此 我们删除这些将多个isoseq基因merge成一个基因的基因 并用原本的isoseq替换回来
同时 还把长度小于200bp的转录本过滤掉
使用自编脚本del_mis_fusion_into_isoseq.ipynb完成

统计来自isoseq和ssRNAseq的转录本

自编脚本 bed12togtf 将bed12 文件转化为转录本的GTF
用TBTOOLS提取转录本序列  feature tag为exon
