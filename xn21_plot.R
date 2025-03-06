## pacbio coverage
library(karyoploteR)
library(GenomicFeatures)
library(magrittr)
library(tidyverse)

setwd("d:/data/XN21/xn1genome/ChrInfo/")

gs <- read.table("XN1.Chr.fasta.fai", header = F)#fasta 索引
gs.gr <- GRanges(seqnames = gs$V1, ranges = IRanges(start=1, end = gs$V2), col="red")
pp = getDefaultPlotParams(plot.type = 2)
pp$leftmargin = 0.05
pp$topmargin = 80
pp$bottommargin = 15
pp$ideogramheight = 0
pp$data1inmargin = 0
pp$data1outmargin = 0

kp = plotKaryotype(genome = gs.gr, plot.params = pp, cex=0.6)


# add track of minimum PacBio coverage in bins

## samtools depth ydj_pbSorted_mapped.bam > ydj_pb.depth
## python3 minBinDepth.py ydj_pb.depth ydj_pb.cov

bincov <- read.table("二代回帖/meanCov.txt.txt")
bincov <- GRanges(seqnames = as.character(bincov$V1),
                  ranges = IRanges(start=bincov$V2, end = bincov$V3),
                  value = as.numeric(bincov$V6))

bincov$value[bincov$value>100] = 100
# kp = kpLines(kp, chr = seqnames(bincov),
#              x = start(bincov) + (end(bincov) - start(bincov)) / 2, 
#              y = bincov$value, 
#              ymin = 0, ymax = 200,
#              col = "darkblue", lwd = 1.5, clipping = T, r0 = 0, r1 = 1)

kp <- kpArea(kp, chr = seqnames(bincov), x = start(bincov) + (end(bincov) - start(bincov)) / 2,
             y = bincov$value, ymin = 0, ymax = 100, col = "lightblue", border = F)

#write.csv(bincov,"D:/data/XN21/pacbio/GAP_fill/bincov2.csv")

telo <- read.table("telo.txt", header = F)[,c(1,3,4)]
# left
kp = kpPoints(kp, chr = telo$V1, 
              x = as.numeric(telo$V3), 
              y = 0, col = "red")
# right
kp = kpPoints(kp, chr = telo$V1, 
              x = as.numeric(telo$V4), 
              y = 0, col = "red")

## add centro
centro<-read.table("D:/data/XN21/pacbio/chr_centro/centro.bed", header = F)
kp=kpRect(kp,chr = centro$V1,x0=centro$V2,x1=centro$V3,y0=0.2,y1=0.8)


###########################
######绘制染色体图像#######
###########################

library(RIdeogram)
library(tidyverse)
setwd("d:/data/XN21/xn1genome/ChrInfo/")
chrinfo<-read.table("XN1.Chr.fasta.fai",header = F)

centro<-read.table("../centro/centromic2.txt",header = F)

chrinfo<-chrinfo%>%mutate(Chr=V1,Start=0,End=V2)%>%select(c(6,7,8))
chrinfo$CE_start<-centro$V2
chrinfo$CE_end<-centro$V3
write.table(chrinfo,"Chrinfo.txt",sep = "\t",quote = F,row.names = F,col.names = T)
gene_density <-GFFex(input = "XN1.typical.genome.gff",karyotype = "Chrinfo.txt",feature = "mRNA",window = 10000)

telo<-read.table("telo.txt",header = F)%>%pivot_longer(cols = c(-1))%>%drop_na()%>%filter(name!="V2")%>%mutate(Chr=V1,Start=value,End=value,Type="Telomere",Shape="circle",color="6a3d9a" )%>%select(c(7,8,4,5,6,9))
gap<-read.table("genomeGap.txt",header = F)%>%mutate(Chr=V1,Start=V2,End=V3,Type="Gap",Shape="box",color="33a02c" )%>%select(c(7,8,4,5,6,9))
#centro<-read.table("quarTeT.best.candidate",header = F)%>%mutate(Chr=V1,Start=V2,End=V3,Type="Centromere",Shape="triangle",color="ff7f00" )%>%select(c(13,14,10,11,12,15))
telogap<-rbind(telo,gap)

ideogram(karyotype = chrinfo, overlaid = gene_density, label = telogap, label_type = "marker")
convertSVG("chromosome.svg", device = "pdf")

#############################################
#######基因转座子串联重复山峦图##############
#############################################
library(RIdeogram)
library(tidyverse)
setwd("d:/data/XN21/xn1genome/centro/")

gene_density <-GFFex(input = "XN1.typical.genome.gff",karyotype = "Chrinfo.txt",feature = "mRNA",window = 10000)

te_density <-GFFex(input = "XN1.Chr.fa.repeat.gff",karyotype = "Chrinfo.txt",feature = "repeat_region",window = 10000)


te_density_Chr1<-te_density%>%filter(Chr=="Chr01")



ggplot(data = te_density_Chr1, mapping = aes(x = Start, y = Value)) + geom_line()



#### xn21 circos ####

library(circlize)
library(Rsamtools)
library(GenomicFeatures)
## prepare data
xn21.fai <- read.table("d:/data/XN21/xn1genome/ChrInfo/XN1.Chr.fasta.fai")
canu.gs =  xn21.fai$V2
#txdb <- makeTxDbFromGFF("d:/data/pb_data/plot/circlize/EVM.all.sort.gff3", format = "gff3")
#ge <- genes(txdb)




chro <- xn21.fai$V1
starts = rep(0, length(canu.gs))
ends = canu.gs
genoCir <- data.frame(chr=chro, start=starts, end=ends)
genoCir$chr <- as.vector(genoCir[,1])

#### start plot######################################### only show fisrt 28 scaffolds of Morchella genome
circos.clear()
circos.par(start.degree = 87, track.height = 0.02, cell.padding = c(0,0,0,0), gap.degree=c(rep(1,14), 5))
circos.genomicInitialize(data = genoCir[1:15,],
                         sector.names = chro,
                         labels.cex = 0.5, track.height = 0.05, plotType = "labels")
### a: ideagram of 16 Chrs
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  idx = which(chro[1:15] == sector.index)
  circos.axis(labels.cex = 0.5,direction = "outside", labels.niceFacing = T, labels = "", minor.ticks = 5, lwd = 0.8, 
              major.at = c(0, canu.gs[idx]), major.tick.length = 0.4)
}, track.height = 0.05, bg.border = NA)

### b: plot gap and telo #绘制gap和端粒的位置
gap.bed<-read.table( "d:/data/XN21/xn1genome/ChrInfo/genomeGap.txt",header = F)#读取gap位置的bed文件
colnames(gap.bed)<-c("chr","start","end")
gap.bed$value<-gap.bed$end-gap.bed$start
#输入telo.bed文件 来自chr.telo.txt 需要awk处理
#awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$4}' chr.telo.txt|awk 'BEGIN{FS=OFS="\t"}{for (f=1; f <= NF; f+=1) {if (f !=1 && $f != "nd") {print $1,$f}}}' > telo.bed
telo.bed<-read.table( "d:/data/XN21/xn1genome/ChrInfo/telo.txt",header = F)
telo.bed1<-telo.bed[c(1,3)]
colnames(telo.bed1)<-c("chr","start")
telo.bed2<-telo.bed[c(1,4)]
colnames(telo.bed2)<-c("chr","start")
telo.bed<-rbind(telo.bed1,telo.bed2)
telo.bed<-drop_na(telo.bed)

telo.bed$end<-telo.bed$start
telo.bed$value<-1

all.bed.list<-list(gap.bed,telo.bed)


circos.genomicTrackPlotRegion(all.bed.list, track.height = 0.03,ylim = c(0,1), bg.border="lightgray",
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                circos.genomicRect(region, value, col = "black", border = "black", lwd=1)}
                                else{
                                  circos.genomicRect(region, value, col = "red", border = "red", lwd=1)
                                }
                                })

### c: reads mean coverage barplot   #这里的reads覆盖度是结合了pacbio和二代数据，求滑窗内平均覆盖度
read.df <- read.table("d:/data/XN21/pacbio/assembleQuality/GAP_fill/combine/minCov.txt.txt",header = F)
read.df<-read.df[,1:4]
colnames(read.df) <- c("chr", "start", "end", "value")
read.df$value <- log2(read.df$value + 1)
circos.genomicTrackPlotRegion(read.df,track.height=0.08, bg.border=NA, ylim = c(min(read.df$value), max(read.df$value)),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ybottom = 0, ytop.column = 1, 
                                                    col = "SlateBlue1", border = NA, lwd = 0.5, area = T)
                              })

### d: repeat heatmap
## repeat 
###awk 'BEGIN{FS = " ";OFS="\t"}{print $5,$6,$7}' ../repeat/result/MC.XN21.genome.fasta.out >repeat.anno.out.txt
repeat.df <- read.delim("d:/data/XN21/xn1genome/repeat/result/repeat.anno.out.txt",header = F,sep = '\t')

repeat.df <- data.frame(chr=repeat.df$V1, start=repeat.df$V2, end=repeat.df$V3)
bed.repeat <- genomicDensity(region = repeat.df, window.size = 10000, overlap = F)###circos.clear() 不然报错
cc2 <- colorRamp2(breaks = c(0,1), colors = c("white", "#a193c6"))
circos.genomicTrackPlotRegion(bed.repeat, track.height = 0.08, bg.border = NA,
                              panel.fun = function(region, value, ...){
                                circos.genomicRect(region, value, col = cc2(value), border = NA, lwd = 0.1)
                              })

####plot gene heatmap
txdb <- makeTxDbFromGFF("d:/data/XN21/xn1genome/ChrInfo/XN1.genome.gff", format = "gff")
ge <- genes(txdb)
ge.df <- data.frame(chr = seqnames(ge), start= start(ge), end=end(ge))%>%filter(chr!="ChrUnanchored")
## calculate gene density in slide windows
bed.gene <- genomicDensity(region = ge.df, window.size = 10000, overlap = F)
cc1 = colorRamp2(breaks = c(0, 1), colors = c("white","#E79397"))
circos.genomicTrackPlotRegion(bed.gene, track.height=0.08, bg.border=NA,
                              panel.fun = function(region, value, ...){
                                circos.genomicRect(region, value, col = cc1(value),border = NA, lwd=0.1)
                              })

####plot LTR
ltr<-read.delim("d:/data/XN21/xn1genome/ltr/XN1.Chr.fasta.out.gff3",header = F)
Copia_LTR_retrotransposon<-filter(ltr,V3=="Copia_LTR_retrotransposon")[c(1,4,5)]
Copia.df <- data.frame(chr=Copia_LTR_retrotransposon$V1, start=Copia_LTR_retrotransposon$V4, end=Copia_LTR_retrotransposon$V5)
bed.Copia <- genomicDensity(region = Copia.df, window.size = 10000, overlap = F)###circos.clear() 不然报错
vvs <- bed.Copia$value

circos.genomicTrackPlotRegion(bed.Copia, track.height = 0.05, bg.border = NA, ylim = c(min(vvs),max(vvs)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){

                              circos.genomicLines(region, value, col = "#f09640", border = NA,lwd = 0.02,
                                                     ybottom = 0, ytop.column = 1)
                              })

Gypsy_LTR_retrotransposon<-filter(ltr,V3=="Gypsy_LTR_retrotransposon")[c(1,4,5)]
Gypsy.df <- data.frame(chr=Gypsy_LTR_retrotransposon$V1, start=Gypsy_LTR_retrotransposon$V4, end=Gypsy_LTR_retrotransposon$V5)
bed.Gypsy <- genomicDensity(region = Gypsy.df, window.size = 10000, overlap = F)###circos.clear() 不然报错
vvs <- bed.Gypsy$value

circos.genomicTrackPlotRegion(bed.Gypsy, track.height = 0.05, bg.border = NA, ylim = c(min(vvs),max(vvs)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){
                                
                                circos.genomicLines(region, value, col = "#23749f", border = NA,lwd = 0.02,
                                                    ybottom = 0, ytop.column = 1)
                              })


### e: plot GC content
## GC content
GC = 0.4389398933180977

gc_content <- read.delim('d:/data/XN21/xn1genome/ChrInfo/binGC.txt',header = F)

bed.gc <- data.frame(chr = gc_content$V1, start = gc_content$V2, 
                     end = gc_content$V3,value=gc_content$V4-GC)


vv <- bed.gc$value
vv[which(vv %in% boxplot.stats(vv)$out)] <- mean(vv)   #检验数据中的异常值，将异常值赋值为平均值
bed.gc$value <- vv
bed.gc.list <- list(bed.gc[bed.gc$value>0, ], 
                    bed.gc[bed.gc$value<0, ])
circos.genomicTrackPlotRegion(bed.gc.list, track.height = 0.08, bg.border = NA, ylim = c(min(vv),max(vv)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                  # circos.genomicLines(region, value, col = "forestgreen", border = NA,
                                  #                     baseline = 0,area = T)
                                  circos.genomicRect(region, value, col = "forestgreen", border = NA,
                                                     ybottom = 0, ytop.column = 1)
                                }else{
                                  # circos.genomicLines(region, value, col = "firebrick1", border = NA,
                                  #                     baseline = 0,area = T)
                                  circos.genomicRect(region, value, col = "firebrick1", border = NA,
                                                     ybottom = 0, ytop.column = 1)
                                }})




### e: full-length LTR 
ltr.df <- read.table( "d:/data/XN21/pacbio/chr_ltr/MC.XN21.genome.chr.gapless.fasta.pass.list")
spl <- strsplit(ltr.df$V1, "[:.]")
chr = unlist(lapply(spl, function(x) x[1] ))
s = as.numeric(unlist(lapply(spl, function(x) x[2] )))
e = as.numeric(unlist(lapply(spl, function(x) x[4] )))
LTR.bed <- data.frame(chr = chr, start = s, end = e, type = ltr.df$V10)
LTR.bed$value <- 1
LTR.bed[LTR.bed$end < LTR.bed$start, ]
LTR.bed.list = list(LTR.bed[LTR.bed$type == "Gypsy", ], 
                    LTR.bed[LTR.bed$type == "Copia",  ],
                    LTR.bed[LTR.bed$type != "Gypsy" & LTR.bed$type != "Copia",  ])
circos.genomicTrackPlotRegion(LTR.bed.list, track.height = 0.05,ylim = c(0,1), bg.border="lightgray",
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                  circos.genomicRect(region, value, col = "orange", border = "orange", lwd = 1)
                                } else if(i == 2){
                                  circos.genomicRect(region, value, col = "red", border = "red", lwd = 1)
                                }else{
                                  circos.genomicRect(region, value, col = "blue", border = "blue", lwd = 1)
                                } })



#####绘制着丝粒位置

library(RIdeogram)
library(tidyverse)
setwd("d:/data/XN21/pacbio/chr_centro/")
data("human_karyotype")
centro<-read.delim("d:/data/XN21/pacbio/chr_centro/centroi.txt",header = T)
centro1<-drop_na(centro)
ideogram(karyotype = centro)
convertSVG("chromosome.svg",device = "png")


##############
#####转座子###
##############

library(ape)
library(Biostrings)
library(GenomicFeatures)
library(Rsamtools)
library(systemPipeR)

ltr1 <- read.table("d:/data/XN21/xn1genome/ltr/XN1.Chr.fa.pass.list")
ltr2 <- read.table("d:/data/XN21/xn1genome/ltr/XN1.Chr.fa.nmtf.pass.list")
ltr<-rbind(ltr1,ltr2)
genome <- "d:/data/XN21/xn1genome/ChrInfo/XN1.Chr.fa"

LTR_insert <- function(ltr, genome){
  ltr.df <- ltr
  a1 <- strsplit(ltr.df$V1, split="[.:]", perl=T)
  chrom <- unlist(lapply(a1, function(x) x[1]))
  s_lLTR <- as.numeric(unlist(lapply(a1, function(x) x[2])))
  e_rLTR <- as.numeric(unlist(lapply(a1, function(x) x[4])))
  
  a2 <- strsplit(ltr.df$V7, split="[.:]", perl=T)
  e_lLTR <- as.numeric(unlist(lapply(a2, function(x) x[2]))) - 1
  s_rLTR <- as.numeric(unlist(lapply(a2, function(x) x[4]))) + 1
  
  genome.fasta <- Rsamtools::FaFile(genome)
  l.ltr <- GRanges(seqnames = chrom, 
                   ranges = IRanges(start = s_lLTR, end = e_lLTR))
  r.ltr <- GRanges(seqnames = chrom, 
                   ranges = IRanges(start = s_rLTR, end = e_rLTR))
  
  l.ltr.seq <- getSeq(genome.fasta, l.ltr)
  r.ltr.seq <- getSeq(genome.fasta, r.ltr)
  
  pair_align <- pairwiseAlignment(l.ltr.seq, r.ltr.seq)
  seq1 <- as.character(pattern(pair_align))
  seq2 <- as.character(subject(pair_align))
  res <- list()
  # res2 <- list()
  for (i in 1:length(seq1)) {
    x <- list(strsplit(seq1, "")[[i]], strsplit(seq2, "")[[i]])
    xx <- as.DNAbin(x)
    d <- dist.dna(xx)[1]
    # res2[[i]] <- d
    # fungal substitution rate  nucleotides per site per year
    r <- 1.05*10e-9 
    divergent_time <- d/(2*r)
    res[[i]] <- divergent_time/1000000
  }
  rr1 <- unlist(res)
  return(rr1)
}
ydj.ltr.res <- density(LTR_insert(ltr, genome))
ydj.ltr.res$x[which.max(ydj.ltr.res$y)]
options(scipen = 200)
p3 <- ~plot(ydj.ltr.res, main = "Time of insertion of LTR",col="purple", 
            lwd=1.5, xlab="Age (MYA)",ylim=c(0,0.5), xlim=c(0,3)
            # xlim=c(0,9), ylim=c(0, 0.3)
)
abline(v=ydj.ltr.res$x[which.max(ydj.ltr.res$y)], col = "red")
cowplot::plot_grid(plotlist = list(p3))


###### repeat pie plot 
retro = 36.8
transpoon = 0.29
unclassi = 5.26
simple = 0.90
pie.data <- data.frame(group = c("Retroelement", "DNA transposon", "Unclassified repeat", "Simple repeat"), 
                       ratio=c(retro, transpoon, unclassi, simple))
aa <- pie.data$ratio/sum(pie.data$ratio)
labs <- paste0(round(aa,3)*100,"%")
p <- ggpubr::ggpie(pie.data, x = "ratio", label = labs[c(1,2,3,4)], fill = "group", lab.pos = "in",
                   color="white", palette = RColorBrewer::brewer.pal(6, "Set2"), legend.title="")
pp <- ggpubr::ggpar(p, legend = "right", tickslab = F)
#### LTR class barplot 
gypsy <- 1609
copia <- 4253
unknow <- 8056-1609-4253
par(mar=c(2, 5, 5, 2))
dd <- c(gypsy, copia, unknow)
names(dd) <- c("Gypsy", "Copia", "Unknow")
p2 <- ~barplot(dd, col = RColorBrewer::brewer.pal(3, "Set1"), ylab = "Counts")

cowplot::plot_grid(plotlist = list(p2,p3))


############WGCN绘图
setwd("d:/data/XN21/xn1genome/DEG/WGCNA/")
###图1  11个模块的代表eigengene
library(ggpubr)
library(tidyverse)


mat<-read.csv("MEs.csv",header = T,row.names = 1)
mat$Group<-c(rep("0h",3),rep("264h",2),rep("12h",3),rep("36h",3),rep("72h",3),rep("192h",2))
mat<-mat%>%pivot_longer(cols=-ncol(mat),names_to = "module")
mat$Group<-factor(mat$Group,levels = c("0h","12h","36h","72h","192h","264h"))
mat$module<-factor(mat$module,levels = c(paste("M",seq(10),sep = "")))
p<-ggline(mat,x="Group",y="value",
          facet.by = "module",
          color = "module",
          size = 1,
          add=c("mean_sd"))
ggpar(p,legend = "none")


####绘制性状和模块相关性图
library(corrplot)
library(tidyverse)
library(ggplotify)
mat<-read.csv("moduleTraitCor.csv",header = T,row.names = 1)
mat<-mat[order(as.integer(substr(row.names(mat),2,3))),]%>%as.matrix()
colnames(mat)<-substr(colnames(mat),2,5)

pmat<-read.csv("moduleTraitPvalue.csv",header = T,row.names = 1)
pmat<-pmat[order(as.integer(substr(row.names(pmat),2,3))),]%>%as.matrix()
colnames(pmat)<-substr(colnames(pmat),2,5)

col <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(50)
p1<-corrplot(mat,p.mat = pmat, insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white",col = col,method = "square",tl.col = "black")


##绘制每个模块hub基因的数量
library(ggprism)
hubdf<-tibble(Module=names(c(M1=1775,M2=762,M3=619,M4=625,M5=329,M6=397,M7=262,M8=262,M9=218,M10=297,M11=347)),Nums=c(M1=1775,M2=762,M3=619,M4=625,M5=329,M6=397,M7=262,M8=262,M9=218,M10=297,M11=347))
p2<-ggplot(hubdf)+
  geom_col(aes(x=Module,y=Nums,fill=Module))+
  scale_x_discrete(limits=c("M1", "M2","M3","M4","M5","M6","M7","M8","M9","M10","M11"))+
  theme_prism(palette = "flames",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16,
              base_line_size = 0.8,
              axis_text_angle = 45)+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0))+
  xlab("")+
  ylab("Number of transcripts")


###绘制hub SP的表达量热图
##library(ggplotify)
##library(pheatmap)
library(ComplexHeatmap)
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(50)
m1mat<-read.csv("hubgene/M1hubSP.tpm.csv",header = T,row.names = 1)
#m1mat<-as.matrix(m1mat)
##p1<-pheatmap(m1mat,scale = "row",color = col1,show_rownames = F,cluster_cols = F, show_colnames = F)
##p1<-as.ggplot(p1)
m2mat<-read.csv("hubgene/M2hubSP.tpm.csv",header = T,row.names = 1)

m3mat<-read.csv("模块性状相关/m3.secreted.csv",header = F,row.names = 1)


m4mat<-read.csv("模块性状相关/m4.secreted.csv",header = F,row.names = 1)

m5mat<-read.csv("模块性状相关/m5.secreted.csv",header = F,row.names = 1)


m6mat<-read.csv("模块性状相关/m6.secreted.csv",header = F,row.names = 1)


m7mat<-read.csv("hubgene/M7hubSP.tpm.csv",header = T,row.names = 1)


m8mat<-read.csv("模块性状相关/m8.secreted.csv",header = F,row.names = 1)

m9mat<-read.csv("hubgene/M9hubSP.tpm.csv",header = T,row.names = 1)
m10mat<-read.csv("hubgene/M10hubSP.tpm.csv",header = T,row.names = 1)
m11mat<-read.csv("hubgene/M11hubSP.tpm.csv",header = T,row.names = 1)

mat<-rbind(m3mat,m4mat,m5mat,m6mat,m8mat)
colnames(mat)<-c("0h","12h","36h","72h","192h","264h")
mat<-as.matrix(mat)
mat=t(scale(t(mat)))

Heatmap(mat, name = "z-score",
        row_split = factor(c(rep("M3", 59),
                             rep("M4", 43),rep("M5", 81),rep("M6", 24),
                             rep("M8", 12)),
                           levels = c("M3","M4","M5","M6","M8")),
        cluster_row_slices = FALSE,
        cluster_columns = F,
        show_row_names =F,
        col = col1
        )



########基因家族扩张收缩分析

##物种超度量树
library(ggtree)
library(ggtreeExtra)
library(tidytree)
library(tidyverse)
setwd("d:/data/XN21/xn1genome/CG/ortho2/")

tree<-read.tree("tree/spp_r8s_ultrametric.txt")
cafenum<-read.table("cafe/k2p/Gamma_clade_results.txt",sep = "\t",header = F)
cafenum$label<-gsub("<.+>","",cafenum$V1)
cafenum$info<-paste("+",cafenum$V2,"/","-",cafenum$V3,sep = "")%>%as.character()
cafe<-select(cafenum,c(4,5))%>%filter(label!="")

tibbletree<-tidytree::as_tibble(tree)%>%left_join(cafe)
tibbletree[is.na(tibbletree$info),][,5]<-""



tree<-as.treedata(tibbletree)

geneInfo<-read.delim("OrthoFinder/statMatrix.txt",header = T,row.names = 1)
geneInfo<-t(geneInfo)%>%as.data.frame()%>% select(c(2,3,9))
geneInfo$Non<-geneInfo[,1]-geneInfo[,3]
geneInfo$species<-row.names(geneInfo)
geneInfo<-geneInfo%>%select(c(2,3,4,5))%>%
  rename_with(~c("Not assigned","Species-species orthogroups","Non species-species orthogroups"),c(1,2,3))%>%
  pivot_longer(cols = c(1,2,3),names_to = "Type")


p<-ggtree(tree,layout = "roundrect")+
  geom_tiplab(align = T, fontface="italic")+
  geom_text(aes(label=tibbletree$info),vjust=-.3,hjust=1)+
  theme_tree2()+
  geom_treescale()

p1<-p+geom_fruit(
    data = geneInfo,
    geom = geom_col,
    mapping = aes(y=species,x=value,fill=Type),
    orientation="y",
    pwidth = 1,
    width=0.3,
    offset = 0.9,
    axis.params = list(
      axis="x",
      line.size=0.5,
      line.color="black",
      text.size=1.8,
      title="No. of Genes"))+ theme(legend.position = "top")+
  scale_x_continuous(breaks = c(0,50,100,150),expand = c(0,0),labels = c("150","100","50","0"))+
  xlab("MYA")



p2<-p1+geom_fruit(
  data =Pergenesinorthogroups,
  geom = geom_bar,
  stat="identity",
  mapping = aes(y=V1,x= as.numeric(V2)),
  fill="#F8AC8C",
  pwidth = 0.4,
  width=0.2,
  offset = 0.4,
  axis.params = list(
    axis="x",
    line.size=0.7,
    line.color="black",
    text.size=1.8,
    title="Genes in orthogroups (%)"
  )
  
)

p3<-p2+geom_fruit(
  data =peruniqgenes,
  geom = geom_bar,
  stat="identity",
  mapping = aes(y=V1,x= as.numeric(V2),fill=V3),
  pwidth = 0.4,
  width=0.2,
  offset = 0.7,
  axis.params = list(
    axis="x",
    line.size=0.7,
    line.color="black",
    text.size=1.8,
    title="Uniqe Genes (%)"
  )
  
)


###发生扩张的基因GO功能富集分析
library(clusterProfiler)
library(tidyverse)


term2gene <-read.delim('d:/data/XN21/xn1genome/function/GO/go2term.txt',header = F)
term2name<-read.delim('d:/data/XN21/xn1genome/function/GO/go2name.txt',header = T)

gene1<-read.delim('d:/data/XN21/xn1genome/CG/ortho2/OrthoFinder/Results_Jan30/Orthogroups/DM.species_specific.gene.list',header = F)
gene1 <- gene1$V1[1:nrow(gene1)]#genelist必须为vector

df <- enricher(gene = gene1, TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 1, qvalueCutoff = 1,pAdjustMethod="fdr")



input<-df@result
input<-separate(input,GeneRatio,c('Study.term','Study.total'),sep = '\\/')
input<-separate(input,BgRatio,c('Pop.term','Pop.total'),sep = '\\/')
#mergedf<-input%>%filter(p.adjust<0.05)


mergedf<-merge(input,term2name,by.x = 'Description',by.y = "GO_name",all.x = T)

#mergedf<-mergedf%>%filter(p.adjust<0.05,GO_ontology=="biological_process")
#write.table(mergedf,"d:/data/XN21/CG/ortho5/cafe/k2p/sample.expanded.og.gene.GO.tsv",row.names = F,quote = F,sep = "\t")

mergedf<-mergedf%>%filter(p.adjust<0.05)

##画柱状图
mergedf <- mergedf[order(mergedf$GO_ontology),]
CPCOLS <- c('#8DA1CB','#FD8D62','#66C3A5')
g <- ggplot(mergedf,aes(Description,Count,fill=GO_ontology))+geom_bar(stat = 'identity',width = 0.6)+coord_flip() +
  scale_fill_manual(values = CPCOLS,name="Ontology",labels=c("BP","CC","MF"))+theme_test()  +
  ylab('Number of gene')+
  scale_x_discrete(limits = factor(mergedf[,1]))+
  xlab("")
#scale_x_discrete(limits = factor(dfs[,13]))
g


##################################
#  比较基因组 dm特异家族进化分析 #
##################################

library(ggtree)
library(ggnewscale)
library(tidyverse)
library(ggtreeExtra)
library(jjPlot)


###整个界
setwd("d:/data/XN21/CG/ortho4/jinhua")

#读取树文件，并根据内节点名称分组
tree<-read.tree("all_class_name.txt.tree")
tree <- groupClade(tree, .node = c("Ascomycota","Basidiomycota"))

#读取物种数目文件
annotation<-read.csv("num.CSV",header = T)

#读取矩阵并处理
mat<-read.delim("utral_stat.tsv",header = T,row.names = 1)
mat<-as.data.frame(t(mat))
colnames(mat)<-c(seq(1:(ncol(mat)-1)),"species_sum")
mat$name<-rownames(mat)

matdf<-pivot_longer(mat,cols = c(-ncol(mat),-(ncol(mat)-1)),names_to = "type")
matdf$Presence<-matdf$value/matdf$species_sum
matdf$Absence<-1-matdf$value/matdf$species_sum
matdf<-matdf[order(matdf$type),]
matdf$seq<-seq(1:nrow(matdf))

matdf<-pivot_longer(matdf,cols = c(5,6),names_to = "ngroup",values_to = "nvalue")
matdf$type<-factor(matdf$type,levels = seq(1:23))

p<-ggtree(tree, aes(color=group), branch.length = 'none', 
          layout = "roundrect",
          ladderize=F,size=0.1) + 
  geom_tiplab(fontface = "bold.italic",size=3.5) + 
  scale_color_manual(values=c("black","orange","darkgreen"))

p<-ggtree::rotate(p,35)
p<-ggtree::rotate(p,36)

p1<-p+
  geom_fruit(
    data=annotation,
    geom=geom_text,
    mapping=aes(y=name, x=num,label=paste("(",num,")",sep = "")), 
    angle=0.1,
    size=3,  ##文字注释大小
    pwidth=0.2,
    offset = 8
    
  )

p2 <- p1 + 
  new_scale_fill() +  #添加新的scale
  geom_fruit(
    data=matdf,
    geom=geom_jjPointPie,  
    mapping=aes(y=name, x=type, group=seq,fill=ngroup, pievar=nvalue),  
    color="white",
    line.size=0.01,
    width=0.25,
    pwidth=16,
    offset = 12
    
  ) + 
  scale_fill_manual(
    values=c(Presence="#E64A19", Absence="grey"),  #设置扇形颜色
    limits=c("Presence","Absence"), #设置图例顺序
    name="" #隐藏图例名称
  )+
  guides(color="none")+ #隐藏进化树分组图例
  theme(
    legend.justification = c("right", "top")  #调整图例位置到右上角
  )


##一个纲
setwd("d:/data/XN21/CG/ortho4/jinhua/genus/")
#读取树文件，并根据内节点名称分组
tree<-read.tree("id.list.tree")  
###注意，读取的树文件属名种名分割不能是空格，空格在phylo对象里会被消除,导致后续的图行名对应不上



#读取矩阵并处理
mat<-read.csv("lmat.CSV",header = T)
mat$species<-gsub(" ","",mat$species)

matdf<-pivot_longer(mat,cols = c(-1))
matdf$value<-as.character(matdf$value)
matdf$name<-factor(matdf$name,level=colnames(mat[2:ncol(mat)]))


p<-ggtree(tree, branch.length = 'none', 
          layout = "roundrect",
          ladderize=F,size=0.1) + 
  geom_tiplab(fontface = "bold.italic",size=3.5) 

p1<-p+
  geom_fruit(
    data=matdf,
    geom=geom_text,
    mapping=aes(y=species, x=name,label=value,color=value), 
    angle=0.1,
    size=3,  ##文字注释大小
    pwidth=16,
    offset = 8
    
  )

p1<-p+
  geom_fruit(
    data=matdf,
    geom=geom_point,
    mapping=aes(y=species, x=name,shape=value,color=value), 
    size=3,  ##文字注释大小
    pwidth=16,
    offset = 4
    
  )+scale_shape_manual(values=c(4,19))+
  scale_color_manual(values = c("red","#008000"))

p2 <- p1 + 
  new_scale_fill() +  #添加新的scale
  geom_fruit(
    data=matdf,
    geom=geom_jjPointPie,  
    mapping=aes(y=name, x=type, group=seq,fill=ngroup, pievar=nvalue),  
    color="white",
    line.size=0.01,
    width=0.25,
    pwidth=16,
    offset = 12
    
  ) + 
  scale_fill_manual(
    values=c(Presence="#E64A19", Absence="grey"),  #设置扇形颜色
    limits=c("Presence","Absence"), #设置图例顺序
    name="" #隐藏图例名称
  )+
  guides(color="none")+ #隐藏进化树分组图例
  theme(
    legend.justification = c("right", "top")  #调整图例位置到右上角
  )

######################################
#########ltr插入基因各阶段表达量######
######################################

library(ComplexHeatmap)
setwd("d:/data/XN21/CG/ortho5/LTR插入相关基因/")
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(50)

mat<-read.csv("插入内部的基因/uniq.ltr.insert.csv",header = F,row.names = 1)

mat<-as.matrix(mat)
mat=t(scale(t(mat)))
colnames(mat)<-c("0h","12h","36h","72h","192h","264h")

Heatmap(mat, name = "z-score",
        cluster_columns = F,
        show_row_names =F,
        col = col1
)





####ABA合成和Chl降解基因
setwd("d:/data/XN21/RNA_seq_data/侵染时间ss/seq_name.outdir/")
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(50)

mat1<-read.csv("ABA.tpmmean.csv",header = T,row.names = 1)
mat2<-read.csv("Chl.tpmmean.csv",header = T,row.names = 1)



mat<-rbind(mat1,mat2)
mat<-as.matrix(mat)
mat=t(scale(t(mat)))
colnames(mat)<-c("Co","12h","36h","72h","192h","264h")

Heatmap(mat, name = "z-score",
        row_split = factor(c(rep("ABA biosynthesis", 7), rep("Chlorophyll degradation", 6)),
                           levels = c("ABA biosynthesis","Chlorophyll degradation")),
        cluster_row_slices = FALSE,
        cluster_columns = F,
        #show_row_names =F,
        col = col1
)


##对应的外泌蛋白表达谱

col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(50)

mat<-read.csv("../fungi.outdir/with264.sp.tpmmean.csv",header = T,row.names = 1)

mat<-as.matrix(mat)
mat=t(scale(t(mat)))
colnames(mat)<-c("Co","12h","36h","72h","192h","264h")

Heatmap(mat, name = "z-score",
        cluster_columns = F,
        show_row_names =T,
        col = col1
)



####不同功能分泌蛋白表达谱 效应 酶  保守 不保守
setwd("d:/data/XN21/CG/ortho5/比较基因组分泌蛋白/")
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(50)

mat1<-read.csv("coreSP.TPM.csv",header = F,row.names = 1)
mat2<-read.csv("shareSP.TPM.csv",header = F,row.names = 1)
mat3<-read.csv("uniqSP.TPM.csv",header = F,row.names = 1)


mat<-rbind(mat1,mat2,mat3)
mat<-as.matrix(mat)
mat=t(scale(t(mat)))
colnames(mat)<-c("0h","12h","36h","72h","192h","264h")
mat[is.na(mat)]<-0


Heatmap(mat, name = "z-score",
        row_split = factor(c(rep("Core", 57), rep("Share", 391),rep("Uniq", 58)),
                           levels = c("Core","Share","Uniq")),
        cluster_row_slices = FALSE,
        cluster_columns = F,
        show_row_names =F,
        #cluster_rows = F,
        col = col1
)

setwd("d:/data/XN21/isoseq/0609rebuild/10function/secreted/")
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(50)

mat1<-read.csv("uniq.Effector.csv",header = F,row.names = 1)
mat2<-read.csv("uniq.Enzyme.csv",header = F,row.names = 1)
mat3<-read.csv("Effector.Enzyme.csv",header = F,row.names = 1)
mat4<-read.csv("others.TPM.csv",header = F,row.names = 1)




mat<-rbind(mat1,mat2,mat3,mat4)
mat<-as.matrix(mat)
mat=t(scale(t(mat)))
colnames(mat)<-c("0h","12h","36h","72h","192h","264h")
mat[is.na(mat)]<--2



Heatmap(mat, name = "z-score",
        row_split = factor(c(rep("Effector", 121), rep("Enzyme", 171),rep("Effector&Enzyme", 30),rep("Others", 184)),
                           levels = c("Effector","Enzyme","Effector&Enzyme","Others")),
        cluster_row_slices = FALSE,
        cluster_columns = F,
        show_row_names =F,
        col = col1
)








#########################################
#######与D. r的染色体共线区域############
#########################################
library(karyoploteR)
library(magrittr)
library(tidyverse)


setwd("d:/data/XN21/xn1genome/CG/species/染色体共线/")


gs <- read.table("XN1.Chr.fa.fai", header = F)%>%mutate(V3=V2,V2=0)%>%select(c(1,2,3))

pp = getDefaultPlotParams(plot.type = 2)
pp$leftmargin = 0.1
pp$topmargin = 80
pp$bottommargin = 15
pp$ideogramheight = 0
pp$data1inmargin = 0
pp$data1outmargin = 0

kp = plotKaryotype(genome = gs, plot.params = pp, cex=0.6)

synteny<-read.table("Chrsynteny.merge.bed",header = F)


kp=kpPlotRegions(kp,synteny,col="#8FB075",avoid.overlapping=F,r0=0, r1=0.6)

ltrgff<-read.table("XN1.Chr.fa.LTR.gff3",header = F)
LTR.gr <- GRanges(seqnames = ltrgff$V1, ranges = IRanges(start=ltrgff$V4, end = ltrgff$V5))

kp <- kpPlotDensity(kp, LTR.gr, window.size = 10000,col="#D6E5FF",r0=0.6, r1=1.2)

specialgene<-read.table("species_specific.gene.gff",header = F)%>%filter(V3=="mRNA")

specialgene.gr <- GRanges(seqnames = specialgene$V1, ranges = IRanges(start=specialgene$V4, end = specialgene$V5))
kp=kpPlotRegions(kp,specialgene.gr,col="#833AB4",avoid.overlapping=F,r0=1.1, r1=1.6)

telo <- read.table("D:/data/XN21/xn1genome/ChrInfo/telo.txt", header = F)[,c(1,3,4)]
# left
kp = kpPoints(kp, chr = telo$V1, 
              x = as.numeric(telo$V3), 
              y = 0, col = "red",cex = 1)
# right
kp = kpPoints(kp, chr = telo$V1, 
              x = as.numeric(telo$V4), 
              y = 0, col = "red",cex = 1)



#####################################
#####Dm侵染阶段早期中期晚期##########
#####################################
library(ComplexHeatmap)
###数量
df=data.frame(nums=c(477,208,561),stage=c("Early","Middle","Late"))

CPCOLS <- c('#F1C89A','#E79397','#A797DA')
ggplot(df,aes(stage,nums,fill=stage))+
  geom_bar(stat = 'identity',width = 0.3,linewidth=1)+
  scale_fill_manual(values = CPCOLS)+
  theme_classic()  + 
  ylab('Number of Genes')+
  labs(fill='category')+
  scale_x_discrete(limits = factor(df[,2]))+
  xlab("")+
  theme(legend.position='none')


####各个阶段核心基因表达量

library(ComplexHeatmap)
setwd("d:/data/XN21/xn1genome/DEG/")


mat_R<-read.csv('d:/data/XN21/xn1genome/DEG/TPM/geneTPMmean2.CSV',header=T,row.names = 1)
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(100)
ha = HeatmapAnnotation(Stage = c("CK",rep("Biotroph",3),rep("Necrotroph",2)),
                       col = list(Stage = c("CK"="white","Biotroph"="#D6AFB9","Necrotroph"="#7E9BB7")))


listde<-read.table("wgcna_deg_trans.list",header = F)
listde<-c(listde[1])
listde<-unlist(listde)


mat1<-mat_R[rownames(mat_R) %in% listde,]
row_indices <- match(listde, rownames(mat1))  
mat1<-mat1[row_indices,]
mat1<-as.matrix(mat1)
mat1=t(scale(t(mat1)))
mat1<-na.omit(mat1)

mat2<-mat_R[rownames(mat_R) %in% listde,]
row_indices <- match(listde, rownames(mat2))  
mat2<-mat2[row_indices,]
mat2<-as.matrix(mat2)
mat2=t(scale(t(mat2)))
mat2<-na.omit(mat2)

mat3<-mat_R[rownames(mat_R) %in% listde,]
row_indices <- match(listde, rownames(mat3))  
mat3<-mat3[row_indices,]
mat3<-as.matrix(mat3)
mat3=t(scale(t(mat3)))
mat3<-na.omit(mat3)

mat<-rbind(mat1,mat2,mat3)

colnames(mat)<-c("0h","12h","36h","72h","192h","264h")


Heatmap(mat, name = "z-score",
        row_split = factor(c(rep("Biotrophic module", 1319), rep("Necrotrophic module", 2621),rep("Transition module", 698))),
        cluster_row_slices = FALSE,
        cluster_columns = F,
        cluster_row = T,
        show_row_names =F,
        col = col1,
        top_annotation = ha
)



###########不同阶段核心基因GO富集
library(clusterProfiler)
library(tidyverse)


term2gene <-read.delim('d:/data/XN21/xn1genome/function/GO/go2term.txt',header = F)
term2name<-read.delim('d:/data/XN21/xn1genome/function/GO/go2name.txt',header = T)

gene1<-read.delim('d:/data/XN21/xn1genome/DEG/wgcna_deg_trans.T0.list',header = F)
gene1 <- gene1$V1[1:nrow(gene1)]#genelist必须为vector

df <- enricher(gene = gene1, TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 1, qvalueCutoff = 1,pAdjustMethod="fdr")
input<-df@result
input<-separate(input,GeneRatio,c('Study.term','Study.total'),sep = '\\/')
input<-separate(input,BgRatio,c('Pop.term','Pop.total'),sep = '\\/')
#mergedf<-input%>%filter(p.adjust<0.05)
mergedf<-merge(input,term2name,by.x = 'Description',by.y = "GO_name",all.x = T)
mergedf<-mergedf%>%filter(p.adjust<0.05)

biodf<-mergedf #mid early
biodf$Group<-"biotrophic"

necdf<-mergedf
necdf$Group<-"necrotrophic"

#transdf<-mergedf
#transdf$Group<-"transition"

alldf<-rbind(biodf,necdf)

ggplot(alldf,aes(Group,Description))+
  geom_point(aes(size=Count,color=GO_ontology))+
  #scale_x_discrete(limits=factor(c("Early−infection","Mid−infection","Late−infection")))+
  scale_y_discrete(limits=factor(alldf$Description))+
  scale_size_continuous(range = c(3,6))+
  scale_color_manual(values = c("#E1C855","#E07B54","#51B1B7"))+ 
  labs(color='Ontology',size="Count",x="",y="",title="")+
  theme_bw()+theme(axis.text.y = element_text(color ="black",size = 12),legend.text = element_text(size = 14),legend.title=element_text(size=14),axis.title.x = element_text(size = 14))




#############################
########绘制寄主通路TPM######
#############################
library(ComplexHeatmap)

setwd("d:/data/XN21/xn1genome/DEG/host/")


mat_R<-read.csv('TPMmean.CSV',header=T,row.names = 1)
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(100)
ha = HeatmapAnnotation(Stage = c("CK",rep("biotrophic",3),rep("necrotrophic",2)),
                       col = list(Stage = c("CK"="white","biotroph"="#D6AFB9","necrotroph"="#7E9BB7")))



listde<-read.table("MdPathway/Md",header = F)
listde<-c(listde[2])
listde<-unlist(listde)
mat<-mat_R[rownames(mat_R) %in% listde,]


mat<-as.matrix(mat)
mat=t(scale(t(mat)))
mat<-na.omit(mat)

colnames(mat)<-c("0h","12h","36h","72h","192h","264h")


Heatmap(mat, name = "z-score",
        cluster_columns = F,
        show_row_names =T,
        col = col1,
        top_annotation = ha
)


#############################
####绘制lnc和gene桑吉图######
#############################
library(ggalluvial)

setwd("d:/data/XN21/xn1genome/lncRNA/hubLnc_target/targetFunction/")
data <- read.table("necLnc2allhub.out",header = F, check.names = F)
df <- to_lodes_form(data[,1:ncol(data)],
                    axes = 1:ncol(data),
                    id = "value")

col <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(102)

ggplot(df, aes(x = x, fill=stratum, label=stratum,
               stratum = stratum, alluvium  = value))+#数据
  geom_flow(width = 0.3,#连线宽度
            curve_type = "sine",#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
            alpha = 0.5,#透明度
            color = 'white',#间隔颜色
            size = 0.1)+#间隔宽度
  geom_stratum(width = 0.28)+#图中方块的宽度
  geom_text(stat = 'stratum', size = 2, color = 'black')+
  scale_fill_manual(values = col)+#自定义颜色
  theme_void()+#主题（无轴及网格线）
  theme(legend.position = 'none')#去除图例


#####LNC靶向CFEM效应蛋白
setwd("d:/data/XN21/xn1genome/lncRNA/hubLnc_target/CFEM/")
mat<-read.csv("cfemTPM.CSV",header = F,row.names = 1)
colnames(mat)<-c("0h","12h","36h","72h","192h","264h")

mat$group<-rownames(mat)
mat_df<-pivot_longer(mat,cols = -group)
group<-factor(mat_df$group,levels = c("lncRNA01282","lncRNA01417","DMXN08729"))

ggplot(mat_df,aes(x=name,y=value,group=group, color=group))+
  geom_point()+
  geom_line(size=1)+
  labs(x="",y="Transcripts Per Kilobase Million (TPM)")+
  scale_x_discrete(limit=c("0h","12h","36h","72h","192h","264h"))+
  theme_classic()


###########################
#######绘制吸器相关########
###########################


###1.表达量
library(ComplexHeatmap)

setwd("d:/data/XN21/xn1genome/DEG/WGCNA/hubGene/M6过渡模块/最终版本/")


mat_R<-read.csv('d:/data/XN21/xn1genome/DEG/TPM/geneTPMmean2.CSV',header=T,row.names = 1)
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(100)
ha = HeatmapAnnotation(Stage = c("CK",rep("biotroph",3),rep("necrotroph",2)),
                       col = list(Stage = c("CK"="white","biotroph"="#D6AFB9","necrotroph"="#7E9BB7")))



listde<-read.table("M6.utral.key.hub.gene.list",header = F)
listde<-c(listde[1])
listde<-unlist(listde)
mat<-mat_R[rownames(mat_R) %in% listde,]


mat<-as.matrix(mat)
mat=t(scale(t(mat)))
mat<-na.omit(mat)

colnames(mat)<-c("0h","12h","36h","72h","192h","264h")


Heatmap(mat, name = "z-score",
        cluster_columns = F,
        show_row_names =T,
        col = col1,
        top_annotation = ha
)



library(clusterProfiler)
library(tidyverse)


term2gene <-read.delim('d:/data/XN21/xn1genome/function/GO/go2term.txt',header = F)
term2name<-read.delim('d:/data/XN21/xn1genome/function/GO/go2name.txt',header = T)

gene1<-read.delim('DM.xiqi.houxuan.v1.list',header = F)
gene1 <- gene1$V1[1:nrow(gene1)]#genelist必须为vector

df <- enricher(gene = gene1, TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 1, qvalueCutoff = 1,pAdjustMethod="fdr")
input<-df@result
input<-separate(input,GeneRatio,c('Study.term','Study.total'),sep = '\\/')
input<-separate(input,BgRatio,c('Pop.term','Pop.total'),sep = '\\/')
#mergedf<-input%>%filter(p.adjust<0.05)
mergedf<-merge(input,term2name,by.x = 'Description',by.y = "GO_name",all.x = T)

#mergedf<-mergedf%>%filter(p.adjust<0.05,GO_ontology=="biological_process")
#write.table(mergedf,"d:/data/XN21/CG/ortho5/cafe/k2p/sample.expanded.og.gene.GO.tsv",row.names = F,quote = F,sep = "\t")

mergedf<-mergedf%>%filter(p.adjust<0.05)

##画柱状图
mergedf <- mergedf[order(mergedf$GO_ontology),]
CPCOLS <- c('#8DA1CB','#FD8D62','#66C3A5')
g <- ggplot(mergedf,aes(Description,Count,fill=GO_ontology))+geom_bar(stat = 'identity',width = 0.6)+coord_flip() +
  scale_fill_manual(values = CPCOLS,name="Type",labels=c("BP","CC","MF"))+theme_test()  +
  ylab('Number of gene')+
  scale_x_discrete(limits = factor(mergedf[,1]))+
  xlab("")
#scale_x_discrete(limits = factor(dfs[,13]))
g



sim_matrix<-read.csv("d:/data/XN21/xn1genome/function/secreted/SP/needPDB/matrix.csv",header = T)
sim_matrix<-as.matrix(sim_matrix)
dist_matrix <- as.dist(1 - sim_matrix)  
fit1<-hclust(dist_matrix,method = "average")

plot(fit1, main = "Hierarchical Clustering Dendrogram", sub = "", xlab = "", cex = 0.9)  

library(factoextra)
library(igraph)
colors = c(  
  "#FF0000",  # 红色  
  "#00FF00",  # 绿色  
  "#0000FF",  # 蓝色  
  "#FFFF00",  # 黄色  
  "#00FFFF",  # 青色  
  "#FF00FF",  # 紫色  
  "#C0C0C0",  # 银色  
  "#808080",  # 灰色  
  "#FFFFFF",  # 白色  
  "#000000",  # 黑色  
  "#FF5733",  # 珊瑚色  
  "#8B4513",  # 深棕色  
  "#2E8B57",  # 海绿色  
  "#9ACD32",  # 黄绿色  
  "#4B0082",  # 靛色  
  "#D2691E",  # 巧克力色  
  "#FF69B4",  # 热粉色  
  "#FFA500",  # 橙色  
  "#7FFFD4",  # 水绿色  
  "#008080",  # 茶绿色  
  "#ADD8E6",  # 淡蓝色  
  "#F08080",  # 浅珊瑚色  
  "#D3D3D3"   # 淡灰色  
)
fviz_dend(fit1,k=23,rect =T,rect_fill = T,type = "circular")

library(ComplexHeatmap)
mat<-read.csv('d:/data/XN21/xn1genome/function/secreted/SP/result/最优轮廓系数簇/cluster_matrix.csv',header=T,row.names = NULL)
mat<-as.matrix(mat)
col1 <- colorRampPalette(c("white", "#1E90FF"))(100)
ha = HeatmapAnnotation(Stage = c("CK",rep("biotroph",3),rep("necrotroph",2)),
                       col = list(Stage = c("CK"="white","biotroph"="#D6AFB9","necrotroph"="#7E9BB7")))






Heatmap(mat, name = " ",
        cluster_rows  = F,
        cluster_columns = F,
        #show_row_names =T,
        show_column_names = F,
        col = col1
        #top_annotation = ha
)


########################################
##########deepredeff预测效应蛋白########
########################################

#install.packages("deepredeff")
library(deepredeff)
library(tidyverse)
#install_tensorflow()
input_aas <- Biostrings::readAAStringSet("D:/data/XN21/xn1genome/function/secreted/secreted.pep.oneTranscript.fa")

pred_result <- predict_effector(
  input = input_aas,
  taxon = "fungi"
)

pred_result<-pred_result%>%filter(prediction=="effector")
write.table(pred_result$name,"D:/data/XN21/xn1genome/function/secreted/deepredeff_effector.list",quote = F,col.names = F,row.names = F)



#################################################
######结构聚类后阶段特异效应蛋白的时间表达模式
##################################################




library(ComplexHeatmap)
setwd("d:/data/XN21/xn1genome/function/secreted/pdb进化/")


mat_R<-read.csv('d:/data/XN21/xn1genome/DEG/TPM/geneTPMmean2.CSV',header=T,row.names = 1)
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(100)
ha = HeatmapAnnotation(Stage = c("CK",rep("Biotroph",3),rep("Necrotroph",2)),
                       col = list(Stage = c("CK"="white","Biotroph"="#D6AFB9","Necrotroph"="#7E9BB7")))

right_annotation=rowAnnotation(cluster=c(rep("cluster321",3),rep("cluster323",4),rep("cluster324",3)),
                               col = list(cluster = c("cluster321"="#A9C4E6","cluster323"="#F1DFA4","cluster324"="#9CD1C8")))

listde<-read.table("DM_uniq_cluster.txt",header = F)
listde<-c(listde[1])
listde<-unlist(listde)
listde<-substr(listde,1,9)
mat<-mat_R[rownames(mat_R) %in% listde,]
row_indices <- match(listde, rownames(mat))  
mat<-mat[row_indices,]
mat<-as.matrix(mat)
mat=t(scale(t(mat)))
mat<-na.omit(mat)

colnames(mat)<-c("0h","12h","36h","72h","192h","264h")


Heatmap(mat, name = "z-score",
        #row_split = factor(c(rep("Biotrophic module", 2), rep("Necrotrophic module", 10),rep("Transition module", 2))),
        cluster_row_slices = FALSE,
        cluster_columns = F,
        cluster_row = F,
        show_row_names =F,
        col = col1,
        top_annotation = ha,
        right_annotation = right_annotation
)


#################################################
######结构聚类后阶段特异效应蛋白的互作网络
##################################################
library(igraph)
setwd("d:/data/XN21/xn1genome/function/secreted/SP/result/最优轮廓系数簇/")
mat<-read.csv("cluster_matrix.csv",header = T)
rownames(mat)<-colnames(mat)

genelist<-read.table("不同侵染阶段wgcna_deg_spClusters/stage_special_cluster_gene.txt",header = F)
genelist<-c(genelist[1])
genelist<-unlist(genelist)
mat<-mat[rownames(mat) %in% genelist,]
mat<-mat[,colnames(mat) %in% genelist]

sim_matrix<-as.matrix(mat)
# 设置阈值  
threshold <- 0.5  

# 将相似性矩阵转换为二进制邻接矩阵  
adj_matrix <- sim_matrix > threshold  
diag(adj_matrix) <- 0  # 如果不需要自环，将对角线设置为0

g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = NULL, diag = FALSE)  

groups <- c(1, 1,1,1,1, 1,2, 2,3,3,1,1,1,1)  
names(groups) <- V(g)$name 

group_colors <- c('#8DA1CB','#FD8D62','#66C3A5')  

# 创建一个颜色向量，根据节点的组属性来决定颜色  
vertex_colors <- group_colors[groups[V(g)$name]] 


plot(g,edge.color = "grey",
     vertex.color = vertex_colors,
     vertex.frame.color="white",
     edge.width=3,
     layout = layout_with_fr)
legend(x=-1.5, y=-1.1, c("other","proto-oncogene"), pch=21,
       col="#777777", pt.bg=c("#80b1d3", "#fb8072"), 
       pt.cex=2, cex=.8, bty="n", ncol=1)


plot(g,   
     #vertex.color = component_colors, # 设置节点的颜色  
     vertex.label = NA, # 去除节点标签  
     edge.color = "grey", # 设置边的颜色  
     layout = layout_with_fr())

########################
###根据DEG进行筛选######
########################

library(tidyverse)
setwd("d:/data/XN21/xn1genome/DEG/deg/")

d12h<-read.csv("12hALLDESeq2.csv",header = T)
d36h<-read.csv("36hALLDESeq2.csv",header = T)
d72h<-read.csv("72hALLDESeq2.csv",header = T)
d192h<-read.csv("192hALLDESeq2.csv",header = T)
d264h<-read.csv("264hALLDESeq2.csv",header = T)


cut_off_fdr<-0.05
cut_off_logFC<-1
d12h$change<-ifelse(d12h$padj < cut_off_fdr & abs(d12h$log2FoldChange) >= cut_off_logFC, 
                    ifelse(d12h$log2FoldChange> cut_off_logFC ,'Up','Down'),
                    'Stable')

d36h$change<-ifelse(d36h$padj < cut_off_fdr & abs(d36h$log2FoldChange) >= cut_off_logFC, 
                    ifelse(d36h$log2FoldChange> cut_off_logFC ,'Up','Down'),
                    'Stable')

d72h$change<-ifelse(d72h$padj < cut_off_fdr & abs(d72h$log2FoldChange) >= cut_off_logFC, 
                    ifelse(d72h$log2FoldChange> cut_off_logFC ,'Up','Down'),
                    'Stable')


d192h$change<-ifelse(d192h$padj < cut_off_fdr & abs(d192h$log2FoldChange) >= cut_off_logFC, 
                    ifelse(d192h$log2FoldChange> cut_off_logFC ,'Up','Down'),
                    'Stable')

d264h$change<-ifelse(d264h$padj < cut_off_fdr & abs(d264h$log2FoldChange) >= cut_off_logFC, 
                    ifelse(d264h$log2FoldChange> cut_off_logFC ,'Up','Down'),
                    'Stable')

d12h<-select(d12h,c(1,8))
d36h<-select(d36h,c(1,8))
d72h<-select(d72h,c(1,8))
d192h<-select(d192h,c(1,8))
d264h<-select(d264h,c(1,8))

df<-d12h%>%full_join(d36h,by="X")%>%full_join(d72h,by="X")%>%full_join(d192h,by="X")%>%full_join(d264h,by="X")


biotro<-df%>%filter(change.x=="Up" |change.y=="Up" | change.x.x=="Up",change.y.y!="Up",change!="Up")%>%select(1)
nectro<-df%>%filter(change.x!="Up",change.y!="Up",change.x.x!="Up",change.y.y=="Up"|change=="Up")%>%select(1)
trans<-df%>%filter(change.x!="Up",change.y!="Up",change.x.x=="Up",change.y.y=="Up",change!="Up")%>%select(1)

write.table(biotro,"bio_deg_gene.list",quote = F,col.names = F,row.names = F)
write.table(nectro,"nec_deg_gene.list",quote = F,col.names = F,row.names = F)
write.table(trans,"trans_deg_gene.list",quote = F,col.names = F,row.names = F)



#########################################
library(ggnewscale)
library(tidyverse)
library(ggtreeExtra)
library(ggtree)

tree <- read.tree("d:/data/XN21/xn1genome/function/secreted/pdb进化/进化绘图/all_class_name.txt.tree")
tree <- groupClade(tree, .node = c("Basidiomycota"))
p0<-ggtree(tree, aes(color=group), branch.length = 'none', layout = "slanted",ladderize=F,size=0.1) + 
  geom_tiplab(fontface = "bold.italic",size=3.5) + scale_color_manual(values=c("black","#000093")) + 
  xlim(c(0, 35)) + theme(legend.position = "none") +geom_text2(aes(label=node),hjust=-.3,color="red")#+ geom_nodelab(node='internal')
p0  ##先看节点号 根据号旋转节点
p1 <- ggtree(tree, aes(color=group), branch.length = 'none', layout = "ellipse",ladderize=F,size=0.1) + 
  geom_tiplab(fontface = "bold.italic",size=3.5) + scale_color_manual(values=c("#24903A","#000093")) + 
  xlim(c(0, 35)) + theme(legend.position = "none")
p2<-ggtree::rotate(p1,18)
p2

geneInfo<-read.csv("d:/data/XN21/xn1genome/function/secreted/pdb进化/jinhua_matrix_.csv",header = T,row.names = 1)

geneInfo$species<-row.names(geneInfo)
geneInfo<-geneInfo%>%
  pivot_longer(cols = c(-ncol(geneInfo)),names_to = "cluster")
geneInfo$classification<-ifelse(geneInfo$value==1,"Presence","Absence")

p3<-p2+geom_fruit(
  data = geneInfo,
  geom = geom_point,
  mapping = aes(y=species,x=cluster,shape=classification),
  pwidth = 1,
  size=3,
  color="#8DA1CB",
  offset = 0.7)+ 
  scale_shape_manual(
    values=c(Presence=16, Absence=4),  #设置扇形颜色
    limits=c("Presence","Absence"), #设置图例顺序
    name="" #隐藏图例名称
  )+
  theme(legend.position = "top")

geneInfo2<-read.csv("d:/data/XN21/xn1genome/function/secreted/pdb进化/jinhua_matrix_len.csv",header = T,row.names = 1)

geneInfo2$species<-row.names(geneInfo2)
geneInfo2<-geneInfo2%>%
  pivot_longer(cols = c(-ncol(geneInfo2)),names_to = "cluster")
geneInfo2$value_category <- ifelse(geneInfo2$value == 0, "Zero",   
                                    ifelse(geneInfo2$value == 1, "One", "GreaterThanOne")) 

p3<-p2+geom_fruit_list(
  geom_fruit(
    data = geneInfo2,
    geom =geom_tile,
    mapping = aes(y=species,x=cluster,fill=value_category),
    offset = 0.7,
    alpha=0.5,
    pwidth = 1
  ),
  scale_fill_manual(values = c("#a6e3e9","#e3fdfd","white")),
  geom_fruit(
  data = geneInfo2,
  geom = geom_text,
  mapping = aes(y=species,x=cluster,label=value),
  pwidth = 1,
  size=3,
  angle=0.1,
  color="black",
  offset = 0.7))+ 
  theme(legend.position = "top")

##########真菌界进化分析

library(ggtree)
library(ggnewscale)
library(tidyverse)
library(ggtreeExtra)
library(jjPlot)

setwd("d:/data/XN21/xn1genome/CG/species/共线/jinhua")

#读取树文件，并根据内节点名称分组
tree<-read.tree("all_class_name.txt.tree")
tree <- groupClade(tree, .node = c("Ascomycota","Basidiomycota"))

#读取物种数目文件
annotation<-read.csv("num.CSV",header = T)

#读取矩阵并处理
mat<-read.csv("pie.CSV",header = T)
#write.csv(mat,"pie.CSV",quote = F,row.names = F)
matdf<-pivot_longer(mat,cols = c(-1,-ncol(mat),-(ncol(mat)-1)),names_to = "type")



p<-ggtree(tree, aes(color=group), branch.length = 'none', 
          layout = "ellipse",
          ladderize=F,size=0.1) + 
  geom_tiplab(fontface = "bold.italic",size=3.5) + 
  scale_color_manual(values=c("black","orange","darkgreen")) + 
  xlim(c(0, 70))#+geom_text2(aes(label=node),hjust=-.3,color="red")

p<-ggtree::rotate(p,14)

p1<-p+
  geom_fruit(
    data=annotation,
    geom=geom_text,
    mapping=aes(y=name, x=num,label=num), 
    angle=0.1,
    size=4,  ##文字注释大小
    pwidth=0.02,
    offset = 1.1
    
  )

p2 <- p1 + 
  new_scale_fill() +  #添加新的scale
  geom_fruit(
    data=matdf,
    geom=geom_jjPointPie,  
    mapping=aes(y=name, x=gene, group=seq, fill=type,pievar=value),  
    color="white",
    line.size=0.01,
    width=0.4,
    pwidth=6,
    offset = 4
    
  ) + 
  scale_fill_manual(
    values=c("#3f72af","grey"),  #设置扇形颜色
    limits=c("presence","absence"), #设置图例顺序
    name="" #隐藏图例名称
  )+
  guides(color="none")+ #隐藏进化树分组图例
  theme(
    legend.justification = c("right", "top")  #调整图例位置到右上角
  )


library(RIdeogram)
setwd("d:/data/XN21/xn1genome/CG/species/共线/")

karyotype_dual_comparison<-read.table("karyotype_dual_comparison.txt",header = T)
synteny_dual_comparison<-read.table("synteny_dual_comparison.txt",header = T)
ideogram(karyotype = karyotype_dual_comparison, synteny = synteny_dual_comparison)
convertSVG("chromosome.svg", device = "pdf")


##############trf 分布
library(RIdeogram)
library(karyoploteR)
library(tidyverse)
setwd("d:/data/XN21/xn1genome/ChrInfo/")


gs <- read.table("XN1.Chr.fasta.fai", header = F)%>%mutate(V3=V2,V2=0)%>%select(c(1,2,3))

pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin <- 20

kp = plotKaryotype(genome = gs, plot.params = pp, cex=0.6)

kp <- plotKaryotype(genome=gs, ideogram.plotter = NULL,
                    labels.plotter = NULL, plot.params = pp)
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp)

trfgff<-read.table("../trf/87bp_trf.gff3",header = F)
trf.gr <- GRanges(seqnames = trfgff$V1, ranges = IRanges(start=trfgff$V4, end = trfgff$V5))


kpPlotDensity(kp, trf.gr, col="#8FB075",window.size = 1000, r0=0, r1=0.4)



centro<-read.table("../centro/centromic2.txt", header = F)
kp=kpRect(kp,chr = centro$V1,x0=centro$V2,x1=centro$V3,y0=0,y1=0.3,col="#E64A19")




#########MUMmer dot-plot
# 加载所需的库  
library(tidyverse) 


# 读取并处理数据  
df <- read.table("d:/data/XN21/xn1genome/panGenome/mummer/dc1.delta.filter.coords", sep = "\t")[,c(1:7,12,13)]  # 从chrdelta.tsv文件中读取前9列数据  
# 重命名列名  
colnames(df) <- c("r_start", "r_end", "q_start", "q_end", "r_len", "q_len",   
                  "identity", "r_info", "q_info")  

# 读取基因组1的信息  
xchr <- read.table("d:/data/XN21/xn1genome/ChrInfo/XN1.Chr.fasta.fai", sep = "\t")  
# 选择并重命名需要的列  
xchrdf <- xchr %>% select(V1, V2) %>% `colnames<-`(c("Chr", "End"))  
# 计算累积长度（在每个染色体前加0，然后计算累积和，去掉最后一个值）  
xchrdf$cumulative_length <- cumsum(c(0, xchrdf$End[-nrow(xchrdf)]))  
# 计算总长度（直接计算累积和）  
xchrdf$sum_length <- cumsum(xchrdf$End)  

# 读取基因组2的信息，处理步骤与基因组1相同  
ychr <- read.table("d:/data/XN21/xn1genome/panGenome/DC1_JKI.genome.fasta.fai", sep = "\t")  
ychrdf <- ychr %>% select(V1, V2) %>% `colnames<-`(c("Chr", "End"))  
ychrdf$cumulative_length <- cumsum(c(0, ychrdf$End[-nrow(ychrdf)]))  
ychrdf$sum_length <- cumsum(ychrdf$End)  

# 将累积位置信息添加到df中  
df <- df %>%  
  # 通过r_info和Chr连接xchrdf的累积长度信息  
  left_join(xchrdf %>% select(Chr, cumulative_length), by = c("r_info" = "Chr")) %>%  
  # 通过q_info和Chr连接ychrdf的累积长度信息  
  left_join(ychrdf %>% select(Chr, cumulative_length), by = c("q_info" = "Chr")) %>%  
  # 根据累积长度计算新的起始和结束位置  
  mutate(r_start_cum = r_start + cumulative_length.x,  
         r_end_cum = r_end + cumulative_length.x,  
         q_start_cum = q_start + cumulative_length.y,  
         q_end_cum = q_end + cumulative_length.y)  

# 创建图形  
ggplot(df) +  
  # 添加线段，表示两个基因组之间的对应关系  
  geom_segment(aes(x = r_start_cum, xend = r_end_cum, y = q_start_cum, yend = q_end_cum,  
                   color = q_start < q_end)) +  
  # 设置颜色  
  scale_color_manual(values = c("red", "blue")) +  
  # 设置x轴（基因组1）  
  scale_x_continuous(name = "Genome1",   
                     breaks = c(0,xchrdf$sum_length),  # 设置断点  
                     labels = c(xchrdf$Chr,""),        # 设置标签  
                     limits = c(min(xchrdf$cumulative_length),max(xchrdf$sum_length)),  # 设置范围  
                     expand = c(0, 0)) +  # 不扩展范围  
  # 设置y轴（基因组2），步骤与x轴相同  
  scale_y_continuous(name = "Genome2",   
                     breaks = c(0,ychrdf$sum_length),   
                     labels = c(ychrdf$Chr,""),  
                     limits = c(min(ychrdf$cumulative_length),max(ychrdf$sum_length)),  
                     expand = c(0, 0)) +  
  # 设置主题  
  theme(    
    plot.background = element_blank(),  # 设置背景为空白    
    panel.background = element_blank(),  # 设置面板背景为空白    
    # 设置主要和次要网格线  
    panel.grid.major = element_line(color = "lightgrey", linewidth = 0.5, linetype = "solid"),    
    panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.25, linetype = "dotted"),    
    # 设置边框  
    panel.border = element_rect(colour = "black", fill = NA, size = 1),   
    # 设置刻度线  
    axis.ticks = element_line(colour = "black"),    
    # 移除图例  
    legend.position = "none"    
  ) +  
  # 设置x轴文本的角度和对齐方式  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+  
  # 设置图形标题  
  #labs(title = "Dot Plot")



##########绘制染色体SV密度
library(RIdeogram)
library(tidyverse)
setwd("d:/data/XN21/xn1genome/panGenome/SV")

chrinfo<-read.table("d:/data/XN21/xn1genome/ChrInfo/Chrinfo.txt",header = T)


gene_density <-GFFex(input = "d:/data/XN21/xn1genome/ChrInfo/XN1.typical.genome.gff",karyotype = "d:/data/XN21/xn1genome/ChrInfo/Chrinfo.txt",feature = "mRNA",window = 10000)

ltr_density<-read.table("d:/data/XN21/xn1genome/ltr/XN1.Chr.fasta.LTR.gff3",header = F)
ltr_density<-ltr_density%>%select(c(1,4,5))%>%`colnames<-`(c("Chr","Start","End"))

sv<-read.table("dc1_sv.clean.vcf",header = F)
sv<-sv%>%select(c(1,2))%>%mutate(V3=sv$V2)%>%`colnames<-`(c("Chr","Start","End"))

library(circlize)
ltr_density<-genomicDensity(region = ltr_density, window.size = 10000)
colnames(ltr_density)<-c("Chr","Start","End","Value")

sv.gene <- genomicDensity(region = sv, window.size = 10000)
sv.gene$color<-"fc8d62"
colnames(sv.gene)<-c("Chr","Start","End","Value","Color")


ideogram(karyotype = chrinfo, overlaid = ltr_density, label = sv.gene, label_type = "line", colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f"))
convertSVG("chromosome.svg", device = "pdf")


############SV非同义突变的基因表达热图
library(ComplexHeatmap)

setwd("d:/data/XN21/xn1genome/panGenome/SV/")


mat_R<-read.csv('SV_spTPM.CSV',header=F,row.names = 1)
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(100)
ha = HeatmapAnnotation(Stage = c("CK",rep("biotrophic",3),rep("necrotrophic",2)),
                       col = list(Stage = c("CK"="white","biotroph"="#D6AFB9","necrotroph"="#7E9BB7")))


right_annotation=rowAnnotation(DEG=c(rep("up",3),rep("other",3),rep("up",2),"other","up",rep("other",3)),
                               col = list(DEG = c("up"="#F1DFA4","other"="#D6AFB9")))


mat<-as.matrix(mat_R)
mat=t(scale(t(mat)))
mat<-na.omit(mat)

colnames(mat)<-c("hyphae","early","late")


Heatmap(mat, name = "z-score",
        cluster_columns = F,
        show_row_names =F,
        col = col1,
        rect_gp = gpar(col = "black", lwd = 1),
        right_annotation = right_annotation
)


####其他菌株二代测序比对到XN1
library(karyoploteR)
library(GenomicFeatures)
library(magrittr)
library(tidyverse)

setwd("d:/data/XN21/xn1genome/panGenome/")

gs <- read.table("d:/data/XN21/xn1genome/ChrInfo/XN1.Chr.fasta.fai", header = F)#fasta 索引
#gs<-gs[15,]
gs.gr <- GRanges(seqnames = gs$V1, ranges = IRanges(start=1, end = gs$V2), col="red")
pp = getDefaultPlotParams(plot.type = 2)
pp$leftmargin = 0.05
pp$topmargin = 80
pp$bottommargin = 15
pp$ideogramheight = 0
pp$data1inmargin = 0
pp$data1outmargin = 0

kp = plotKaryotype(genome = gs.gr, plot.params = pp, cex=0.6)


# add track of minimum PacBio coverage in bins

## samtools depth ydj_pbSorted_mapped.bam > ydj_pb.depth
## python3 minBinDepth.py ydj_pb.depth ydj_pb.cov

bincov <- read.table("phgb_second/meanCov.txt.txt")


####平均覆盖值
bincov <- GRanges(seqnames = as.character(bincov$V1),
                  ranges = IRanges(start=bincov$V2, end = bincov$V3),
                  value = as.numeric(bincov$V4))

bincov$value[bincov$value>100] = 100
# kp = kpLines(kp, chr = seqnames(bincov),
#              x = start(bincov) + (end(bincov) - start(bincov)) / 2, 
#              y = bincov$value, 
#              ymin = 0, ymax = 200,
#              col = "darkblue", lwd = 1.5, clipping = T, r0 = 0, r1 = 1)

kp <- kpArea(kp, chr = seqnames(bincov), x = start(bincov) + (end(bincov) - start(bincov)) / 2,
             y = bincov$value, ymin = 0, ymax = 100, col = "lightblue", border = F,r0=0, r1=0.3)

bincov2 <- read.table("dc1_second/meanCov.txt.txt")
####最小覆盖值
bincov2 <- GRanges(seqnames = as.character(bincov2$V1),
                  ranges = IRanges(start=bincov2$V2, end = bincov2$V3),
                  value = as.numeric(bincov2$V4))
bincov2$value[bincov2$value>100] = 100
kp <- kpArea(kp, chr = seqnames(bincov2), x = start(bincov2) + (end(bincov2) - start(bincov2)) / 2,
             y = bincov2$value, ymin = 0, ymax = 100, col="#ddaacc", border = F,r0=0.35, r1=0.65)

bincov3 <- read.table("dc1_ont/meanCov.txt.txt")
####最小覆盖值
bincov3 <- GRanges(seqnames = as.character(bincov3$V1),
                   ranges = IRanges(start=bincov3$V2, end = bincov3$V3),
                   value = as.numeric(bincov3$V4))
bincov3$value[bincov3$value>100] = 100
kp <- kpArea(kp, chr = seqnames(bincov3), x = start(bincov3) + (end(bincov3) - start(bincov3)) / 2,
             y = bincov3$value, ymin = 0, ymax = 100, col="grey", border = F,r0=0.7, r1=1)


telo <- read.table("d:/data/XN21/xn1genome/ChrInfo/telo.txt", header = F)[,c(1,3,4)]
# left
kp = kpPoints(kp, chr = telo$V1, 
              x = as.numeric(telo$V3), 
              y = 0, col = "red")
# right
kp = kpPoints(kp, chr = telo$V1, 
              x = as.numeric(telo$V4), 
              y = 0, col = "red")


###chr15上LTR分布
library(karyoploteR)
library(GenomicFeatures)
library(magrittr)
library(tidyverse)

gs <- read.table("d:/data/XN21/xn1genome/ChrInfo/XN1.Chr.fasta.fai", header = F)#fasta 索引
gs<-gs[15,]
gs.gr <- GRanges(seqnames = gs$V1, ranges = IRanges(start=1, end = gs$V2), col="red")
pp = getDefaultPlotParams(plot.type = 2)
pp$leftmargin = 0.05
pp$topmargin = 10
pp$bottommargin = 0
pp$ideogramheight = 15

kp = plotKaryotype(genome = gs.gr, plot.params = pp, cex=0.6)

Copiagff<-read.table("d:/data/XN21/xn1genome/ltr/LTR.Copia.gff3",header = F)
Copia.gr <- GRanges(seqnames = Copiagff$V1, ranges = IRanges(start=Copiagff$V4, end = Copiagff$V5))

kp <- kpPlotDensity(kp, Copia.gr, window.size = 10000,col="#D6E5FF",r0=0.42, r1=0.82)

Gypsygff<-read.table("d:/data/XN21/xn1genome/ltr/LTR.Gypsy.gff3",header = F)
Gypsy.gr <- GRanges(seqnames = Gypsygff$V1, ranges = IRanges(start=Gypsygff$V4, end = Gypsygff$V5))

kp <- kpPlotDensity(kp, Gypsy.gr, window.size = 10000,col="#8FB075",r0=0, r1=0.4)



####Chr15上基因 的表达量
############SV非同义突变的基因表达热图
library(ComplexHeatmap)

setwd("d:/data/XN21/xn1genome/Chr15")


mat_R<-read.csv('chr15_geneTPM2.CSV',header=F,row.names = 1)
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(100)

ha = HeatmapAnnotation(Stage = c("CK",rep("biotroph",3),rep("necrotroph",2)),
                       col = list(Stage = c("CK"="white","biotroph"="#D6AFB9","necrotroph"="#7E9BB7")))


right_annotation=rowAnnotation(DEG=c("conserved","specific","HGT","specific",rep("specific in Leotiomycetes",2),"specific","HGT",rep("specific",3),rep("conserved",3)),
                               col = list(DEG = c("conserved"="#e6875c","specific"="#a393c5","HGT"="#1bafa1","specific in Leotiomycetes"="#65ade0")))


mat<-as.matrix(mat_R)
mat=t(scale(t(mat)))
mat<-na.omit(mat)

#colnames(mat)<-c("hyphae","early","late")
colnames(mat)<-c("0h","12h","36h","72h","192h","264h")

Heatmap(mat, name = "z-score",
        cluster_columns = F,
        show_row_names =F,
        col = col1,
        rect_gp = gpar(col = "black", lwd = 1),
        top_annotation = ha,
        right_annotation = right_annotation
)



#########Fig S8 AS统计
library(tidyverse)  

# 从指定路径读取一个以制表符分隔的文件，并设置header为True表示文件的第一行是列名  
df <- read.delim("d:/data/XN21/xn1genome/IsoseqA/AS鉴定/gencode.DMXN.all.events.ioe", header = T)  

# 使用管道操作符对df进行处理  
# 1. 使用separate函数将event_id列拆分为gID和event_message两列，分隔符是";"  
# 2. 使用mutate函数和正则表达式从event_message中提取事件类别到新的列event_class  
# 3. 使用select函数选择gene_id, event_class, alternative_transcripts三列  
df <- df %>% separate(col = event_id, into = c("gID", "event_message"), sep = ";") %>%  
  mutate(event_class = str_extract(event_message, "^[^:]+")) %>%  
  select(gene_id, event_class, alternative_transcripts)  

# 打印出发生可变剪接的基因数量  
print(paste("Number of AS genes is", length(unique(df$gene_id))))  

# 将df中的alternative_transcripts列转换为数据框，并重命名列为asline  
# 然后使用separate_rows函数将asline列中的每一行按","分隔成多行  
df2 <- as.data.frame(df$alternative_transcripts) %>%   
  `colnames<-`("asline") %>%  
  separate_rows(asline, sep = ",")  

# 打印出经过可变剪接的转录本数量  
print(paste("Number of AS transcripts is", length(unique(df2$asline))))  

# 对df进行处理，按gene_id分组，然后统计每个基因的不同事件类别数量  
# 过滤出那些具有至少两种不同事件类别的基因  
multi_group_genes <- df %>%    
  group_by(gene_id) %>%    
  summarise(group_count = n_distinct(event_class)) %>%    
  filter(group_count >= 2)    

# 打印出具有多种组合可变剪接事件的基因数量  
print(paste("genes had isoforms resulting from multiple combinatory AS events: ", nrow(multi_group_genes)))  

# 统计df中event_class的频数，并转换为数据框  
# 然后添加Freq的整数版本和百分比版本  
classdf <- as.data.frame(table(df$event_class))  
classdf <- classdf %>% mutate(Freq = as.integer(Freq)) %>% mutate(percent = Freq / sum(classdf$Freq))  

# 打印出事件类别的频数和百分比  
print(classdf)  

# 使用ggplot2绘制事件类别的柱状图  
ggplot(classdf, aes(Var1, Freq / 1000, fill = Var1)) +  
  geom_bar(stat = 'identity', width = 0.8) +  
  scale_fill_manual(values = c("#88d2d2", "#bae6fa", "#4c98ce", "#c4c4d3", "#ffef96", "#eb9d9e", "#fbcac5")) +  
  geom_text(aes(label = Freq), vjust = -0.5, size = 5, color = "black") + 
  ylab(expression("No. of AS events (10"^"4)")) +  
  xlab("") +  
  labs(fill="")+
  theme_classic()

##########谱系特异性区域绘制
library(RIdeogram)
library(tidyverse)
setwd("d:/data/XN21/xn1genome/panGenome/mummer")

chrinfo<-read.table("d:/data/XN21/xn1genome/ChrInfo/Chrinfo.txt",header = T)
#bedtools complement -i specific.bed -g ../../ChrInfo/XN1.Chr.fasta.fai > non_covered_regions.bed
specific.bed<-read.table("specific.bed",header = F)
shared.bed<-read.table("shared.bed",header = F)
conserved.bed<-read.table("conserved.bed",header = F)

colnames(specific.bed)<-c("Chr","Start","End")
specific.bed$Value=1
colnames(shared.bed)<-c("Chr","Start","End")
shared.bed$Value=0
colnames(conserved.bed)<-c("Chr","Start","End")
conserved.bed$Value=-1


gene_density <-rbind(specific.bed,shared.bed,conserved.bed)

ltr_density<-read.table("d:/data/XN21/xn1genome/ltr/XN1.Chr.fasta.LTR.gff3",header = F)
ltr_density<-ltr_density%>%select(c(1,4,5))%>%`colnames<-`(c("Chr","Start","End"))



library(circlize)
ltr_density<-genomicDensity(region = ltr_density, window.size = 10000)
ltr_density$color<-"ff7f00"
colnames(ltr_density)<-c("Chr","Start","End","Value","Color")




ideogram(karyotype = chrinfo, overlaid = gene_density, label = ltr_density, label_type = "polygon", colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f"))
convertSVG("chromosome.svg", device = "pdf")






###################双速基因组鉴定

library(mixtools)
library(depmixS4)
state2col <- function(x) {
  col <- rep("black", length(x))
  for (i in 1:length(x)) {
    if (x[i]==1) {
      col[i] <- c("fill_color=yellow")
    } else {
      col[i] <- c("fill_color=purple")
    }
  }
  return(col)
}

# fit 2 speed model
## with EM
vcf = read.table("d:/data/XN21/xn1genome/panGenome/twospeed/vcfdensity.txt")
vcfem = normalmixEM(vcf$V4/10,k=2)
#### model assessment
mu <- mean(vcf$V4/10)
sigma <- sd(vcf$V4/10)
# D'Test and KS test check normality
ks.test(vcf$V4/10, "pnorm",mu, sigma)
library(fBasics)
dagoTest(x = vcf$V4/10)
## log likehood function of Gaussian
# methods 1
normal <- function(theta,x){
  mu <- theta[1]
  sigma <- theta[2]
  n <- length(x)
  logL <- -0.5*n*log(2*pi)-0.5*n*log(sigma^2)-(1/(2*sigma^2))*sum((x-mu)**2)
  return (-logL)
}
res <- optim(c(mu,sigma),normal,x=vcf$V4/10)  # optimization
# method 2
LL <- function(mu, sigma) {
  R = dnorm(vcf$V4/10, mu, sigma)
  -sum(log(R))
}
library(stats4)
fit <- mle(LL, start = list(mu = mu, sigma = sigma))
# summary(fit)
# logLik(fit)
aic.gaussian <- 2*2 - (-2*res$value)
aic.mixGaussian <- 2*5 - 2*vcfem$loglik
# Likelihood Rate Test judge whether normal distribution can substitude mixGaussian model
LR = 2*((vcfem$loglik)- (-res$value))
dchisq(x = LR, df = 2)
# expand y lim for plot.mixEM
plot.em <- function(x,
                    whichplots = 1,
                    loglik = 1 %in% whichplots,
                    density = 2 %in%
                      whichplots,
                    xlab1 = "Iteration",
                    ylab1 = "Log-Likelihood",
                    main1 = "Observed Data Log-Likelihood",
                    col1 = 1,
                    lwd1 = 2,
                    xlab2 = NULL,
                    ylab2 = NULL,
                    main2 = NULL,
                    col2 = NULL,
                    lwd2 = 2,
                    alpha = 0.05,
                    marginal = FALSE,
                    ...) {
  mix.object <- x
  k <- ncol(mix.object$posterior)
  x <- sort(mix.object$x)
  a <- hist(x, plot = FALSE)
  maxy <-
    max(max(a$density), 0.3989 * mix.object$lambda / mix.object$sigma) + 0.01
  if (is.null(main2)) {
    main2 <- "Density Curves"
  }
  if (is.null(xlab2)) {
    xlab2 <- "Data"
  }
  if (is.null(col2)) {
    col2 <- 2:(k + 1)
  }
  hist(
    x,
    prob = TRUE,
    main = main2,
    xlab = xlab2,
    ylim = c(0, maxy),
    ...
  )
  if (length(mix.object$mu) == 1) {
    arbvar <- TRUE
    mix.object$sigma <- mix.object$scale * mix.object$sigma
    arbmean <- FALSE
  }
  if (length(mix.object$mu) == k && length(mix.object$sigma) ==
      1) {
    arbmean <- TRUE
    arbvar <- FALSE
  }
  if (length(mix.object$sigma) == k && length(mix.object$mu) == k) {
    arbmean <- TRUE
    arbvar <- TRUE
  }
  for (i in 1:k) {
    lines(
      x,
      mix.object$lambda[i] * dnorm(x, mean = mix.object$mu[i *arbmean + (1 - arbmean)], sd = mix.object$sigma[i *arbvar + (1 - arbvar)]),
      col = col2[i],
      lwd = lwd2
    )
  }
}
par(mar=c(4.5,4.5,1,1))
plot.em(vcfem,breaks=50,lwd2 = 2, w=1.1, xlab2="Variants/kb",main2="",
        col2=c("#1B9E77","#D95F02"))
# str(vcfem)
library(shape)
Arrows(5.5, 0.05, 4.5, 0.057, arr.length = 0.2, segment = T, code = 1, arr.adj = 0.5, col="purple")
Arrows(18, 0.055, 19, 0.062, arr.length = 0.2, segment = T, code = 1, arr.adj = 0.5, col="orange" )
text(4.4, 0.060, expression(paste("slow: ", mu, " = 8.35,  ", sigma," = 6.19")))
text(21, 0.066, expression(paste("fast: ", mu, " = 24.04,  ", sigma," = 12.72")))

# mtext("A", adj=0.01, line=-2, outer=T)

## decode with veterbi
msp <- depmix(V4/10~1,nstates=2,data=vcf)
set.seed(1)
fmsp <- fit(msp)
#plot(ts(posterior(fmsp)))  
state = posterior(fmsp)$state

# check accuracy
ss <- vcfem$posterior
aa <- ifelse(ss[,1]>ss[,2],1,2)
xx <- aa == state
length(xx[which(xx==TRUE)])/length(xx)
# write.table(paste(vcf$V1, vcf$V2, vcf$V3, state),file="/Users/alexwang/Downloads/density/state.tsv", sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(paste(vcf$V1, vcf$V2, vcf$V3, state),file="d:/data/XN21/xn1genome/panGenome/twospeed/twospeed_state3.txt", sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

###############
library(circlize)
library(Rsamtools)
library(GenomicFeatures)
library(tidyverse)
## prepare data
xn21.fai <- read.table("d:/data/XN21/xn1genome/ChrInfo/XN1.Chr.fasta.fai")
canu.gs =  xn21.fai$V2
#txdb <- makeTxDbFromGFF("d:/data/pb_data/plot/circlize/EVM.all.sort.gff3", format = "gff3")
#ge <- genes(txdb)




chro <- xn21.fai$V1
starts = rep(0, length(canu.gs))
ends = canu.gs
genoCir <- data.frame(chr=chro, start=starts, end=ends)
genoCir$chr <- as.vector(genoCir[,1])

#### start plot######################################### only show fisrt 28 scaffolds of Morchella genome
circos.clear()
circos.par(start.degree = 87, track.height = 0.02, cell.padding = c(0,0,0,0), gap.degree=c(rep(1,14), 5))
circos.genomicInitialize(data = genoCir[1:15,],
                         sector.names = chro,
                         labels.cex = 0.5, track.height = 0.05, plotType = "labels")
### a: ideagram of 16 Chrs
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  idx = which(chro[1:15] == sector.index)
  circos.axis(labels.cex = 0.5,direction = "outside", labels.niceFacing = T, labels = "", minor.ticks = 5, lwd = 0.8, 
              major.at = c(0, canu.gs[idx]), major.tick.length = 0.4)
}, track.height = 0.05, bg.border = NA)

### b: plot gap and telo #绘制gap和端粒的位置


twospeed.bed<-read.table("d:/data/XN21/xn1genome/panGenome/twospeed/twospeed_state.txt",header = F)

all.bed.list<-list(twospeed.bed[twospeed.bed$V4=="fill_color=yellow", ], 
                   twospeed.bed[twospeed.bed$V4=="fill_color=purple", ])


circos.genomicTrackPlotRegion(all.bed.list, track.height = 0.1,ylim = c(0,1), bg.border="lightgray",
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                  circos.genomicRect(region, value, col = "#ff8002", border = "#ff8002", lwd=1)}
                                else{
                                  circos.genomicRect(region, value, col = "#08baff", border = "#08baff", lwd=1)
                                }
                              })

### c: reads mean coverage barplot   #这里的reads覆盖度是结合了pacbio和二代数据，求滑窗内平均覆盖度
read.df <- read.table("d:/data/XN21/xn1genome/panGenome/twospeed/vcfdensity.txt",header = F)

colnames(read.df) <- c("chr", "start", "end", "value")
read.df$value <- log2(read.df$value + 1)
circos.genomicTrackPlotRegion(read.df,track.height=0.06, bg.border=NA, ylim = c(min(read.df$value), max(read.df$value)),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ybottom = 0, ytop.column = 1, 
                                                    col = "forestgreen", border = NA, lwd = 0.1, area = T)
                              })



####plot LTR
ltr<-read.delim("d:/data/XN21/xn1genome/ltr/XN1.Chr.fasta.out.gff3",header = F)
Copia_LTR_retrotransposon<-filter(ltr,V3=="Copia_LTR_retrotransposon")[c(1,4,5)]
Copia.df <- data.frame(chr=Copia_LTR_retrotransposon$V1, start=Copia_LTR_retrotransposon$V4, end=Copia_LTR_retrotransposon$V5)
bed.Copia <- genomicDensity(region = Copia.df, window.size = 10000, overlap = F)###circos.clear() 不然报错
vvs <- bed.Copia$value

circos.genomicTrackPlotRegion(bed.Copia, track.height = 0.05, bg.border = NA, ylim = c(min(vvs),max(vvs)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){
                                
                                circos.genomicLines(region, value, col = "#DE582B", border = NA,lwd = 0.06,
                                                    ybottom = 0, ytop.column = 1)
                              })

Gypsy_LTR_retrotransposon<-filter(ltr,V3=="Gypsy_LTR_retrotransposon")[c(1,4,5)]
Gypsy.df <- data.frame(chr=Gypsy_LTR_retrotransposon$V1, start=Gypsy_LTR_retrotransposon$V4, end=Gypsy_LTR_retrotransposon$V5)
bed.Gypsy <- genomicDensity(region = Gypsy.df, window.size = 10000, overlap = F)###circos.clear() 不然报错
vvs <- bed.Gypsy$value

circos.genomicTrackPlotRegion(bed.Gypsy, track.height = 0.05, bg.border = NA, ylim = c(min(vvs),max(vvs)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){
                                
                                circos.genomicLines(region, value, col = "#808080", border = NA,lwd = 0.06,
                                                    ybottom = 0, ytop.column = 1)
                              })


circos.genomicTrackPlotRegion(all.bed.list, track.height = 0.02,ylim = c(0,1), bg.border="lightgray",
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                  circos.genomicRect(region, value, col = "#ff8002", border = "#ff8002", lwd=1)}
                                else{
                                  circos.genomicRect(region, value, col = "#08baff", border = "#08baff", lwd=1)
                                }
                              })



#############################################
## repeat region associate accessory chromosome and 2-speed regions

cc <- c(rgb(8,186,255,maxColorValue = 255), rgb(255,128,2, maxColorValue = 255))
library(regioneR)
library(GenomicFeatures)

repeat.gr <- toGRanges(read.table("d:/data/XN21/xn1genome/ltr/LTR.bed", header = F))


df.speed <- read.table("d:/data/XN21/xn1genome/panGenome/twospeed/twospeed_state3.txt")
# df.speed <- read.table(paste0("/Users/alexwang/0data/0mango/genome/2speed/", strain,"_EMstate.tsv"))

speed.gr <- toGRanges(df.speed)
fast.gr <- speed.gr[speed.gr$V4 == 2,]
slow.gr <- speed.gr[speed.gr$V4 == 1,]

# calculate overlapped base of repeat && fast/slow regions
ovBase.fast <- overlapRegions(repeat.gr, fast.gr, get.bases = T)$ov.bases
# (sum(ovBase.fast)*10000)/sum(width(fast.gr))
m.fast <- mean(ovBase.fast)
ovBase.slow <- overlapRegions(repeat.gr, slow.gr, get.bases = T)$ov.bases
# (sum(ovBase.slow)*10000)/sum(width(slow.gr))
m.slow <- mean(ovBase.slow)

# chisq.test
p <- c(1/2, 1/2)
x <- c(m.fast, m.slow)
chisq.test(x = x, p = p,correct =F)
## plot bar
cc <- c(rgb(8,186,255,maxColorValue = 255), rgb(255,128,2, maxColorValue = 255))
par(mar=c(1.5,2,5,2))
barplot(c(m.fast/1000, m.slow/1000), horiz = T, space = 0.5,col = cc, las=1, xaxt="n")
legend(x = 0.4, y = 0.5, legend = c("Fast", "Slow"), 
       col = cc, xpd = T, pch=15, bty = "n", ncol = 2)
axis(side = 3, cex=0.5)
#text(x = 0.4, y=1.7, "*** X-squared p<2.2e-16")
title("Number of repeats / 10 kb")



#####gene


###
gene.gr <-toGRanges(read.table("d:/data/XN21/xn1genome/panGenome/twospeed/result/XN1.Chr.gene.bed", header = F))

# calculate overlapped base of repeat && fast/slow regions
ovBase.fast <- overlapRegions(gene.gr, fast.gr, get.bases = T)$ov.bases
# (sum(ovBase.fast)*10000)/sum(width(fast.gr))
m.fast <- mean(ovBase.fast)
ovBase.slow <- overlapRegions(gene.gr, slow.gr, get.bases = T)$ov.bases
# (sum(ovBase.slow)*10000)/sum(width(slow.gr))
m.slow <- mean(ovBase.slow)

# chisq.test
p <- c(1/2, 1/2)
x <- c(m.fast, m.slow)
chisq.test(x = x, p = p,correct =F)
## plot bar
cc <- c(rgb(8,186,255,maxColorValue = 255), rgb(255,128,2, maxColorValue = 255))
par(mar=c(1.5,2,5,2))
barplot(c(m.fast/1000, m.slow/1000), horiz = T, space = 0.5,col = cc, las=1, xaxt="n")
legend(x = 0.4, y = 0.5, legend = c("Fast", "Slow"), 
       col = cc, xpd = T, pch=15, bty = "n", ncol = 2)
axis(side = 3, cex=0.5)
#text(x = 0.4, y=1.7, "*** X-squared p<2.2e-16")
title("Number of genes / 10 kb")


### permutation test of the location correlation secretome to closest repeat
library(GenometriCorr)
library(regioneR)
library(tidyverse)

repeat.gr <- toGRanges(read.table("d:/data/XN21/xn1genome/ltr/LTR.bed", header = F))

fast.repeat.gr<-overlapRegions(repeat.gr, fast.gr, get.bases = T)%>%
  select(1,2,3)%>%distinct()%>%toGRanges()

slow.repeat.gr<-overlapRegions(repeat.gr, slow.gr, get.bases = T)%>%
  select(1,2,3)%>%distinct()%>%toGRanges()

gene1.gr <-toGRanges(read.table("d:/data/XN21/xn1genome/panGenome/twospeed/result/fast.speed.gene.bed", header = F))

gs <- read.table("d:/data/XN21/xn1genome/centro/chrom.sizes", header = F)
gs.len <- gs$V2
#names(gs.len) <- paste0(strain,"_",gs$chr)

### 统计检验

geCor <- GenometriCorr::GenometriCorrelation(query = gene1.gr, reference = fast.repeat.gr,
                                             permut.number = 100, keep.distributions = TRUE,
                                             chromosomes.to.proceed = gs$V1, chromosomes.length = gs.len)
GenometriCorr::graphical.report(geCor, pdffile = "d:/data/XN21/xn1genome/panGenome/twospeed/result/repeatPermutationTest_fast.SM.pdf")

########################################################################## plot LTR percent

library(ggpubr)


set.seed(123)  # 设置随机种子以确保结果的可重复性
group1 <- read.table("d:/data/XN21/xn1genome/panGenome/jcvi/fast.speed.gene.kaks.txt",header = F)
group2 <- read.table("d:/data/XN21/xn1genome/panGenome/jcvi/slow.speed.gene.kaks.txt",header = F)

data <- data.frame(
  value = c(group1$V1, group2$V1),
  group = c(rep("Fast", length(group1$V1)), rep("Slow", length(group2$V1)))
)


ggboxplot(data, x = "group", y = "value", color="black",
          fill = "group", palette = c("#08baff", "#ff8002"),
          #order = c("Group1", "Group2"),
          ylab = "Ration of non-syn to syn", xlab = "",ylim = c(0, 2))
t.test(value ~ group, data = data)




###################特异和扩张基因家族基因与LTR的距离

setwd("d:/data/XN21/xn1genome/CG/ortho2/analysis/")


repeat.gr <- toGRanges(read.table("DM.expanded.family.gene.bed", header = F))


df.speed <- read.table("d:/data/XN21/xn1genome/panGenome/twospeed/twospeed_state2.txt")
# df.speed <- read.table(paste0("/Users/alexwang/0data/0mango/genome/2speed/", strain,"_EMstate.tsv"))

speed.gr <- toGRanges(df.speed)
fast.gr <- speed.gr[speed.gr$V4 == 2,]
slow.gr <- speed.gr[speed.gr$V4 == 1,]

# calculate overlapped base of repeat && fast/slow regions
ovBase.fast <- overlapRegions(repeat.gr, fast.gr, get.bases = T)$ov.bases
# (sum(ovBase.fast)*10000)/sum(width(fast.gr))
m.fast <- mean(ovBase.fast)
ovBase.slow <- overlapRegions(repeat.gr, slow.gr, get.bases = T)$ov.bases
# (sum(ovBase.slow)*10000)/sum(width(slow.gr))
m.slow <- mean(ovBase.slow)

# chisq.test
p <- c(1/2, 1/2)
x <- c(m.fast, m.slow)
chisq.test(x = x, p = p,correct =F)
## plot bar
cc <- c(rgb(8,186,255,maxColorValue = 255), rgb(255,128,2, maxColorValue = 255))
par(mar=c(1.5,2,5,2))
barplot(c(m.fast/1000, m.slow/1000), horiz = T, space = 0.5,col = cc, las=1, xaxt="n")
legend(x = 0.4, y = 0.5, legend = c("Fast", "Slow"), 
       col = cc, xpd = T, pch=15, bty = "n", ncol = 2)
axis(side = 3, cex=0.5)
text(x = 0.4, y=1.7, "*** X-squared p<2.2e-16")
title("Number of repeats / 10 kb")

#####特异和扩张基因家族在双速基因组和谱系特异性区域的分布

library(RIdeogram)
library(tidyverse)


chrinfo<-read.table("d:/data/XN21/xn1genome/ChrInfo/Chrinfo.txt",header = T)
#bedtools complement -i specific.bed -g ../../ChrInfo/XN1.Chr.fasta.fai > non_covered_regions.bed
specific.bed<-read.table("d:/data/XN21/xn1genome/panGenome/mummer/specific.bed",header = F)
shared.bed<-read.table("d:/data/XN21/xn1genome/panGenome/mummer/shared.bed",header = F)
conserved.bed<-read.table("d:/data/XN21/xn1genome/panGenome/mummer/conserved.bed",header = F)

colnames(specific.bed)<-c("Chr","Start","End")
specific.bed$Value=1
colnames(shared.bed)<-c("Chr","Start","End")
shared.bed$Value=0
colnames(conserved.bed)<-c("Chr","Start","End")
conserved.bed$Value=-1


gene_density <-rbind(specific.bed,shared.bed,conserved.bed)

geneGff<-read.table("d:/data/XN21/xn1genome/ChrInfo/XN1.genome.gff",header = F)%>%filter(V3=="gene")%>%
  separate(V9,c("V10","V11"),sep = "=")%>%select(V1,V4,V5,V11)


listde<-read.table("d:/data/XN21/xn1genome/CG/ortho2/analysis/DM.species_specific.gene.list",header = F)
listde<-c(listde[1])
listde<-unlist(listde)


specificGene<-geneGff[geneGff$V11 %in% listde,]%>%
  `colnames<-`(c("Chr","Start","End","V11"))%>%
  mutate(Type="specific",Shape="circle",color="6a3d9a")%>%
  select(Type,Shape,Chr,Start,End,color)



listde<-read.table("d:/data/XN21/xn1genome/CG/ortho2/analysis/sample.expanded.gene.list",header = F)
listde<-c(listde[1])
listde<-unlist(listde)


expandedGene<-geneGff[geneGff$V11 %in% listde,]%>%
  `colnames<-`(c("Chr","Start","End","V11"))%>%
  mutate(Type="expanded",Shape="box",color="ff7f00")%>%
  select(Type,Shape,Chr,Start,End,color)


Random_Genes<-rbind(specificGene,expandedGene)





ideogram(karyotype = chrinfo, overlaid = gene_density, label = specificGene, label_type = "marker",colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f"))
convertSVG("chromosome.svg", device = "pdf")


speed_density <-read.table("d:/data/XN21/xn1genome/panGenome/twospeed/twospeed_state2.txt",header = F)
colnames(speed_density)<-c("Chr","Start","End","Value")
ideogram(karyotype = chrinfo, overlaid = speed_density, label = Random_Genes, label_type = "marker",colorset1 = c("#ff8002","#08baff"))
convertSVG("chromosome.svg", device = "pdf")



library(GenometriCorr)
library(regioneR)

repeat.gr <- toGRanges(read.table("d:/data/XN21/xn1genome/ltr/LTR.bed", header = F))

gene1<-read.table("d:/data/XN21/xn1genome/CG/ortho2/analysis/DM.uniq.family.gene.bed", header = F)
gene2<-read.table("d:/data/XN21/xn1genome/CG/ortho2/analysis/DM.expanded.family.gene.bed", header = F)
gene<-rbind(gene1,gene2)


gene1.gr <-toGRanges(gene)

gs <- read.table("d:/data/XN21/xn1genome/centro/chrom.sizes", header = F)
gs.len <- gs$V2
#names(gs.len) <- paste0(strain,"_",gs$chr)

### 统计检验

geCor <- GenometriCorr::GenometriCorrelation(query = gene1.gr, reference = repeat.gr,
                                             permut.number = 100, keep.distributions = TRUE,
                                             chromosomes.to.proceed = gs$V1, chromosomes.length = gs.len)
GenometriCorr::graphical.report(geCor, pdffile = "repeatPermutationTest_uniq.SM.pdf")





#####################
setwd("d:/data/XN21/xn1genome/DEG/deg/")

library(tidyverse)

csv_12<-read.csv("12hALLDESeq2.csv",header = T)%>%select(1,3,7)
csv_36<-read.csv("36hALLDESeq2.csv",header = T)%>%select(1,3,7)
csv_72<-read.csv("72hALLDESeq2.csv",header = T)%>%select(1,3,7)
csv_192<-read.csv("192hALLDESeq2.csv",header = T)%>%select(1,3,7)
csv_264<-read.csv("264hALLDESeq2.csv",header = T)%>%select(1,3,7)


merge_csv<-full_join(csv_12,csv_36,by="X")%>%full_join(csv_72,by="X")%>%full_join(csv_192,by="X")%>%full_join(csv_264,by="X")


listde<-read.table("d:/data/XN21/xn1genome/panGenome/graphgenome/gaf/coreGene/Dc_Dr_core_genes_ID.list",header = F)

listde<-c(listde[1])
listde<-unlist(listde)
listde<-substr(listde,1,9)
mat<-merge_csv[merge_csv$X %in% listde,]

mat_filter<-filter(mat,(log2FoldChange.x>1 & padj.x <0.05) | (log2FoldChange.y>1 & padj.y <0.05) | (log2FoldChange.x.x>1 & padj.x.x <0.05))

#mat_filter<-filter(mat,((log2FoldChange.x>1 & padj.x <0.05) | (log2FoldChange.y>1 & padj.y <0.05) | (log2FoldChange.x.x>1 & padj.x.x <0.05)) & log2FoldChange.y.y<=1 & log2FoldChange<=1 )



mat<-select(mat_filter,c(2,4,6,8,10))

rownames(mat)<-mat_filter$X
write.table(mat_filter,"d:/data/XN21/xn1genome/panGenome/graphgenome/gaf/coreGene/xiqi_gene_fc.txt",sep = "\t",col.names = T,quote = F)



colnames(mat)<-c("12h","36h","72h","192h","264h")


library(ComplexHeatmap)
col1 <- colorRampPalette(c("#3CB371", "white", "#1E90FF"))(100)

split = rownames(mat)

ha = HeatmapAnnotation(Stage = c(rep("biotroph",3),rep("necrotroph",2)),
                       col = list(Stage = c("biotroph"="#D6AFB9","necrotroph"="#7E9BB7")))


mat<-as.matrix(mat)
Heatmap(mat, name = "Log2FC",
        cluster_columns = F,
        row_split = split,
        row_title = NULL,
        show_row_names =F,
        col = col1,
        top_annotation = ha,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        }
)



################绘制物种分歧进化树
library(ggtree)

library(ggnewscale)
library(tidytree)
library(tidyverse)
setwd("d:/data/XN21/xn1genome/panGenome/xiqiCG/ortho_o/tree/")

tree<-read.tree("spp_r8s_ultrametric.txt")

tibbletree<-as_tibble(tree)

calculate_divergence_time<-function(treesdata){
  dftip <- data.frame(node = numeric(), divergence_time = numeric())
  internode<-data.frame(parent=numeric(),node = numeric(), divergence_time = numeric())
  for (i in 1:nrow(treesdata)) {
    if (isTip(treesdata,treesdata$node[i])){
      newtip<-data.frame(node = treesdata$parent[i],divergence_time =treesdata$branch.length[i])
      dftip<-rbind(dftip,newtip)
    }
    else{
      if (treesdata$node[i] != treesdata$parent[i]){
        newdata<-data.frame(parent=treesdata$parent[i],node = treesdata$node[i],divergence_time =treesdata$branch.length[i])
        internode<-rbind(internode,newdata)
        
      }
    }
  }
  #print(internode)
  unfindnode<-data.frame(parent=numeric(),node = numeric(), divergence_time = numeric())
  for (j in 1:nrow(internode)) {
    if (internode$node[j] %in% dftip$node){
      newtip<-data.frame(node = internode$parent[j],divergence_time =internode$divergence_time[j]+dftip[which(dftip$node==internode$node[j]),2])
      dftip<-rbind(dftip,newtip)
    }else{
      newnode<-data.frame(parent=internode$parent[j],node = internode$node[j],divergence_time =internode$divergence_time[j])
      unfindnode<-rbind(unfindnode,newnode)
    }
  }
  print(unfindnode)
  for (s in 1:nrow(unfindnode)){
    if (unfindnode$node[s] %in% dftip$node){
      newtip<-data.frame(node = unfindnode$parent[s],divergence_time =unfindnode$divergence_time[s]+dftip[which(dftip$node==unfindnode$node[s]),2])
      dftip<-rbind(dftip,newtip)
    }
  }
  dftip$divergence_time<-round(dftip$divergence_time, digits = 1)
  dftip<-dftip%>% distinct()
  return(dftip)
}

myainfo<-calculate_divergence_time(tibbletree)

tibbletree<-left_join(tibbletree,myainfo,by="node")

tree<-as.treedata(tibbletree)

p<-ggtree(tree,layout = "ellipse")+
  geom_tiplab(align = T, fontface="italic")+
  #geom_text(aes(label=tibbletree$branch.length),vjust=-.3,hjust=1)+
  geom_point(aes(shape=isTip, color=isTip))+
  geom_text2(aes(label=divergence_time),vjust=-.3,hjust=1)+
  geom_treescale()+
  guides(shape="none",color="none")



# 计算 x 轴的范围（从右开始）
#x_range <- range(tree_data$xend, na.rm = TRUE)
start <- max(tibbletree$branch.length,na.rm = T)
end <- 0




sequence <- seq(from = start, to = end, by = -50)
label_seq<-seq(from = end, to = start, by = 50)

p1<-p +
  scale_x_continuous(breaks = seq(from = start, to = end, by = -150),
                     labels = seq(from = end, to = start, by = 150),
                     expand = expansion(mult=c(0.15,0.7)))+
  theme(axis.line.x = element_line(colour = "black"), # 添加 x 轴线
        axis.ticks.x = element_line(colour = "black"), # 添加 x 轴刻度线
        axis.text.x = element_text(colour = "black")) # 添加 x 轴刻度文本

p2<-p1+
  new_scale_fill()+
  geom_fruit(
    data=tibbletree,
    geom=geom_text,
    mapping=aes(y=label, x=label,label=label), 
    angle=0.1,
    size=3,  ##文字注释大小
    pwidth=0.2,
    offset = 1
    
  )

#############泛基因组统计绘图
###############
library(circlize)
library(Rsamtools)
library(GenomicFeatures)
library(tidyverse)
## prepare data
xn21.fai <- read.table("d:/data/XN21/xn1genome/ChrInfo/XN1.Chr.fasta.fai")
canu.gs =  xn21.fai$V2
#txdb <- makeTxDbFromGFF("d:/data/pb_data/plot/circlize/EVM.all.sort.gff3", format = "gff3")
#ge <- genes(txdb)




chro <- xn21.fai$V1
starts = rep(0, length(canu.gs))
ends = canu.gs
genoCir <- data.frame(chr=chro, start=starts, end=ends)
genoCir$chr <- as.vector(genoCir[,1])

#### start plot######################################### only show fisrt 28 scaffolds of Morchella genome
circos.clear()
circos.par(start.degree = 87, track.height = 0.02, cell.padding = c(0,0,0,0), gap.degree=c(rep(1,14), 5))
circos.genomicInitialize(data = genoCir[1:15,],
                         sector.names = chro,
                         labels.cex = 0.5, track.height = 0.05, plotType = "labels")
### a: ideagram of 16 Chrs
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  idx = which(chro[1:15] == sector.index)
  circos.axis(labels.cex = 0.5,direction = "outside", labels.niceFacing = T, labels = "", minor.ticks = 5, lwd = 0.8, 
              major.at = c(0, canu.gs[idx]), major.tick.length = 0.4)
}, track.height = 0.05, bg.border = NA)

### b: plot gap and telo #绘制gap和端粒的位置


twospeed.bed<-read.table("d:/data/XN21/xn1genome/panGenome/twospeed/twospeed_state3.txt",header = F)

all.bed.list<-list(twospeed.bed[twospeed.bed$V4=="1", ], 
                   twospeed.bed[twospeed.bed$V4=="2", ])


circos.genomicTrackPlotRegion(all.bed.list, track.height = 0.1,ylim = c(0,1), bg.border="lightgray",
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                  circos.genomicRect(region, value, col = "#ff8002", border = "#ff8002", lwd=NULL)}
                                else{
                                  circos.genomicRect(region, value, col = "#08baff", border = "#08baff", lwd=NULL)
                                }
                              })

### c: reads mean coverage barplot   #这里的reads覆盖度是结合了pacbio和二代数据，求滑窗内平均覆盖度
core.bed<-read.table("d:/data/XN21/xn1genome/panGenome/graphgenome/gaf/coreGene/core_graph_len.bed",header = F)
dis.bed<-read.table("d:/data/XN21/xn1genome/panGenome/graphgenome/gaf/coreGene/Dispensable_graph_len.bed",header = F)

all.bed.list<-list(core.bed, 
                   dis.bed)

circos.genomicTrackPlotRegion(all.bed.list, track.height = 0.1,ylim = c(0,1), bg.border="lightgray",
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                  circos.genomicRect(region, value, col = "#4da0a0", border = NA, lwd=NULL)}
                                else{
                                  circos.genomicRect(region, value, col = "#9b3a74", border = NA, lwd=NULL)
                                }
                              })


### b: plot gap and telo #绘制gap和端粒的位置
sp.bed<-read.table( "d:/data/XN21/xn1genome/CG/ortho2/analysis/DM.uniq.family.gene.bed",header = F)[1:3]#读取gap位置的bed文件
colnames(sp.bed)<-c("chr","start","end")
sp.bed$value<-1
#输入telo.bed文件 来自chr.telo.txt 需要awk处理
#awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$4}' chr.telo.txt|awk 'BEGIN{FS=OFS="\t"}{for (f=1; f <= NF; f+=1) {if (f !=1 && $f != "nd") {print $1,$f}}}' > telo.bed
exp.bed<-read.table( "d:/data/XN21/xn1genome/CG/ortho2/analysis/DM.expanded.family.gene.bed",header = F)[1:3]
colnames(exp.bed)<-c("chr","start","end")
exp.bed$value<-1

all.bed.list<-list(sp.bed,exp.bed)


circos.genomicTrackPlotRegion(all.bed.list, track.height = 0.06,ylim = c(0,1), bg.border="lightgray",
                              panel.fun = function(region, value, ...){
                                i=getI(...)
                                if(i == 1){
                                  circos.genomicRect(region, value, col = "black", border = NA, lwd=0.01)}
                                else{
                                  circos.genomicRect(region, value, col = "red", border = NA, lwd=0.01)
                                }
                              })


####plot LTR
ltr<-read.delim("d:/data/XN21/xn1genome/ltr/XN1.Chr.fasta.out.gff3",header = F)
Copia_LTR_retrotransposon<-filter(ltr,V3=="Copia_LTR_retrotransposon")[c(1,4,5)]
Copia.df <- data.frame(chr=Copia_LTR_retrotransposon$V1, start=Copia_LTR_retrotransposon$V4, end=Copia_LTR_retrotransposon$V5)
bed.Copia <- genomicDensity(region = Copia.df, window.size = 10000, overlap = F)###circos.clear() 不然报错
vvs <- bed.Copia$value

circos.genomicTrackPlotRegion(bed.Copia, track.height = 0.05, bg.border = NA, ylim = c(min(vvs),max(vvs)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){
                                
                                circos.genomicLines(region, value, col = "#DE582B", border = NA,lwd = 0.06,
                                                    ybottom = 0, ytop.column = 1)
                              })

Gypsy_LTR_retrotransposon<-filter(ltr,V3=="Gypsy_LTR_retrotransposon")[c(1,4,5)]
Gypsy.df <- data.frame(chr=Gypsy_LTR_retrotransposon$V1, start=Gypsy_LTR_retrotransposon$V4, end=Gypsy_LTR_retrotransposon$V5)
bed.Gypsy <- genomicDensity(region = Gypsy.df, window.size = 10000, overlap = F)###circos.clear() 不然报错
vvs <- bed.Gypsy$value

circos.genomicTrackPlotRegion(bed.Gypsy, track.height = 0.05, bg.border = NA, ylim = c(min(vvs),max(vvs)),
                              # ylim=c(min(bed.gc$value), max(bed.gc$value)),
                              panel.fun = function(region, value, ...){
                                
                                circos.genomicLines(region, value, col = "#808080", border = NA,lwd = 0.06,
                                                    ybottom = 0, ytop.column = 1)
                              })



